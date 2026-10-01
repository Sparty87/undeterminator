#!/usr/bin/env python3
"""
Benchmark undeterminator on simulated data (see simulate.py).

Experiments
  A  recovery and precision over insert-size means x error rates (read length 301)
  B  recovery as a function of insert length (uniform inserts 50-400 bp)
  C  ablation: full method vs exact matching vs i7-only sample sheet
     (results/C2_tolerance_1pct_error.tsv: anchor/index tolerance sweep on the same dataset)
  D  automatic detection of index orientation (i5 as-is vs reverse complement)
  E  throughput (read pairs per second) by number of worker processes

Results go to results/*.tsv; figure.py draws Figure 1 from them.
"""
import csv
import gzip
import os
import platform
import re
import shutil
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import undeterminator as rd  # noqa: E402
from simulate import simulate  # noqa: E402

SHEET = HERE / "example_samplesheet_24plex.csv"
SHEET_I7 = HERE / "example_samplesheet_24plex_i7only.csv"
WORK = HERE / "work"
RES = HERE / "results"
NAME = re.compile(r"^@sim(\d+)_s(\d+)_L(\d+)")
READ_LEN = 301
PAIRS = int(os.environ.get("BENCH_PAIRS", "100000"))
THREADS = max(1, (os.cpu_count() or 2))
quiet = lambda *a, **k: None  # noqa: E731


def truth_from_fastq(prefix):
    t = {}
    with gzip.open(f"{prefix}_R1.fastq.gz", "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                m = NAME.match(line)
                t[int(m.group(1))] = (int(m.group(2)), int(m.group(3)))
    return t


def evaluate(outdir, samples):
    """Return {read_id: assigned sample index} from the per-sample output files."""
    got = {}
    for i, s in enumerate(samples):
        p = Path(outdir) / f"{s['id']}_R1.fastq.gz"
        if not p.exists():
            continue
        with gzip.open(p, "rt") as f:
            for j, line in enumerate(f):
                if j % 4 == 0:
                    got[int(NAME.match(line).group(1))] = i
    return got


def demux(prefix, outdir, sheet=SHEET, **kw):
    if Path(outdir).exists():
        shutil.rmtree(outdir)
    t0 = time.time()
    st = rd.run_demux(f"{prefix}_R1.fastq.gz", f"{prefix}_R2.fastq.gz", sheet, outdir,
                      write_undet=False, threads=kw.pop("threads", THREADS), log=quiet, **kw)
    return st, time.time() - t0


def recoverable(L, samples):
    """Both anchors + both indices fit inside the read."""
    l7 = max(len(s["i7"]) for s in samples)
    return L + len(rd.PRESETS[rd.DEFAULT_PRESET]["r1"]) + l7 <= READ_LEN


def scores(truth, got, samples):
    n = len(truth)
    assigned = len(got)
    correct = sum(1 for k, v in got.items() if truth[k][0] == v)
    rec_ids = [k for k, (s, L) in truth.items() if recoverable(L, samples)]
    rec_ok = sum(1 for k in rec_ids if got.get(k) == truth[k][0])
    return {"pairs": n, "assigned": assigned, "correct": correct,
            "recovery_pct": 100 * assigned / n,
            "precision_pct": 100 * correct / assigned if assigned else float("nan"),
            "theoretical_pct": 100 * len(rec_ids) / n,
            "sensitivity_recoverable_pct": 100 * rec_ok / len(rec_ids) if rec_ids else float("nan")}


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in r.items()})


def main():
    samples = rd.read_samplesheet(SHEET)
    WORK.mkdir(exist_ok=True)

    # A — insert-size means x error rates
    rows = []
    for mean in (150, 200, 250, 300, 350):
        for err in (0.001, 0.005, 0.01):
            pre = WORK / f"A_m{mean}_e{err}"
            simulate(samples, PAIRS, mean, 50, err, READ_LEN, "rc", pre, seed=mean + int(err * 1e4))
            demux(pre, WORK / "outA")
            r = scores(truth_from_fastq(pre), evaluate(WORK / "outA", samples), samples)
            rows.append({"insert_mean": mean, "insert_sd": 50, "error_rate": err, **r})
            print("A", mean, err, {k: round(v, 2) if isinstance(v, float) else v for k, v in r.items()}, flush=True)
            for f in WORK.glob(f"{pre.name}_R*.fastq.gz"):
                f.unlink()
    write_tsv(RES / "A_insert_error_grid.tsv", rows)

    # B — recovery by insert length (uniform 50–400 bp approximated by a wide normal, clipped)
    pre = WORK / "B_uniform"
    simulate_uniform(samples, PAIRS * 2, 50, 400, 0.005, pre)
    demux(pre, WORK / "outB")
    truth = truth_from_fastq(pre)
    got = evaluate(WORK / "outB", samples)
    bins = defaultdict(lambda: [0, 0, 0])
    for k, (s, L) in truth.items():
        b = (L // 10) * 10
        bins[b][0] += 1
        if k in got:
            bins[b][1] += 1
            bins[b][2] += got[k] == s
    write_tsv(RES / "B_recovery_by_insert.tsv",
              [{"insert_bin_bp": b, "pairs": v[0], "assigned_pct": 100 * v[1] / v[0],
                "correct_pct": 100 * v[2] / v[0]} for b, v in sorted(bins.items())])
    print("B done", flush=True)

    # C — ablation on one dataset (mean 200, error 1%)
    pre = WORK / "C_m200_e0.01"
    simulate(samples, PAIRS, 200, 50, 0.01, READ_LEN, "rc", pre, seed=7)
    truth = truth_from_fastq(pre)
    rows = []
    for label, kw, sheet in (("full (anchor<=3 mm, index<=1 mm, dual index)", {}, SHEET),
                             ("exact anchor and index", {"anchor_mm": 0, "index_mm": 0}, SHEET),
                             ("i7 only (no i5)", {}, SHEET_I7)):
        demux(pre, WORK / "outC", sheet=sheet, **kw)
        smp = rd.read_samplesheet(sheet)
        r = scores(truth, evaluate(WORK / "outC", smp), samples)
        rows.append({"configuration": label, **r})
        print("C", label, round(r["recovery_pct"], 2), round(r["precision_pct"], 2), flush=True)
    write_tsv(RES / "C_ablation.tsv", rows)

    # D — orientation detection
    rows = []
    for orient in ("fw", "rc"):
        pre = WORK / f"D_{orient}"
        simulate(samples, 50000, 200, 50, 0.005, READ_LEN, orient, pre, seed=11)
        base = {"a1": rd.PRESETS[rd.DEFAULT_PRESET]["r1"], "a2": rd.PRESETS[rd.DEFAULT_PRESET]["r2"],
                "amm": 3, "l7": 8, "l5": 8}
        base["s1"], base["s2"] = rd.make_seeds(base["a1"], 3), rd.make_seeds(base["a2"], 3)
        d7, d5 = rd.detect_orientation(f"{pre}_R1.fastq.gz", f"{pre}_R2.fastq.gz", samples, base, log=quiet)
        rows.append({"simulated_i5": orient, "detected_i7": d7, "detected_i5": d5,
                     "correct": d5 == orient and d7 == "fw"})
    write_tsv(RES / "D_orientation.tsv", rows)
    print("D", rows, flush=True)

    # E — throughput
    rows = []
    pre = WORK / "C_m200_e0.01"
    for t in sorted({1, THREADS}):
        _, dt = demux(pre, WORK / "outE", threads=t)
        rows.append({"processes": t, "pairs": PAIRS, "seconds": dt, "pairs_per_second": PAIRS / dt,
                     "cpu": platform.processor() or platform.machine(), "python": platform.python_version()})
    write_tsv(RES / "E_throughput.tsv", rows)
    print("E", rows, flush=True)


def simulate_uniform(samples, pairs, lo, hi, err, out, seed=3):
    """Uniform insert lengths: reuse simulate() one length at a time is slow, so inline it."""
    import random
    from simulate import R1_ADAPTER, R1_TAIL, R2_ADAPTER, R2_TAIL, POLY, mutate
    rng = random.Random(seed)
    qual = "I" * READ_LEN
    with gzip.open(f"{out}_R1.fastq.gz", "wt", compresslevel=1) as f1, \
            gzip.open(f"{out}_R2.fastq.gz", "wt", compresslevel=1) as f2:
        for k in range(pairs):
            si = rng.randrange(len(samples))
            s = samples[si]
            L = rng.randint(lo, hi)
            ins = "".join(rng.choice("ACGT") for _ in range(L))
            r1 = (ins + R1_ADAPTER + s["i7"] + R1_TAIL + POLY)[:READ_LEN]
            r2 = (rd.revcomp(ins) + R2_ADAPTER + rd.revcomp(s["i5"]) + R2_TAIL + POLY)[:READ_LEN]
            f1.write(f"@sim{k}_s{si}_L{L} 1:N:0:1\n{mutate(r1, err, rng)}\n+\n{qual}\n")
            f2.write(f"@sim{k}_s{si}_L{L} 2:N:0:1\n{mutate(r2, err, rng)}\n+\n{qual}\n")


if __name__ == "__main__":
    main()
