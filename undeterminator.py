#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Undeterminator — recover sample identity from a run without index reads
===========================================================================
When the index cycles were not sequenced, reads from fragments shorter than
the read length run past the insert into the adapter and the index:

  R1:  insert … CTGTCTCTTATACACATCT CCGAGCCCACGAGAC [i7] ATCTCGTATGCC…   (Nextera)
  R2:  insert … CTGTCTCTTATACACATCT GACGCTGCCGACGA  [i5] GTGTAGATCTCG…

The script finds the anchor (adapter + spacer) allowing mismatches, reads the
index that follows it, assigns the read pair to a sample of the sample sheet
and writes one R1/R2 FASTQ pair per sample (adapter trimmed).

Only fragments shorter than about (read length − 42) bp can be recovered, so
each sample's output is partial and enriched in short inserts.

IMPORTANT: the pooled input FASTQ files must be produced WITHOUT adapter
trimming (remove the Adapter lines from the sample sheet used for conversion);
otherwise the part of the read that holds the index has already been cut off.

Usage:
  GUI:  python3 undeterminator.py
  CLI:  python3 undeterminator.py --r1 Pool_R1.fastq.gz --r2 Pool_R2.fastq.gz \
            --samplesheet SampleSheet.csv --out demux_out [--threads 8]
Requires only Python >= 3.8 (tkinter for the GUI). pigz, if installed, speeds up decompression.
"""

import argparse
import csv
import gzip
import io
import itertools
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import sys
import threading
import queue
import time
from collections import Counter
from pathlib import Path

PRESETS = {
    "Nextera / Nextera XT / Illumina DNA Prep": {
        "r1": "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
        "r2": "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
    },
    "TruSeq": {
        "r1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        "r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    },
}
DEFAULT_PRESET = "Nextera / Nextera XT / Illumina DNA Prep"
ADAPTER_CORE = 19          # length of the adapter part of the anchor (trim point)
BATCH = 20000
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


# ── samplesheet ────────────────────────────────────────────────────────────

def read_samplesheet(path):
    raw = Path(path).read_bytes()
    for enc in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    section, header, samples = None, None, []
    for row in csv.reader(io.StringIO(text)):
        cells = [c.strip() for c in row]
        while cells and cells[-1] == "":
            cells.pop()
        if not cells:
            continue
        m = re.match(r"^\[(.+?)\]$", cells[0])
        if m:
            section = m.group(1).lower()
            header = None
            continue
        if section not in ("data", "bclconvert_data"):
            continue
        if header is None:
            header = [h.lower() for h in cells]
            continue
        d = dict(zip(header, cells + [""] * (len(header) - len(cells))))
        sid = d.get("sample_id", "").strip()
        if not sid:
            continue
        sid = re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9_-]", "_", sid)).strip("_")
        samples.append({"id": sid,
                        "i7": re.sub(r"\s", "", d.get("index", "")).upper(),
                        "i5": re.sub(r"\s", "", d.get("index2", "")).upper()})
    if not samples:
        raise RuntimeError("No samples found in the [Data] section of the sample sheet.")
    if not any(s["i7"] for s in samples):
        raise RuntimeError("The sample sheet contains no indices (index column).")
    seen = {}
    for s in samples:            # repeated Sample_IDs -> unique names
        if s["id"] in seen:
            seen[s["id"]] += 1
            s["id"] = f"{s['id']}_{seen[s['id']]}"
        else:
            seen[s["id"]] = 1
    return samples


# ── FASTQ reading / writing ───────────────────────────────────────────────

def open_fastq(path):
    path = str(path)
    if path.endswith(".gz"):
        pigz = shutil.which("pigz")
        if pigz:
            p = subprocess.Popen([pigz, "-dc", "-p", "2", path], stdout=subprocess.PIPE)
            return io.TextIOWrapper(p.stdout, encoding="ascii", errors="replace")
        return gzip.open(path, "rt", encoding="ascii", errors="replace")
    return open(path, "rt", encoding="ascii", errors="replace")


def read_pairs(r1, r2):
    f1, f2 = open_fastq(r1), open_fastq(r2)
    it1, it2 = iter(f1), iter(f2)
    try:
        while True:
            a = list(itertools.islice(it1, 4))
            b = list(itertools.islice(it2, 4))
            if len(a) < 4 or len(b) < 4:
                if a or b:
                    if len(a) != len(b):
                        raise RuntimeError("R1 and R2 contain a different number of reads.")
                return
            yield (a[0].rstrip("\n"), a[1].rstrip("\n"), a[3].rstrip("\n"),
                   b[0].rstrip("\n"), b[1].rstrip("\n"), b[3].rstrip("\n"))
    finally:
        f1.close()
        f2.close()


def batches(gen, n=BATCH):
    while True:
        b = list(itertools.islice(gen, n))
        if not b:
            return
        yield b


# ── anchor search and assignment ──────────────────────────────────────────

def make_seeds(anchor, maxmm):
    k = maxmm + 1                           # pigeonhole principle: one seed is exact
    size = len(anchor) // k
    return [(i * size, anchor[i * size:(i + 1) * size]) for i in range(k)]


def find_anchor(seq, anchor, seeds, maxmm, need):
    p = seq.find(anchor)
    if p >= 0:
        return p if p + len(anchor) + need <= len(seq) else -1
    if maxmm == 0:
        return -1
    L, tried = len(anchor), set()
    for off, sd in seeds:
        start = 0
        while True:
            q = seq.find(sd, start)
            if q < 0:
                break
            start = q + 1
            p = q - off
            if p < 0 or p in tried or p + L + need > len(seq):
                continue
            tried.add(p)
            if hamming(seq[p:p + L], anchor) <= maxmm:
                return p
    return -1


def extract(seq, cfg, which):
    anchor = cfg["a1"] if which == 1 else cfg["a2"]
    seeds = cfg["s1"] if which == 1 else cfg["s2"]
    L = cfg["l7"] if which == 1 else cfg["l5"]
    if L == 0:
        return -1, ""
    p = find_anchor(seq, anchor, seeds, cfg["amm"], L)
    if p < 0:
        return -1, ""
    s = p + len(anchor)
    return p, seq[s:s + L]


def match(obs, table, mm):
    """table: {index_sequence: set(samples)} -> samples at the minimum mismatch count (<= mm)."""
    if not obs:
        return None
    hit = table.get(obs)
    if hit is not None:
        return hit
    best, res = mm + 1, set()
    for seq, ids in table.items():
        d = hamming(obs, seq)
        if d < best:
            best, res = d, set(ids)
        elif d == best:
            res |= ids
    return res if best <= mm else set()


_CFG = None


def _init_worker(cfg):
    global _CFG
    _CFG = cfg


def process_batch(batch):
    cfg = _CFG
    out = {}
    stats = Counter()
    conflicts = Counter()
    for h1, s1, q1, h2, s2, q2 in batch:
        p1, o7 = extract(s1, cfg, 1)
        p2, o5 = extract(s2, cfg, 2)
        c7 = match(o7, cfg["t7"], cfg["imm"]) if o7 else None
        c5 = match(o5, cfg["t5"], cfg["imm"]) if o5 else None
        dest, how = -1, "no_readthrough"
        if c7 is None and c5 is None:
            pass
        elif c7 is not None and c5 is not None:
            both = c7 & c5
            if len(both) == 1:
                dest, how = next(iter(both)), "both"
            elif len(both) > 1:
                how = "ambiguous"
            else:
                how = "conflict"
                conflicts[(o7, o5)] += 1
        else:
            c = c7 if c7 is not None else c5
            if cfg["single"] and len(c) == 1:
                dest, how = next(iter(c)), "single"
            else:
                how = "ambiguous" if c and len(c) > 1 else "unknown_index"
                if not c:
                    conflicts[(o7 or "-", o5 or "-")] += 1
        stats[how] += 1
        if dest >= 0:
            stats[("sample", dest, how)] += 1
        if dest < 0 and not cfg["write_undet"]:
            continue
        if cfg["trim"] and dest >= 0:
            if p1 >= 0:
                s1, q1 = s1[:p1], q1[:p1]
            if p2 >= 0:
                s2, q2 = s2[:p2], q2[:p2]
            if not s1:
                s1, q1 = "N", "#"
            if not s2:
                s2, q2 = "N", "#"
        r1, r2 = out.setdefault(dest, ([], []))
        r1.append(f"{h1}\n{s1}\n+\n{q1}\n")
        r2.append(f"{h2}\n{s2}\n+\n{q2}\n")
    return {k: ("".join(v[0]), "".join(v[1])) for k, v in out.items()}, stats, conflicts


# ── index orientation ─────────────────────────────────────────────────────

def detect_orientation(r1, r2, samples, cfg_base, n=200000, log=print):
    fw7 = {s["i7"] for s in samples if s["i7"]}
    fw5 = {s["i5"] for s in samples if s["i5"]}
    rc7, rc5 = {revcomp(x) for x in fw7}, {revcomp(x) for x in fw5}
    cnt = Counter()
    pairs = 0
    lens = Counter()
    for h1, s1, q1, h2, s2, q2 in itertools.islice(read_pairs(r1, r2), n):
        pairs += 1
        lens[len(s1)] += 1
        _, o7 = extract(s1, cfg_base, 1)
        _, o5 = extract(s2, cfg_base, 2)
        if o7:
            cnt["seen7"] += 1
            cnt["fw7"] += o7 in fw7
            cnt["rc7"] += o7 in rc7
        if o5:
            cnt["seen5"] += 1
            cnt["fw5"] += o5 in fw5
            cnt["rc5"] += o5 in rc5
    if not pairs:
        raise RuntimeError("The FASTQ files are empty.")
    maxlen = max(lens)
    short = sum(v for k, v in lens.items() if k < maxlen) / pairs
    log(f"Sample of {pairs:,} pairs: anchor found in R1 {cnt['seen7']:,} "
        f"({100 * cnt['seen7'] / pairs:.1f}%), in R2 {cnt['seen5']:,} ({100 * cnt['seen5'] / pairs:.1f}%).")
    if cnt["seen7"] + cnt["seen5"] < 0.001 * pairs:
        if short > 0.2:
            log("⚠ Almost no anchors found and many reads are shortened: the FASTQ files look "
                "adapter-trimmed. Regenerate them without the Adapter lines in the sample sheet.", "warn")
        else:
            log("⚠ Almost no anchors found: fragments too long or wrong adapter preset.", "warn")
    o7 = "rc" if cnt["rc7"] > cnt["fw7"] else "fw"
    o5 = "rc" if cnt["rc5"] > cnt["fw5"] else "fw"
    log(f"i7: as in sample sheet {cnt['fw7']:,} vs reverse complement {cnt['rc7']:,} -> "
        f"{'reverse complement' if o7 == 'rc' else 'as in sample sheet'}")
    if fw5:
        log(f"i5: as in sample sheet {cnt['fw5']:,} vs reverse complement {cnt['rc5']:,} -> "
            f"{'reverse complement' if o5 == 'rc' else 'as in sample sheet'}")
    return o7, o5


# ── main engine ───────────────────────────────────────────────────────────

def run_demux(r1, r2, sheet, outdir, preset=DEFAULT_PRESET, anchor1=None, anchor2=None,
              anchor_mm=3, index_mm=1, single=True, trim=True, write_undet=True,
              threads=4, orient7="auto", orient5="auto", log=print, stop=None):
    t0 = time.time()
    samples = read_samplesheet(sheet)
    a1 = (anchor1 or PRESETS[preset]["r1"]).upper()
    a2 = (anchor2 or PRESETS[preset]["r2"]).upper()
    anchor_mm = max(0, min(int(anchor_mm), 5))
    l7 = max(len(s["i7"]) for s in samples)
    l5 = max([len(s["i5"]) for s in samples] + [0])
    log(f"{len(samples)} samples; i7 index {l7} bp, i5 index {l5} bp.")

    base = {"a1": a1, "a2": a2, "s1": make_seeds(a1, anchor_mm), "s2": make_seeds(a2, anchor_mm),
            "amm": anchor_mm, "l7": l7, "l5": l5}
    if orient7 == "auto" or orient5 == "auto":
        log("Detecting index orientation…")
        d7, d5 = detect_orientation(r1, r2, samples, base, log=log)
        orient7 = d7 if orient7 == "auto" else orient7
        orient5 = d5 if orient5 == "auto" else orient5

    t7, t5 = {}, {}
    for i, s in enumerate(samples):
        if s["i7"]:
            t7.setdefault(revcomp(s["i7"]) if orient7 == "rc" else s["i7"], set()).add(i)
        if s["i5"]:
            t5.setdefault(revcomp(s["i5"]) if orient5 == "rc" else s["i5"], set()).add(i)
    if not t5:
        base["l5"] = 0
    # warn about index pairs that collide with the chosen mismatch tolerance
    combos = [(s["i7"], s["i5"]) for s in samples]
    for a in range(len(combos)):
        for b in range(a + 1, len(combos)):
            d7 = hamming(combos[a][0], combos[b][0])
            d5 = hamming(combos[a][1], combos[b][1]) if combos[a][1] else 0
            if d7 <= 2 * index_mm and (not t5 or d5 <= 2 * index_mm) and (d7 + d5) > 0 and index_mm:
                log(f"⚠ {samples[a]['id']} and {samples[b]['id']} are close at {index_mm} mismatch(es): "
                    "equidistant pairs will be discarded as ambiguous.", "warn")

    cfg = dict(base, t7=t7, t5=t5, imm=int(index_mm), single=bool(single), trim=bool(trim),
               write_undet=bool(write_undet))
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    names = {i: s["id"] for i, s in enumerate(samples)}
    names[-1] = "Undetermined"
    handles = {}

    def fh(dest):
        if dest not in handles:
            n = names[dest]
            handles[dest] = (gzip.open(outdir / f"{n}_R1.fastq.gz", "wt", compresslevel=3),
                             gzip.open(outdir / f"{n}_R2.fastq.gz", "wt", compresslevel=3))
        return handles[dest]

    stats, conflicts = Counter(), Counter()
    total = 0
    last = time.time()
    threads = max(1, int(threads))
    log(f"Demultiplexing with {threads} process(es)…")
    pool = mp.Pool(threads, initializer=_init_worker, initargs=(cfg,)) if threads > 1 else None
    if pool is None:
        _init_worker(cfg)
    try:
        it = batches(read_pairs(r1, r2))
        results = pool.imap(process_batch, it) if pool else map(process_batch, it)
        for out, st, cf in results:
            if stop is not None and stop.is_set():
                raise RuntimeError("Stopped by the user.")
            for dest, (t1, t2) in out.items():
                f1, f2 = fh(dest)
                f1.write(t1)
                f2.write(t2)
            stats.update(st)
            conflicts.update(cf)
            total += sum(v for k, v in st.items() if isinstance(k, str))
            if time.time() - last > 5:
                last = time.time()
                ok = stats["both"] + stats["single"]
                log(f"  {total:,} pairs — assigned {ok:,} ({100 * ok / total:.2f}%) — "
                    f"{total / (time.time() - t0):,.0f} pairs/s")
    finally:
        if pool:
            pool.terminate()
        for f1, f2 in handles.values():
            f1.close()
            f2.close()

    # summary
    ok = stats["both"] + stats["single"]
    rep = outdir / "demux_stats.tsv"
    with rep.open("w") as f:
        f.write("Sample_ID\tindex\tindex2\tpairs_dual_index\tpairs_single_index\ttotal\tpct_of_all_pairs\n")
        for i, s in enumerate(samples):
            b, sg = stats[("sample", i, "both")], stats[("sample", i, "single")]
            f.write(f"{s['id']}\t{s['i7']}\t{s['i5']}\t{b}\t{sg}\t{b + sg}\t{100 * (b + sg) / max(total, 1):.3f}\n")
        f.write("\ncategory\tpairs\n")
        for k in ("both", "single", "ambiguous", "conflict", "unknown_index", "no_readthrough"):
            f.write(f"{k}\t{stats[k]}\n")
        f.write(f"total\t{total}\n")
        if conflicts:
            f.write("\nobserved_combination_not_in_sample_sheet (i7+i5)\tpairs\n")
            for (x, y), n in conflicts.most_common(30):
                f.write(f"{x}+{y}\t{n}\n")

    log("")
    log("── Summary ──", "head")
    w = max(len(s["id"]) for s in samples)
    for i, s in enumerate(samples):
        n = stats[("sample", i, "both")] + stats[("sample", i, "single")]
        log(f"  {s['id']:<{w}}  {n:>10,}", "warn" if n == 0 else None)
    labels = {"both": "assigned (i7+i5)", "single": "assigned (single index)",
              "ambiguous": "ambiguous", "conflict": "incompatible i7+i5",
              "unknown_index": "unknown index", "no_readthrough": "no read-through"}
    for k, lab in labels.items():
        log(f"  {lab:<28} {stats[k]:>12,}  {100 * stats[k] / max(total, 1):6.2f}%")
    log(f"  {'total':<28} {total:>12,}")
    if conflicts:
        log("Most frequent combinations not in the sample sheet:")
        for (x, y), n in conflicts.most_common(5):
            log(f"  {x}+{y}  {n:,}")
    log(f"✔ Done in {time.time() - t0:.0f} s. Recovered {ok:,} of {total:,} pairs "
        f"({100 * ok / max(total, 1):.2f}%). Output and statistics in {outdir}", "ok")
    return stats


# ── CLI ───────────────────────────────────────────────────────────────────

def cli(argv):
    ap = argparse.ArgumentParser(description="Demultiplexing from adapter read-through "
                                             "(runs without index reads).")
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--preset", default="nextera", choices=["nextera", "truseq"])
    ap.add_argument("--anchor1", help="custom R1 anchor (adapter + spacer immediately before i7)")
    ap.add_argument("--anchor2", help="custom R2 anchor (adapter + spacer immediately before i5)")
    ap.add_argument("--anchor-mm", type=int, default=3)
    ap.add_argument("--index-mm", type=int, default=1)
    ap.add_argument("--no-single", action="store_true", help="always require both indices")
    ap.add_argument("--no-trim", action="store_true", help="do not trim adapter and index from assigned reads")
    ap.add_argument("--no-undetermined", action="store_true")
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--orient-i7", default="auto", choices=["auto", "fw", "rc"])
    ap.add_argument("--orient-i5", default="auto", choices=["auto", "fw", "rc"])
    a = ap.parse_args(argv)
    preset = DEFAULT_PRESET if a.preset == "nextera" else "TruSeq"
    run_demux(a.r1, a.r2, a.samplesheet, a.out, preset, a.anchor1, a.anchor2, a.anchor_mm,
              a.index_mm, not a.no_single, not a.no_trim, not a.no_undetermined, a.threads,
              a.orient_i7, a.orient_i5, log=lambda m, t=None: print(m, flush=True))


# ── GUI ───────────────────────────────────────────────────────────────────

def gui():
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    from tkinter.scrolledtext import ScrolledText

    root = tk.Tk()
    root.title("Undeterminator — recover samples from runs without index reads")
    root.geometry("980x760")
    try:
        ttk.Style().theme_use("clam")
    except tk.TclError:
        pass
    q = queue.Queue()
    stop = threading.Event()
    state = {"busy": False}

    v = {k: tk.StringVar() for k in ("r1", "r2", "sheet", "out")}
    preset = tk.StringVar(value=DEFAULT_PRESET)
    a1 = tk.StringVar(value=PRESETS[DEFAULT_PRESET]["r1"])
    a2 = tk.StringVar(value=PRESETS[DEFAULT_PRESET]["r2"])
    amm, imm = tk.IntVar(value=3), tk.IntVar(value=1)
    thr = tk.IntVar(value=max(1, (os.cpu_count() or 2) - 1))
    single, trim, undet = tk.BooleanVar(value=True), tk.BooleanVar(value=True), tk.BooleanVar(value=True)
    o7, o5 = tk.StringVar(value="auto"), tk.StringVar(value="auto")

    main = ttk.Frame(root, padding=8)
    main.pack(fill="both", expand=True)
    ttk.Label(main, wraplength=940, foreground="#b26a00", text=(
        "The pooled FASTQ files must be produced WITHOUT adapter trimming (remove the 'Adapter' lines "
        "from the sample sheet used for conversion). Only short fragments can be recovered.")).pack(
        fill="x", pady=(0, 6))

    io_f = ttk.LabelFrame(main, text="Files", padding=6)
    io_f.pack(fill="x", pady=3)

    def pick_r1():
        f = filedialog.askopenfilename(title="Pooled FASTQ R1",
                                       filetypes=[("FASTQ", "*.fastq.gz *.fq.gz *.fastq *.fq"), ("All files", "*")])
        if not f:
            return
        v["r1"].set(f)
        cand = re.sub(r"_R1([_.])", r"_R2\1", f, count=1)
        if cand != f and Path(cand).exists():
            v["r2"].set(cand)
        if not v["out"].get():
            v["out"].set(str(Path(f).parent / "undeterminator_out"))

    def pick_file(key, title, types):
        f = filedialog.askopenfilename(title=title, filetypes=types)
        if f:
            v[key].set(f)

    def pick_out():
        d = filedialog.askdirectory(title="Output folder")
        if d:
            v["out"].set(d)

    rows = [("FASTQ R1:", "r1", pick_r1),
            ("FASTQ R2:", "r2", lambda: pick_file("r2", "Pooled FASTQ R2",
                                                  [("FASTQ", "*.fastq.gz *.fq.gz *.fastq *.fq"), ("All files", "*")])),
            ("Samplesheet:", "sheet", lambda: pick_file("sheet", "Samplesheet", [("CSV", "*.csv"), ("All files", "*")])),
            ("Output:", "out", pick_out)]
    for r, (lab, key, cmd) in enumerate(rows):
        ttk.Label(io_f, text=lab, width=13).grid(row=r, column=0, sticky="w")
        ttk.Entry(io_f, textvariable=v[key]).grid(row=r, column=1, sticky="we", padx=4, pady=2)
        ttk.Button(io_f, text="Browse…", command=cmd).grid(row=r, column=2)
    io_f.columnconfigure(1, weight=1)

    of = ttk.LabelFrame(main, text="Parameters", padding=6)
    of.pack(fill="x", pady=3)
    ttk.Label(of, text="Adapter chemistry:").grid(row=0, column=0, sticky="w")
    cb = ttk.Combobox(of, textvariable=preset, values=list(PRESETS), state="readonly", width=42)
    cb.grid(row=0, column=1, columnspan=3, sticky="w")

    def on_preset(*_):
        a1.set(PRESETS[preset.get()]["r1"])
        a2.set(PRESETS[preset.get()]["r2"])
    cb.bind("<<ComboboxSelected>>", on_preset)
    ttk.Label(of, text="R1 anchor (→ i7):").grid(row=1, column=0, sticky="w")
    ttk.Entry(of, textvariable=a1, font=("Monospace", 9)).grid(row=1, column=1, columnspan=3, sticky="we")
    ttk.Label(of, text="R2 anchor (→ i5):").grid(row=2, column=0, sticky="w")
    ttk.Entry(of, textvariable=a2, font=("Monospace", 9)).grid(row=2, column=1, columnspan=3, sticky="we")
    ttk.Label(of, text="Anchor mismatches:").grid(row=3, column=0, sticky="w")
    ttk.Spinbox(of, from_=0, to=5, width=5, textvariable=amm, state="readonly").grid(row=3, column=1, sticky="w")
    ttk.Label(of, text="Index mismatches:").grid(row=3, column=2, sticky="e")
    ttk.Spinbox(of, from_=0, to=2, width=5, textvariable=imm, state="readonly").grid(row=3, column=3, sticky="w")
    ttk.Label(of, text="i7 orientation:").grid(row=4, column=0, sticky="w")
    ttk.Combobox(of, textvariable=o7, values=["auto", "fw", "rc"], state="readonly", width=6).grid(
        row=4, column=1, sticky="w")
    ttk.Label(of, text="i5 orientation:").grid(row=4, column=2, sticky="e")
    ttk.Combobox(of, textvariable=o5, values=["auto", "fw", "rc"], state="readonly", width=6).grid(
        row=4, column=3, sticky="w")
    ttk.Label(of, text="Processes:").grid(row=5, column=0, sticky="w")
    ttk.Spinbox(of, from_=1, to=128, width=5, textvariable=thr).grid(row=5, column=1, sticky="w")
    ttk.Checkbutton(of, text="Assign with a single index when it is unique", variable=single).grid(
        row=6, column=0, columnspan=2, sticky="w")
    ttk.Checkbutton(of, text="Trim adapter from assigned reads", variable=trim).grid(
        row=6, column=2, columnspan=2, sticky="w")
    ttk.Checkbutton(of, text="Write unassigned reads (Undetermined)", variable=undet).grid(
        row=7, column=0, columnspan=2, sticky="w")
    of.columnconfigure(1, weight=1)
    of.columnconfigure(3, weight=1)

    bf = ttk.Frame(main)
    bf.pack(fill="x", pady=4)
    btn_run = ttk.Button(bf, text="Run")
    btn_stop = ttk.Button(bf, text="Stop", state="disabled", command=lambda: stop.set())
    btn_run.pack(side="left", padx=4)
    btn_stop.pack(side="left", padx=4)
    pb = ttk.Progressbar(main, mode="indeterminate")
    pb.pack(fill="x", pady=3)
    txt = ScrolledText(main, font=("Monospace", 9), wrap="word")
    txt.pack(fill="both", expand=True)
    for tag, col in (("err", "#c62828"), ("warn", "#b26a00"), ("ok", "#2e7d32"), ("head", "#000")):
        txt.tag_configure(tag, foreground=col)

    def log(m, tag=None):
        q.put(("log", m, tag))

    def poll():
        try:
            while True:
                item = q.get_nowait()
                if item[0] == "log":
                    txt.insert("end", item[1] + "\n", item[2] or ())
                    txt.see("end")
                elif item[0] == "done":
                    state["busy"] = False
                    pb.stop()
                    btn_run.configure(state="normal")
                    btn_stop.configure(state="disabled")
                    if item[1]:
                        messagebox.showinfo("Undeterminator", "Done: the summary is in the log.")
        except queue.Empty:
            pass
        root.after(150, poll)

    def start():
        for k, lab in (("r1", "FASTQ R1"), ("r2", "FASTQ R2"), ("sheet", "samplesheet")):
            if not Path(v[k].get()).is_file():
                messagebox.showerror("Undeterminator", f"Select the {lab} file.")
                return
        if not v["out"].get().strip():
            messagebox.showerror("Undeterminator", "Choose an output folder.")
            return
        stop.clear()
        state["busy"] = True
        btn_run.configure(state="disabled")
        btn_stop.configure(state="normal")
        pb.start(12)
        args = dict(r1=v["r1"].get(), r2=v["r2"].get(), sheet=v["sheet"].get(), outdir=v["out"].get(),
                    preset=preset.get(), anchor1=a1.get().strip(), anchor2=a2.get().strip(),
                    anchor_mm=amm.get(), index_mm=imm.get(), single=single.get(), trim=trim.get(),
                    write_undet=undet.get(), threads=thr.get(), orient7=o7.get(), orient5=o5.get(),
                    log=log, stop=stop)

        def worker():
            ok = False
            try:
                run_demux(**args)
                ok = True
            except Exception as e:
                log(f"ERROR: {e}", "err")
            q.put(("done", ok))
        threading.Thread(target=worker, daemon=True).start()

    btn_run.configure(command=start)
    root.after(150, poll)
    root.mainloop()


if __name__ == "__main__":
    mp.freeze_support()
    if len(sys.argv) > 1:
        cli(sys.argv[1:])
    else:
        try:
            gui()
        except ImportError:
            print("tkinter is not available: use the command-line mode (--help).")
        except Exception as e:
            if e.__class__.__name__ == "TclError":
                print(f"Cannot open the window ({e}). Use 'ssh -X' or the command-line mode (--help).")
            else:
                raise
