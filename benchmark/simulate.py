#!/usr/bin/env python3
"""
Simulate paired-end Illumina reads from a pooled, dual-indexed Nextera library
whose index reads were NOT sequenced.

Each read pair is drawn from a sample of the sample sheet; its insert length is
drawn from a normal distribution. Reads longer than the insert run through into
the adapter, the i7 (R1) or i5 (R2) index and the flow-cell adapter. Substitution
errors are added with a rate that rises linearly along the read (x1 at the start,
x4 at the end), as in MiSeq 2x300 data.

The true sample and insert length are written into each read name:
    @sim<k>_s<sample index>_L<insert length>

Usage:
    python3 simulate.py --samplesheet SampleSheet.csv --pairs 100000 \
        --insert-mean 200 --insert-sd 50 --error 0.005 --read-len 301 \
        --i5-orientation rc --out sim
"""
import argparse
import gzip
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import undeterminator as rd  # noqa: E402

R1_ADAPTER = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
R1_TAIL = "ATCTCGTATGCCGTCTTCTGCTTG"
R2_ADAPTER = "CTGTCTCTTATACACATCTGACGCTGCCGACGA"
R2_TAIL = "GTGTAGATCTCGGTGGTCGCCGTATCATT"
POLY = "G" * 400  # dark cycles after the end of the molecule (2-colour chemistry reads G)


def mutate(seq, base_rate, rng):
    n = len(seq)
    out = list(seq)
    for i in range(n):
        if rng.random() < base_rate * (1 + 3 * i / max(n - 1, 1)):
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def simulate(samples, pairs, mean, sd, err, read_len, i5_orient, out, seed=1):
    rng = random.Random(seed)
    qual = "I" * read_len
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(f"{out}_R1.fastq.gz", "wt", compresslevel=1) as f1, \
            gzip.open(f"{out}_R2.fastq.gz", "wt", compresslevel=1) as f2:
        for k in range(pairs):
            si = rng.randrange(len(samples))
            s = samples[si]
            L = max(30, int(round(rng.gauss(mean, sd))))
            insert = "".join(rng.choice("ACGT") for _ in range(L))
            i5 = rd.revcomp(s["i5"]) if i5_orient == "rc" else s["i5"]
            r1 = (insert + R1_ADAPTER + s["i7"] + R1_TAIL + POLY)[:read_len]
            r2 = (rd.revcomp(insert) + R2_ADAPTER + i5 + R2_TAIL + POLY)[:read_len]
            name = f"@sim{k}_s{si}_L{L}"
            f1.write(f"{name} 1:N:0:1\n{mutate(r1, err, rng)}\n+\n{qual}\n")
            f2.write(f"{name} 2:N:0:1\n{mutate(r2, err, rng)}\n+\n{qual}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--pairs", type=int, default=100000)
    ap.add_argument("--insert-mean", type=float, default=200)
    ap.add_argument("--insert-sd", type=float, default=50)
    ap.add_argument("--error", type=float, default=0.005, help="substitution rate at read start")
    ap.add_argument("--read-len", type=int, default=301)
    ap.add_argument("--i5-orientation", choices=["fw", "rc"], default="rc")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True, help="output prefix")
    a = ap.parse_args()
    samples = rd.read_samplesheet(a.samplesheet)
    simulate(samples, a.pairs, a.insert_mean, a.insert_sd, a.error, a.read_len,
             a.i5_orientation, a.out, a.seed)


if __name__ == "__main__":
    main()
