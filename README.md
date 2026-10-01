# Undeterminator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23083065.svg)](https://doi.org/10.5281/zenodo.23083065)

Undeterminator turns `Undetermined` reads back into samples. It recovers sample identity from **adapter read-through** when an Illumina run was
sequenced **without index reads** (for example after a misparsed sample sheet).

When an insert is shorter than the read, Read 1 continues into the adapter and the
**i7** index, and Read 2 into the adapter and the **i5** index. Undeterminator finds
the adapter anchor with an error-tolerant search, reads the index that follows it,
assigns each read pair to the sample of the sample sheet, trims adapter and index,
and writes one pair of FASTQ files per sample.

```
R1: insert … CTGTCTCTTATACACATCT CCGAGCCCACGAGAC [i7] ATCTCGTATGCC…
R2: insert … CTGTCTCTTATACACATCT GACGCTGCCGACGA  [i5] GTGTAGATCTCG…
```

Only fragments shorter than about *read length − 42 bp* can be recovered
(≈ 259 bp for 2 × 301 cycles, ≈ 109 bp for 2 × 151 cycles with 8-nt indices).

## Features

- Paired-end, dual-index assignment (resolves combinatorial designs where samples share i7 or i5)
- Anchor search tolerant to substitutions (pigeonhole seeds + Hamming verification)
- Automatic detection of index orientation (i5 as written or reverse complement)
- Nextera and TruSeq presets, custom anchors
- Report of every read pair: assigned (dual / single index), ambiguous, conflicting, unknown, no read-through
- Standard library only (Python ≥ 3.8); `pigz` used if present; GUI (Tkinter) and CLI

## Installation

```bash
git clone https://github.com/<user>/undeterminator.git
cd undeterminator
python3 undeterminator.py --help
```

No dependencies beyond Python. For the GUI, Tkinter must be available
(`python3-tk` on Debian/Ubuntu, `python3-tkinter` on RHEL/Fedora).

## Usage

### 1. Make pooled FASTQ files without adapter trimming

Convert the run with bcl2fastq or BCL Convert using a sample sheet that has **one
sample, no index columns and no `Adapter` lines**:

```
[Header]
IEMFileVersion,4
[Reads]
301
301
[Data]
Sample_ID,Sample_Name
Pool,Pool
```

If adapters are trimmed, the index is trimmed with them and nothing can be recovered.

### 2. Demultiplex from read-through

```bash
python3 undeterminator.py \
    --r1 Pool_S1_R1_001.fastq.gz --r2 Pool_S1_R2_001.fastq.gz \
    --samplesheet SampleSheet_original.csv \
    --out undeterminator_out --threads 8
```

or launch `python3 undeterminator.py` without arguments for the graphical interface.

| Option | Default | Meaning |
|---|---|---|
| `--preset` | `nextera` | adapter chemistry (`nextera`, `truseq`) |
| `--anchor1`, `--anchor2` | preset | custom anchor ending right before i7 / i5 |
| `--anchor-mm` | 3 | substitutions allowed in the anchor |
| `--index-mm` | 1 | substitutions allowed in each index |
| `--orient-i7`, `--orient-i5` | `auto` | `fw` = as in sheet, `rc` = reverse complement |
| `--no-single` | off | require both indices for assignment |
| `--no-trim` | off | keep adapter and index in assigned reads |
| `--no-undetermined` | off | do not write unassigned pairs |

### Output

- `<Sample_ID>_R1.fastq.gz`, `<Sample_ID>_R2.fastq.gz` — assigned pairs, adapter trimmed
- `Undetermined_R1/R2.fastq.gz` — pairs not assigned
- `demux_stats.tsv` — pairs per sample, category totals, most frequent unmatched index combinations

## Benchmark

`benchmark/` simulates 2 × 301 Nextera reads for a 24-plex combinatorial dual-index
design and measures recovery, precision, the effect of insert size and error rate,
an ablation and throughput:

```bash
cd benchmark
python3 run_benchmark.py      # results/*.tsv (BENCH_PAIRS=100000 by default)
python3 figure.py             # figure1.png / figure1.pdf
```

## Tests

```bash
python3 -m pytest tests
```

## Citation

See `CITATION.cff`. Manuscript in preparation.

## License

MIT — see `LICENSE`.
