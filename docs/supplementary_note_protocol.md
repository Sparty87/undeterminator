# Supplementary Note — From a run without index cycles to per-sample FASTQ files

## 1. Confirm that the index reads are missing

Open `RunInfo.xml` in the run folder. A run without index cycles lists only the
template reads, for example:

```xml
<Read Number="1" NumCycles="301" IsIndexedRead="N" />
<Read Number="2" NumCycles="301" IsIndexedRead="N" />
```

Count the cycle folders to check that no index cycles exist on disk:

```bash
ls -d <run>/Data/Intensities/BaseCalls/L001/C*.1 | wc -l
```

For a 2 × 301 run the count is 602 without index cycles and 618 with two 8-nt
index reads. If index cycles are on disk but missing from `RunInfo.xml`, add the
index reads to `RunInfo.xml` and demultiplex normally; Undeterminator is not needed.

## 2. Convert the run to pooled FASTQ files without adapter trimming

Use a sample sheet with one sample, no index columns and no `Adapter` settings:

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

```bash
bcl2fastq --runfolder-dir <run> --output-dir pool_fastq \
          --sample-sheet SampleSheet_pool.csv --no-lane-splitting
```

## 3. Check that read-through is present

```bash
zcat pool_fastq/Pool_S1_R1_001.fastq.gz | head -400000 | awk 'NR%4==2{print length($0)}' \
    | sort -n | uniq -c | sort -k1,1nr | head        # all reads full length = not trimmed
zcat pool_fastq/Pool_S1_R1_001.fastq.gz | head -4000000 | awk 'NR%4==2' \
    | grep -c CTGTCTCTTATAC                          # reads containing the Nextera adapter
```

If fewer than about 1 read in 1,000 contains the adapter, the inserts are longer
than the reads and the run cannot be rescued; resequence with index reads enabled.

## 4. Demultiplex from read-through

```bash
python3 undeterminator.py --r1 pool_fastq/Pool_S1_R1_001.fastq.gz \
    --r2 pool_fastq/Pool_S1_R2_001.fastq.gz --samplesheet SampleSheet_original.csv \
    --out undeterminator_out --threads 8
```

## 5. Read the report

`undeterminator_out/demux_stats.tsv` gives pairs per sample and the category of every
pair. A high `conflict` count, or unmatched combinations that differ from sheet
entries by a reverse complement, point to an orientation or sample-sheet error:
re-run with `--orient-i5 rc` (or `fw`) as indicated.

## 6. Use the rescued data

Assigned reads are trimmed of adapter and index and can be aligned directly.
Expect lower and uneven coverage, enriched for short fragments; check coverage at
the loci of interest before reporting results.
