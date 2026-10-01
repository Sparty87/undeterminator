"""Unit tests for undeterminator (run with: python3 -m pytest tests)."""
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import undeterminator as rd  # noqa: E402

A1 = rd.PRESETS[rd.DEFAULT_PRESET]["r1"]
A2 = rd.PRESETS[rd.DEFAULT_PRESET]["r2"]
SHEET = """[Header]
IEMFileVersion,4
[Data]
Sample_ID,Sample_Name,index,index2
S1,S1,GTAGAGGA,TATCCTCT
S2,S2,GTAGAGGA,GTAAGGAG
S3,S3,GGACTCCT,GTAAGGAG
"""


def cfg(samples, **kw):
    t7, t5 = {}, {}
    for i, s in enumerate(samples):
        t7.setdefault(s["i7"], set()).add(i)
        t5.setdefault(rd.revcomp(s["i5"]), set()).add(i)
    c = {"a1": A1, "a2": A2, "s1": rd.make_seeds(A1, 3), "s2": rd.make_seeds(A2, 3), "amm": 3,
         "l7": 8, "l5": 8, "t7": t7, "t5": t5, "imm": 1, "single": True, "trim": True,
         "write_undet": True}
    c.update(kw)
    return c


def pair(ins, i7, i5, n=301):
    r1 = (ins + A1 + i7 + "ATCTCGTATGCC" + "G" * 400)[:n]
    r2 = (rd.revcomp(ins) + A2 + rd.revcomp(i5) + "GTGTAGATCTCG" + "G" * 400)[:n]
    return ("@a", r1, "I" * n, "@a", r2, "I" * n)


def samples(tmp_path):
    p = tmp_path / "ss.csv"
    p.write_text(SHEET)
    return rd.read_samplesheet(p)


def test_revcomp():
    assert rd.revcomp("ACGTN") == "NACGT"


def test_samplesheet_v2_and_bom(tmp_path):
    p = tmp_path / "v2.csv"
    p.write_bytes("﻿[Header]\r\nFileFormatVersion,2\r\n[BCLConvert_Data]\r\nSample_ID,Index,Index2\r\n"
                  "A b,acgtacgt,TTTTAAAA\r\n".encode())
    s = rd.read_samplesheet(p)
    assert s == [{"id": "A_b", "i7": "ACGTACGT", "i5": "TTTTAAAA"}]


def test_anchor_with_mismatches():
    seq = "N" * 50 + A1[:5] + "A" + A1[6:20] + "T" + A1[21:] + "ACGTACGT" + "G" * 20
    seeds = rd.make_seeds(A1, 3)
    assert rd.find_anchor(seq, A1, seeds, 3, 8) == 50
    assert rd.find_anchor(seq, A1, seeds, 0, 8) == -1


def test_anchor_needs_room_for_index():
    seq = "N" * 50 + A1 + "ACG"
    assert rd.find_anchor(seq, A1, rd.make_seeds(A1, 3), 3, 8) == -1


def test_assignment_dual_index(tmp_path):
    smp = samples(tmp_path)
    rd._init_worker(cfg(smp))
    out, st, _ = rd.process_batch([pair("ACGT" * 30, "GTAGAGGA", "GTAAGGAG")])
    assert st["both"] == 1 and 1 in out            # S2, not S1 (shared i7) nor S3 (shared i5)
    r1 = out[1][0].split("\n")[1]
    assert r1 == "ACGT" * 30                        # adapter and index trimmed


def test_one_mismatch_in_index(tmp_path):
    smp = samples(tmp_path)
    rd._init_worker(cfg(smp))
    out, st, _ = rd.process_batch([pair("ACGT" * 30, "GTAGAGGT", "GTAAGGAG")])
    assert 1 in out


def test_long_insert_not_assigned(tmp_path):
    smp = samples(tmp_path)
    rd._init_worker(cfg(smp))
    out, st, _ = rd.process_batch([pair("ACGT" * 80, "GTAGAGGA", "GTAAGGAG")])
    assert st["no_readthrough"] == 1 and -1 in out


def test_conflicting_combination(tmp_path):
    smp = samples(tmp_path)
    rd._init_worker(cfg(smp))
    _, st, cf = rd.process_batch([pair("ACGT" * 30, "GGACTCCT", "TATCCTCT")])
    assert st["conflict"] == 1 and cf


def test_end_to_end(tmp_path):
    ss = tmp_path / "ss.csv"
    ss.write_text(SHEET)
    smp = rd.read_samplesheet(ss)
    with gzip.open(tmp_path / "p_R1.fastq.gz", "wt") as f1, gzip.open(tmp_path / "p_R2.fastq.gz", "wt") as f2:
        for k in range(300):
            s = smp[k % 3]
            h1, r1, q1, h2, r2, q2 = pair("ACGT" * 25 + "T" * (k % 7), s["i7"], s["i5"])
            f1.write(f"@r{k}\n{r1}\n+\n{q1}\n")
            f2.write(f"@r{k}\n{r2}\n+\n{q2}\n")
    st = rd.run_demux(tmp_path / "p_R1.fastq.gz", tmp_path / "p_R2.fastq.gz", ss, tmp_path / "out",
                      threads=1, log=lambda *a, **k: None)
    assert st["both"] == 300
    assert (tmp_path / "out" / "demux_stats.tsv").exists()
