#!/usr/bin/env python3
"""Draw Figure 1 (schematic + recovery by insert length) from results/B_recovery_by_insert.tsv."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
ACCENT, QUIET, LIGHT = "#2a78d6", "#8c8b86", "#e4e3df"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2})


def segment(ax, x, w, y, label, fc, tc=INK, h=0.34):
    ax.add_patch(FancyBboxPatch((x, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.06",
                                fc=fc, ec="white", lw=1.5))
    ax.text(x + w / 2, y, label, ha="center", va="center", color=tc, fontsize=7)


def panel_a(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.5)
    ax.axis("off")
    rows = [(2.55, "Read 1", "CTGTCTCTTATACACATCT", "CCGAGCCCACGAGAC", "i7", "P7"),
            (0.95, "Read 2", "CTGTCTCTTATACACATCT", "GACGCTGCCGACGA", "i5", "P5")]
    for y, name, me, spacer, idx, flow in rows:
        ax.text(0.0, y, name, ha="left", va="center", fontsize=7.5, fontweight="bold", color=INK)
        segment(ax, 0.9, 3.4, y, "insert (≤ 259 bp)", LIGHT)
        segment(ax, 4.3, 1.9, y, f"adapter  {me[:7]}…", "#c9c8c3")
        segment(ax, 6.2, 1.5, y, f"spacer {len(spacer)} nt", "#c9c8c3")
        segment(ax, 7.7, 0.8, y, f"{idx} 8 nt", ACCENT, tc="white")
        segment(ax, 8.5, 1.4, y, f"{flow} adapter", "#f1f0ed", tc=QUIET)
        # read extent
        ax.annotate("", xy=(8.5, y + 0.36), xytext=(0.9, y + 0.36),
                    arrowprops=dict(arrowstyle="->", color=INK2, lw=0.9))
        ax.text(4.7, y + 0.5, "301 cycles", ha="center", va="bottom", fontsize=6.5, color=INK2)
        ax.plot([4.3, 7.7], [y - 0.3, y - 0.3], color=INK2, lw=0.8)
        ax.text(6.0, y - 0.42, f"anchor {len(me) + len(spacer)} nt (≤ 3 mismatches)", ha="center",
                va="top", fontsize=6.5, color=INK2)
    ax.text(-0.05, 3.45, "A", fontsize=11, fontweight="bold", va="top")


def panel_b(ax):
    rows = list(csv.DictReader(open(HERE / "results" / "B_recovery_by_insert.tsv"), delimiter="\t"))
    x = [int(r["insert_bin_bp"]) + 5 for r in rows]
    y = [float(r["assigned_pct"]) for r in rows]
    ok = [float(r["correct_pct"]) for r in rows]
    ax.axvline(259, color=QUIET, lw=0.9, ls=(0, (3, 3)))
    ax.text(268, 55, "limit for 2 × 301 cycles:\n301 − 34 − 8 = 259 bp", fontsize=6.5, color=INK2, va="center")
    ax.step(x, y, where="mid", color=ACCENT, lw=2)
    ax.set_xlim(50, 400)
    ax.set_ylim(0, 102)
    ax.set_xlabel("Insert length (bp)")
    ax.set_ylabel("Read pairs assigned (%)")
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    wrong = sum(a - b for a, b in zip(y, ok))
    msg = "no pair assigned to the wrong sample" if wrong == 0 else f"misassigned: {wrong:.2f}%"
    ax.text(60, 8, msg, fontsize=6.5, color=INK2)
    ax.text(-0.13, 1.06, "B", transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


def main():
    fig = plt.figure(figsize=(6.9, 5.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.2, 1.3], hspace=0.22)
    panel_a(fig.add_subplot(gs[0]))
    panel_b(fig.add_subplot(gs[1]))
    for ext in ("png", "pdf", "tiff"):
        fig.savefig(HERE / f"figure1.{ext}", dpi=600 if ext != "pdf" else None, bbox_inches="tight",
                    facecolor="white")


if __name__ == "__main__":
    main()
