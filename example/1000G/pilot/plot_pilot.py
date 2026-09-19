#!/usr/bin/env python3
"""Figure for the pilot report (needs matplotlib): python example/1000G/pilot/plot_pilot.py

Four panels, one question each:
  a  does the library GC curve differ between samples and between library generations?
  b  what does the 45S unit look like after the GC model (the window efficiencies)?
  c  do two independent libraries of the same DNA source give the same copy number?
  d  does sequence of known copy number come out at its known copy number?
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ngsdose import cohort, gcmodel, io, resources, estimate  # noqa: E402
from evaluate_pilot import SEX, run  # noqa: E402

# categorical slots 1-3 of the reference palette (validated all-pairs), chart chrome and ink
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"

plt.rcParams.update({
    "font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"], "font.size": 9, "text.color": INK,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.spines.top": False, "axes.spines.right": False, "axes.facecolor": SURFACE,
    "figure.facecolor": SURFACE, "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2,
    "ytick.labelcolor": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False,
})


def curves(counts_dir):
    out = []
    for f in sorted((HERE / counts_dir).glob("*.json.gz")):
        c = io.load_counts(f)
        t = estimate.nearest_table(c, None)
        out.append(gcmodel.fit_gc_curve(t["n"], t["o"], t["l"]).relative())
    return out


def main():
    bundle = resources.Bundle(ROOT / "resources" / "GRCh38")
    panel, units, feats = io.load_panel(bundle.panel), bundle.units(), bundle.features()
    ny, hg = run("counts_nygc", bundle, panel, units, feats), run("counts_replicates", bundle, panel, units, feats)
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.6), constrained_layout=True)
    g = np.arange(101)

    # ---- a: GC curves
    a = ax[0, 0]
    for col, d in ((BLUE, "counts_nygc"), (ORANGE, "counts_replicates")):
        for y in curves(d):
            a.plot(g, y, color=col, lw=1.1, alpha=0.75)
    a.axhline(1, color=AXIS, lw=0.8)
    a.set_xlim(15, 85)
    a.set_ylim(0.6, 1.5)
    a.set_xlabel("fragment GC (%)")
    a.set_ylabel("5′-end rate relative to genome mean")
    a.set_title("a  Library GC behaviour, fitted per sample on single-copy controls")
    for yy, col, lab in ((1.46, BLUE, f"NYGC, NovaSeq 2×150, 2019 (n={len(ny)})"),
                         (1.405, ORANGE, f"older libraries, HiSeq 2×100 / 2×126 (n={len(hg)})")):
        a.plot([30, 32.5], [yy, yy], color=col, lw=2)
        a.text(33.5, yy, lab, color=INK, fontsize=8.5, va="center")

    # ---- b: window efficiencies along the unit
    b = ax[0, 1]
    cal = cohort.calibrate(list(ny.values()), "rDNA45S")
    x = (cal.window_start + 125) / 1000
    eff = np.exp(cal.a)
    b.axhline(1, color=AXIS, lw=0.8)
    b.vlines(x, 1, eff, color=BLUE, lw=1.0, alpha=0.55)
    b.scatter(x[cal.anchor], eff[cal.anchor], s=9, color=BLUE, zorder=3, label="anchor window (fragment GC 40–60%)")
    b.scatter(x[~cal.anchor], eff[~cal.anchor], s=9, facecolor=SURFACE, edgecolor=BLUE, linewidth=0.9, zorder=3, label="other window")
    for name, s0, e0 in feats["rDNA45S"]:
        if name in ("18S", "28S", "5.8S"):
            b.plot([s0 / 1000, e0 / 1000], [1.52, 1.52], color=INK2, lw=3, solid_capstyle="butt")
            if name != "5.8S":                                   # 157 bp: drawn, too narrow to label
                b.text((s0 + e0) / 2000, 1.56, name, ha="center", fontsize=8, color=INK2)
    b.text(29, 1.56, "intergenic spacer", ha="center", fontsize=8, color=INK2)
    b.set_ylim(0.4, 1.65)
    b.set_xlim(0, 44.838)
    b.set_xlabel("position in the 45S unit (kb, KY962518.1)")
    b.set_ylabel("window efficiency  exp(a_w)")
    b.set_title("b  The same unit, 250-bp windows, after the GC model")
    b.legend(loc="lower right", fontsize=8, handletextpad=0.3)

    # ---- c: replicate concordance (out-of-sample values written by evaluate_pilot.py)
    c = ax[1, 0]
    held = HERE / "pilot_heldout.tsv"
    if held.exists():
        import csv
        rows = list(csv.DictReader(open(held), delimiter="\t"))
        xs1, ys1 = np.array([float(r["nygc_18S_flat"]) for r in rows]), np.array([float(r["replicate_18S_flat"]) for r in rows])
        xs2, ys2 = np.array([float(r["nygc"]) for r in rows]), np.array([float(r["replicate"]) for r in rows])
        lo, hi = 0.85 * min(xs1.min(), xs2.min(), ys1.min(), ys2.min()), 1.08 * max(xs1.max(), xs2.max(), ys1.max(), ys2.max())
        c.plot([lo, hi], [lo, hi], color=AXIS, lw=0.9)
        for xs, ys, col, lab in ((xs1, ys1, ORANGE, "18S read-depth ratio, no GC model"),
                                 (xs2, ys2, BLUE, "GC model + consensus anchor windows (held out)")):
            lr = np.log(ys / xs)
            c.scatter(xs, ys, s=34, color=col, edgecolor=SURFACE, linewidth=1.2, zorder=3,
                      label=f"{lab}: offset {100 * (np.exp(lr.mean()) - 1):+.0f}%, pair SD {100 * lr.std(ddof=1):.1f}%")
        c.set_xlim(lo, hi)
        c.set_ylim(lo, hi)
        c.legend(loc="upper left", fontsize=8, handletextpad=0.3)
    c.set_xlabel("45S copies, NYGC library")
    c.set_ylabel("45S copies, older library of the same sample")
    c.set_title("c  Two independent libraries per sample")

    # ---- d: known-truth controls
    d = ax[1, 1]
    rows = [("held-out autosomal\n(truth 2)", lambda r, s: r["truth_regions"]["auto"]["cn"] / 2),
            ("chrX\n(truth 1 male, 2 female)", lambda r, s: r["truth_regions"]["chrX"]["cn"] / (1 if SEX.get(s) == "M" else 2)),
            ("distal junction\n(truth 10)", lambda r, s: r["classes"]["DJ"]["cn"] / 10)]
    rng = np.random.default_rng(0)
    for i, (label, fn) in enumerate(rows):
        for off, col, src in ((0.13, BLUE, ny), (-0.13, ORANGE, hg)):
            v = np.array([fn(r, s) for s, r in src.items()])
            d.scatter(v, np.full(len(v), i + off) + rng.uniform(-0.05, 0.05, len(v)), s=22, color=col, edgecolor=SURFACE, linewidth=0.9, zorder=3)
    d.axvline(1, color=AXIS, lw=0.9)
    d.set_yticks(range(len(rows)))
    d.set_yticklabels([r[0] for r in rows])
    d.set_ylim(-0.6, len(rows) - 0.25)
    d.set_xlim(0.78, 1.06)
    d.annotate("HG00732: part of the culture\nhas lost an X (both libraries)", xy=(0.808, 1.13), xytext=(0.785, 0.42),
               fontsize=8, color=INK2, arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.8))
    d.grid(axis="y", visible=False)
    d.set_xlabel("estimate / truth")
    d.set_title("d  Sequence of known copy number, every sample")
    for yy, col, lab in ((2.5, BLUE, "NYGC library"), (2.33, ORANGE, "older library of the same sample")):
        d.scatter([0.79], [yy], s=22, color=col, edgecolor=SURFACE, linewidth=0.9, zorder=3, clip_on=False)
        d.text(0.797, yy, lab, color=INK, fontsize=8.5, va="center")

    out = HERE / "pilot_figure.png"
    fig.savefig(out, dpi=170)
    print("wrote", out)


if __name__ == "__main__":
    main()
