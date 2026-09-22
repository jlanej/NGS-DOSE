#!/usr/bin/env python3
"""The evidence figure: eight panels, each one thing a skeptic would ask for, drawn from the
cohort page's data and the pilot's tables. Everything is read from files that `ngsdose report`
and `evaluate_pilot.py` write; nothing is typed in.

    evidence_figure.py --report NGS-DOSE-1000G/docs/report.json --pilot example/1000G/pilot -o docs/evidence.png
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

BLUE, ORANGE, AQUA, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#898781", "#e1e0d9"


def rows_of(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(row, col):
    v = row.get(col)
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True, help="report.json of the cohort page")
    ap.add_argument("--pilot", required=True, help="the pilot directory (pilot_heldout.tsv)")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    d = json.loads(Path(a.report).read_text())
    S = d["samples"]
    heldout = rows_of(Path(a.pilot) / "pilot_heldout.tsv")
    rng = np.random.default_rng(1)
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.titleweight": "bold", "axes.titlesize": 10})
    fig, axes = plt.subplots(2, 4, figsize=(16, 7.6))
    ax = axes.ravel()
    n = d["meta"]["n"]

    # a. known copy numbers, every sample -------------------------------------------------------
    A = ax[0]
    items = [("autosomal\n(2)", [s["truth.auto"] / 2 for s in S if s.get("truth.auto")], BLUE),
             ("chrX\nmen\n(1)", [s["truth.chrX"] for s in S if s.get("sex_inferred") == "M" and s.get("truth.chrX")], BLUE),
             ("chrX\nwomen\n(2)", [s["truth.chrX"] / 2 for s in S if s.get("sex_inferred") == "F" and s.get("truth.chrX")], ORANGE),
             ("chrY\nmen\n(1)", [s["truth.chrY"] for s in S if s.get("sex_inferred") == "M" and s.get("truth.chrY") is not None], BLUE),
             ("distal\njunction\n(10)", [s[d["known_truth"]["DJ_col"]] / 10 for s in S if s.get(d["known_truth"]["DJ_col"])], BLUE)]
    for i, (lab, v, col) in enumerate(items):
        v = np.array(v)
        A.scatter(i + rng.uniform(-0.28, 0.28, len(v)), v, s=7, color=col, alpha=0.55, linewidths=0)
        A.hlines(np.median(v), i - 0.35, i + 0.35, color=INK, lw=2)
    A.axhline(1, color=GRID, lw=1, zorder=0)
    A.set_xticks(range(len(items)))
    A.set_xticklabels([x[0] for x in items], fontsize=7.5)
    A.set_xlim(-0.6, len(items) - 0.4)
    A.set_ylim(0.6, 1.12)
    A.set_ylabel("estimate / known copy number")
    A.set_title(f"a  Known copy numbers in {n:,} genomes")
    A.text(0.02, 0.04, "below 0.95: cultures that lost part of an X or Y; distal-junction steps (panel b)", transform=A.transAxes, fontsize=7, color=MUTED)

    # b. distal-junction steps ------------------------------------------------------------------
    B = ax[1]
    dj = d["known_truth"]["DJ_steps"]
    steps = np.array([s["DJ.step"] for s in S if s.get("DJ.step") is not None])
    B.hist(steps, bins=np.arange(-2.35, 1.75, 0.1), color=BLUE, edgecolor="white", linewidth=0.5)
    for k in (-2, -1, 0, 1):
        B.axvline(k, color=GRID, lw=1, zorder=0)
    B.set_yscale("log")
    B.set_xlabel("distal-junction copies relative to the cohort's level")
    B.set_ylabel("samples (log)")
    B.set_title("b  A ten-copy paralog steps in whole copies")
    tot = dj["transmitted"] + dj["not_transmitted"]
    B.text(0.02, 0.96, f"−2: {dj['near'].get('-2', 0)}   −1: {dj['near'].get('-1', 0)}   0: {dj['near'].get('0', 0)}   +1: {dj['near'].get('1', 0)}\nparent→child transmitted {dj['transmitted']} of {tot}; de novo {len(dj['de_novo'])}",
           transform=B.transAxes, fontsize=8, va="top", bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2))
    B.set_ylim(0.8, B.get_ylim()[1] * 4)

    # c. child against midparent ---------------------------------------------------------------
    C = ax[2]
    sc = d["trios"]["scatter"]
    if sc:
        mid, ch = np.array([p["mid"] for p in sc]), np.array([p["c"] for p in sc])
        lo, hi = min(mid.min(), ch.min()) * 0.95, max(mid.max(), ch.max()) * 1.05
        C.plot([lo, hi], [lo, hi], color=GRID, lw=1, zorder=0)
        C.scatter(mid, ch, s=18, color=BLUE, edgecolor="white", linewidth=0.8)
        C.set_ylim(lo, hi * 1.08)
        b, a0 = np.polyfit(mid, ch, 1)
        C.plot([lo, hi], [a0 + b * lo, a0 + b * hi], color=INK, lw=1.5, alpha=0.7)
        t = next((t for t in d["trios"]["table"] if t["column"] == d["trios"]["scatter_column"]), None)
        if t:
            ci = f" ({t['R_lo']:.2f}–{t['R_hi']:.2f})" if "R_lo" in t else ""
            C.text(0.03, 0.96, f"{len(sc)} trios: slope {t['slope']:.2f} ± {t['slope_se']:.2f}\nreliability {min(t['R'], 1):.2f}{ci}", transform=C.transAxes, fontsize=8, va="top", bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2))
    C.set_xlabel("midparent 45S copies")
    C.set_ylabel("child 45S copies")
    C.set_title("c  Inherited: child against midparent")

    # d. the same person, two sequencing technologies (pilot) ---------------------------------
    D = ax[3]
    x1, y1 = np.array([num(r, "nygc") for r in heldout]), np.array([num(r, "replicate") for r in heldout])
    x2, y2 = np.array([num(r, "nygc_18S_flat") for r in heldout]), np.array([num(r, "replicate_18S_flat") for r in heldout])
    lo, hi = 0.9 * min(x1.min(), y1.min(), x2.min(), y2.min()), 1.05 * max(x1.max(), y1.max(), x2.max(), y2.max())
    D.plot([lo, hi], [lo, hi], color=GRID, lw=1, zorder=0)
    D.scatter(x2, y2, s=22, color=ORANGE, edgecolor="white", linewidth=0.8, label=f"18S depth ratio: {np.exp(np.mean(np.log(y2 / x2))) - 1:+.0%}, pair SD {np.std(np.log(y2 / x2), ddof=1):.0%}")
    D.scatter(x1, y1, s=22, color=BLUE, edgecolor="white", linewidth=0.8, label=f"calibrated, anchors held out: {np.exp(np.mean(np.log(y1 / x1))) - 1:+.0%}, pair SD {np.std(np.log(y1 / x1), ddof=1):.1%}")
    D.legend(fontsize=7, frameon=False, loc="upper left")
    D.set_xlabel("45S copies, NovaSeq 2×150 (2019)")
    D.set_ylabel("45S copies, HiSeq 2×100 / 2×126 (2012–15)")
    D.set_title("d  Same person, two technologies")

    # e. what the GC model removes --------------------------------------------------------------
    E = ax[4]
    g = np.array([s["gc_rel_65"] for s in S if s.get("gc_rel_65") and s.get("rDNA45S.18S.flat_over_cn")])
    f = np.array([s["rDNA45S.18S.flat_over_cn"] for s in S if s.get("gc_rel_65") and s.get("rDNA45S.18S.flat_over_cn")])
    m = np.array([s["rDNA45S.18S_over_cn"] for s in S if s.get("gc_rel_65") and s.get("rDNA45S.18S.flat_over_cn")])
    gb = d["biology"]["gc_bias"]
    E.scatter(g, f, s=9, color=ORANGE, alpha=0.6, linewidths=0, label=f"18S depth ratio, no GC model: r = {gb['flat_vs_gc']['r']:.2f}")
    E.scatter(g, m, s=9, color=BLUE, alpha=0.6, linewidths=0, label=f"18S under the fragment-GC model: r = {gb['modelled_vs_gc']['r']:.2f}")
    E.legend(fontsize=7, frameon=False, loc="upper left")
    E.set_xlabel("the library's rate at 65% GC, relative to its mean")
    E.set_ylabel("18S estimate / calibrated 45S estimate")
    E.set_title("e  The model removes the library's GC bias")

    # f. another pipeline, the same files -------------------------------------------------------
    F = ax[5]
    h = d.get("hall") or {}
    if h.get("points"):
        xs, ys = np.array([p["flat"] for p in h["points"]]), np.array([p["theirs"] for p in h["points"]])
        lo, hi = 0.9 * min(xs.min(), ys.min()), 1.05 * max(xs.max(), ys.max())
        F.plot([lo, hi], [lo, hi], color=GRID, lw=1, zorder=0)
        F.scatter(xs, ys, s=9, color=BLUE, alpha=0.6, linewidths=0)
        F.text(0.03, 0.96, f"{h['n']:,} shared samples: r = {h['flat']['r']:.3f}\ntheirs / ours {h['flat_ratio']:.3f}; {h['dup_corrected_ratio']:.3f} after their duplicate exclusion", transform=F.transAxes, fontsize=8, va="top", bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2))
    F.set_xlabel("NGS-DOSE, 18S depth ratio / 2")
    F.set_ylabel("Hall et al. 2021, 18S")
    F.set_title("f  Another pipeline, the same CRAMs")

    # g. against assemblies ---------------------------------------------------------------------
    G = ax[6]
    hp = (d["satellites"].get("hprc") or {})
    if hp.get("rows"):
        pts = [r for r in hp["rows"] if r["assembly_Mb"] > 0 and r["ngsdose_Mb"] > 0 and r["assembly_gapped_Mb"] <= 0.02 * (r["assembly_Mb"] + r["assembly_gapped_Mb"])]
        xs, ys = np.array([r["assembly_Mb"] for r in pts]), np.array([r["ngsdose_Mb"] for r in pts])
        G.plot([0.1, 300], [0.1, 300], color=GRID, lw=1, zorder=0)
        G.scatter(xs, ys, s=14, color=BLUE, edgecolor="white", linewidth=0.6)
        G.set_xscale("log")
        G.set_yscale("log")
        good = [c for c, st in hp["stats"].items() if st.get("n", 0) >= 4 and st.get("pearson", 0) >= 0.95]
        G.text(0.97, 0.04, f"{hp['n_samples']} people, {len(hp['stats'])} satellite families\n{len(good)} families track the assembly with r ≥ 0.95", transform=G.transAxes, fontsize=8, va="bottom", ha="right")
    G.set_xlabel("HPRC assembly, Mb (both haplotypes)")
    G.set_ylabel("NGS-DOSE, Mb")
    G.set_title("g  Against long-read assemblies")

    # h. fetch = scan ---------------------------------------------------------------------------
    H = ax[7]
    r45 = np.array([s["fetch_ratio.rDNA45S"] for s in S if s.get("fetch_ratio.rDNA45S")])
    if len(r45):
        H.hist(r45, bins=30, color=BLUE, edgecolor="white", linewidth=0.5)
        H.set_ylim(0, H.get_ylim()[1] * 1.35)
        H.axvline(1, color=GRID, lw=1, zorder=0)
        H.text(0.03, 0.96, f"{len(r45):,} samples counted both ways\nmedian {np.median(r45):.4f}, range {r45.min():.4f}–{r45.max():.4f}\n0.5 GB read instead of 15 GB", transform=H.transAxes, fontsize=8, va="top", bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2))
    H.set_xlabel("45S copies: targeted fetch / whole-file scan")
    H.set_ylabel("samples")
    H.set_title("h  Fetch equals scan")

    for x in ax:
        x.grid(axis="y", color=GRID, lw=0.6)
        x.set_axisbelow(True)
    fig.suptitle(f"NGS-DOSE on the 1000 Genomes 30× cohort: the evidence, as of {d['meta']['as_of']} ({n:,} of {d['meta']['total']:,} genomes)", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(a.out, dpi=160)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
