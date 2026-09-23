#!/usr/bin/env python3
"""A printable assessment of NGS-DOSE on the cohort, written for a sceptic and focused on the trios:
the methods in brief, a heatmap of transmission statistics for every metric (what is claimed to be
inherited, what must be, what cannot be), a child-against-midparent scatter for each metric, the
same test run on the fetch counts alone, and the satellite arrays against long-read assemblies.
Everything is read from what `ngsdose report` writes; nothing is typed in.

    trio_report.py --report docs/report.json --data docs/data -o docs/trio_report.pdf [--png-dir DIR]
"""
import argparse
import csv
import json
import textwrap
from pathlib import Path

import numpy as np

BLUE, ORANGE, AQUA, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#898781", "#e1e0d9"
GROUP_COLOR = {"rDNA": BLUE, "satellites": AQUA, "truth": MUTED, "culture": ORANGE}
GROUP_TITLE = {"rDNA": "rDNA copy number (the claim)", "satellites": "satellite arrays, mass (genomic: positive controls)",
               "truth": "known copy number (nothing to inherit except the DJ's whole-copy steps)", "culture": "culture and library (not in the nuclear genome)"}
GROUP_ORDER = ["rDNA", "satellites", "truth", "culture"]


def rows_of(path):
    p = Path(path)
    if not p.exists():
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def f(x, nd=2):
    return "–" if x is None or not np.isfinite(x) else f"{x:,.{nd}f}"


def wrap(s, width=124):
    return "\n".join(textwrap.fill(par, width) for par in s.split("\n\n"))


class Doc:
    """Letter pages; text pages laid out top-down in a monospace-free sans font."""

    def __init__(self, out, png_dir=None):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        # no creation date: the same report.json gives the same bytes, so a regeneration that changes nothing commits nothing
        self.plt, self.pdf, self.n, self.png_dir = plt, PdfPages(out, metadata={"CreationDate": None}), 0, Path(png_dir) if png_dir else None
        plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                             "axes.spines.top": False, "axes.spines.right": False, "axes.titleweight": "bold", "axes.titlesize": 9.5, "pdf.fonttype": 42})

    def page(self, title=None):
        fig = self.plt.figure(figsize=(8.5, 11))
        self.n += 1
        if title:
            fig.text(0.07, 0.955, title, fontsize=13 if len(title) < 70 else 11.5, fontweight="bold", va="top")
        fig.text(0.93, 0.03, f"NGS-DOSE trio assessment · page {self.n}", fontsize=7.5, color=MUTED, ha="right")
        return fig

    def close(self, fig):
        self.pdf.savefig(fig)
        if self.png_dir:
            self.png_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(self.png_dir / f"page{self.n:02d}.png", dpi=110)
        self.plt.close(fig)

    def text_page(self, title, blocks, subtitle=None):
        """blocks: list of (heading or None, text) laid out top to bottom; overflow starts a new page."""
        fig, y = self.page(title), 0.92
        if subtitle:
            fig.text(0.07, 0.93, subtitle, fontsize=8.5, color=MUTED, va="top"); y = 0.90
        for head, body in blocks:
            lines = wrap(body).count("\n") + 1
            need = 0.0136 * lines + (0.028 if head else 0.010)
            if y - need < 0.06:
                self.close(fig); fig, y = self.page(title + " (continued)"), 0.92
            if head:
                fig.text(0.07, y, head, fontsize=10, fontweight="bold", va="top"); y -= 0.022
            fig.text(0.07, y, wrap(body), fontsize=8.1, va="top", linespacing=1.3); y -= need - (0.022 if head else 0)
        self.close(fig)

    def done(self):
        self.pdf.close()


def load(a):
    d = json.loads(Path(a.report).read_text())
    D = Path(a.data)
    tr = [r for r in rows_of(D / "transmission.tsv") if not r["column"].endswith(".adj")]
    return d, tr, rows_of(D / "trios.tsv"), rows_of(D / "fetch_check.tsv"), rows_of(D / "satellites_hprc.tsv")


def summary_blocks(d, tr, fc):
    m, kt, md, t = d["meta"], d["known_truth"], d["modes"], d["trios"]
    a, sx, DJ = kt["auto"], kt["sex"], kt["DJ"]
    by = {r["column"]: r for r in tr}
    g = lambda col, k: num(by[col][k]) if col in by else float("nan")
    ok = lambda col: col in by and np.isfinite(g(col, "R_lo"))
    sat_R = [min(num(r["R"]), 1) for r in tr if r["group"] == "satellites" and not r["column"].startswith(("HSat1B", "TEL")) and np.isfinite(num(r["R"]))]
    tel_R = next((num(r["R"]) for r in tr if r["column"] == "TEL.mass_Mb" and np.isfinite(num(r["R"]))), None)
    cul_R = [num(r["R"]) for r in tr if r["group"] == "culture" and np.isfinite(num(r["R"]))]
    tru_R = [num(r["R"]) for r in tr if r["group"] == "truth" and np.isfinite(num(r["R"]))]
    fR = {r["metric"]: num(r["R_fetch"]) for r in fc} if fc else {}
    q = ("Can the copy number of a sequence that a reference genome collapses — the ribosomal DNA arrays, several hundred copies per genome — "
         "be measured from ordinary short-read whole-genome sequencing? No assay of rDNA copy number exists for these samples, so the answer has to "
         "come from things that are known: sequence of known copy number in every genome, Mendelian transmission in trios, the same genome counted "
         "two ways, and long-read assemblies of the same people for the satellite arrays that the same machinery measures.")
    meth = ("Reads are assigned to a class (45S rDNA, 5S rDNA, the distal junction, ten satellite families, the telomeric repeat) by 31-mers "
            "that occur in the class and nowhere else in GRCh38 or CHM13; where the aligner put a read is not used. Fragment 5′ ends are counted "
            "on the unit and in 800 single-copy control regions. A per-library Poisson spline of fragment-end density on fragment GC content, fitted "
            "on the controls, gives the expected count of any sequence; copy number is 2 × observed / expected per 250-bp window, and the windows of "
            "the 45S unit are calibrated across the cohort with the scale set by anchor windows on which three Illumina chemistries agreed. "
            "Two counting modes: scan reads the whole CRAM, fetch retrieves only the controls and the intervals where the aligner places class "
            "reads. Cohort: the 1000 Genomes 30× CRAMs (NYGC, NovaSeq 2×150, GRCh38), all lymphoblastoid cell lines.\n\n"
            "Three 45S estimators travel through every table, two of them NGS-DOSE's. '45S, NGS-DOSE calibrated': the cohort model log C(i,w) = c(i) + a(w) + e(i,w) over every "
            "retained 250-bp window w of the unit, fitted by median polish across samples i with the window efficiencies a(w) pinned to a median of "
            "zero over the anchor windows; the estimate is exp(c(i)). '45S, NGS-DOSE single-sample': 2 × observed / expected fragment ends over the "
            "anchor windows alone, under the sample's own GC model, with no information from any other sample. '45S, 18S depth ratio (published)': "
            "2 × fragment ends in the 18S gene / (positions × the control regions' mean rate), no GC model, no calibration — the estimator of published "
            "studies, computed from the same reads for comparison; it is not NGS-DOSE's estimate.")
    test = ("In a trio, a child's copy number is the mean of the parents' plus segregation; measurement error is not inherited. Values are centred "
            "within population; the correlations are Pearson's. The child–midparent correlation cannot reach 1 even for a perfectly measured heritable "
            "trait, because half of a child's variance is segregation, which the midparent does not predict (r is bounded by about √((1 + ρ)/2), 0.71 "
            "at ρ = 0). The least-squares slope b of child on midparent has no such ceiling, so the slope, corrected for the spousal correlation ρ, is "
            "the statistic that answers what share of the measured variance is real: reliability R = b − ρ(1 − b), 1 for a perfectly measured "
            "heritable trait and 0 for pure error, with family-bootstrap 95% intervals and a one-sided p-value for the slope from 1,000 permutations "
            "of children among families. The test is applied to every metric "
            "in four groups. The rDNA estimators are the claim. The satellite arrays are positive controls: their mass is a property of the genome, "
            "measured by the same k-mer machinery, and must be inherited. Known copy numbers (held-out autosomal sequence, the distal junction) have "
            "nothing to inherit. The culture's and the library's properties (mitochondrial and EBV content, depth, duplicate rate, GC bias, insert "
            "size) are not in the nuclear genome. A method that measured library artefacts would light up the last group; one that measured "
            "nothing would light up none. The same test is then run on the fetch counts alone, and the satellite masses are held against HPRC "
            "assemblies of the same people.")
    res = [f"{m['n']:,} of {m['total']:,} genomes, {t['n_complete']:,} complete trios, as of {m['as_of']}."]
    if a.get("n"):
        res.append(f"Known copy numbers: held-out autosomal sequence {f(a['mean'], 3)} ± {f(a.get('sd'), 3)} (expected 2); the distal junction {f(DJ.get('mean'), 2)} ± {f(DJ.get('sd'), 2)} (expected 10); "
                   f"sex from the reads agrees with the pedigree in {sx['n_inferred'] - len(sx['mismatch']):,} of {sx['n_inferred']:,}.")
    if ok("rDNA45S.cn"):
        res.append(f"Transmission, 45S calibrated estimate: reliability {f(min(g('rDNA45S.cn', 'R'), 1))} ({f(g('rDNA45S.cn', 'R_lo'))}–{f(g('rDNA45S.cn', 'R_hi'))}), child–midparent r {f(g('rDNA45S.cn', 'r_mid'))}, spousal r {f(g('rDNA45S.cn', 'spousal_r'))}."
                   + (f" 5S: {f(min(g('rDNA5S.cn', 'R'), 1))} ({f(g('rDNA5S.cn', 'R_lo'))}–{f(g('rDNA5S.cn', 'R_hi'))})." if ok("rDNA5S.cn") else ""))
    if sat_R:
        res.append(f"Positive controls, {len(sat_R)} satellite families: reliability {f(min(sat_R))} to {f(max(sat_R))} (HSat1B, on Yq, excluded: father-to-son only)"
                   + (f"; the telomeric repeat, which changes with age and in culture, {f(tel_R)}." if tel_R is not None else "."))
    if cul_R or tru_R:
        res.append(f"Negative controls: culture and library {f(min(cul_R))} to {f(max(cul_R))}; known copy numbers {f(min(tru_R))} to {f(max(tru_R))}." if cul_R and tru_R else "")
    if "rDNA45S.cn" in fR and np.isfinite(fR["rDNA45S.cn"]):
        res.append(f"Fetch counts alone: 45S reliability {f(min(fR['rDNA45S.cn'], 1))}; fetch/scan agreement of the 45S estimate across genomes r {f(next((num(r['r']) for r in fc if r['metric'] == 'rDNA45S.cn'), float('nan')), 4)}.")
    if md.get("n_both"):
        c = md["columns"].get("rDNA45S.cn_single", {})
        res.append(f"Targeted fetch against whole-file scan, 45S: median {f(c.get('median'), 4)} ({f(c.get('min'), 4)}–{f(c.get('max'), 4)}) over {md['n_both']:,} genomes.")
    bt = t.get("batches") or {}
    if bt.get("n"):
        (kb, kn), (pb, pn) = (max(bt[w].items(), key=lambda kv: kv[1]) for w in ("child", "parent"))
        lab = lambda b_: f"{int(b_):,}" if str(b_).isdigit() else str(b_)
        if kb != pb:
            # genomic metrics whose children spread differently from their parents: the batch's scale, which moves R and not R rescaled
            moved = [r for r in tr if r["group"] in ("rDNA", "satellites") and np.isfinite(num(r.get("sd_ratio_lo")))
                     and (num(r["sd_ratio_lo"]) > 1 or num(r["sd_ratio_hi"]) < 1)]
            res.append(f"Parents and children were sequenced apart: {kn} of {bt['n']} children in release batch {lab(kb)}, {pn} of {2 * bt['n']} parents in "
                       f"{lab(pb)}. No batch is shared within a family, but one that reads a generation on another scale moves R with it"
                       + (": " + "; ".join(f"{r['label']} spreads {num(r['sd_ratio']):.2f}× as much in the children as in their parents and reads R "
                                           f"{num(r['R']):.2f}, or {num(r['R_rescaled']):.2f} with the children on their parents' scale" for r in moved) if moved else "")
                       + ".")
    hp = (d["satellites"].get("hprc") or {})
    if hp:
        good = [c for c, s_ in hp["stats"].items() if s_.get("n", 0) >= 4 and s_.get("pearson", 0) >= 0.95]
        res.append(f"Assemblies: {len(good)} of {len(hp['stats'])} satellite families track HPRC assemblies of {hp['n_samples']} of these genomes with r ≥ 0.95 ({', '.join(good)}).")
    caveat = ("What this does not show: an absolute calibration of rDNA copy number (none exists for these samples; the scale rests on unit windows "
              "on which three chemistries agree); anything about DNA that is not from a lymphoblastoid cell line or a chemistry other than NovaSeq "
              "2×150; and, because people differ in 45S copy number far more than any competent estimator errs, the trios show that the measured "
              "variation is inherited, not which rDNA estimator is best — that ranking rests on the same people sequenced on two technologies "
              "(the pilot), where the calibrated estimate reproduced and the depth ratio did not.")
    return [("The question", q), ("What is measured", meth), ("How it is tested", test), ("Results", " ".join(x for x in res if x)), ("Limits", caveat)]


def heatmap_page(doc, tr, fc):
    by = {r["column"]: r for r in tr}
    fR = {r["metric"]: (num(r["R_fetch"]), num(r["R_fetch_lo"]), num(r["R_fetch_hi"])) for r in fc} if fc else {}
    order = []
    try:
        from ngsdose.report import TRIO_GROUPS
        for key, _, cols in TRIO_GROUPS:
            order += [(key, col, label) for col, label in cols if col in by]
    except Exception:
        for key in GROUP_ORDER:
            order += [(key, r["column"], r["label"]) for r in tr if r["group"] == key]
    if not order:
        return
    stats = [("r_father", "child–father r"), ("r_mother", "child–mother r"), ("r_mid", "child–midparent r"), ("spousal_r", "spousal r"), ("slope", "midparent slope"), ("R", "reliability R")]
    cols = [l for _, l in stats] + (["R from fetch"] if fR else [])
    M = np.array([[num(by[col][k]) for k, _ in stats] + ([fR.get(col, (np.nan,))[0]] if fR else []) for _, col, _ in order])
    fig = doc.page("Transmission of every metric, at a glance")
    y0, h0 = 0.19, 0.69
    ax = fig.add_axes([0.42, y0, 0.54, h0])
    ax.imshow(np.clip(np.nan_to_num(M, nan=-1), 0, 1), cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7.5, color="white" if v > 0.55 else INK)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=30, ha="left", fontsize=8); ax.xaxis.tick_top()
    ax.set_yticks(range(len(order))); ax.set_yticklabels([lab for _, _, lab in order], fontsize=8)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    last, start = None, 0
    for i, (key, _, _) in enumerate(order + [(None, None, None)]):
        if key != last:
            if last is not None:
                ax.axhline(i - 0.5, color="white", lw=3)
                fig.text(0.05, y0 + h0 * (1 - (start + i) / 2 / len(order)), textwrap.fill(GROUP_TITLE[last], 20), fontsize=7.5, fontweight="bold", color=GROUP_COLOR[last], va="center")
            last, start = key, i
    fig.text(0.07, 0.15, wrap("Colour and number: the statistic, from 0 (white) to 1 (blue); values are within population, over the complete trios. "
                               "A metric that is inherited reads near 1 in the correlation columns and in R; a property of the culture or the library reads "
                               "near 0; sequence of known copy number has no variation to inherit and reads near 0 too. HSat1B lives mostly on Yq and passes "
                               "from father to son only, which dilutes its midparent statistics by design. The last column repeats R with the cohort layer and "
                               "the trio test run on the fetch counts alone.", 120), fontsize=7.8, va="top", color=INK)
    doc.close(fig)


def scatter_pages(doc, tr, trios):
    by = {r["column"]: r for r in tr}
    items = [(r["group"], r["column"], r["label"]) for r in tr if r["column"] in by]
    per = 16
    for page in range(0, len(items), per):
        fig = doc.page("Child against midparent, every metric" + (" (continued)" if page else ""))
        chunk = items[page:page + per]
        for k, (grp, col, label) in enumerate(chunk):
            ax = fig.add_axes([0.08 + (k % 4) * 0.235, 0.75 - (k // 4) * 0.205, 0.185, 0.15])
            c = np.array([num(t.get(f"{col}.child")) for t in trios]); fa = np.array([num(t.get(f"{col}.father")) for t in trios]); mo = np.array([num(t.get(f"{col}.mother")) for t in trios])
            ok = np.isfinite(c) & np.isfinite(fa) & np.isfinite(mo)
            mid, ch = (fa[ok] + mo[ok]) / 2, c[ok]
            if len(ch) >= 3:
                lo, hi = min(mid.min(), ch.min()), max(mid.max(), ch.max()); pad = (hi - lo or 1) * 0.06
                ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=GRID, lw=1, zorder=0)
                ax.scatter(mid, ch, s=9, color=GROUP_COLOR.get(grp, BLUE), alpha=0.7, linewidths=0)
                b, a0 = np.polyfit(mid, ch, 1)
                ax.plot([lo, hi], [a0 + b * lo, a0 + b * hi], color=INK, lw=1.2, alpha=0.7)
                ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
                t = by[col]
                pp = num(t.get("perm_p"))
                ax.text(0.03, 0.97, f"n {int(num(t['n_trios']))}  r {num(t['r_mid']):.2f}" + (f"  p {'<0.001' if pp < 0.001 else f'{pp:.3f}'}" if np.isfinite(pp) else "") + f"\nR {min(num(t['R']), 1):.2f}" + (f" ({num(t['R_lo']):.2f}–{num(t['R_hi']):.2f})" if np.isfinite(num(t.get("R_lo"))) else ""),
                        transform=ax.transAxes, fontsize=6.5, va="top", bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5))
            ax.set_title(textwrap.shorten(label, 34, placeholder="…"), fontsize=7.5, color=GROUP_COLOR.get(grp, INK))
            ax.tick_params(labelsize=6)
            if k // 4 == min(3, (len(chunk) - 1) // 4):
                ax.set_xlabel("midparent", fontsize=6.5)
            if k % 4 == 0:
                ax.set_ylabel("child", fontsize=6.5)
        fig.text(0.07, 0.085, wrap("Raw values (not centred); grey diagonal: child equals midparent; black line: least squares. Printed: the report's statistics, computed within "
                                   "population (n trios, Pearson child–midparent r, the permutation p of the midparent slope, reliability R with its family-bootstrap 95% interval). Colours: blue rDNA estimators, green satellite "
                                   "arrays (positive controls), grey known copy numbers, orange culture and library (negative controls).", 140),
                 fontsize=7.4, va="top")
        doc.close(fig)


def fetch_page(doc, d, fc):
    if not fc:
        return
    fig = doc.page("The same test on the fetch counts alone")
    rows = [r for r in fc if r.get("identical") not in ("True", "true", "1")]
    ax = fig.add_axes([0.10, 0.55, 0.38, 0.33])
    xs, ys, labs = [], [], []
    for r in rows:
        a_, b_ = num(r.get("R_scan")), num(r.get("R_fetch"))
        if np.isfinite(a_) and np.isfinite(b_):
            xs.append(min(a_, 1.05)); ys.append(min(b_, 1.05)); labs.append(r["label"])
    if xs:
        ax.plot([-0.3, 1.1], [-0.3, 1.1], color=GRID, lw=1, zorder=0)
        ax.scatter(xs, ys, s=22, color=BLUE, edgecolor="white", linewidth=0.8)
        ax.set_xlabel("reliability from the scans"); ax.set_ylabel("reliability from the fetches"); ax.set_title("Transmission: one point per metric (see the table)")
    ax2 = fig.add_axes([0.58, 0.55, 0.38, 0.33])
    med = [num(r["ratio_median"]) for r in rows]; lab2 = [textwrap.shorten(r["label"], 26, placeholder="…") for r in rows]
    y = np.arange(len(rows))
    ax2.barh(y, [m_ - 1 for m_ in med], left=1, color=BLUE, height=0.6)
    for i, r in enumerate(rows):
        ax2.plot([num(r["q10"]), num(r["q90"])], [i, i], color=INK, lw=1)
    ax2.axvline(1, color=GRID, lw=1); ax2.set_yticks(y); ax2.set_yticklabels(lab2, fontsize=6.5); ax2.invert_yaxis()
    ax2.set_xlabel("fetch / scan estimate per genome: median, 10–90%"); ax2.set_title("Agreement across genomes")
    cell = [[r["label"], r["n"], f"{num(r['ratio_median']):.4f}", f"{num(r.get('r')):.4f}", f(min(num(r.get('R_scan')), 1)) + (f" ({f(num(r.get('R_scan_lo')))}–{f(num(r.get('R_scan_hi')))})" if np.isfinite(num(r.get("R_scan_lo"))) else ""),
             f(min(num(r.get('R_fetch')), 1)) + (f" ({f(num(r.get('R_fetch_lo')))}–{f(num(r.get('R_fetch_hi')))})" if np.isfinite(num(r.get("R_fetch_lo"))) else "")] for r in rows]
    axt = fig.add_axes([0.07, 0.17, 0.89, 0.30]); axt.axis("off")
    tb = axt.table(cellText=cell, colLabels=["metric", "genomes", "fetch / scan", "r across genomes", "R from scan (95% CI)", "R from fetch (95% CI)"], loc="upper center", cellLoc="left", colWidths=[0.30, 0.09, 0.12, 0.14, 0.18, 0.18])
    tb.auto_set_font_size(False); tb.set_fontsize(7); tb.scale(1, 1.25)
    for (i, j), c in tb.get_celld().items():
        c.set_edgecolor(GRID)
        if i == 0:
            c.set_text_props(fontweight="bold")
    ident = [r["label"] for r in fc if r.get("identical") in ("True", "true", "1")]
    fig.text(0.07, 0.11, wrap(f"The cohort layer (window calibration, control-region PCs) and the trio test were run on the fetch counts of {d['fetch_check']['n']:,} genomes without "
                              "reference to the scans. Columns made from the same reads in both modes come out identical and are omitted"
                              + (f" ({', '.join(ident)})" if ident else "") + ". The classes differ by what the sinks miss and by a calibration learned twice.", 120), fontsize=7.8, va="top")
    doc.close(fig)


def assembly_page(doc, d, hp_rows):
    hp = d["satellites"].get("hprc") or {}
    if not hp or not hp_rows:
        return
    fig = doc.page("Satellite arrays against long-read assemblies of the same people")
    st = hp["stats"]
    classes = [c for c, s_ in st.items() if s_.get("n", 0) >= 3]
    for k, cls in enumerate(classes[:12]):
        ax = fig.add_axes([0.08 + (k % 4) * 0.235, 0.66 - (k // 4) * 0.22, 0.185, 0.16])
        pts = [(num(r["assembly_Mb"]), num(r["ngsdose_Mb"])) for r in hp_rows if r["cls"] == cls and num(r["assembly_gapped_Mb"]) <= 0.02 * (num(r["assembly_Mb"]) + num(r["assembly_gapped_Mb"])) and num(r["assembly_Mb"]) > 0]
        if pts:
            x, y = np.array(pts).T
            lo, hi = 0.9 * min(x.min(), y.min()), 1.1 * max(x.max(), y.max())
            ax.plot([lo, hi], [lo, hi], color=GRID, lw=1, zorder=0)
            ax.scatter(x, y, s=16, color=AQUA, edgecolor="white", linewidth=0.8)
            ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        s_ = st[cls]
        ax.set_title(f"{cls}: r {s_.get('pearson', float('nan')):.2f}, ratio {s_.get('ratio_median', float('nan')):.2f}, n {s_.get('n', 0)}", fontsize=7.5)
        ax.tick_params(labelsize=6)
        if k % 4 == 0:
            ax.set_ylabel("NGS-DOSE, Mb", fontsize=6.5)
        ax.set_xlabel("assembly, Mb", fontsize=6.5)
    good = [c for c, s_ in st.items() if s_.get("n", 0) >= 4 and s_.get("pearson", 0) >= 0.95]
    fig.text(0.07, 0.15, wrap(f"{hp['n_samples']} of the genomes have HPRC release-2 assemblies; the CenSat annotation of both haplotypes gives the size of every satellite array, "
                              "the same kind of sequence the k-mer machinery measures. A genome is compared in a class only when the arrays its assembly did not close "
                              f"are immaterial. Grey diagonal: equality. The relative measures (β-satellite, CER, ACRO) sit below it by a constant factor, their k-mer recall. "
                              f"Tracking the assemblies across people with r ≥ 0.95: {', '.join(good) if good else 'none yet'}. Assemblies collapse the rDNA itself and are no truth for it.", 120),
             fontsize=7.8, va="top")
    doc.close(fig)


def truth_page(doc, d):
    kt, md = d["known_truth"], d["modes"]
    a, X, Y, DJ, sx = kt["auto"], kt["chrX"], kt["chrY"], kt["DJ"], kt["sex"]
    rows = [["held-out autosomal sequence (80 regions)", "2", f"{f(a.get('mean'), 3)} ± {f(a.get('sd'), 3)}", f"{a.get('n', 0):,}"],
            ["chrX in men (60 regions)", "1", f"{f(X['M'].get('mean'), 3)} ± {f(X['M'].get('sd'), 3)}", f"{X['M'].get('n', 0):,}"],
            ["chrX in women with an intact culture", "2", f"{f(sx.get('women_intact', {}).get('mean'), 3)} ± {f(sx.get('women_intact', {}).get('sd'), 3)}", f"{sx.get('women_intact', {}).get('n', 0):,}"],
            ["chrY in men with an intact Y (40 regions)", "1", f"{f(sx.get('men_intact_Y', {}).get('mean'), 3)} ± {f(sx.get('men_intact_Y', {}).get('sd'), 3)}", f"{sx.get('men_intact_Y', {}).get('n', 0):,}"],
            ["chrY in women", "0", f"at most {f(Y['F'].get('max'), 4)}", f"{Y['F'].get('n', 0):,}"],
            ["distal junction (one per acrocentric short arm)", "10", f"{f(DJ.get('mean'), 2)} ± {f(DJ.get('sd'), 2)}", f"{DJ.get('n', 0):,}"],
            ["sex inferred from the reads", "pedigree", f"agrees in {sx['n_inferred'] - len(sx['mismatch']):,} of {sx['n_inferred']:,}", ""]]
    for col, dd in md.get("columns", {}).items():
        if dd.get("n"):
            rows.append([f"fetch / scan, {dd['label']}", "1", f"{f(dd['median'], 4)} ({f(dd['min'], 4)}–{f(dd['max'], 4)})", f"{dd['n']:,}"])
    fig = doc.page("Sequence of known copy number, and fetch against scan")
    axt = fig.add_axes([0.07, 0.50, 0.89, 0.40]); axt.axis("off")
    tb = axt.table(cellText=rows, colLabels=["what", "expected", "measured", "n"], loc="upper center", cellLoc="left", colWidths=[0.42, 0.10, 0.34, 0.10])
    tb.auto_set_font_size(False); tb.set_fontsize(7.5); tb.scale(1, 1.3)
    for (i, j), c in tb.get_celld().items():
        c.set_edgecolor(GRID)
        if i == 0:
            c.set_text_props(fontweight="bold")
    dj = kt.get("DJ_steps") or {}
    txt = ("Every genome carries sequence whose copy number is not in question, measured by exactly the code that measures the rDNA: held-out autosomal "
           "regions, the X and Y, and the distal junction, a 400-kb sequence present once on each of the ten acrocentric short arms — multi-copy, paralogous "
           "and acrocentric, the kind of sequence the rDNA is. Women read the X and every genome reads the distal junction a few percent below expectation; "
           "both are late-replicating sequence, which DNA from a growing culture under-represents.")
    if dj.get("carriers") is not None:
        nr = lambda k: dj["near"].get(k, dj["near"].get(str(k), 0))
        tot = dj["transmitted"] + dj["not_transmitted"]
        txt += (f"\n\nRelative to the cohort's level the distal junction sits at whole numbers: {nr(-2)} genomes at −2, {nr(-1)} at −1, {nr(0):,} at 0, {nr(1)} at +1 "
                f"(robust SD {f(dj['spread'])} copies). A step is a structural variant of an acrocentric short arm; where a carrier parent and a child were both counted it was "
                f"transmitted in {dj['transmitted']} of {tot}, and {len(dj['de_novo'])} child(ren) carry a step neither counted parent has.")
    fig.text(0.07, 0.44, wrap(txt, 120), fontsize=8.2, va="top", linespacing=1.35)
    doc.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--png-dir", help="also write each page as a PNG here")
    a = ap.parse_args()
    d, tr, trios, fc, hp_rows = load(a)
    doc = Doc(a.out, a.png_dir)
    doc.text_page("NGS-DOSE on the 1000 Genomes 30× cohort: a trio-focused assessment", summary_blocks(d, tr, fc), subtitle=f"as of {d['meta']['as_of']}, from the tables ngsdose report writes; nothing is typed in")
    if tr:
        heatmap_page(doc, tr, fc)
        scatter_pages(doc, tr, trios)
    fetch_page(doc, d, fc)
    assembly_page(doc, d, hp_rows)
    truth_page(doc, d)
    doc.done()
    print(f"wrote {a.out} ({doc.n} pages)")


if __name__ == "__main__":
    main()
