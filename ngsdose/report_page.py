"""The report page: prose and layout around the numbers `report.build` computed.

Every number in the text comes from `data`; nothing is typed in. Sections that have nothing to
show yet say so and say what would make them appear, so that the page can be published at any
stage of a cohort run.
"""
from __future__ import annotations

import html
import json
import math

from .report import ASSETS, fmt

esc = lambda x: html.escape(str(x))
SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
SUPERPOP_NAMES = {"AFR": "African", "AMR": "American", "EAS": "East Asian", "EUR": "European", "SAS": "South Asian"}


def pm(d: dict, nd=3, key="mean") -> str:
    """'mean ± SD (n)' from a describe() dict."""
    if not d or not d.get("n"):
        return "–"
    s = fmt(d[key], nd)
    if d.get("sd") is not None:
        s += f" ± {fmt(d['sd'], nd)}"
    return s


class Page:
    def __init__(self, data):
        self.d = data
        self.parts: list[str] = []
        self.charts: dict = {}
        self.toc: list[tuple[str, str]] = []

    def h(self, s: str):
        self.parts.append(s)

    def section(self, id_: str, title: str, short: str | None = None):
        self.toc.append((id_, short or title))
        self.h(f'<section id="{id_}"><h2>{esc(title)}</h2>')

    def end(self):
        self.h("</section>")

    def chart(self, id_: str, spec: dict, title: str, caption: str = ""):
        self.charts[id_] = spec
        self.h(f'<figure><div class="title">{esc(title)}</div><div class="chart" id="chart-{id_}"></div>'
               + (f"<figcaption>{caption}</figcaption>" if caption else "") + "</figure>")

    def tiles(self, items: list[tuple[str, str, str]]):
        self.h('<div class="tiles">' + "".join(f'<div class="tile"><div class="label">{esc(a)}</div><div class="value">{b}</div><div class="note">{esc(c)}</div></div>'
                                             for a, b, c in items) + "</div>")

    def table(self, rows: list[list], header: list[str], numeric: set[int] = frozenset(), flagged=None, filter_box=False, wrap=False):
        head = "".join(f'<th class="{"num" if i in numeric else ""}">{esc(h)}</th>' for i, h in enumerate(header))
        body = []
        for r in rows:
            cls = ' class="flagged"' if flagged and flagged(r) else ""
            body.append(f"<tr{cls}>" + "".join(f'<td class="{"num" if i in numeric else ""}">{c if isinstance(c, str) and c.startswith("<") else esc(c)}</td>' for i, c in enumerate(r)) + "</tr>")
        t = f'<table class="data"><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table>'
        if filter_box:
            self.h('<input class="filter" type="search" placeholder="filter rows…" aria-label="filter rows">')
        self.h(f'<div class="tablewrap">{t}</div>' if wrap else t)


def page(data: dict, rows: list[dict]) -> str:
    m, kt, md, tr, bio, pcs, sat, qc = (data["meta"], data["known_truth"], data["modes"], data["trios"], data["biology"], data["pcs"],
                                         data["satellites"], data["qc"])
    n, total = m["n"], m["total"]
    P = Page(data)
    sex_levels = [["M", "male"], ["F", "female"]]
    superpop_order = [s for s in SUPERPOPS if any(r.get("superpop") == s for r in rows)]

    # ---------------------------------------------------------------- masthead
    P.h(f'''<header class="mast"><h1>{esc(m["title"])}</h1>
<p class="sub">Ribosomal DNA copy number, and other sequence that reference genomes collapse, measured from
short-read genomes and checked against what is known — on the 1000 Genomes cohort, as the run proceeds.</p>
<div class="hero"><div class="n">{n:,}</div><div class="of">of {total:,} genomes {"scanned" if m["primary_mode"] == "scan" else "counted"}
&middot; {m["n_fetch"]:,} also fetched &middot; {tr["n_complete"]:,} of {tr["n_total"]:,} trios complete</div></div>
<div class="progress"><div style="width:{100 * n / max(total, 1):.1f}%"></div></div>
<div class="stamp">As of {esc(m["as_of"])}. Every number and figure on this page is recomputed from the counts files in this repository by
<code>ngsdose report</code> ({esc(m["generator"])}); nothing is typed in. Partial results are published as they stand.</div></header>''')

    # ---------------------------------------------------------------- the case, in numbers
    a, X, Y, DJ, sx = kt["auto"], kt["chrX"], kt["chrY"], kt["DJ"], kt["sex"]
    dj = kt.get("DJ_steps") or {}
    t45 = next((t for t in tr["table"] if t["column"] == "rDNA45S.cn"), None)
    gcb = bio.get("gc_bias") or {}
    hall = data.get("hall") or {}
    hp = (sat.get("hprc") or {}).get("stats") or {}
    tracking = [c for c, st_ in hp.items() if st_.get("n", 0) >= 4 and st_.get("pearson", 0) >= 0.95]
    case = [("Known copy numbers, every sample", f'{fmt(a.get("mean"), 3)} ± {fmt(a.get("sd"), 3)}', f"held-out autosomal sequence, truth 2, n = {a.get('n', 0):,}", "#truth")]
    if sx.get("women_intact", {}).get("n"):
        case.append(("Sex from the reads", f'{sx["n_inferred"] - len(sx["mismatch"]):,} of {sx["n_inferred"]:,}',
                     f"chrX {fmt(sx['chrX_men_max'], 2)} at most in men, {fmt(sx['women_intact']['min'], 2)} at least in women", "#truth"))
    if dj.get("carriers") is not None:
        tot = dj["transmitted"] + dj["not_transmitted"]
        case.append(("A ten-copy paralog", f'{fmt(DJ.get("median"), 2)} ± {fmt(dj.get("spread"), 2)}', f"distal junction; {sum(dj['near'].get(k, 0) for k in ('-2', '-1', '1', '2', -2, -1, 1, 2))} people one or two copies off, steps transmitted {dj['transmitted']} of {tot}", "#truth"))
    if md["n_both"]:
        c45 = md["columns"].get("rDNA45S.cn_single", {})
        case.append(("Fetch = scan", fmt(c45.get("median"), 4), f"45S, {md['n_both']:,} samples both ways, range {fmt(c45.get('min'), 4)}–{fmt(c45.get('max'), 4)}", "#modes"))
    if t45 and "R_lo" in t45:
        case.append(("Inherited", fmt(min(t45["R"], 1.0), 2), f"45S transmission reliability ({fmt(t45['R_lo'], 2)}–{fmt(t45['R_hi'], 2)}), {t45['n_trios']} trios", "#trios"))
    if gcb.get("flat_vs_gc", {}).get("n"):
        case.append(("What the GC model removes", f'r {fmt(gcb["flat_vs_gc"].get("r"), 2)} → {fmt(gcb["modelled_vs_gc"].get("r"), 2)}', "how much the estimate follows the library's GC bias, before and after", "#gcmodel"))
    if hall.get("n", 0) >= 3:
        case.append(("Another pipeline, same files", f'r = {fmt(hall["flat"].get("r"), 3)}', f"Hall et al. 2021, {hall['n']:,} shared samples; offset explained by the duplicate flag", "#published"))
    if tracking:
        case.append(("Against assemblies", f"{len(tracking)} of {len(hp)} families", f"track HPRC assemblies with r ≥ 0.95 in {(sat.get('hprc') or {}).get('n_samples', 0)} people", "#satellites"))
    P.h('<div class="tiles case">' + "".join(f'<a class="tile" href="{h}"><div class="label">{esc(l)}</div><div class="value">{v}</div><div class="note">{esc(n)}</div></a>' for l, v, n, h in case) + "</div>")

    # ---------------------------------------------------------------- 1. what and why
    P.section("what", "What is being measured, and why", "What and why")
    P.h('''<p class="lede">Every human genome carries a few hundred copies of the ribosomal DNA unit — the 43-kb sequence that
encodes the ribosome's RNA — in tandem arrays on the short arms of the five acrocentric chromosomes, next to the
5S arrays and a shared block of sequence called the distal junction. How many copies a person has varies
severalfold, is inherited, and has been linked to phenotypes; it is also invisible to ordinary genome analysis,
because a read from a sequence that occurs four hundred times cannot be placed, and reference genomes collapse the
arrays into a token copy.</p>
<p>Hundreds of thousands of short-read genomes already exist, in biobanks and consortia, and almost all of them have
been aligned to GRCh38 and stored as CRAM. If the dosage of multi-copy sequence can be measured accurately from those
files — without re-sequencing and, ideally, without reading whole files — it becomes a trait that can be studied at
that scale. That is what NGS-DOSE is for.</p>
<p><strong>How it measures.</strong> Reads are assigned to a class (45S rDNA, 5S rDNA, the distal junction, and
several satellite families) by 31-mers that occur in that class and nowhere else in two reference genomes, so that
where the aligner put a read does not matter. Fragment ends are counted; a per-sample model of how coverage depends on
fragment GC content is fitted on 800 single-copy control regions; the copy number of each 250-bp window of the unit
is the ratio of what was seen to what that model expects, and the windows are calibrated across the cohort so that
sequence-specific dropout does not masquerade as dosage. Two ways to count: a <em>scan</em> reads the whole file
(about 14 CPU-minutes for a 30× genome), a <em>fetch</em> retrieves only the control regions and the few places
the aligner puts class reads (about half a gigabyte and a minute), and the two must agree.</p>
<p><strong>What "working" means here.</strong> There is no assay of rDNA copy number for these samples to compare
with. What there is, in every sample, is sequence whose copy number <em>is</em> known, measured by exactly the code
that measures the classes: held-out autosomal regions (two copies), the X and Y chromosomes (one or two, one or
none, by sex), and the distal junction, present once on each of the ten acrocentric short arms. There are 602 trios,
in which a child's dosage must be the average of the parents' plus segregation, so the share of an estimate's
variance that is transmitted can be measured — and compared between estimators. There are two independent counting
modes of the same files, and a published table of values for the same CRAMs. This page shows each of those as the
cohort accumulates. The method and its own audit are described in
<a href="https://github.com/jlanej/NGS-DOSE/blob/main/docs/DESIGN.md">DESIGN.md</a>.</p>''')
    P.end()

    # ---------------------------------------------------------------- 2. the run
    P.section("run", "The run so far", "The run")
    eng = ", ".join(f"{esc(k)} ({v:,})" for k, v in sorted(m["engines"].items(), key=lambda kv: -kv[1]))
    P.tiles([("Genomes scanned", f"{m['n_scan']:,}", f"of {total:,}"), ("Fetched as well", f"{m['n_fetch']:,}", "targeted mode, same files"),
             ("Complete trios", f"{tr['n_complete']:,}", f"of {tr['n_total']:,}"),
             ("Median depth", fmt(qc["depth"].get("median"), 1) + "×", f"{fmt(qc['depth'].get('q10'), 1)}–{fmt(qc['depth'].get('q90'), 1)} (10–90%)"),
             ("Insert size", fmt(qc["insert"].get("median"), 0) + " bp", "median of medians"),
             ("Duplicate-flagged", fmt(qc["dup"].get("median"), 3, pct=True), "of control reads, median"),
             ("Scan time", fmt(qc["elapsed"].get("median") / 60 if qc["elapsed"].get("n") else None, 1) + " min", "per genome, median" if qc["elapsed"].get("n") else "no scans yet")])
    if m.get("by_superpop"):
        sp = m["by_superpop"]
        cats = [s for s in SUPERPOPS if s in sp]
        P.chart("progress", dict(type="meters", categories=[f"{c} · {SUPERPOP_NAMES[c]}" for c in cats], values=[sp[c]["done"] for c in cats], totals=[sp[c]["total"] for c in cats]),
                "Samples counted, by super-population", "The cohort's 26 populations in five groups; the pale track is the whole cohort.")
    prov = [f"Engine builds: {eng}.", f"Resource bundle: {esc(m.get('bundle'))}; panel(s) {', '.join(esc(x) for x in m['panel_sha']) or '–'}; controls {', '.join(esc(x) for x in m['controls_sha']) or '–'}"
            + (f"; sinks {', '.join(esc(x) for x in m['sinks_sha'])}" if m["sinks_sha"] else "") + " (SHA-256 prefixes; one of each means one cohort)."]
    if len(m["panel_sha"]) > 1 or len(m["controls_sha"]) > 1:
        prov.append('<strong>More than one resource set is in play: these samples should not be analysed as one cohort until that is resolved.</strong>')
    prov.append(f"Read length {', '.join(str(x) for x in qc['read_length'])}; placement grid {', '.join(qc['placement_bins'])} bp. "
                + (f"<strong>{len(m['eof_absent'])} file(s) lacked an end-of-file marker</strong>: {', '.join(esc(s) for s in m['eof_absent'])}." if m["eof_absent"] else "Every input carried its end-of-file marker."))
    P.h("<p class=\"small\">" + " ".join(prov) + "</p>")
    if data["flags"]:
        P.h(f'<p class="small">{len(data["flags"])} sample(s) carry a flag (aneuploid chromosome, sex mismatch, mosaic X loss, unusual DJ, low depth, another engine build); they are marked in the <a href="#samples">sample table</a> and listed in <code>data/flags.tsv</code>. A flag is a thing to look at, not a verdict.</p>')
    P.end()

    # ---------------------------------------------------------------- 3. known truth
    P.section("truth", "Evidence 1 — sequence of known copy number, in every sample", "Known truth")
    P.h(f'''<p class="lede">If the method is right, held-out autosomal sequence reads 2, the X reads 1 in men and 2 in women, the
Y reads 1 and 0, and the distal junction reads 10 — in every sample, by the same code that measures the rDNA.</p>
<p>Across {a.get("n", 0):,} samples the held-out autosomal regions read <strong>{pm(a)}</strong> copies (expected 2).
chrX: <strong>{pm(X["M"])}</strong> in {X["M"].get("n", 0):,} men and a median of <strong>{fmt(X["F"].get("median"), 3)}</strong> in {X["F"].get("n", 0):,} women
(mean {pm(X["F"])}: the spread is a handful of cultures that have lost an X in part of their cells, listed below);
chrY: <strong>{pm(Y["M"])}</strong> in men (a few cultures have lost the Y in part of their cells too) and <strong>{pm(Y["F"], 4)}</strong> in women. The distal junction — ten copies, on five different chromosomes, measured
by the k-mer path that the rDNA uses — reads <strong>{pm(DJ)}</strong>{" (cohort-calibrated)" if kt["DJ_col"] == "DJ.cn" else ""}.
{"Women read the X a little below 2 and every sample reads the distal junction a little below 10: both are late-replicating sequence, and DNA from a growing cell culture under-represents it (see the cell-line section)." if X["F"].get("median", 2) < 1.98 else ""}</p>''')
    P.h('<div class="grid2">')
    P.chart("auto", dict(type="hist", col="truth.auto", xlabel="copies", ref=[dict(x=2, label="truth: 2")], xfmt=3), "Held-out autosomal regions (80 regions, 0.8 Mb)")
    P.chart("chrX", dict(type="hist", col="truth.chrX", group=dict(col="sex_inferred", levels=sex_levels), xlabel="copies", ref=[dict(x=1, label="1"), dict(x=2, label="2")], xfmt=2), "chrX (60 regions), by sex")
    P.chart("chrY", dict(type="hist", col="truth.chrY", group=dict(col="sex_inferred", levels=sex_levels), xlabel="copies", ref=[dict(x=0, label="0"), dict(x=1, label="1")], xfmt=2), "chrY (40 X-degenerate regions), by sex")
    P.chart("dj", dict(type="hist", col=kt["DJ_col"], xlabel="copies", ref=[dict(x=10, label="truth: 10")], xfmt=2), "Distal junction (one per acrocentric short arm)")
    P.h("</div>")
    sx = kt["sex"]
    if sx["n_pedigree"]:
        P.h(f'<p>Sex read from the X and Y agrees with the pedigree in {sx["n_inferred"] - len(sx["mismatch"]):,} of {sx["n_inferred"]:,} samples'
            + (f'; it does not in <strong>{", ".join(esc(s) for s in sx["mismatch"][:20])}{" and " + str(len(sx["mismatch"]) - 20) + " more" if len(sx["mismatch"]) > 20 else ""}</strong> — worth a look at those files (a swapped sample, or a line that has lost its Y).' if sx["mismatch"] else ".") + "</p>")
    outl = [(s, f) for s, f in data["flags"] if "chrX" in f or "chrY" in f or "autosomal" in f]
    if outl:
        P.h(f'<details><summary>{len(outl)} sample(s) off the expected value</summary>')
        P.table([[s, f] for s, f in outl], ["sample", "what"])
        P.h("<p class=\"small\">A woman whose X reads well below 2, or a man whose Y does, has lost that chromosome in part of the cell culture — the known behaviour of lymphoblastoid lines, and a reason the controls are measured in every sample.</p></details>")
    dj = kt.get("DJ_steps") or {}
    if dj.get("carriers") is not None:
        near = dj["near"]
        P.h(f'''<h3>Steps of one copy in the distal junction</h3>
<p>Everyone has ten distal junctions, one per acrocentric short arm; a person with a rearranged short arm has nine, and a
Robertsonian translocation — the commonest structural rearrangement in humans, about one person in a thousand — fuses two
acrocentrics and loses both their short arms: eight. Relative to the cohort's level ({fmt(dj["median"], 2)}), the copy number of
each sample should therefore sit near a whole number, and a step, being a structural variant, should pass to half of a
carrier's children and arise de novo in almost none. Within ±0.3 of a step: <strong>{near.get(-2, 0)}</strong> samples at −2,
<strong>{near.get(-1, 0)}</strong> at −1, {near.get(0, 0):,} at 0, <strong>{near.get(1, 0)}</strong> at +1; the main mode has a robust SD of
{fmt(dj["spread"], 2)} copies, and {dj["between"]} samples sit between steps.</p>''')
        if dj["carriers"]:
            rows_c = []
            for c in dj["carriers"]:
                rel = "; ".join(f"{r['who']} {esc(r['sample'])} {r['step']:+.2f}" for r in c["relatives"]) or "none counted"
                rows_c.append([c["sample"], c.get("pop") or "", c.get("sex") or "", f"{c['step']:+.2f}", rel])
            P.table(rows_c, ["sample", "population", "sex", "step (copies)", "relatives counted, and their step"], numeric={3})
            two = [c for c in dj["carriers"] if c["step"] <= -1.5 and c.get("arm_content")]
            if two and dj.get("arm_ref"):
                arm = [cls for cls in ("ACRO", "SST1", "bSat", "HSat3", "CER", "HSat1A", "aSatHOR") if cls in dj["arm_ref"]]
                P.h("<p>A lost short arm takes its satellite arrays with it. The satellite families of the acrocentric short arms, in the two-copy carriers, as a fraction of the cohort's median (the pan-centromeric α-satellite, which every chromosome carries, is the control):</p>")
                P.table([[c["sample"], f"{c['step']:+.2f}"] + [fmt(c["arm_content"].get(cls), 2) for cls in arm] for c in two]
                        + [["cohort SD", ""] + [fmt(dj["arm_ref"][cls]["sd_rel"], 2) for cls in arm]], ["sample", "DJ step"] + arm, numeric=set(range(1, len(arm) + 2)))
            tot = dj["transmitted"] + dj["not_transmitted"]
            P.h(f'''<p>Where a carrier parent and a child were both counted: the step was transmitted in <strong>{dj["transmitted"]} of {tot}</strong>
(half is the expectation for a heterozygous variant){"; " + ", ".join(esc(x) for x in dj["de_novo"]) + " carr" + ("ies" if len(dj["de_novo"]) == 1 else "y") + " a step that neither counted parent has" if dj["de_novo"] else "; no child carries a step that neither parent has"}.
This is the sharpest test the known truths offer: not that the average is right, but that a single-copy change in a
ten-copy paralogous sequence is seen in one person and then again in their child.</p>''')
    P.end()

    # ---------------------------------------------------------------- 4. modes
    P.section("modes", "Evidence 2 — the one-minute fetch reproduces the whole-file scan", "Fetch vs scan")
    if md["n_both"]:
        c45 = md["columns"].get("rDNA45S.cn_single", {})
        P.h(f'''<p class="lede">A biobank will not read 15 GB per genome. Fetch mode reads about half a gigabyte — the control regions and
the few places the aligner puts class reads — and for the {md["n_both"]:,} samples counted both ways it returns
{fmt(c45.get("median"), 4)} of the scan's 45S estimate (range {fmt(c45.get("min"), 4)}–{fmt(c45.get("max"), 4)}).</p>''')
        P.table([[d["label"], d["n"], fmt(d["median"], 4), fmt(d["min"], 4), fmt(d["max"], 4), fmt(d.get("sd_log", 0), 5)] for col, d in md["columns"].items() if d.get("n")],
                ["estimate", "n", "median fetch / scan", "min", "max", "SD of log ratio"], numeric={1, 2, 3, 4, 5})
        P.chart("fetch45", dict(type="hist", col="fetch_ratio.rDNA45S", xlabel="fetch / scan, 45S copies", ref=[dict(x=1, label="1")], xfmt=4), "45S: targeted fetch against the whole-file scan, per sample")
        P.h("<p>The known-truth and dosage columns are made from the same reads in both modes and agree exactly; the classes differ by what the sinks miss.</p>")
    else:
        P.h("<p>No sample has been counted in both modes yet; this section fills in when the per-sample jobs, which do both, land.</p>")
    cap = md["capture"]
    if any(c.get("n") for c in cap.values()):
        P.h("<p>The <em>sinks</em> — the intervals fetch mode retrieves — were learned from two whole-file scans. In every scan of this run the share of each class's reads that fell inside them is:</p>")
        P.table([[cls, c["n"], fmt(c["median"], 5), fmt(c["min"], 5), ", ".join(c["below_99"][:10]) + (" …" if len(c["below_99"]) > 10 else "") or "none"] for cls, c in cap.items() if c.get("n")],
                ["class", "scans", "median capture", "minimum", "below 99%"], numeric={1, 2, 3})
        P.h("<p class=\"small\">A sample below 99% would mean its aligner put class reads somewhere the sinks do not cover; the cohort re-learns the sinks from all its scans at the end (<code>03_compare_modes.sh</code>).</p>")
    P.end()

    # ---------------------------------------------------------------- 5. trios
    P.section("trios", "Evidence 3 — inheritance", "Inheritance")
    P.h('''<p class="lede">A child's rDNA dosage is the average of the parents' plus whatever segregation adds; measurement error is not
inherited. The slope of child on midparent therefore measures the share of an estimate's variance that is real —
its <em>reliability</em> — and does so separately for every estimator, so it can decide between them. Two traits
of the same genomes whose values are <em>not</em> transmitted are run beside them as controls: held-out autosomal
sequence (no true variance, only error) and the mitochondrial and EBV content of the cell culture (large true
variance, none of it in the nuclear genome).</p>''')
    if tr["n_complete"] >= 3 and tr["table"]:
        ci = tr["n_complete"] >= 20
        P.h(f'<p><strong>{tr["n_complete"]:,} complete trios</strong> so far.' + ("" if ci else f" Bootstrap intervals and paired comparisons appear at 20 trios; slopes and reliabilities at n = {tr['n_complete']} are indicative only.") + "</p>")
        P.chart("trio", dict(type="scatter", points=[dict(x=p["mid"], y=p["c"], label=p["child"], extra=[f"father {fmt(p['f'], 0)}, mother {fmt(p['m'], 0)}", p["pop"]]) for p in tr["scatter"]],
                             xlabel="midparent copies", ylabel="child copies", identity=True, fit=True), f"Child against midparent: {esc(tr['scatter_column'])}",
                "Grey diagonal: child equals midparent. Fitted line: the midparent slope, which is the reliability times a correction for spousal correlation.")
        rows_t = []
        for t in tr["table"]:
            r_ci = f" ({fmt(t['R_lo'], 2)} to {fmt(t['R_hi'], 2)})" if "R_lo" in t else ""
            err = ("≤ " + fmt(t["error_cv_max"], 1, pct=True)) if "error_cv_max" in t else fmt(t["error_cv"], 1, pct=True)
            rows_t.append([t["label"], t["n_trios"], fmt(t["slope"], 3) + " ± " + fmt(t["slope_se"], 3), fmt(t["spousal_r"], 3), fmt(min(t["R"], 1.0), 3) + r_ci,
                           fmt(min(t["R_single"], 1.0), 3), fmt(min(t["R_mendel"], 1.0), 3), err])
        P.table(rows_t, ["estimator", "trios", "midparent slope", "spousal r", "reliability (95% CI)", "single-parent R", "Mendelian R", "error CV the interval allows"], numeric={1, 2, 3, 4, 5, 6, 7})
        P.h('<p class="small">Reliability = slope − ρ(1 − slope), with ρ the spousal correlation after centring within population, capped at 1 (a slope above 1 is noise around 1; the interval says how much). The three estimators agree when error is independent between family members. A spousal correlation far from zero means members of a family share something other than DNA (a batch), and every reliability in the table is inflated by about as much. The negative-control rows should read near zero.</p>')
        t45 = next((t for t in tr["table"] if t["column"] == "rDNA45S.cn"), None) or next((t for t in tr["table"] if t["column"].startswith("rDNA45S")), None)
        t5 = next((t for t in tr["table"] if t["column"] == "rDNA5S.cn"), None)
        if t45 and ci:
            cv = bio.get("cn45_cv")
            P.h(f'''<p><strong>What this says.</strong> The 45S copy number differs between people by a CV of {fmt(cv, 2, pct=True) if cv else "about 20%"} — far more
than the few percent by which two libraries of the same person disagree — so any competent estimator has a reliability near 1 within one
cohort and one pipeline, and the trios cannot rank estimators whose errors are all small next to that spread: the paired differences below
are the test, and they are small. What the trios establish is that the variation being measured is inherited: the 45S reliability is
{fmt(min(t45["R"], 1.0), 2)} ({fmt(t45["R_lo"], 2)} to {fmt(t45["R_hi"], 2)}), and the interval allows a measurement error of at most
{fmt(t45.get("error_cv_max"), 1, pct=True)} of a person's value. The case for the calibrated estimator over the 18S depth ratio is not
made here; it was made across library chemistries, where the ratio moved by 27% and the calibrated estimate by 2% (the pilot).</p>''')
        if t5 and ci:
            P.h(f'''<p>The 5S array reads {fmt(min(t5["R"], 1.0), 2)} ({fmt(t5["R_lo"], 2)} to {fmt(t5["R_hi"], 2)}) at {t5["n_trios"]} trios: an interval too wide to say
whether its copy number is transmitted like the 45S's or not, which is worth watching — its spread between people is as large as the 45S's,
its estimate is as precise, and it reproduced across libraries in the pilot. A tandem array whose copy number did not pass from parent to
child would be remarkable; the full cohort decides.</p>''')
        if tr["compare"]:
            P.h("<p>Paired family bootstrap of the reliability difference against the estimator used in the literature (the 18S read-depth ratio, no GC model, no calibration):</p>")
            P.table([[c["label"], fmt(c["delta"], 3), f"{fmt(c['lo'], 3)} to {fmt(c['hi'], 3)}", fmt(c["p_better"], 3)] for c in tr["compare"]],
                    ["estimator", "ΔR vs 18S ratio", "95% CI", "P(better)"], numeric={1, 2, 3})
        if data.get("trios_adjusted", {}).get("table"):
            ta = data["trios_adjusted"]
            P.h(f'<details><summary>The same, after regressing out {pcs["adjusted"]["k"]} control-region PCs</summary>')
            P.table([[t["label"], t["n_trios"], fmt(t["slope"], 3), fmt(t["spousal_r"], 3), fmt(t["R"], 3) + (f" ({fmt(t['R_lo'], 2)} to {fmt(t['R_hi'], 2)})" if "R_lo" in t else "")] for t in ta["table"]],
                    ["estimator", "trios", "midparent slope", "spousal r", "reliability (95% CI)"], numeric={1, 2, 3, 4})
            P.h("</details>")
    else:
        P.h(f'<p>{tr["n_complete"]} complete trio(s) among the samples counted so far (the cohort has {tr["n_total"]:,}); the transmission analysis appears at three, its confidence intervals at twenty.</p>')
    P.end()

    # ---------------------------------------------------------------- 5b. what the model removes
    P.section("gcmodel", "Evidence 4 — what the model removes is the library, not the person", "GC model")
    if gcb.get("flat_vs_gc", {}).get("n", 0) >= 10:
        fv, mv, g65 = gcb["flat_vs_gc"], gcb["modelled_vs_gc"], gcb["gc65"]
        P.h(f'''<p class="lede">Every library has its own GC bias, even within one chemistry: here the rate at which 65%-GC fragments were
sequenced, relative to the library's mean, runs from {fmt(g65.get("q10"), 2)} to {fmt(g65.get("q90"), 2)} across the middle 80% of samples. The rDNA
is GC-rich. A depth ratio that ignores this reads a library property as copy number.</p>
<p>The 18S read-depth ratio as the literature computes it, divided by the calibrated estimate of the same sample, follows the library's GC bias with
<strong>r = {fmt(fv.get("r"), 2)}</strong> ({fmt(fv.get("r_lo"), 2)} to {fmt(fv.get("r_hi"), 2)}; n = {fv["n"]:,}). The same 18S region under the fragment-GC
model: r = {fmt(mv.get("r"), 2)} ({fmt(mv.get("r_lo"), 2)} to {fmt(mv.get("r_hi"), 2)}). Within one cohort the effect is a few percent (SD of the log ratio
{fmt(gcb["flat_sd_log"], 3)}), small next to the {fmt(bio.get("cn45_cv"), 0, pct=True)} by which people differ — which is why every estimator has the same
transmission reliability here — but it is systematic, and between sequencing technologies it is not small: in the pilot the same twelve people, sequenced years
apart on different instruments, differed by 27% on the depth ratio and by 2% on the calibrated estimate.</p>''')
        P.h('<div class="grid2">')
        P.chart("gcflat", dict(type="scatter", x="gc_rel_65", y="rDNA45S.18S.flat_over_cn", xlabel="library: rate at 65% GC relative to its mean", ylabel="18S depth ratio / calibrated estimate", fit=True),
                "Uncorrected 18S ratio against the library's GC bias", f"r = {fmt(fv.get('r'), 2)}: the ratio rises with the library's appetite for GC-rich fragments.")
        P.chart("gcmodelled", dict(type="scatter", x="gc_rel_65", y="rDNA45S.18S_over_cn", xlabel="library: rate at 65% GC relative to its mean", ylabel="18S under the GC model / calibrated estimate", fit=True),
                "The same region under the fragment-GC model", f"r = {fmt(mv.get('r'), 2)}: what is left does not follow the library.")
        P.h("</div>")
    else:
        P.h("<p>Appears at ten samples.</p>")
    P.end()

    # ---------------------------------------------------------------- 6. rDNA across people
    P.section("rdna", "rDNA dosage across people", "rDNA")
    rd = data["rdna"]
    c45 = rd["rDNA45S.cn"] if rd["rDNA45S.cn"].get("n") else rd["rDNA45S.cn_single"]
    col45 = "rDNA45S.cn" if rd["rDNA45S.cn"].get("n") else "rDNA45S.cn_single"
    P.h(f'''<p class="lede">The 45S array holds <strong>{fmt(c45.get("median"), 0)}</strong> copies per diploid genome in the median person
(10–90%: {fmt(c45.get("q10"), 0)}–{fmt(c45.get("q90"), 0)}; range {fmt(c45.get("min"), 0)}–{fmt(c45.get("max"), 0)}; n = {c45.get("n", 0):,}),
the 5S array {fmt(rd["rDNA5S.cn"].get("median"), 0)} ({fmt(rd["rDNA5S.cn"].get("q10"), 0)}–{fmt(rd["rDNA5S.cn"].get("q90"), 0)}).</p>''')
    P.h('<div class="grid2">')
    P.chart("cn45", dict(type="hist", col=col45, xlabel="45S copies per diploid genome", xfmt=0), "45S rDNA copy number" + (" (cohort-calibrated)" if col45 == "rDNA45S.cn" else ""))
    P.chart("cn5", dict(type="hist", col="rDNA5S.cn" if rd["rDNA5S.cn"].get("n") else "rDNA5S.cn_single", xlabel="5S copies per diploid genome", xfmt=0), "5S rDNA copy number")
    P.h("</div>")
    if superpop_order:
        P.chart("pop", dict(type="strip", col=col45, by="superpop", order=superpop_order, labels=SUPERPOP_NAMES, ylabel="45S copies"), "45S copy number by super-population",
                "Each dot one sample; the bar is the median. Population differences in the mean are expected from earlier work; whether they survive adjustment for technical structure is a question the cohort answers.")
        if len(rd["by_pop"]) > 1:
            P.h("<details><summary>By population</summary>")
            P.table([[g["group"], g["n"], fmt(g["median"], 0), fmt(g.get("q10"), 0) + "–" + fmt(g.get("q90"), 0)] for g in rd["by_pop"]], ["population", "n", "median 45S", "10–90%"], numeric={1, 2, 3})
            P.h("</details>")
    cx = bio["cn45_vs_5S"]
    if cx.get("n", 0) >= 3:
        P.chart("c45v5", dict(type="scatter", x=col45, y="rDNA5S.cn" if rd["rDNA5S.cn"].get("n") else "rDNA5S.cn_single", xlabel="45S copies", ylabel="5S copies", fit=True),
                "45S against 5S", f"n = {cx['n']:,}: Pearson r = {fmt(cx.get('r'), 2)} ({fmt(cx.get('r_lo'), 2)} to {fmt(cx.get('r_hi'), 2)}), Spearman {fmt(cx.get('spearman'), 2)}. "
                "Gibbons et al. (2015) reported the two arrays' copy numbers to be correlated; Hall et al. (2021) did not see it in these same genomes. This is the test, with the two classes measured the same way and no shared denominator artefact.")
    P.end()

    # ---------------------------------------------------------------- 7. the cell line
    P.section("cells", "The cell line the DNA came from", "Cell line")
    P.h(f'''<p class="lede">These genomes were sequenced from lymphoblastoid cell lines, and a culture has a state: how many mitochondrial
genomes per cell, how many copies of the Epstein–Barr virus that immortalised it, how large a fraction of its cells
are replicating. All three leave marks on coverage, and the same counts measure two of them directly.</p>
<p>Mitochondrial genomes per cell: median <strong>{fmt(bio["chrM"].get("median"), 0)}</strong> (10–90%: {fmt(bio["chrM"].get("q10"), 0)}–{fmt(bio["chrM"].get("q90"), 0)}).
EBV episomes per cell: median <strong>{fmt(bio["chrEBV"].get("median"), 1)}</strong> ({fmt(bio["chrEBV"].get("q10"), 1)}–{fmt(bio["chrEBV"].get("q90"), 1)}).
Both vary far more between people than the rDNA does, and neither is inherited through the nuclear genome, which is why they sit in the transmission table as controls.</p>''')
    P.h('<div class="grid2">')
    P.chart("chrM", dict(type="hist", col="chrM.copies", xlabel="mitochondrial genomes per cell", xfmt=0), "Mitochondrial DNA content")
    P.chart("ebv", dict(type="hist", col="chrEBV.copies", xlabel="EBV episomes per cell", xfmt=0), "EBV load")
    P.h("</div>")
    cm, sph = bio["log_cn45_vs_log_chrM"], bio["DJ_vs_chrX_female"]
    P.h('<div class="grid2">')
    if cm.get("n", 0) >= 3:
        P.chart("c45vM", dict(type="scatter", x="chrM.copies", y=col45, xlabel="mitochondrial genomes per cell", ylabel="45S copies", log="xy", fit=True),
                "45S copy number against mitochondrial content", f"On log scales, n = {cm['n']:,}: r = {fmt(cm.get('r'), 2)} ({fmt(cm.get('r_lo'), 2)} to {fmt(cm.get('r_hi'), 2)}). Gibbons et al. (2014) reported rDNA copy number to be coupled with mitochondrial DNA abundance in lymphoblastoid lines; here the two are measured from the same reads, with the rDNA under a GC model and the mitochondrial genome as copies per cell.")
    if sph.get("n", 0) >= 3:
        P.chart("sphase", dict(type="scatter", x="truth.chrX", y=kt["DJ_col"], where=dict(sex_inferred="F"), xlabel="chrX copies (women)", ylabel="distal junction copies", fit=False, xref=2, yref=10),
                "Two late-replicating controls, in women", f"Women whose culture has not lost an X (chrX between 1.85 and 2.15), n = {sph['n']:,}: r = {fmt(sph.get('r'), 2)} ({fmt(sph.get('r_lo'), 2)} to {fmt(sph.get('r_hi'), 2)}). The inactive X and the acrocentric short arms both replicate late; DNA from a culture with more cells in S phase should under-represent both together. If they moved together across people, that would be the S-phase fraction of the culture showing — and a covariate that rDNA, also late-replicating, would need. An r near zero says the two deficits are not one thing.")
    P.h("</div>")
    du = bio["dup"]
    P.chart("dup", dict(type="scatter", x="ctrl_dup_frac", y="rDNA45S.dup_flag_frac", xlabel="duplicate-flagged, control reads", ylabel="duplicate-flagged, 45S reads", identity=True, xfmt=3),
            "The duplicate flag inside and outside the rDNA", f"Median {fmt(du['control'].get('median'), 3, pct=True)} of control reads carry the duplicate flag, {fmt(du['rDNA'].get('median'), 3, pct=True)} of 45S reads (ratio {fmt(du['ratio'].get('median'), 2)}, range {fmt(du['ratio'].get('min'), 2)}–{fmt(du['ratio'].get('max'), 2)}). A pipeline that drops flagged reads under-reads the rDNA by a different amount in every sample; NGS-DOSE counts all primary reads.")
    P.end()

    # ---------------------------------------------------------------- 8. published values
    hall = data.get("hall")
    P.section("published", "Against published values for the same files", "Published values")
    if hall and hall.get("n", 0) >= 3:
        P.h(f'''<p class="lede">Hall, Turner &amp; Queitsch (2021) estimated rDNA copy number for 2,419 of these genomes from the same CRAMs,
as the read depth of the 18S gene relative to chromosome 1, duplicates excluded, per haploid genome.</p>
<p>On the {hall["n"]:,} samples shared so far, their 18S value against our 18S depth ratio with no GC model (halved to their scale):
r = <strong>{fmt(hall["flat"].get("r"), 3)}</strong>, their values {fmt(hall["flat_ratio"], 3)}× ours. Re-applying the duplicate-flag exclusion to
our counts brings the ratio to {fmt(hall["dup_corrected_ratio"], 3)} (SD {fmt(hall["dup_corrected_ratio_sd"], 3)}): most of the offset between the two pipelines
is the duplicate flag. Against the calibrated estimate ({esc(hall["calibrated_column"])}/2): r = {fmt(hall["calibrated"].get("r"), 3)}, ratio {fmt(hall["calibrated_ratio"], 3)}.</p>''')
        P.chart("hall", dict(type="scatter", points=[dict(x=p["cal"], y=p["theirs"], label=p["sample"]) for p in hall["points"]], xlabel="NGS-DOSE, calibrated 45S / 2", ylabel="Hall et al. 2021, 18S", identity=True, fit=True),
                "Same CRAMs, two pipelines", "Per haploid genome, as they report it.")
    else:
        P.h("<p>Appears when samples in Hall et al.'s table (their Supplementary Data 1) have been counted; the table covers the 2,504 unrelated samples.</p>")
    P.end()

    # ---------------------------------------------------------------- 9. dispersed sequence
    P.section("satellites", "Dispersed sequence: satellite families and the telomeric repeat", "Satellites")
    if sat["classes"]:
        P.h('''<p class="lede">Satellite arrays are spread over the alignment, so only a whole-file scan measures them — which is why the scan carries
experimental panels for ten families and the telomeric repeat. What each panel can see was measured on the genome it was
built from; here is what it sees in people.</p>''')
        rows_s = [[cls, d["n"], fmt(d["median"], 1), fmt(d.get("q10"), 1) + "–" + fmt(d.get("q90"), 1), fmt(d["by_sex"]["M"].get("median"), 1), fmt(d["by_sex"]["F"].get("median"), 1)] for cls, d in sat["classes"].items()]
        P.table(rows_s, ["class", "n", "median Mb (diploid)", "10–90%", "men", "women"], numeric={1, 2, 3, 4, 5})
        P.h('<div class="grid2">')
        P.chart("hsat3", dict(type="hist", col="HSat3.mass_Mb", xlabel="Mb per diploid genome", xfmt=0), "HSat3")
        P.chart("hsat1b", dict(type="hist", col="HSat1B.mass_Mb", group=dict(col="sex_inferred", levels=sex_levels), xlabel="Mb per diploid genome", xfmt=1), "HSat1B, by sex (the family lives mostly on Yq)")
        P.chart("ahor", dict(type="hist", col="aSatHOR.mass_Mb", xlabel="Mb per diploid genome", xfmt=0), "α-satellite higher-order repeats")
        P.chart("tel", dict(type="hist", col="TEL.mass_Mb", xlabel="Mb of (TTAGGG)n-bearing reads, diploid", xfmt=2), "Telomeric repeat (a relative measure)")
        P.h("</div>")
        hp = sat.get("hprc")
        if hp:
            P.h(f'<p>{hp["n_samples"]} of these samples have an HPRC release-2 assembly. An assembly is a truth only for the arrays it spans, so a sample is compared in a class only when the arrays it did not close are immaterial:</p>')
            st = hp["stats"]
            P.table([[cls, s["n"], s.get("n_gapped", 0), fmt(s.get("ratio_median"), 2), fmt(s.get("sd_log"), 3), fmt(s.get("pearson"), 2), fmt(s.get("spearman"), 2)] for cls, s in st.items() if s.get("n")],
                    ["class", "samples", "left out (gaps)", "median estimate / assembly", "SD of log ratio", "Pearson r", "Spearman"], numeric={1, 2, 3, 4, 5, 6})
            pts = [dict(x=r["assembly_Mb"], y=r["ngsdose_Mb"], label=r["sample"], extra=[r["cls"]]) for r in hp["rows"] if r["assembly_gapped_Mb"] <= 0.02 * (r["assembly_Mb"] + r["assembly_gapped_Mb"]) and r["assembly_Mb"] > 0 and r["ngsdose_Mb"] > 0]
            good = [cls for cls, st_ in st.items() if st_.get("n", 0) >= 4 and st_.get("pearson", 0) >= 0.95]
            P.chart("hprc", dict(type="scatter", points=pts, xlabel="assembly, Mb (both haplotypes)", ylabel="NGS-DOSE, Mb", identity=True, log="xy"), "Every class, every sample with an assembly",
                    "Log scales; the diagonal is agreement. The relative measures (β-satellite, CER, ACRO) sit below it by a constant factor, their k-mer recall."
                    + (f" What matters for association work is whether the estimates track the assemblies across people: {', '.join(good)} do (r ≥ 0.95)." if good else ""))
        else:
            P.h('<p class="small">The comparison with HPRC assemblies appears when their CenSat annotations are given (<code>--censat</code>); 200 samples of the cohort have one.</p>')
    else:
        P.h("<p>No scan-mode counts with the experimental panels yet.</p>")
    P.end()

    # ---------------------------------------------------------------- 10. coverage PCs
    P.section("pcs", "Technical structure: coverage PCs and how many to remove", "Coverage PCs")
    ctrl = pcs.get("control") or {}
    if ctrl.get("describe"):
        P.h(f'''<p class="lede">The 800 control regions' residual depth, after the GC model, carries whatever library and sample structure is
left; its principal components are technical covariates computed on sequence disjoint from every class. How many to regress
out is decided at the edge of the noise bulk of their spectrum — and then checked against the known truths and the trios.</p>
<p>{esc(ctrl["describe"])}.</p>''')
        var = ctrl.get("variance", [])
        if var:
            P.chart("scree", dict(type="lines", series=[dict(name="variance explained", x=list(range(1, len(var) + 1)), y=var)], xlabel="component", ylabel="fraction of variance"), "Control-region PCs: variance explained")
        adj = pcs.get("adjusted")
        if adj:
            P.h(f"<p>Regressing out the {adj['k']} components above the edge removes this much of each estimate's variance (log scale), against what {adj['k']} random regressors would remove by chance:</p>")
            P.table([[c, fmt(v["r2"], 3, pct=True), fmt(v["chance"], 3, pct=True), fmt(v["r2_adj"], 3, pct=True)] for c, v in adj["columns"].items()], ["column", "variance removed", "expected by chance", "adjusted R²"], numeric={1, 2, 3})
        sw = pcs.get("sweep")
        if sw:
            rec = sw["recommend"]
            P.h(f"<p>The sweep: 0 to {sw['max_pc']} control PCs regressed out, cross-validated. The known truths say when adjustment stops removing noise; transmission ({sw['n_trios']} trios) says when it starts removing signal. By the one-standard-error rule the picks are: "
                + ", ".join(f"<strong>{esc(c)}: {v['pick']}</strong>" for c, v in rec.items()) + ".</p>")
            series = []
            for col in ("truth.auto", "truth.chrX", "truth.chrY", kt["DJ_col"]):
                rr = [r for r in sw["rows"] if r["column"] == col and "sd_log_robust" in r]
                if rr:
                    base = rr[0]["sd_log_robust"] or 1
                    series.append(dict(name=col, x=[r["n_pc"] for r in rr], y=[r["sd_log_robust"] / base for r in rr]))
            if series:
                P.chart("sweep_truth", dict(type="lines", series=series[:3], xlabel="control PCs regressed out", ylabel="error relative to none", ref=1), "Known truths: cross-validated error against the number of PCs", "Below 1 the PCs remove error; a curve that keeps falling is what in-sample fitting would show, and does not happen out of fold.")
            series = []
            for col in ("rDNA45S.cn", "rDNA45S.cn_single", "rDNA45S.18S.flat"):
                rr = [r for r in sw["rows"] if r["column"] == col and "R_midparent" in r]
                if rr:
                    series.append(dict(name=col, x=[r["n_pc"] for r in rr], y=[r["R_midparent"] for r in rr], lo=[r.get("R_lo", r["R_midparent"]) for r in rr], hi=[r.get("R_hi", r["R_midparent"]) for r in rr]))
            if series:
                P.chart("sweep_R", dict(type="lines", series=series, xlabel="control PCs regressed out", ylabel="transmission reliability"), "Transmission reliability against the number of PCs", "Bands: family-bootstrap 95% intervals.")
        ng = pcs.get("ngspca")
        if ng:
            P.h(f'<p class="small">NGS-PCA coverage PCs: {esc(ng["describe"])}; {ng["n_with_pcs"]:,} of the samples have them.</p>')
    else:
        P.h("<p>Appears at ten samples (the components need a cohort); the sweep at sixty.</p>")
    P.end()

    # ---------------------------------------------------------------- 11. every sample
    P.section("samples", "Every sample", "Samples")
    cols = [("sample", "sample"), ("pop", "pop"), ("sex", "sex (ped)"), ("sex_inferred", "sex (reads)"), ("depth", "depth"), ("truth.auto", "auto"), ("truth.chrX", "chrX"),
            ("truth.chrY", "chrY"), (kt["DJ_col"], "DJ"), (col45, "45S"), ("rDNA45S.cn_single", "45S single"), ("rDNA45S.18S.flat", "18S flat"),
            ("rDNA5S.cn" if rd["rDNA5S.cn"].get("n") else "rDNA5S.cn_single", "5S"), ("chrM.copies", "chrM"), ("chrEBV.copies", "EBV"), ("fetch_ratio.rDNA45S", "fetch/scan 45S"),
            ("HSat3.mass_Mb", "HSat3 Mb"), ("aSatHOR.mass_Mb", "αSat Mb"), ("TEL.mass_Mb", "TEL Mb"), ("flags", "flags")]
    cols = [(c, h) for c, h in cols if any(r.get(c) not in (None, "") for r in rows)]
    nd = {"depth": 1, "truth.auto": 3, "truth.chrX": 3, "truth.chrY": 3, kt["DJ_col"]: 2, col45: 0, "rDNA45S.cn_single": 0, "rDNA45S.18S.flat": 0,
          "rDNA5S.cn": 0, "rDNA5S.cn_single": 0, "chrM.copies": 0, "chrEBV.copies": 1, "fetch_ratio.rDNA45S": 4, "HSat3.mass_Mb": 1, "aSatHOR.mass_Mb": 1, "TEL.mass_Mb": 3}
    body = []
    for r in sorted(rows, key=lambda r: r["sample"]):
        body.append([("–" if r.get(c) in (None, "") or (isinstance(r.get(c), float) and not math.isfinite(r[c])) else (fmt(r[c], nd[c]) if c in nd else r[c])) for c, _ in cols])
    P.h(f"<p>{len(rows):,} samples; click a heading to sort, type to filter. Flagged rows are shaded. The full table with every column is <code>data/cohort.tsv</code>.</p>")
    P.table(body, [h for _, h in cols], numeric={i for i, (c, _) in enumerate(cols) if c in nd}, flagged=lambda r: bool(r[-1] and r[-1] != "–") if cols[-1][0] == "flags" else None, filter_box=True, wrap=True)
    P.end()

    # ---------------------------------------------------------------- 12. caveats
    P.section("caveats", "What this page does not show, and what to hold against it", "Caveats")
    P.h(f'''<ul>
<li><strong>No absolute calibration.</strong> No orthogonal assay of rDNA copy number exists for these samples. The absolute level
rests on windows of the unit on which three Illumina chemistries were shown to agree; the known truths test the model, not the
absolute scale of the rDNA specifically.</li>
<li><strong>Cell-line DNA.</strong> Every sample here is a lymphoblastoid line, and its culture's replication state, EBV load and
mitochondrial content are measured but not removed. Blood-derived biobank genomes will not carry the first of these.</li>
<li><strong>Trios bound reliability from above</strong> where members of a family were prepared together; the spousal correlation
is the check, and the negative-control rows are the calibration of the table itself. And because the rDNA varies so much more between
people than any estimator errs, trios say <em>that</em> the measured variation is real, not which estimator measures it best.</li>
<li><strong>The satellite panels</strong> were built from one genome (CHM13). Four of the ten families are relative measures, under-read
by a known factor; HSat2 is unjudged until assemblies without gaps in it have been compared.</li>
<li><strong>The telomere class</strong> is a relative measure of (TTAGGG)n content, not a telomere length.</li>
<li><strong>Partial cohort.</strong> {n:,} of {total:,}: population comparisons and the number of complete trios depend on which samples have
landed, not on the cohort.</li>
</ul>''')
    P.end()

    # ---------------------------------------------------------------- 13. reproduce
    P.section("reproduce", "Reproduce")
    P.h(f'''<p>The counts files under <code>counts_scan/</code> and <code>counts_fetch/</code> are the primary data: ~240 kB per whole-file scan, ~70 kB per fetch,
no reads, no genotypes. From them, everything above:</p>
<pre>pip install ngsdose        # or: apptainer pull ngs-dose.sif docker://ghcr.io/jlanej/ngs-dose:latest
ngsdose report --scan counts_scan/ --fetch counts_fetch/ -p pedigree.txt --hall hall2021_MOESM1.txt -o docs/</pre>
<p>Tables behind every figure: <code>data/cohort.tsv</code> (one row per sample, every column), <code>data/modes.tsv</code>, <code>data/transmission.tsv</code>,
<code>data/pcsweep.tsv</code>, <code>data/satellites_hprc.tsv</code>, <code>data/flags.tsv</code>; the numbers in the prose: <code>report.json</code>.
The counts were made by <code>ngs-dose count</code> ({eng}) from the 1000 Genomes 30× CRAMs (Byrska-Bishop et al., <em>Cell</em> 2022; AWS Open Data),
with the {esc(m.get("bundle"))} resource bundle. Method: <a href="https://github.com/jlanej/NGS-DOSE">github.com/jlanej/NGS-DOSE</a>.</p>''')
    P.end()

    toc = "".join(f'<a href="#{i}">{esc(t)}</a>' for i, t in P.toc)
    nav = f'<nav class="toc">{toc}<button id="theme" type="button" title="light / dark">theme</button></nav>'
    payload = json.dumps(dict(charts=P.charts, samples=data["samples"]), separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    css = (ASSETS / "report.css").read_text()
    js = (ASSETS / "report.js").read_text()
    body = "".join(P.parts)
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(m["title"])}</title>
<meta name="description" content="rDNA copy number and other multi-copy sequence measured from short-read genomes on the 1000 Genomes cohort, with the evidence that it works, recomputed as the run proceeds.">
<style>{css}</style></head>
<body><main>{body[:body.index("<section")]}{nav}{body[body.index("<section"):]}
<footer>Generated {esc(m["as_of"])} by {esc(m["generator"])}. NGS-DOSE was developed by Claude (Anthropic) with @jlanej; the data are 1000 Genomes open-access.</footer>
</main>
<script id="report-data" type="application/json">{payload}</script>
<script>{js}</script>
</body></html>'''
