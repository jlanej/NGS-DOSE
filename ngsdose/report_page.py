"""The report page: prose and layout around the numbers `report.build` computed.

Every number in the text comes from `data`; nothing is typed in. Sections that have nothing to
show yet say so and say what would make them appear, so that the page can be published at any
stage of a cohort run. The order is that of a paper: summary, rationale, methods, the validation
results in the order a sceptic would ask for them, the descriptive results, limitations, data.
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
    """'mean ± SD' from a describe() dict."""
    if not d or not d.get("n"):
        return "–"
    s = fmt(d[key], nd)
    if d.get("sd") is not None:
        s += f" ± {fmt(d['sd'], nd)}"
    return s


def ci(d: dict, nd=2) -> str:
    """'r (lo to hi)' from a corr() dict."""
    if not d or d.get("r") is None:
        return "–"
    return f"{fmt(d['r'], nd)} ({fmt(d.get('r_lo'), nd)} to {fmt(d.get('r_hi'), nd)})"


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
    rep, hall, rd, nq = data.get("replicates") or {}, data.get("hall") or {}, data["rdna"], data.get("ngspca_qc") or {}
    qc_ok = nq.get("n", 0) >= 3 and nq.get("mtdna", {}).get("r") is not None
    n, total = m["n"], m["total"]
    P = Page(data)
    sex_levels = [["M", "male"], ["F", "female"]]
    superpop_order = [s for s in SUPERPOPS if any(r.get("superpop") == s for r in rows)]
    a, X, Y, DJ, sx = kt["auto"], kt["chrX"], kt["chrY"], kt["DJ"], kt["sex"]
    dj = kt.get("DJ_steps") or {}
    near = dj.get("near") or {}
    nr = lambda k: near.get(k, near.get(str(k), 0))
    n_off = nr(-2) + nr(-1) + nr(1) + nr(2)
    t45 = next((t for t in tr["table"] if t["column"] == "rDNA45S.cn"), None) or next((t for t in tr["table"] if t["column"].startswith("rDNA45S")), None)
    t5 = next((t for t in tr["table"] if t["column"] == "rDNA5S.cn"), None)
    have_ci = tr["n_complete"] >= 20
    gcb = bio.get("gc_bias") or {}
    hpst = (sat.get("hprc") or {}).get("stats") or {}
    tracking = [c for c, st_ in hpst.items() if st_.get("n", 0) >= 4 and st_.get("pearson", 0) >= 0.95]
    col45 = "rDNA45S.cn" if rd["rDNA45S.cn"].get("n") else "rDNA45S.cn_single"
    col5 = "rDNA5S.cn" if rd["rDNA5S.cn"].get("n") else "rDNA5S.cn_single"
    c45 = rd[col45]
    modes45 = md["columns"].get("rDNA45S.cn_single", {}) if md["n_both"] else {}
    rt = rep.get("table") or {}
    rc, rf, rfc = rt.get("calibrated"), rt.get("flat"), rt.get("flat_centred")
    women_ok = sx.get("women_intact", {}).get("n")
    eng = ", ".join(f"{esc(k)} ({v:,})" for k, v in sorted(m["engines"].items(), key=lambda kv: -kv[1]))

    # ---------------------------------------------------------------- masthead
    P.h(f'''<header class="mast"><h1>{esc(m["title"])}</h1>
<p class="sub">Ribosomal DNA copy number from short-read whole-genome sequencing: the method and its validation on the 1000 Genomes
30× cohort, updated as the run proceeds.</p>
<div class="hero"><div class="n">{n:,}</div><div class="of">of {total:,} genomes {"scanned" if m["primary_mode"] == "scan" else "counted"}
&middot; {m["n_fetch"]:,} also fetched &middot; {tr["n_complete"]:,} of {tr["n_total"]:,} trios complete</div></div>
<div class="progress"><div style="width:{100 * n / max(total, 1):.1f}%"></div></div>
<div class="stamp">As of {esc(m["as_of"])}. Every number and figure on this page is recomputed from the counts files in this repository by
<code>ngsdose report</code> ({esc(m["generator"])}); nothing is typed in. Partial results are published as they stand.</div></header>''')

    # ---------------------------------------------------------------- summary
    P.section("summary", "Summary")
    s = ['<div class="summary"><p>NGS-DOSE measures the copy number of the 45S and 5S ribosomal DNA arrays, and of other sequence that reference genomes '
         'collapse, from aligned short-read genomes. Reads are assigned to a sequence class by class-specific 31-mers and counted as fragment ends; a '
         'per-library model of fragment-GC bias, fitted on 800 single-copy control regions, gives the expected count of any sequence, and the windows '
         'of the rDNA unit are calibrated across the cohort. No orthogonal assay of rDNA copy number exists for these samples, so the measurement is '
         'validated against sequence of known copy number in every genome, inheritance in trios, the same individuals sequenced on two technologies, '
         'and independent measurements of the same files.</p><p>']
    r_ = []
    if a.get("n"):
        r_.append(f'In {n:,} genomes, held-out autosomal sequence reads {pm(a)} copies (expected 2)')
        if women_ok:
            r_[-1] += (f'; chrX reads {pm(X["M"])} in men and {pm(sx["women_intact"])} in women with an intact culture (expected 1 and 2); '
                       f'chrY reads {pm(sx["men_intact_Y"])} in men and at most {fmt(Y["F"].get("max"), 3)} in women (expected 1 and 0)')
        r_[-1] += "."
    if sx.get("n_pedigree"):
        r_.append(f'Sex inferred from the reads agrees with the pedigree in {sx["n_inferred"] - len(sx["mismatch"]):,} of {sx["n_inferred"]:,}.')
    if DJ.get("n"):
        t = f'The distal junction, present once on each acrocentric short arm, reads {pm(DJ)} (expected 10)'
        if dj.get("carriers") is not None and n_off:
            tot = dj["transmitted"] + dj["not_transmitted"]
            t += (f'; departures from the cohort\'s level are whole copies ({n_off} carriers)'
                  + (f', transmitted in {dj["transmitted"]} of {tot} carrier-parent–child pairs and de novo in {len(dj["de_novo"])}' if tot else ""))
        r_.append(t + ".")
    if md["n_both"]:
        r_.append(f'The targeted fetch returns {fmt(modes45.get("median"), 4)} of the whole-file scan\'s 45S estimate ({md["n_both"]:,} genomes; range {fmt(modes45.get("min"), 4)}–{fmt(modes45.get("max"), 4)}).')
    r_.append("</p><p>")
    if t45 and have_ci:
        neg = {t["column"]: t for t in tr["table"]}
        negs = ", ".join(f"{lab} {fmt(neg[c]['R'], 2)}" for c, lab in (("truth.auto", "held-out autosomal sequence"), ("chrM.copies", "mitochondrial content")) if c in neg)
        r_.append(f'In {t45["n_trios"]} trios the calibrated 45S estimate has a transmission reliability of {fmt(min(t45["R"], 1.0), 2)} ({fmt(t45["R_lo"], 2)}–{fmt(t45["R_hi"], 2)})'
                  + (f'; traits that are not transmitted read near zero ({negs})' if negs else "") + ".")
    if rc and rf:
        r_.append(f'In {rep["n"]} individuals sequenced on two technologies, the test–retest intraclass correlation is {fmt(rc["icc"], 2)} for the calibrated estimate and {fmt(rf["icc"], 2)} for the 18S depth ratio ({fmt(rfc["icc"], 2)} after removing its offset between technologies).')
    if gcb.get("flat_vs_gc", {}).get("n", 0) >= 10:
        r_.append(f'Within one chemistry the depth ratio follows each library\'s GC bias (r = {fmt(gcb["flat_vs_gc"]["r"], 2)}); under the fragment-GC model it does not (r = {fmt(gcb["modelled_vs_gc"]["r"], 2)}).')
    if hall.get("n", 0) >= 3:
        r_.append(f'Against the published estimates of Hall et al. (2021) for the same files, r = {fmt(hall["flat"].get("r"), 3)} on {hall["n"]:,} shared samples; their exclusion of duplicate-flagged reads accounts for the offset.')
    if qc_ok:
        r_.append(f'NGS-PCA\'s coverage-based mitochondrial copy number for the same files agrees with ours at r = {fmt(nq["mtdna"]["r"], 3)} (theirs {fmt(nq["mtdna"]["ratio"]["median"], 2)}× ours, the duplicate flag again), its chrX ratio at r = {fmt(nq["chrX"]["r"], 4)}.')
    if tracking:
        r_.append(f'{len(tracking)} of {len(hpst)} satellite families track HPRC assemblies of {(sat.get("hprc") or {}).get("n_samples", 0)} of these individuals with r ≥ 0.95.')
    s.append(" ".join(r_) + "</p>")
    s.append("<p>What is not shown: the absolute scale of the rDNA rests on unit windows on which three Illumina chemistries agree, not on an assay; "
             + ("5S transmission is undecided at this number of trios; " if t5 and have_ci and t5["R_lo"] < 0.5 else "")
             + "every sample is a lymphoblastoid cell line sequenced with one chemistry and one pipeline.</p></div>")
    P.h("".join(s))
    case = [("Known copy numbers, every genome", f'{fmt(a.get("mean"), 3)} ± {fmt(a.get("sd"), 3)}', f"held-out autosomal sequence, expected 2, n = {a.get('n', 0):,}", "#truth")]
    if women_ok:
        case.append(("Sex from the reads", f'{sx["n_inferred"] - len(sx["mismatch"]):,} of {sx["n_inferred"]:,}',
                     f"chrX {fmt(sx['chrX_men_max'], 2)} at most in men, {fmt(sx['women_intact']['min'], 2)} at least in women", "#truth"))
    if dj.get("carriers") is not None:
        tot = dj["transmitted"] + dj["not_transmitted"]
        case.append(("A ten-copy paralog", f'{fmt(DJ.get("median"), 2)} ± {fmt(dj.get("spread"), 2)}', f"distal junction; {n_off} people one or two copies off, steps transmitted {dj['transmitted']} of {tot}", "#djsteps"))
    if md["n_both"]:
        case.append(("Fetch = scan", fmt(modes45.get("median"), 4), f"45S, {md['n_both']:,} genomes both ways, range {fmt(modes45.get('min'), 4)}–{fmt(modes45.get('max'), 4)}", "#modes"))
    if rc and rf:
        case.append(("Two technologies", f'ICC {fmt(rc["icc"], 2)} vs {fmt(rf["icc"], 2)}', f"calibrated estimate vs 18S depth ratio, {rep['n']} people sequenced twice", "#replicates"))
    if t45 and "R_lo" in t45:
        case.append(("Inherited", fmt(min(t45["R"], 1.0), 2), f"45S transmission reliability ({fmt(t45['R_lo'], 2)}–{fmt(t45['R_hi'], 2)}), {t45['n_trios']} trios", "#trios"))
    if gcb.get("flat_vs_gc", {}).get("n"):
        case.append(("What the GC model removes", f'r {fmt(gcb["flat_vs_gc"].get("r"), 2)} → {fmt(gcb["modelled_vs_gc"].get("r"), 2)}', "how much the estimate follows the library's GC bias, before and after", "#gcmodel"))
    if hall.get("n", 0) >= 3:
        case.append(("Another pipeline, same files", f'r = {fmt(hall["flat"].get("r"), 3)}', f"Hall et al. 2021, {hall['n']:,} shared samples; offset explained by the duplicate flag", "#published"))
    if qc_ok:
        case.append(("Coverage QC, same files", f'r = {fmt(nq["mtdna"]["r"], 3)}', f"mitochondrial copies per cell vs NGS-PCA (mosdepth), {nq['n']:,} genomes; chrX r = {fmt(nq['chrX']['r'], 4)}", "#published"))
    if tracking:
        case.append(("Against assemblies", f"{len(tracking)} of {len(hpst)} families", f"track HPRC assemblies with r ≥ 0.95 in {(sat.get('hprc') or {}).get('n_samples', 0)} people", "#assemblies"))
    P.h('<div class="tiles case">' + "".join(f'<a class="tile" href="{h}"><div class="label">{esc(l)}</div><div class="value">{v}</div><div class="note">{esc(t)}</div></a>' for l, v, t, h in case) + "</div>")
    P.end()

    # ---------------------------------------------------------------- 1. rationale
    P.section("rationale", "1. Rationale", "Rationale")
    P.h('''<p>Each human genome carries several hundred copies of the 43-kb 45S ribosomal DNA unit, in tandem arrays on the short arms of
the five acrocentric chromosomes, and a tandem array of the 5S unit on chromosome 1. Copy number varies severalfold between individuals
and is heritable. It is not measured by standard genome analysis: a read from a sequence present in hundreds of copies has no unique
alignment, the reference genome holds a single token copy, and variant and copy-number callers mask these regions.</p>
<p>Hundreds of thousands of short-read genomes aligned to GRCh38 exist in biobanks and consortia. A measurement of multi-copy dosage
that is accurate, cheap enough to run on those files and does not require reading them in full would make rDNA copy number a trait
available at that scale.</p>
<p>Whether the measurement is correct cannot be settled by comparison with an assay, because none exists for these samples. It can be
settled by comparison with what is known: (i) sequence of known copy number in every sample, measured by the same code;
(ii) Mendelian transmission in the cohort's ''' + f"{tr['n_total']:,}" + ''' trios; (iii) the same individuals sequenced on different technologies;
(iv) an independent estimate from the same files; (v) long-read assemblies, for the satellite arrays measured by the same k-mer method.
Each comparison excludes a different failure. The results are presented in that order.</p>''')
    P.end()

    # ---------------------------------------------------------------- 2. methods
    P.section("methods", "2. Methods", "Methods")
    P.h(f'''<dl class="methods">
<dt>Samples</dt><dd>The expanded 1000 Genomes cohort: {total:,} individuals in 26 populations, including {tr["n_total"]:,} trios, sequenced by the New York
Genome Center to about 30× (Illumina NovaSeq, 2×150 bp, PCR-free) and aligned to GRCh38 with decoy and HLA contigs (Byrska-Bishop et al.,
<em>Cell</em> 2022). All DNA is from lymphoblastoid cell lines. Counting runs on the public CRAMs; {n:,} genomes have been counted so far,
{m["n_both"]:,} in both modes, with {tr["n_complete"]:,} of {tr["n_total"]:,} trios complete.</dd>
<dt>Class assignment</dt><dd>Every 31-mer of every read is tested against a panel of class-diagnostic k-mers: k-mers of the class's unit
sequence (45S, KY962518.1; 5S, X12811.1; distal junction, 400 kb of CHM13 chr21) that occur nowhere in GRCh38 or T2T-CHM13 outside the
class's own loci. A read with at least four panel k-mers is assigned to the class and placed on the unit by its hits; its alignment position
is not used. Experimental panels built from the CHM13 CenSat annotation add ten satellite families and the telomeric repeat.</dd>
<dt>Counting</dt><dd>Fragment 5′ ends are counted: per 50-bp bin and strand of the unit for class reads, and per position in 800 single-copy
control regions (10.1 Mb) for the library model. <em>Scan</em> mode reads the whole CRAM; <em>fetch</em> mode retrieves only the control
regions and 80 sink intervals (3.3 Mb) where the NYGC pipeline places class reads, learned from whole-file scans of two genomes.</dd>
<dt>Library model</dt><dd>A Poisson spline of fragment-end density on fragment GC content is fitted per sample on the control regions. The
expected count of any sequence follows from its fragment-GC composition; the copy number of each 250-bp window of a unit is
2 × observed / expected.</dd>
<dt>Calibration</dt><dd>Windows of the 45S unit drop out beyond what the GC curve predicts, by amounts that depend on the sequencing
chemistry. Per-window efficiencies are learned across the cohort by median polish; the absolute scale is set by anchor windows on which
three Illumina chemistries agreed in the pilot. Two comparators are carried alongside: a single-sample estimate from the anchor windows
alone, and the 18S read-depth ratio used in the literature, with no GC model and no calibration.</dd>
<dt>Known-truth controls</dt><dd>80 held-out autosomal regions (two copies), 60 chrX regions (one in men, two in women) and 40 X-degenerate
chrY regions (one, none), measured by the alignment-position path; and the distal junction, present once on each of the ten acrocentric
short arms, measured by the k-mer path with its class restricted to k-mers that occur exactly once in each of the five CHM13 junctions.
Mitochondrial genomes and EBV episomes per cell are measured as covariates of the culture.</dd>
<dt>Transmission</dt><dd>For each estimator, the regression of child on midparent in complete trios. Reliability R = b − ρ(1 − b), with
b the midparent slope and ρ the spousal correlation after centring within population; 95% intervals by family bootstrap at 20 trios or
more; paired bootstrap for differences between estimators. Held-out autosomal sequence (no true variance) and mitochondrial and EBV
content (not in the nuclear genome) are the negative controls.</dd>
<dt>Technical structure</dt><dd>Principal components of the control regions' residual depth after the GC model. The number retained is
chosen at the Marchenko–Pastur edge of the noise bulk and checked by a cross-validated sweep against the known truths and the trios.
NGS-PCA's genome-wide coverage PCs are applied where available.</dd>
<dt>External comparisons</dt><dd>The pilot's twelve genomes have an independent older library of the same cell line (HGSVC, HiSeq 2500
2×126, 2015; Illumina Platinum, HiSeq 2000 2×100, 2012–13), compared with anchor windows chosen with the family held out. Hall, Turner
&amp; Queitsch (<em>Sci Rep</em> 2021) published 18S copy number for 2,419 of these CRAMs as read depth relative to chromosome 1 with
duplicate-flagged reads excluded. HPRC release-2 assemblies of cohort samples give the size of every satellite array (CenSat annotation,
both haplotypes); a class is compared only where arrays containing gaps are immaterial. NGS-PCA's per-sample QC for the same cohort
(mosdepth, 1-kb bins, duplicate-flagged reads excluded) gives mitochondrial copies per cell, the X and Y coverage ratios and the autosomal
depth by a coverage route.</dd>
<dt>Provenance</dt><dd>Engine builds: {eng}. Resource bundle {esc(m.get("bundle"))}; panel {", ".join(esc(x) for x in m["panel_sha"]) or "–"},
controls {", ".join(esc(x) for x in m["controls_sha"]) or "–"}{(", sinks " + ", ".join(esc(x) for x in m["sinks_sha"])) if m["sinks_sha"] else ""}
(SHA-256 prefixes; one of each means one cohort). Read length {", ".join(str(x) for x in qc["read_length"])}; placement grid
{", ".join(qc["placement_bins"])} bp. {(f'<strong>{len(m["eof_absent"])} file(s) lacked an end-of-file marker</strong>: ' + ", ".join(esc(x) for x in m["eof_absent"]) + ".") if m["eof_absent"] else "Every input carried its end-of-file marker."}
{'<strong>More than one resource set is in play: these samples should not be analysed as one cohort until that is resolved.</strong>' if len(m["panel_sha"]) > 1 or len(m["controls_sha"]) > 1 else ""}</dd>
</dl>''')
    P.h("<h3>The run so far</h3>")
    P.tiles([("Genomes scanned", f"{m['n_scan']:,}", f"of {total:,}"), ("Fetched as well", f"{m['n_fetch']:,}", "targeted mode, same files"),
             ("Complete trios", f"{tr['n_complete']:,}", f"of {tr['n_total']:,}"),
             ("Median depth", fmt(qc["depth"].get("median"), 1) + "×", f"{fmt(qc['depth'].get('q10'), 1)}–{fmt(qc['depth'].get('q90'), 1)} (10–90%)"),
             ("Insert size", fmt(qc["insert"].get("median"), 0) + " bp", "median of medians"),
             ("Duplicate-flagged", fmt(qc["dup"].get("median"), 3, pct=True), "of control reads, median"),
             ("Scan time", fmt(qc["elapsed"].get("median") / 60 if qc["elapsed"].get("n") else None, 1) + " min", "per genome, median" if qc["elapsed"].get("n") else "no scans yet"),
             ("Fetch time", fmt(qc["elapsed_fetch"].get("median") / 60 if qc["elapsed_fetch"].get("n") else None, 1) + " min", "per genome, median" if qc["elapsed_fetch"].get("n") else "no fetches yet")])
    if m.get("by_superpop"):
        sp = m["by_superpop"]
        cats = [c for c in SUPERPOPS if c in sp]
        P.chart("progress", dict(type="meters", categories=[f"{c} · {SUPERPOP_NAMES[c]}" for c in cats], values=[sp[c]["done"] for c in cats], totals=[sp[c]["total"] for c in cats]),
                "Genomes counted, by super-population", "The pale track is the whole cohort.")
    if data["flags"]:
        P.h(f'<p class="small">{len(data["flags"])} sample(s) carry a flag (aneuploid chromosome, sex mismatch, mosaic loss of an X or Y, a distal-junction step, low depth, another engine build); they are marked in the <a href="#samples">sample table</a> and listed in <code>data/flags.tsv</code>. A flag marks something to examine, not a verdict.</p>')
    P.end()

    # ---------------------------------------------------------------- 3.1 known truth
    P.section("truth", "3.1 Sequence of known copy number, in every genome", "Known truth")
    P.h(f'''<p>If the model is right, held-out autosomal sequence reads 2, chrX reads 1 in men and 2 in women, chrY reads 1 and 0, and the
distal junction reads 10, in every genome. Held-out autosomal sequence reads <strong>{pm(a)}</strong> copies (n = {a.get("n", 0):,}).
chrX reads <strong>{pm(X["M"])}</strong> in {X["M"].get("n", 0):,} men'''
        + (f''' and <strong>{pm(sx["women_intact"])}</strong> in the {sx["women_intact"]["n"]:,} women whose culture has kept both X chromosomes (at least 1.85 copies); {sx["n_mosaic_X"]} women read below that.
chrY reads <strong>{pm(sx["men_intact_Y"])}</strong> in men with an intact Y and {fmt(Y["F"].get("mean"), 4)} in women (maximum {fmt(Y["F"].get("max"), 4)}); {sx["n_mosaic_Y"]} men read below 0.85.''' if women_ok
           else f''' and {pm(X["F"])} in {X["F"].get("n", 0):,} women; chrY {pm(Y["M"])} in men and {pm(Y["F"], 4)} in women.''')
        + f''' The distal junction reads <strong>{pm(DJ)}</strong>{" (cohort-calibrated)" if kt["DJ_col"] == "DJ.cn" else ""}.'''
        + (" Women read the X and every genome reads the distal junction a few percent below expectation; both are late-replicating sequence, which DNA from a growing culture under-represents (section 4)." if X["F"].get("median", 2) < 1.98 else "") + "</p>")
    P.h('<div class="grid2">')
    P.chart("auto", dict(type="hist", col="truth.auto", xlabel="copies", ref=[dict(x=2, label="expected 2")], xfmt=3), "Held-out autosomal regions (80 regions)")
    P.chart("chrX", dict(type="hist", col="truth.chrX", group=dict(col="sex_inferred", levels=sex_levels), xlabel="copies", ref=[dict(x=1, label="1"), dict(x=2, label="2")], xfmt=2), "chrX (60 regions), by sex")
    P.chart("chrY", dict(type="hist", col="truth.chrY", group=dict(col="sex_inferred", levels=sex_levels), xlabel="copies", ref=[dict(x=0, label="0"), dict(x=1, label="1")], xfmt=2), "chrY (40 X-degenerate regions), by sex")
    P.chart("dj", dict(type="hist", col=kt["DJ_col"], xlabel="copies", ref=[dict(x=10, label="expected 10")], xfmt=2), "Distal junction (one per acrocentric short arm)")
    P.h("</div>")
    if sx["n_pedigree"]:
        P.h(f'<p>Sex inferred from the reads (a Y above 0.1 copies) agrees with the pedigree in {sx["n_inferred"] - len(sx["mismatch"]):,} of {sx["n_inferred"]:,} samples'
            + (f'; it does not in <strong>{", ".join(esc(x) for x in sx["mismatch"][:20])}{" and " + str(len(sx["mismatch"]) - 20) + " more" if len(sx["mismatch"]) > 20 else ""}</strong> (a swapped sample, or a line that has lost its Y).' if sx["mismatch"] else ".") + "</p>")
    outl = [(s_, f) for s_, f in data["flags"] if "chrX" in f or "chrY" in f or "autosomal" in f]
    if outl:
        P.h(f'<details><summary>{len(outl)} sample(s) off the expected value</summary>')
        P.table([[s_, f] for s_, f in outl], ["sample", "what"])
        P.h("<p class=\"small\">A woman whose X reads well below 2, or a man whose Y does, has lost that chromosome in part of the cell culture: the known behaviour of lymphoblastoid lines, and a reason the controls are measured in every sample.</p></details>")
    P.end()

    # ---------------------------------------------------------------- 3.2 DJ steps
    P.section("djsteps", "3.2 The distal junction changes in whole copies, and the changes are inherited", "DJ steps")
    if dj.get("carriers") is not None:
        P.h(f'''<p>Ten distal junctions is the norm; a rearranged acrocentric short arm leaves nine, and a Robertsonian translocation, which fuses two
acrocentrics and loses both short arms, leaves eight. Copy number relative to the cohort's level ({fmt(dj["median"], 2)}) should therefore sit
near a whole number, and a step, being a structural variant, should be transmitted to half of a carrier's children and arise de novo in almost
none. Within ±0.3 of a step: <strong>{nr(-2)}</strong> genomes at −2, <strong>{nr(-1)}</strong> at −1, {nr(0):,} at 0, <strong>{nr(1)}</strong> at +1;
the main mode has a robust SD of {fmt(dj["spread"], 2)} copies and {dj["between"]} genomes sit between steps.</p>''')
        P.chart("djstep", dict(type="hist", col="DJ.step", xlabel="distal-junction copies relative to the cohort's level", ref=[dict(x=k, label=str(k)) for k in (-2, -1, 0, 1)], xfmt=1, bins=40),
                "Distal-junction copy number relative to the cohort", "Reference lines at whole copies.")
        if dj["carriers"]:
            rows_c = []
            for c in dj["carriers"]:
                rel = "; ".join(f"{r['who']} {esc(r['sample'])} {r['step']:+.2f}" for r in c["relatives"]) or "none counted"
                rows_c.append([c["sample"], c.get("pop") or "", c.get("sex") or "", f"{c['step']:+.2f}", rel])
            P.table(rows_c, ["sample", "population", "sex", "step (copies)", "relatives counted, and their step"], numeric={3})
            two = [c for c in dj["carriers"] if c["step"] <= -1.5 and c.get("arm_content")]
            if two and dj.get("arm_ref"):
                arm = [cls for cls in ("ACRO", "SST1", "bSat", "HSat3", "CER", "HSat1A", "aSatHOR") if cls in dj["arm_ref"]]
                P.h("<p>A lost short arm takes its satellite arrays with it. Satellite families of the acrocentric short arms in the two-copy carriers, as a fraction of the cohort's median; the pan-centromeric α-satellite (aSatHOR), which every chromosome carries, is the control:</p>")
                P.table([[c["sample"], f"{c['step']:+.2f}"] + [fmt(c["arm_content"].get(cls), 2) for cls in arm] for c in two]
                        + [["cohort SD", ""] + [fmt(dj["arm_ref"][cls]["sd_rel"], 2) for cls in arm]], ["sample", "DJ step"] + arm, numeric=set(range(1, len(arm) + 2)))
            tot = dj["transmitted"] + dj["not_transmitted"]
            P.h(f'''<p>Where a carrier parent and a child were both counted, the step was transmitted in <strong>{dj["transmitted"]} of {tot}</strong>
(the expectation for a heterozygous variant is one half){"; " + ", ".join(esc(x) for x in dj["de_novo"]) + " carr" + ("ies" if len(dj["de_novo"]) == 1 else "y") + " a step that neither counted parent has" if dj["de_novo"] else "; no child carries a step that neither parent has"}.
The distal junction is measured by the same k-mer path as the rDNA. A change of one copy in ten, seen in a parent and again in the child,
shows that the path resolves multi-copy acrocentric sequence to a single copy.</p>''')
    else:
        P.h("<p>Appears once the distal junction has been measured.</p>")
    P.end()

    # ---------------------------------------------------------------- 3.3 modes
    P.section("modes", "3.3 The targeted fetch against the whole-file scan", "Fetch vs scan")
    if md["n_both"]:
        P.h(f'''<p>Fetch mode reads about 0.5 GB of a 15-GB CRAM. For the {md["n_both"]:,} genomes counted both ways it returns
<strong>{fmt(modes45.get("median"), 4)}</strong> of the scan's 45S estimate (range {fmt(modes45.get("min"), 4)}–{fmt(modes45.get("max"), 4)}). The
known-truth and dosage columns are made from the same reads in both modes and agree exactly; the classes differ by what the sinks miss.</p>''')
        P.table([[d["label"], d["n"], fmt(d["median"], 4), fmt(d["min"], 4), fmt(d["max"], 4), fmt(d.get("sd_log", 0), 5)] for col, d in md["columns"].items() if d.get("n")],
                ["estimate", "n", "median fetch / scan", "min", "max", "SD of log ratio"], numeric={1, 2, 3, 4, 5})
        P.chart("fetch45", dict(type="hist", col="fetch_ratio.rDNA45S", xlabel="fetch / scan, 45S copies", ref=[dict(x=1, label="1")], xfmt=4), "45S: targeted fetch against the whole-file scan, per genome")
    else:
        P.h("<p>No genome has been counted in both modes yet; this section fills in when the per-sample jobs, which do both, land.</p>")
    cap = md["capture"]
    if any(c.get("n") for c in cap.values()):
        P.h("<p>Share of each class's reads that fell inside the sink intervals, in every whole-file scan of this run:</p>")
        P.table([[cls, c["n"], fmt(c["median"], 5), fmt(c["min"], 5), ", ".join(c["below_99"][:10]) + (" …" if len(c["below_99"]) > 10 else "") or "none"] for cls, c in cap.items() if c.get("n")],
                ["class", "scans", "median capture", "minimum", "below 99%"], numeric={1, 2, 3})
        P.h("<p class=\"small\">A genome below 99% would mean its aligner put class reads where the sinks do not reach; the sinks are re-learned from all scans at the end of the run.</p>")
    P.end()

    # ---------------------------------------------------------------- 3.4 two technologies
    P.section("replicates", "3.4 The same individuals on two sequencing technologies", "Two technologies")
    if rc and rf:
        P.h(f'''<p>The pilot's {rep["n"]} genomes were also sequenced years earlier on a different instrument, with different read length, insert size,
depth, alignment pipeline and a GC response the reverse of NovaSeq's. Agreement between the two libraries of a person, per estimator (the
calibrated estimate with anchor windows chosen with the person's family held out, so the cross-technology level is out of sample):</p>''')
        order = [k for k in ("calibrated", "flat", "flat_centred", "rDNA5S", "DJ") if k in rt]
        P.table([[rt[k]["label"], rt[k]["n"], fmt(rt[k]["offset"], 1, pct=True), fmt(rt[k]["sd_log_ratio"], 3), fmt(rt[k]["r"], 3), fmt(rt[k]["icc"], 3), fmt(rt[k]["within_cv"], 1, pct=True)] for k in order],
                ["estimator", "pairs", "offset, older / NYGC", "pair SD of log ratio", "Pearson r", "ICC", "within-person CV"], numeric={1, 2, 3, 4, 5, 6})
        notes = []
        if "rDNA5S" in rt:
            notes.append(f'The 5S unit is 68% GC throughout and has no anchor windows, so its estimate rests on the GC model alone and carries an offset of {fmt(abs(rt["rDNA5S"]["offset"]), 0, pct=True)} between technologies at a pair SD of {fmt(rt["rDNA5S"]["sd_log_ratio"], 3)}.')
        if "DJ" in rt:
            notes.append(f'The distal junction does not vary between people, so its intraclass correlation is near zero by construction; its within-person CV of {fmt(rt["DJ"]["within_cv"], 1, pct=True)} is the figure that matters.')
        if notes:
            P.h('<p class="small">' + " ".join(notes) + "</p>")
        P.chart("replicates", dict(type="scatter", points=[dict(x=p["x"], y=p["y"], label=p["sample"], si=p["si"]) for p in rep["points"]],
                                   legend=["calibrated, anchors held out", "18S depth ratio, no GC model"], xlabel="45S copies, NovaSeq 2×150 (2019)", ylabel="45S copies, HiSeq 2×100 / 2×126 (2012–15)", identity=True),
                "The same person, two technologies", "Diagonal: agreement.")
        P.h(f'''<p>The intraclass correlation measures absolute agreement, so an offset between technologies counts against it: the depth ratio's
{fmt(rf["offset"], 0, pct=True)} offset is a property of the library, not of the person, and removing it as a batch effect leaves an ICC of
{fmt(rfc["icc"], 2)} against {fmt(rc["icc"], 2)} for the calibrated estimate, whose offset is {fmt(rc["offset"], 0, pct=True)}. The two DNA batches were
drawn from different cultures of each line, so these figures bound measurement error from above. This is the comparison that distinguishes the
calibrated estimate from a depth ratio; within one chemistry (3.6) it cannot be seen.</p>''')
    else:
        P.h("<p>Appears when the pilot's replicate tables are given (<code>--pilot</code>).</p>")
    P.end()

    # ---------------------------------------------------------------- 3.5 GC model
    P.section("gcmodel", "3.5 What the fragment-GC model removes", "GC model")
    if gcb.get("flat_vs_gc", {}).get("n", 0) >= 10:
        fv, mv, g65 = gcb["flat_vs_gc"], gcb["modelled_vs_gc"], gcb["gc65"]
        P.h(f'''<p>Libraries differ in GC bias even within one chemistry: the rate at which 65%-GC fragments were sequenced, relative to each library's
mean, runs from {fmt(g65.get("q10"), 2)} to {fmt(g65.get("q90"), 2)} across the middle 80% of genomes, and the rDNA is GC-rich. The 18S depth ratio
divided by the calibrated estimate of the same genome follows that bias with <strong>r = {ci(fv)}</strong> (n = {fv["n"]:,}); the same 18S
region under the fragment-GC model, r = {ci(mv)}. Within the cohort the effect is a few percent (SD of the log ratio {fmt(gcb["flat_sd_log"], 3)}),
small next to the {fmt(bio.get("cn45_cv"), 0, pct=True)} by which people differ, which is why it does not show in the trios; between technologies
(3.4) it is {fmt(abs(rf["offset"]), 0, pct=True) if rf else "large"}.</p>''')
        P.h('<div class="grid2">')
        P.chart("gcflat", dict(type="scatter", x="gc_rel_65", y="rDNA45S.18S.flat_over_cn", xlabel="library: rate at 65% GC relative to its mean", ylabel="18S depth ratio / calibrated estimate", fit=True),
                "18S depth ratio against the library's GC bias", f"r = {fmt(fv.get('r'), 2)}.")
        P.chart("gcmodelled", dict(type="scatter", x="gc_rel_65", y="rDNA45S.18S_over_cn", xlabel="library: rate at 65% GC relative to its mean", ylabel="18S under the GC model / calibrated estimate", fit=True),
                "The same region under the fragment-GC model", f"r = {fmt(mv.get('r'), 2)}.")
        P.h("</div>")
    else:
        P.h("<p>Appears at ten genomes.</p>")
    P.end()

    # ---------------------------------------------------------------- 3.6 trios
    P.section("trios", "3.6 Transmission in trios", "Inheritance")
    P.h('''<p>A child's dosage is the mean of the parents' plus segregation; measurement error is not inherited. The midparent slope therefore
measures the share of an estimate's variance that is real, its reliability, for each estimator separately. Two traits of the same genomes
that are not transmitted run beside them as negative controls: held-out autosomal sequence, which has error but no true variance, and the
mitochondrial and EBV content of the culture, which vary but are not in the nuclear genome.</p>''')
    if tr["n_complete"] >= 3 and tr["table"]:
        P.h(f'<p><strong>{tr["n_complete"]:,} complete trios.</strong>' + ("" if have_ci else f" Bootstrap intervals and paired comparisons appear at 20 trios; slopes at n = {tr['n_complete']} are indicative only.") + "</p>")
        P.chart("trio", dict(type="scatter", points=[dict(x=p["mid"], y=p["c"], label=p["child"], extra=[f"father {fmt(p['f'], 0)}, mother {fmt(p['m'], 0)}", p["pop"]]) for p in tr["scatter"]],
                             xlabel="midparent copies", ylabel="child copies", identity=True, fit=True), f"Child against midparent: {esc(tr['scatter_column'])}",
                "Grey diagonal: child equals midparent. Fitted line: the midparent slope.")
        rows_t = []
        for t in tr["table"]:
            r_ci = f" ({fmt(t['R_lo'], 2)} to {fmt(t['R_hi'], 2)})" if "R_lo" in t else ""
            err = ("≤ " + fmt(t["error_cv_max"], 1, pct=True)) if "error_cv_max" in t else fmt(t["error_cv"], 1, pct=True)
            rows_t.append([t["label"], t["n_trios"], fmt(t["slope"], 3) + " ± " + fmt(t["slope_se"], 3), fmt(t["spousal_r"], 3), fmt(min(t["R"], 1.0), 3) + r_ci,
                           fmt(min(t["R_single"], 1.0), 3), fmt(min(t["R_mendel"], 1.0), 3), err])
        P.table(rows_t, ["estimator", "trios", "midparent slope", "spousal r", "reliability (95% CI)", "single-parent R", "Mendelian R", "error CV the interval allows"], numeric={1, 2, 3, 4, 5, 6, 7})
        P.h('<p class="small">Reliability is capped at 1; a slope above 1 is noise around 1, and the interval says how much. A spousal correlation far from zero means members of a family share something other than DNA (a batch), and every reliability in the table is inflated by about as much. The negative-control rows should read near zero.</p>')
        if t45 and have_ci:
            cv = bio.get("cn45_cv")
            P.h(f'''<p>The 45S reliability is <strong>{fmt(min(t45["R"], 1.0), 2)}</strong> ({fmt(t45["R_lo"], 2)} to {fmt(t45["R_hi"], 2)}); the interval allows a
measurement error of at most {fmt(t45.get("error_cv_max"), 1, pct=True)} of a person's value. The variation measured is therefore inherited. Because
people differ in 45S copy number by a CV of {fmt(cv, 0, pct=True) if cv else "about 20%"} and every estimator errs by a few percent, every estimator
has a reliability near 1 within one cohort and one pipeline, and the trios cannot rank them: the paired differences below are the test, and they
are small. Ranking rests on the comparison across technologies (3.4).</p>''')
        if t5 and have_ci:
            P.h(f'''<p>The 5S array reads {fmt(min(t5["R"], 1.0), 2)} ({fmt(t5["R_lo"], 2)} to {fmt(t5["R_hi"], 2)}) at {t5["n_trios"]} trios: an interval too wide to
say whether its copy number is transmitted like the 45S's. Its spread between people is as large, its estimate as precise, and it reproduced
across technologies in the pilot (3.4); the full cohort decides.</p>''')
        if tr["compare"]:
            P.h("<p>Paired family bootstrap of the reliability difference against the 18S depth ratio:</p>")
            P.table([[c["label"], fmt(c["delta"], 3), f"{fmt(c['lo'], 3)} to {fmt(c['hi'], 3)}", fmt(c["p_better"], 3)] for c in tr["compare"]],
                    ["estimator", "ΔR vs 18S ratio", "95% CI", "P(better)"], numeric={1, 2, 3})
        if data.get("trios_adjusted", {}).get("table"):
            ta = data["trios_adjusted"]
            P.h(f'<details><summary>The same, after regressing out {pcs["adjusted"]["k"]} control-region PCs</summary>')
            P.table([[t["label"], t["n_trios"], fmt(t["slope"], 3), fmt(t["spousal_r"], 3), fmt(t["R"], 3) + (f" ({fmt(t['R_lo'], 2)} to {fmt(t['R_hi'], 2)})" if "R_lo" in t else "")] for t in ta["table"]],
                    ["estimator", "trios", "midparent slope", "spousal r", "reliability (95% CI)"], numeric={1, 2, 3, 4})
            P.h("</details>")
    else:
        P.h(f'<p>{tr["n_complete"]} complete trio(s) among the genomes counted so far (the cohort has {tr["n_total"]:,}); the analysis appears at three, its confidence intervals at twenty.</p>')
    P.end()

    # ---------------------------------------------------------------- 3.7 published values
    P.section("published", "3.7 An independent pipeline on the same files", "Published values")
    if hall.get("n", 0) >= 3:
        P.h(f'''<p>On the {hall["n"]:,} genomes shared so far with Hall, Turner &amp; Queitsch (2021), their 18S value against our 18S depth ratio with no GC
model, halved to their per-haploid scale: r = <strong>{fmt(hall["flat"].get("r"), 3)}</strong>, their values {fmt(hall["flat_ratio"], 3)}× ours. Re-applying
their exclusion of duplicate-flagged reads to our counts brings the ratio to {fmt(hall["dup_corrected_ratio"], 3)} (SD {fmt(hall["dup_corrected_ratio_sd"], 3)}):
the offset is the duplicate flag, which marks fewer reads inside the collapsed rDNA than outside it (section 4). Against the calibrated
estimate ({esc(hall["calibrated_column"])}/2): r = {fmt(hall["calibrated"].get("r"), 3)}, ratio {fmt(hall["calibrated_ratio"], 3)}.</p>''')
        P.chart("hall", dict(type="scatter", points=[dict(x=p["cal"], y=p["theirs"], label=p["sample"]) for p in hall["points"]], xlabel="NGS-DOSE, calibrated 45S / 2", ylabel="Hall et al. 2021, 18S", identity=True, fit=True),
                "Same CRAMs, two pipelines", "Per haploid genome, as they report it.")
    else:
        P.h("<p>Appears when genomes in Hall et al.'s table (their Supplementary Data 1, the 2,504 unrelated samples) have been counted.</p>")
    if qc_ok:
        mt, xx, yy, dp = nq["mtdna"], nq["chrX"], nq["chrY_men"], nq["depth"]
        mos = ""
        if nq["mosaic_X"] or nq["mosaic_Y"]:
            mos = (f' The {len(nq["mosaic_X"])} women and {len(nq["mosaic_Y"])} men whose cultures have lost part of an X or a Y read the same by both routes'
                   + (f' (for example {esc(nq["mosaic_X"][0]["sample"])}: {fmt(nq["mosaic_X"][0]["ours"], 2)} here, {fmt(nq["mosaic_X"][0]["theirs"], 2)} there).' if nq["mosaic_X"] else "."))
        P.h(f'''<h3>NGS-PCA's coverage QC on the same files</h3>
<p><a href="https://github.com/jlanej/NGS-PCA">NGS-PCA</a> computes, from mosdepth coverage of the same CRAMs in 1-kb bins with duplicate-flagged reads
excluded, mitochondrial copies per cell as twice the chrM mean coverage over the median autosomal coverage, and the X and Y coverage ratios. On the
{nq["n"]:,} shared genomes, mitochondrial copies per cell agree at r = <strong>{ci(mt, 3)}</strong>, theirs {fmt(mt["ratio"]["median"], 3)}× ours
(10–90% {fmt(mt["ratio"]["q10"], 3)}–{fmt(mt["ratio"]["q90"], 3)}; SD of the log ratio {fmt(mt["ratio"]["sd_log"], 3)}); chrX at r = {fmt(xx["r"], 4)}, ratio
{fmt(xx["ratio"]["median"], 3)}; autosomal depth at r = {fmt(dp["r"], 3)}, theirs {fmt(dp["ratio"]["median"], 2)}× ours. Sex inferred by the two agrees in
{nq["sex_agree"]:,} of {nq["sex_n"]:,}.{mos} The chrY ratio in men is compressed and noisier by the coverage route (median
{fmt(nq["chrY_intact_men"].get("median"), 2)} in men with an intact Y, r = {fmt(yy["r"], 2)} against our 40 X-degenerate regions), since a whole-chromosome
mean includes sequence that maps poorly. The mitochondrial offset lies in the direction of the duplicate flag, which mosdepth honours and NGS-DOSE
does not: a 16.6-kb genome at several thousand-fold depth saturates the positions a duplicate marker can distinguish, and the ratio falls with
depth (r = {fmt(nq["mtdna_ratio_vs_depth"].get("r"), 2)}).</p>''')
        P.h('<div class="grid2">')
        P.chart("qc_mtdna", dict(type="scatter", x="chrM.copies", y="ngspca.MTDNA_CN", xlabel="NGS-DOSE, mitochondrial genomes per cell", ylabel="NGS-PCA, mtDNA copy number", identity=True, fit=True),
                "Mitochondrial copies per cell, two routes", f"r = {fmt(mt['r'], 3)}; the diagonal is equality.")
        P.chart("qc_chrX", dict(type="scatter", x="truth.chrX", y="ngspca.chrX", group=dict(col="sex_inferred", levels=sex_levels), xlabel="NGS-DOSE, chrX copies (60 regions)", ylabel="NGS-PCA, 2 × chrX coverage ratio", identity=True),
                "chrX copies, two routes", f"r = {fmt(xx['r'], 4)}; the cultures that have lost part of an X sit off the clusters in both.")
        P.h("</div>")
    P.end()

    # ---------------------------------------------------------------- 3.8 assemblies
    P.section("assemblies", "3.8 Satellite arrays against long-read assemblies", "Assemblies")
    hp = sat.get("hprc")
    if hp:
        st = hp["stats"]
        good = [cls for cls, st_ in st.items() if st_.get("n", 0) >= 4 and st_.get("pearson", 0) >= 0.95]
        P.h(f'''<p>Assemblies collapse the rDNA and are no truth for it, but {hp["n_samples"]} of these genomes have HPRC release-2 assemblies whose CenSat
annotation gives the size of every satellite array: the same kind of sequence, measured by the same k-mer machinery. A genome is compared in a
class only when the arrays its assembly did not close are immaterial.{" " + ", ".join(good) + " track the assemblies across people with r ≥ 0.95." if good else ""}</p>''')
        P.table([[cls, s_["n"], s_.get("n_gapped", 0), fmt(s_.get("ratio_median"), 2), fmt(s_.get("sd_log"), 3), fmt(s_.get("pearson"), 2), fmt(s_.get("spearman"), 2)] for cls, s_ in st.items() if s_.get("n")],
                ["class", "genomes", "left out (gaps)", "median estimate / assembly", "SD of log ratio", "Pearson r", "Spearman"], numeric={1, 2, 3, 4, 5, 6})
        pts = [dict(x=r["assembly_Mb"], y=r["ngsdose_Mb"], label=r["sample"], extra=[r["cls"]]) for r in hp["rows"] if r["assembly_gapped_Mb"] <= 0.02 * (r["assembly_Mb"] + r["assembly_gapped_Mb"]) and r["assembly_Mb"] > 0 and r["ngsdose_Mb"] > 0]
        P.chart("hprc", dict(type="scatter", points=pts, xlabel="assembly, Mb (both haplotypes)", ylabel="NGS-DOSE, Mb", identity=True, log="xy"), "Every class, every genome with an assembly",
                "Log scales; the diagonal is agreement. The relative measures (β-satellite, CER, ACRO) sit below it by a constant factor, their k-mer recall.")
    else:
        P.h('<p class="small">Appears when HPRC CenSat annotations are given (<code>--censat</code>); 200 genomes of the cohort have an assembly.</p>')
    P.end()

    # ---------------------------------------------------------------- 3.9 coverage PCs
    P.section("pcs", "3.9 Technical structure: coverage PCs", "Coverage PCs")
    ctrl = pcs.get("control") or {}
    if ctrl.get("describe"):
        P.h(f'''<p>The residual depth of the 800 control regions after the GC model carries whatever library and sample structure is left; its principal
components are technical covariates computed on sequence disjoint from every class. {esc(ctrl["describe"])}.</p>''')
        var = ctrl.get("variance", [])
        if var:
            P.chart("scree", dict(type="lines", series=[dict(name="variance explained", x=list(range(1, len(var) + 1)), y=var)], xlabel="component", ylabel="fraction of variance"), "Control-region PCs: variance explained")
        adj = pcs.get("adjusted")
        if adj:
            P.h(f"<p>Regressing out the {adj['k']} components above the edge removes this share of each estimate's variance (log scale), against what {adj['k']} random regressors would remove by chance:</p>")
            P.table([[c, fmt(v["r2"], 3, pct=True), fmt(v["chance"], 3, pct=True), fmt(v["r2_adj"], 3, pct=True)] for c, v in adj["columns"].items()], ["column", "variance removed", "expected by chance", "adjusted R²"], numeric={1, 2, 3})
            r45 = adj["columns"].get("rDNA45S.cn") or adj["columns"].get("rDNA45S.cn_single")
            au = adj["columns"].get("truth.auto")
            if r45 and au:
                P.h(f'''<p>The components find technical variance where it exists: {fmt(au["r2"], 0, pct=True)} of the held-out autosomal estimate's variance, which is
nothing but error, against {fmt(r45["r2"], 0, pct=True)} of the 45S estimate's ({fmt(r45["chance"], 1, pct=True)} expected by chance). With a measurement error of a
few percent and a between-person CV of {fmt(bio.get("cn45_cv"), 0, pct=True)}, the technical share of the 45S variance is small; the GC model and the
calibration do the work, and the PCs are insurance.</p>''')
        sw = pcs.get("sweep")
        if sw:
            rec = sw["recommend"]
            P.h(f"<p>The sweep, 0 to {sw['max_pc']} control PCs regressed out, cross-validated: the known truths say when adjustment stops removing noise; transmission ({sw['n_trios']} trios) says when it starts removing signal. One-standard-error picks: "
                + ", ".join(f"<strong>{esc(c)}: {v['pick']}</strong>" for c, v in rec.items()) + ".</p>")
            series = []
            for col in ("truth.auto", "truth.chrX", "truth.chrY", kt["DJ_col"]):
                rr = [r for r in sw["rows"] if r["column"] == col and "sd_log_robust" in r]
                if rr:
                    base = rr[0]["sd_log_robust"] or 1
                    series.append(dict(name=col, x=[r["n_pc"] for r in rr], y=[r["sd_log_robust"] / base for r in rr]))
            if series:
                P.chart("sweep_truth", dict(type="lines", series=series[:3], xlabel="control PCs regressed out", ylabel="error relative to none", ref=1), "Known truths: cross-validated error against the number of PCs", "Below 1 the PCs remove error.")
            series = []
            for col in ("rDNA45S.cn", "rDNA45S.cn_single", "rDNA45S.18S.flat"):
                rr = [r for r in sw["rows"] if r["column"] == col and "R_midparent" in r]
                if rr:
                    series.append(dict(name=col, x=[r["n_pc"] for r in rr], y=[r["R_midparent"] for r in rr], lo=[r.get("R_lo", r["R_midparent"]) for r in rr], hi=[r.get("R_hi", r["R_midparent"]) for r in rr]))
            if series:
                P.chart("sweep_R", dict(type="lines", series=series, xlabel="control PCs regressed out", ylabel="transmission reliability"), "Transmission reliability against the number of PCs", "Bands: family-bootstrap 95% intervals.")
        ng = pcs.get("ngspca")
        if ng:
            P.h(f'<p class="small">NGS-PCA coverage PCs: {esc(ng["describe"])}; {ng["n_with_pcs"]:,} of the genomes have them.</p>')
    else:
        P.h("<p>Appears at ten genomes; the sweep at sixty.</p>")
    P.end()

    # ---------------------------------------------------------------- 4. descriptive results
    P.section("rdna", "4. The measurements: rDNA copy number, the cell line, satellite arrays", "rDNA")
    P.h(f'''<p>The 45S array holds <strong>{fmt(c45.get("median"), 0)}</strong> copies per diploid genome in the median person (10–90%:
{fmt(c45.get("q10"), 0)}–{fmt(c45.get("q90"), 0)}; range {fmt(c45.get("min"), 0)}–{fmt(c45.get("max"), 0)}; n = {c45.get("n", 0):,}), the 5S array
{fmt(rd[col5].get("median"), 0)} ({fmt(rd[col5].get("q10"), 0)}–{fmt(rd[col5].get("q90"), 0)}).</p>''')
    P.h('<div class="grid2">')
    P.chart("cn45", dict(type="hist", col=col45, xlabel="45S copies per diploid genome", xfmt=0), "45S rDNA copy number" + (" (cohort-calibrated)" if col45 == "rDNA45S.cn" else ""))
    P.chart("cn5", dict(type="hist", col=col5, xlabel="5S copies per diploid genome", xfmt=0), "5S rDNA copy number")
    P.h("</div>")
    if superpop_order:
        P.chart("pop", dict(type="strip", col=col45, by="superpop", order=superpop_order, labels=SUPERPOP_NAMES, ylabel="45S copies"), "45S copy number by super-population",
                "Each dot one genome; the bar is the median. Whether population differences survive adjustment for technical structure is a question for the complete cohort.")
        if len(rd["by_pop"]) > 1:
            P.h("<details><summary>By population</summary>")
            P.table([[g["group"], g["n"], fmt(g["median"], 0), fmt(g.get("q10"), 0) + "–" + fmt(g.get("q90"), 0)] for g in rd["by_pop"]], ["population", "n", "median 45S", "10–90%"], numeric={1, 2, 3})
            P.h("</details>")
    cx = bio["cn45_vs_5S"]
    if cx.get("n", 0) >= 3:
        P.chart("c45v5", dict(type="scatter", x=col45, y=col5, xlabel="45S copies", ylabel="5S copies", fit=True),
                "45S against 5S", f"n = {cx['n']:,}: Pearson r = {ci(cx)}, Spearman {fmt(cx.get('spearman'), 2)}. Gibbons et al. (2015) reported the two arrays' copy numbers to be correlated; Hall et al. (2021) did not see it in these genomes. Here the two classes are measured the same way, with no shared denominator.")
    P.h("<h3>The cell line</h3>")
    P.h(f'''<p>Every genome was sequenced from a lymphoblastoid cell line, whose state leaves marks on coverage. Mitochondrial genomes per cell: median
<strong>{fmt(bio["chrM"].get("median"), 0)}</strong> (10–90%: {fmt(bio["chrM"].get("q10"), 0)}–{fmt(bio["chrM"].get("q90"), 0)}); EBV episomes per cell: median
<strong>{fmt(bio["chrEBV"].get("median"), 1)}</strong> ({fmt(bio["chrEBV"].get("q10"), 1)}–{fmt(bio["chrEBV"].get("q90"), 1)}). Both vary far more between people
than the rDNA does and neither is inherited through the nuclear genome, which is why they serve as negative controls in 3.6.</p>''')
    P.h('<div class="grid2">')
    P.chart("chrM", dict(type="hist", col="chrM.copies", xlabel="mitochondrial genomes per cell", xfmt=0), "Mitochondrial DNA content")
    P.chart("ebv", dict(type="hist", col="chrEBV.copies", xlabel="EBV episomes per cell", xfmt=0), "EBV load")
    P.h("</div>")
    cm, sph = bio["log_cn45_vs_log_chrM"], bio["DJ_vs_chrX_female"]
    P.h('<div class="grid2">')
    if cm.get("n", 0) >= 3:
        P.chart("c45vM", dict(type="scatter", x="chrM.copies", y=col45, xlabel="mitochondrial genomes per cell", ylabel="45S copies", log="xy", fit=True),
                "45S copy number against mitochondrial content", f"Log scales, n = {cm['n']:,}: r = {ci(cm)}. Gibbons et al. (2014) reported rDNA copy number to be coupled with mitochondrial DNA abundance in lymphoblastoid lines.")
    if sph.get("n", 0) >= 3:
        P.chart("sphase", dict(type="scatter", x="truth.chrX", y=kt["DJ_col"], where=dict(sex_inferred="F"), xlabel="chrX copies (women)", ylabel="distal junction copies", fit=False, xref=2, yref=10),
                "Two late-replicating controls, in women", f"Women whose culture has kept both X chromosomes (chrX 1.85–2.15), n = {sph['n']:,}: r = {ci(sph)}. The inactive X and the acrocentric short arms replicate late; a culture with more cells in S phase should under-represent both together. An r near zero says the two deficits are not one thing.")
    P.h("</div>")
    du = bio["dup"]
    P.chart("dup", dict(type="scatter", x="ctrl_dup_frac", y="rDNA45S.dup_flag_frac", xlabel="duplicate-flagged, control reads", ylabel="duplicate-flagged, 45S reads", identity=True, xfmt=3),
            "The duplicate flag inside and outside the rDNA", f"Median {fmt(du['control'].get('median'), 3, pct=True)} of control reads carry the duplicate flag, {fmt(du['rDNA'].get('median'), 3, pct=True)} of 45S reads (ratio {fmt(du['ratio'].get('median'), 2)}, range {fmt(du['ratio'].get('min'), 2)}–{fmt(du['ratio'].get('max'), 2)}): the collapsed rDNA hides duplicates from the marker. A pipeline that drops flagged reads under-reads the rDNA by a different amount in every genome; NGS-DOSE counts all primary reads.")
    if sat["classes"]:
        P.h("<h3>Satellite arrays and the telomeric repeat</h3>")
        P.h('''<p>Satellite arrays are dispersed over the alignment, so only a whole-file scan measures them. Diploid mass per family; what each
panel can see was measured on the genome it was built from (four families are relative measures, under-read by their k-mer recall).</p>''')
        rows_s = [[cls, d["n"], fmt(d["median"], 1), fmt(d.get("q10"), 1) + "–" + fmt(d.get("q90"), 1), fmt(d["by_sex"]["M"].get("median"), 1), fmt(d["by_sex"]["F"].get("median"), 1)] for cls, d in sat["classes"].items()]
        P.table(rows_s, ["class", "n", "median Mb (diploid)", "10–90%", "men", "women"], numeric={1, 2, 3, 4, 5})
        P.h('<div class="grid2">')
        P.chart("hsat3", dict(type="hist", col="HSat3.mass_Mb", xlabel="Mb per diploid genome", xfmt=0), "HSat3")
        P.chart("hsat1b", dict(type="hist", col="HSat1B.mass_Mb", group=dict(col="sex_inferred", levels=sex_levels), xlabel="Mb per diploid genome", xfmt=1), "HSat1B, by sex (the family lives mostly on Yq)")
        P.chart("ahor", dict(type="hist", col="aSatHOR.mass_Mb", xlabel="Mb per diploid genome", xfmt=0), "α-satellite higher-order repeats")
        P.chart("tel", dict(type="hist", col="TEL.mass_Mb", xlabel="Mb of (TTAGGG)n-bearing reads, diploid", xfmt=2), "Telomeric repeat (a relative measure)")
        P.h("</div>")
    P.end()

    # ---------------------------------------------------------------- 5. limitations
    P.section("limitations", "5. Limitations", "Limitations")
    P.h(f'''<ul>
<li><strong>No absolute calibration.</strong> No orthogonal assay of rDNA copy number exists for these samples. The absolute level rests on unit
windows on which three Illumina chemistries agree; the known truths test the model and the k-mer path, not the absolute scale of the rDNA.</li>
<li><strong>Cell-line DNA.</strong> Every sample is a lymphoblastoid line; its replication state, EBV load and mitochondrial content are measured
but not removed. Blood-derived genomes will not carry the first of these.</li>
<li><strong>Trios bound reliability from above</strong> where members of a family were prepared together; the spousal correlation is the check.
Because the rDNA varies far more between people than any estimator errs, trios show that the measured variation is real, not which estimator
measures it best.</li>
<li><strong>One chemistry, one pipeline.</strong> The cohort is NovaSeq 2×150 aligned by one pipeline; the cross-technology evidence is twelve
genomes. DRAGEN alignments and other chemistries are untested.</li>
<li><strong>The satellite panels</strong> were built from one genome (CHM13). Four of the ten families are relative measures; HSat2 is unjudged
until assemblies without gaps in it have been compared. The telomere class is a relative measure of (TTAGGG)n content, not a telomere length.</li>
<li><strong>Partial cohort.</strong> {n:,} of {total:,}: population comparisons and the number of complete trios depend on which genomes have
landed.</li>
</ul>''')
    P.end()

    # ---------------------------------------------------------------- 6. data
    P.section("reproduce", "6. Data and reproducibility", "Data")
    P.h(f'''<p>The counts files under <code>counts_scan/</code> and <code>counts_fetch/</code> are the primary data: about 240 kB per whole-file scan and
70 kB per fetch, no reads, no genotypes. Everything above is computed from them:</p>
<pre>pip install ngsdose        # or: apptainer pull ngs-dose.sif docker://ghcr.io/jlanej/ngs-dose:latest
ngsdose report --scan counts_scan/ --fetch counts_fetch/ -p pedigree.txt --hall hall2021_MOESM1.txt --pilot pilot/ --qc ngspca_sample_qc.tsv -o docs/</pre>
<p>Tables behind every figure: <code>data/cohort.tsv</code> (one row per genome, every column), <code>data/modes.tsv</code>,
<code>data/transmission.tsv</code>, <code>data/pcsweep.tsv</code>, <code>data/satellites_hprc.tsv</code>, <code>data/flags.tsv</code>; the numbers
in the prose, <code>report.json</code>. The counts were made by <code>ngs-dose count</code> ({eng}) from the 1000 Genomes 30× CRAMs
(Byrska-Bishop et al., <em>Cell</em> 2022; AWS Open Data) with the {esc(m.get("bundle"))} resource bundle. Method, design document and
audit: <a href="https://github.com/jlanej/NGS-DOSE">github.com/jlanej/NGS-DOSE</a>.</p>''')
    P.end()

    # ---------------------------------------------------------------- 7. every sample
    P.section("samples", "7. Every genome", "Samples")
    cols = [("sample", "sample"), ("pop", "pop"), ("sex", "sex (ped)"), ("sex_inferred", "sex (reads)"), ("depth", "depth"), ("truth.auto", "auto"), ("truth.chrX", "chrX"),
            ("truth.chrY", "chrY"), (kt["DJ_col"], "DJ"), (col45, "45S"), ("rDNA45S.cn_single", "45S single"), ("rDNA45S.18S.flat", "18S flat"),
            (col5, "5S"), ("chrM.copies", "chrM"), ("chrEBV.copies", "EBV"), ("fetch_ratio.rDNA45S", "fetch/scan 45S"),
            ("HSat3.mass_Mb", "HSat3 Mb"), ("aSatHOR.mass_Mb", "αSat Mb"), ("TEL.mass_Mb", "TEL Mb"), ("flags", "flags")]
    cols = [(c, h) for c, h in cols if any(r.get(c) not in (None, "") for r in rows)]
    nd = {"depth": 1, "truth.auto": 3, "truth.chrX": 3, "truth.chrY": 3, kt["DJ_col"]: 2, col45: 0, "rDNA45S.cn_single": 0, "rDNA45S.18S.flat": 0,
          "rDNA5S.cn": 0, "rDNA5S.cn_single": 0, "chrM.copies": 0, "chrEBV.copies": 1, "fetch_ratio.rDNA45S": 4, "HSat3.mass_Mb": 1, "aSatHOR.mass_Mb": 1, "TEL.mass_Mb": 3}
    body = []
    for r in sorted(rows, key=lambda r: r["sample"]):
        body.append([("–" if r.get(c) in (None, "") or (isinstance(r.get(c), float) and not math.isfinite(r[c])) else (fmt(r[c], nd[c]) if c in nd else r[c])) for c, _ in cols])
    P.h(f"<p>{len(rows):,} genomes; click a heading to sort, type to filter. Flagged rows are shaded. The full table with every column is <code>data/cohort.tsv</code>.</p>")
    P.table(body, [h for _, h in cols], numeric={i for i, (c, _) in enumerate(cols) if c in nd}, flagged=lambda r: bool(r[-1] and r[-1] != "–") if cols[-1][0] == "flags" else None, filter_box=True, wrap=True)
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
<meta name="description" content="rDNA copy number from short-read genomes: the method and its validation on the 1000 Genomes cohort, recomputed as the run proceeds.">
<style>{css}</style></head>
<body><main>{body[:body.index("<section")]}{nav}{body[body.index("<section"):]}
<footer>Generated {esc(m["as_of"])} by {esc(m["generator"])}. NGS-DOSE was developed by Claude (Anthropic) with @jlanej; the data are 1000 Genomes open-access.</footer>
</main>
<script id="report-data" type="application/json">{payload}</script>
<script>{js}</script>
</body></html>'''
