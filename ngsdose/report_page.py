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
import re

from .hprc import MAX_GAPPED
from .report import ASSETS, fmt
from .tables import num as num_

esc = lambda x: html.escape(str(x))
SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
SUPERPOP_NAMES = {"AFR": "African", "AMR": "American", "EAS": "East Asian", "EUR": "European", "SAS": "South Asian"}
# the satellite families of the experimental panel: what each is, and the panel's recall on CHM13, the genome it was built from
# (resources/experimental/README.md: the share of 150-bp reads from the family's own CHM13 arrays that the panel can assign)
SAT_FAMILY = {"HSat1A": ("human satellite 1A", 0.998), "HSat1B": ("human satellite 1B", 0.992), "HSat2": ("human satellite 2", 0.990),
              "HSat3": ("human satellite 3", 0.973), "aSatHOR": ("α-satellite higher-order repeats", 0.999), "bSat": ("β-satellite", 0.69),
              "ACRO": ("ACRO1 composites of the acrocentric short arms", 0.90), "SST1": ("SST1 arrays", 0.61), "CER": ("centromeric repeat", 0.42),
              "SATR": ("SATR1/2 satellite", 0.50)}
SAT_NOTE = {"HSat1B": " Most HSat1B lies on the long arm of chrY, so men form the upper cluster.",
            "SST1": " The HPRC annotation labels several times more sequence as SST1 than the CHM13 annotation the panel was built from, so the ratio is not a recall.",
            "SATR": " The HPRC annotation labels several times more sequence as SATR than the CHM13 annotation the panel was built from, so the ratio is not a recall."}


def pm(d: dict, nd=3, key="mean") -> str:
    """'mean ± SD' from a describe() dict."""
    if not d or not d.get("n"):
        return "–"
    s = fmt(d[key], nd)
    if d.get("sd") is not None:
        s += f" ± {fmt(d['sd'], nd)}"
    return s


def doi(text: str, d: str) -> str:
    """A citation that links to its DOI."""
    return f'<a href="https://doi.org/{d}">{text}</a>'


def sentence_case(s: str) -> str:
    """A label as the start of a title: an all-lowercase first word is capitalised; a name such as bSat is left alone."""
    first = s.split(" ", 1)[0]
    return s[:1].upper() + s[1:] if first.islower() else s


def unit_of(col: str) -> str:
    """The unit a trio metric is plotted in."""
    if col.endswith(".mass_Mb"):
        return "Mb"
    return {"chrM.copies": "genomes per cell", "chrEBV.copies": "episomes per cell", "depth": "×", "ctrl_dup_frac": "fraction flagged",
            "gc_rel_65": "rate at 65% GC, relative", "insert_median": "bp"}.get(col, "copies")


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
        self.n_fig = self.n_supp = 0

    def h(self, s: str):
        self.parts.append(s)

    def section(self, id_: str, title: str, short: str | None = None):
        self.toc.append((id_, short or title))
        self.h(f'<section id="{id_}"><h2>{esc(title)}</h2>')

    def end(self):
        self.h("</section>")

    def chart(self, id_: str, spec: dict, title: str, caption: str, supp: bool = False):
        """A numbered figure that stands alone, so that a screenshot of it can be shared: the title says what is plotted, the
        caption which genomes, what every colour and line is and what the figure shows, and the last line where it comes from.
        Supplementary figures are numbered apart (S1, S2, ...)."""
        self.charts[id_] = spec
        if supp:
            self.n_supp += 1
            number = f"S{self.n_supp}"
        else:
            self.n_fig += 1
            number = str(self.n_fig)
        m = self.d["meta"]
        src = f'{esc(m["title"])} · {esc(m["as_of"])} · github.com/jlanej/NGS-DOSE'
        self.h(f'<figure id="fig-{id_}"><div class="title">Figure {number}. {esc(title)}</div><div class="chart" id="chart-{id_}"></div>'
               f'<figcaption>{caption}<span class="src">{src}</span></figcaption></figure>')

    def tiles(self, items: list[tuple[str, str, str]]):
        self.h('<div class="tiles">' + "".join(f'<div class="tile"><div class="label">{esc(a)}</div><div class="value">{b}</div><div class="note">{esc(c)}</div></div>'
                                             for a, b, c in items) + "</div>")

    def table(self, rows: list[list], header: list[str], numeric: set[int] = frozenset(), flagged=None, filter_box=False, wrap=False, cls=""):
        head = "".join(f'<th class="{"num" if i in numeric else ""}">{esc(h)}</th>' for i, h in enumerate(header))
        body = []
        for r in rows:
            row_cls = ' class="flagged"' if flagged and flagged(r) else ""
            body.append(f"<tr{row_cls}>" + "".join(f'<td class="{"num" if i in numeric else ""}">{c if isinstance(c, str) and c.startswith("<") else esc(c)}</td>' for i, c in enumerate(r)) + "</tr>")
        t = f'<table class="data{" " + cls if cls else ""}"><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table>'
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
         'collapse, from aligned short-read genomes, which exist for hundreds of thousands of people. Standard analysis of those files does not report it, '
         'long-read assemblies leave the arrays in pieces, and read-depth ratios measure the library along with the person (<a href="#why">section 1</a>). '
         'Reads are assigned to a sequence class by class-specific 31-mers and counted as fragment ends; a '
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
    fc = data.get("fetch_check") or {}
    fR = {t["column"]: t for t in (fc.get("trios") or {}).get("table", [])}
    if fc and t45 and "rDNA45S.cn" in fR and "R_lo" in fR["rDNA45S.cn"]:
        r_.append(f'Run on the fetch counts alone, the cohort layer and the trio test give the same answer: 45S reliability {fmt(min(fR["rDNA45S.cn"]["R"], 1.0), 2)} ({fmt(fR["rDNA45S.cn"]["R_lo"], 2)}–{fmt(fR["rDNA45S.cn"]["R_hi"], 2)}) against {fmt(min(t45["R"], 1.0), 2)} from the scans.')
    r_.append("</p><p>")
    if t45 and have_ci:
        neg = {t["column"]: t for t in tr["table"]}
        negs = ", ".join(f"{lab} {fmt(neg[c]['R'], 2)}" for c, lab in (("truth.auto", "held-out autosomal sequence"), ("chrM.copies", "mitochondrial content")) if c in neg)
        r_.append(f'In {t45["n_trios"]} trios the calibrated 45S estimate has a transmission reliability of {fmt(min(t45["R"], 1.0), 2)} ({fmt(t45["R_lo"], 2)}–{fmt(t45["R_hi"], 2)})'
                  + (f'; traits that are not transmitted read near zero ({negs})' if negs else "") + ".")
    if rc and rf:
        r_.append(f'In {rep["n"]} individuals sequenced on two technologies, the test–retest intraclass correlation is {fmt(rc["icc"], 2)} for NGS-DOSE\'s calibrated estimate '
                  f'and {fmt(rf["icc"], 2)} for the 18S depth ratio of published studies, computed from the same reads for comparison ({fmt(rfc["icc"], 2)} after removing its offset between technologies).')
    if gcb.get("flat_vs_gc", {}).get("n", 0) >= 10:
        r_.append(f'Within one chemistry the published depth ratio follows each library\'s GC bias (r = {fmt(gcb["flat_vs_gc"]["r"], 2)}); under NGS-DOSE\'s fragment-GC model it does not (r = {fmt(gcb["modelled_vs_gc"]["r"], 2)}).')
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
        case.append(("Two technologies", f'ICC {fmt(rc["icc"], 2)} vs {fmt(rf["icc"], 2)}', f"NGS-DOSE vs the published 18S depth ratio, {rep['n']} people sequenced twice", "#replicates"))
    if t45 and "R_lo" in t45:
        case.append(("Inherited", fmt(min(t45["R"], 1.0), 2), f"45S transmission reliability ({fmt(t45['R_lo'], 2)}–{fmt(t45['R_hi'], 2)}), {t45['n_trios']} trios", "#trios"))
    if gcb.get("flat_vs_gc", {}).get("n"):
        case.append(("What the GC model removes", f'r {fmt(gcb["flat_vs_gc"].get("r"), 2)} → {fmt(gcb["modelled_vs_gc"].get("r"), 2)}', "how much the 18S estimate follows the library's GC bias: published depth ratio → NGS-DOSE's GC model", "#gcmodel"))
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
    P.h(f"""<p>Each human genome carries several hundred copies of the 45S ribosomal DNA unit, about 45 kb long, in tandem arrays on the short
arms of the five acrocentric chromosomes, and a tandem array of the 5S unit on chromosome 1. A single array spans from 50 kb to more than
6 Mb, the two homologues of a chromosome almost never carry arrays of one size, and an array is rearranged in more than one meiosis in ten
({doi("Stults et al., <em>Genome Res</em> 2008", "10.1101/gr.6858507")}). Total copy number varies severalfold between individuals and is heritable.</p>""")
    # why a method is needed: what each technology gives now, with this cohort's numbers wherever it has them
    ukb = 490_640                                               # genomes sequenced in UK Biobank (Nature 2025)
    p01 = math.erfc(0.01 * math.sqrt((ukb - 2) / (1 - 0.01 ** 2)) / math.sqrt(2))      # a correlation of 0.01 at that size
    e10 = math.floor(math.log10(p01))
    du = bio["dup"]
    faults = ["the rDNA is GC-rich", "parts of the unit drop out by amounts that depend on the chemistry"]
    if du["ratio"].get("n"):
        faults.append(f'the duplicate flag marks {fmt(du["rDNA"].get("median"), 1, pct=True)} of 45S reads against {fmt(du["control"].get("median"), 1, pct=True)} '
                      f'of single-copy reads, so dropping flagged reads moves the ratio by a different amount in every genome (the two rates\' ratio runs from '
                      f'{fmt(du["ratio"].get("min"), 2)} to {fmt(du["ratio"].get("max"), 2)})')
    here = []
    if gcb.get("flat_vs_gc", {}).get("n", 0) >= 10:
        here.append(f'within one chemistry the 18S ratio\'s departure from the calibrated estimate follows each library\'s GC bias (r = {fmt(gcb["flat_vs_gc"]["r"], 2)}; <a href="#gcmodel">3.5</a>)')
    if rc and rf:
        here.append(f'the same people read {fmt(abs(rf["offset"]), 0, pct=True)} {"lower" if rf["offset"] < 0 else "higher"} by it on the older of two technologies '
                    f'(intraclass correlation {fmt(rf["icc"], 2)}, against {fmt(rc["icc"], 2)} for the calibrated estimate; <a href="#replicates">3.4</a>)')
    ra = rd.get("assemblies") or {}
    asm = ""
    if ra.get("n"):
        fr = ra["fraction"]
        asm = (f' Of the genomes counted here, {ra["n"]} ha{"ve" if ra["n"] > 1 else "s"} an HPRC release-2 assembly; '
               + ("these hold " + f'{fmt(fr["median"], 0, pct=True)} ({100 * fr["min"]:.0f}–{fmt(fr["max"], 0, pct=True)})' if ra["n"] > 1 else "it holds " + fmt(fr["median"], 0, pct=True))
               + f' of the rDNA the measured copy number implies, in stretches no longer than {fmt(ra["longest_stretch_Mb"]["max"], 2)} Mb, where the average'
               f' array would be {fmt(ra["mean_array_Mb"]["median"], 1)} Mb (<code>data/rdna_hprc.tsv</code>).')
    # where NGS-DOSE and the published ratio can be told apart, and whether NGS-DOSE does better there: computed, not asserted
    where_better = []
    if rc and rf and rc["icc"] > rf["icc"] and abs(rc["offset"]) < abs(rf["offset"]) and rc["sd_log_ratio"] < rf["sd_log_ratio"]:
        where_better.append('across technologies (<a href="#replicates">3.4</a>)')
    if gcb.get("flat_vs_gc", {}).get("n", 0) >= 10 and abs(gcb["modelled_vs_gc"]["r"]) < abs(gcb["flat_vs_gc"]["r"]):
        where_better.append('across libraries\' GC bias (<a href="#gcmodel">3.5</a>)')
    better = (", and wherever the two can be told apart, " + " and ".join(where_better) + ", NGS-DOSE does better") if where_better else ""
    rak = "10.1016/j.xgen.2024.100562"
    P.h(f"""<h3 id="why">Why a method is needed</h3>
<p class="callout">Short-read genomes are the only genomes that exist at population scale: {ukb:,} in UK Biobank alone
({doi("UK Biobank Whole-Genome Sequencing Consortium, <em>Nature</em> 2025", "10.1038/s41586-025-09272-9")}), where rDNA copy number
estimated from them has been associated with blood-cell counts and kidney function ({doi("Rodriguez-Algarra, Evans &amp; Rakyan, <em>Cell Genomics</em> 2024", rak)})
and with metabolic disease ({doi("Raj et al., medRxiv 2026", "10.64898/2026.01.09.26343685")}, a preprint). At that size a correlation of 0.01 between an
estimate and a trait has p ≈ {p01 / 10 ** e10:.0f} × 10<sup>{str(e10).replace("-", "−")}</sup>, so whatever an estimate carries of the library, rather than
the person, becomes a finding wherever the library also tracks the trait. Standard analysis does not give the number, and read-depth ratios,
long reads and laboratory assays each fall short on scale or on the library (below). What is needed is a measurement that finds every rDNA
read wherever the aligner put it (<a href="#modes">3.3</a>), predicts each library's count of GC-rich sequence from single-copy sequence in the same genome (<a href="#gcmodel">3.5</a>), holds its scale
across chemistries (<a href="#replicates">3.4</a>), does not depend on the duplicate flag (<a href="#published">3.7</a>), reads a small part of
each file (<a href="#modes">3.3</a>), and shows in every genome that it reads known copy numbers correctly (<a href="#truth">3.1</a>,
<a href="#djsteps">3.2</a>). NGS-DOSE is one way to meet these, and the tests on this page apply to any other. The 18S depth ratio of
published studies is not NGS-DOSE's estimate: it is computed from the same reads and carried beside it throughout as the comparator{better}.</p>
<dl class="methods">
<dt>Standard analysis</dt><dd>Reports no copy number. GRCh38 carries fewer than ten copies of the 45S transcribed region, on chr21 and two
small contigs, where a genome has hundreds; the reads of every unit fall on those copies with no unique alignment, and variant and
copy-number callers mask the loci.</dd>
<dt>Read-depth ratio</dt><dd>Depth on the reference's rDNA relative to depth elsewhere in the genome, as in published studies
({doi("Gibbons et al., <em>Nat Commun</em> 2014", "10.1038/ncomms5850")}; {doi("Hall et al., <em>Sci Rep</em> 2021", "10.1038/s41598-020-80049-y")};
{doi("Rodriguez-Algarra et al. 2024", rak)}), runs on existing files. Without a model of the library it measures the library along with the
person: {", ".join(faults[:-1])}{"," if len(faults) > 2 else ""} and {faults[-1]}. In UK Biobank its mean differs between the two sequencing centres
({doi("Rodriguez-Algarra et al. 2024", rak)}).{(" Here, " + "; ".join(here) + ".") if here else ""}</dd>
<dt>Long reads</dt><dd>The units are near-identical, which defeats assembly. Even in T2T-CHM13, three of the five arrays are model sequences,
because ultra-long nanopore reads were not long enough to order their units ({doi("Nurk et al., <em>Science</em> 2022", "10.1126/science.abj6987")});
in 156 acrocentric short arms assembled from HiFi, ultra-long nanopore and Hi-C data, the rDNA array collapsed in every one
({doi("Lin et al., <em>Cell</em> 2026", "10.1016/j.cell.2026.05.035")}). Long reads, imaging and methylation have given the size and activity of
each array in fifteen genomes ({doi("Potapova et al., <em>Cell Genomics</em> 2025", "10.1016/j.xgen.2025.101031")}).{asm}</dd>
<dt>Molecular assays</dt><dd>Droplet digital PCR is precise (CHM13: 409 ± 9 copies; Nurk et al. 2022) and pulsed-field gels size single
arrays (Stults et al. 2008), but each needs the DNA and a laboratory assay per person: in a biobank, every participant's stored DNA.</dd>
</dl>
<p class="small">The case is strongest for the rDNA. Long-read assemblies close most satellite arrays and are their truth here (<a href="#assemblies">3.8</a>),
so for the satellites it is one of scale alone; telomeric content already has short-read estimators
({doi("Ding et al., <em>Nucleic Acids Res</em> 2014", "10.1093/nar/gku181")}), and the class measured here is a relative one.</p>""")
    P.h('''<h3>How the measurement is judged</h3>
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
2 × observed / expected (<a href="#fig-m_gc">figure</a>).</dd>
<dt>Calibration</dt><dd>Windows of the 45S unit drop out beyond what the GC curve predicts, by amounts that depend on the sequencing
chemistry. Per-window efficiencies are learned across the cohort by median polish; the absolute scale is set by anchor windows on which
three Illumina chemistries agreed in the pilot. Beside it the page carries NGS-DOSE's single-sample variant, from the anchor windows
alone, and, as the comparator, the 18S read-depth ratio of published studies, with no GC model and no calibration: not NGS-DOSE's estimate
(<a href="#fig-m_unit">figure</a>).</dd>
<dt>Estimators</dt><dd>Three 45S estimates travel through every table, two of them NGS-DOSE's. <em>45S, NGS-DOSE calibrated</em> (<code>rDNA45S.cn</code>): the cohort
model log C<sub>iw</sub> = c<sub>i</sub> + a<sub>w</sub> + e<sub>iw</sub> over every retained 250-bp window w of the unit, fitted by
median polish across samples i, with the window efficiencies a<sub>w</sub> pinned to a median of zero over the anchor windows; the estimate
is exp(c<sub>i</sub>), so every window contributes precision and the anchors set the level. <em>45S, NGS-DOSE single-sample</em>
(<code>rDNA45S.cn_single</code>): 2 × observed / expected fragment ends summed over the anchor windows alone, under the sample's own
fragment-GC model, with no information from any other sample. <em>45S, 18S depth ratio (published)</em> (<code>rDNA45S.18S.flat</code>): 2 × fragment
ends in the 18S gene / (positions × the control regions' mean rate), with no GC model and no calibration: the read-depth ratio of published
studies, computed from the same reads as the comparator (<a href="#fig-m_est">figure</a>). The 5S and distal-junction estimates are calibrated the same way as the 45S; the satellite masses are diploid megabases from the
class's read count under the GC model.</dd>
<dt>Known-truth controls</dt><dd>80 held-out autosomal regions (two copies), 60 chrX regions (one in men, two in women) and 40 X-degenerate
chrY regions (one, none), measured by the alignment-position path; and the distal junction, present once on each of the ten acrocentric
short arms, measured by the k-mer path with its class restricted to k-mers that occur exactly once in each of the five CHM13 junctions.
Mitochondrial genomes and EBV episomes per cell are measured as covariates of the culture.</dd>
<dt>Transmission</dt><dd>For each metric, values on the natural scale with the population mean subtracted; then, over the complete trios,
the Pearson correlations of child with father, mother and midparent, the spousal (father–mother) correlation ρ, and the least-squares slope b
of child on midparent. Reliability R = b − ρ(1 − b) estimates the share of the metric's variance that is transmitted (σ²<sub>T</sub> /
(σ²<sub>T</sub> + σ²<sub>e</sub>)): 1 for a perfectly measured heritable trait, 0 for pure error. 95% intervals by family bootstrap at 20
trios or more; a one-sided p-value for the slope from 1,000 permutations of children among families; paired bootstrap for differences
between estimators. Held-out autosomal sequence (no true variance) and mitochondrial and EBV content (not in the nuclear genome) are the
negative controls.</dd>
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
    # ---- the method, step by step, on this cohort's data
    mp = data.get("methods") or {}
    if mp:
        g, wg, u = mp["gc"], mp["windows_gc"], mp["unit"]
        i65 = g["x"].index(65) if 65 in g["x"] else None
        P.h("<h3>The method, step by step, on this cohort's data</h3>")
        P.chart("m_gc", dict(type="lines", series=[dict(name="median library, with the middle 80%", x=g["x"], y=g["median"], lo=g["q10"], hi=g["q90"], ci=0)],
                             xlabel="fragment GC, %", ylabel="rate relative to the library's mean", ref=1,
                             bands=[dict(x0=100 * wg["q10"], x1=100 * wg["q90"], label="45S unit")]),
                "Library model: how each library's sequencing rate depends on fragment GC",
                f"The per-library Poisson spline of fragment-end density on fragment GC content, fitted on the 800 single-copy control regions, in each "
                f"of {mp['n_gc']:,} genomes of the 1000 Genomes 30× cohort: the rate at each fragment GC relative to the library's mean. Line: the "
                f"cohort's median library; shaded: the middle 80% of libraries; grey band: the fragment GC of the middle 80% of the 45S unit's windows "
                f"({fmt(100 * wg['q10'], 0)}–{fmt(100 * wg['q90'], 0)}%)."
                + (f" These libraries sequence 65%-GC fragments at {fmt(g['median'][i65], 2)}× their mean rate (middle 80% {fmt(g['q10'][i65], 2)}–"
                   f"{fmt(g['q90'][i65], 2)}×), so a count of GC-rich sequence taken without this model reads high, by a different amount in each "
                   "library." if i65 is not None else ""))
        rr = [v for v in u["raw"]["median"] if v == v]
        feats = [f_ for f_ in mp.get("features", []) if f_["name"] in ("18S", "5.8S", "28S")]
        P.chart("m_unit", dict(type="lines", series=[dict(name="GC model alone", x=u["mid_kb"], y=u["raw"]["median"], lo=u["raw"]["q10"], hi=u["raw"]["q90"], ci=2),
                                                     dict(name="after the window efficiencies", x=u["mid_kb"], y=u["cal"]["median"], lo=u["cal"]["q10"], hi=u["cal"]["q90"], ci=0)],
                               xlabel="position in the 45S unit, kb", ylabel="window copy number / the genome's estimate", ref=1,
                               bands=[dict(x0=f_["start"] / 1000, x1=f_["end"] / 1000, label=f_["name"]) for f_ in feats],
                               ticks=dict(x=u["anchor_kb"], label="anchor windows")),
                "Calibration: the windows of the 45S unit before and after the window efficiencies",
                f"Each of the {u['n_retained']} retained 250-bp windows of the 45S unit (KY962518.1; {u['n_windows'] - u['n_retained']} masked windows left "
                f"out), in each of {mp['n_unit']:,} genomes: the window's copy number under the library's GC model divided by the genome's calibrated "
                f"estimate. Lines: the cohort's median; shaded: the middle 80%. Aqua, the GC model alone: windows still read from "
                f"{fmt(min(rr), 2)}× to {fmt(max(rr), 2)}× the genome's level, by amounts that depend on the sequencing chemistry, which the GC curve "
                f"does not predict. Blue, after dividing by the window efficiencies learned across the cohort by median polish: the median is 1 at "
                f"every window by construction, and the band is what remains, each genome's scatter about its own level. Grey bands: the {', '.join(f_['name'] for f_ in feats)} genes; blue ticks: the {len(u['anchor_kb'])} anchor windows, "
                f"on which three Illumina chemistries agreed in the pilot, which set the absolute scale.")
        pts_e = [dict(x=num_(r, "rDNA45S.cn"), y=num_(r, k), label=r["sample"], si=si) for si, k in enumerate(("rDNA45S.cn_single", "rDNA45S.18S.flat"))
                 for r in rows if num_(r, "rDNA45S.cn") > 0 and num_(r, k) > 0]
        est = {}
        for si, k in enumerate(("rDNA45S.cn_single", "rDNA45S.18S.flat")):
            xy = [(p_["x"], p_["y"]) for p_ in pts_e if p_["si"] == si]
            if len(xy) >= 3:
                xv, yv = [x for x, _ in xy], [y for _, y in xy]
                mx, my = sum(xv) / len(xv), sum(yv) / len(yv)
                sxy = sum((x - mx) * (y - my) for x, y in xy)
                est[k] = dict(n=len(xy), ratio=sorted(y / x for x, y in xy)[len(xy) // 2],
                              r=sxy / math.sqrt(sum((x - mx) ** 2 for x in xv) * sum((y - my) ** 2 for y in yv)))
        if len(est) == 2:
            e1, e2 = est["rDNA45S.cn_single"], est["rDNA45S.18S.flat"]
            P.chart("m_est", dict(type="scatter", points=pts_e, legend=["45S, NGS-DOSE single-sample", "45S, 18S depth ratio (published)"],
                                  xlabel="45S, NGS-DOSE calibrated, copies", ylabel="the other estimate, copies", identity=True, fit=True,
                                  fit_labels=["NGS-DOSE single-sample", "published 18S ratio"]),
                    "Estimators: the three 45S estimates of every genome",
                    f"Each of {e1['n']:,} genomes of the 1000 Genomes 30× cohort twice, against its NGS-DOSE calibrated estimate (x). Blue: NGS-DOSE's "
                    f"single-sample estimate, from the anchor windows alone under the genome's own GC model (median ratio {fmt(e1['ratio'], 2)}, "
                    f"r = {fmt(e1['r'], 3)}). Orange: the 18S read-depth ratio of published studies, with no GC model and no calibration, computed from "
                    f"the same reads (median ratio {fmt(e2['ratio'], 2)}, r = {fmt(e2['r'], 3)}). Grey diagonal: equality; lines: least-squares fits."
                    + (" Within this one chemistry the three agree on who carries more rDNA and differ in level; across technologies only NGS-DOSE's "
                       "holds its level (3.4)." if min(e1["r"], e2["r"]) > 0.95 else ""))
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
                "Genomes counted so far, by super-population",
                f"The 1000 Genomes 30× cohort (New York Genome Center; NovaSeq 2×150, PCR-free; aligned to GRCh38) by super-population: the dark bar "
                f"is the genomes counted so far, the pale track all of them. {n:,} of {total:,} genomes counted; {tr['n_complete']:,} of {tr['n_total']:,} "
                f"trios complete.")
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
    cohort_ = "genomes of the 1000 Genomes 30× cohort"
    P.chart("auto", dict(type="hist", col="truth.auto", xlabel="copies", ref=[dict(x=2, label="expected 2")], xfmt=3), "Known copy number: held-out autosomal sequence (expected 2)",
            f"80 autosomal regions present in two copies in every genome, kept out of the library model and measured by the same code as the rDNA, "
            f"in each of {a.get('n', 0):,} {cohort_}. Bars count genomes; the line marks 2. Mean ± SD {pm(a)}.")
    P.chart("chrX", dict(type="hist", col="truth.chrX", group=dict(col="sex_inferred", levels=sex_levels), xlabel="copies", ref=[dict(x=1, label="1"), dict(x=2, label="2")], xfmt=2),
            "Known copy number: chrX, by sex (expected 1 in men, 2 in women)",
            f"60 chrX regions in each of {X['M'].get('n', 0) + X['F'].get('n', 0):,} {cohort_}: men blue, women orange, by the sex the reads show; lines at 1 and 2. Men {pm(X['M'])}"
            + (f"; women whose cell line has kept both X chromosomes {pm(sx['women_intact'])}. The {sx['n_mosaic_X']} women below 1.85 copies have lost an X in part of their cell line."
               if women_ok else f"; women {pm(X['F'])}."))
    P.chart("chrY", dict(type="hist", col="truth.chrY", group=dict(col="sex_inferred", levels=sex_levels), xlabel="copies", ref=[dict(x=0, label="0"), dict(x=1, label="1")], xfmt=2),
            "Known copy number: chrY, by sex (expected 1 in men, 0 in women)",
            f"40 X-degenerate chrY regions in each of {Y['M'].get('n', 0) + Y['F'].get('n', 0):,} {cohort_}: men blue, women orange; lines at 0 and 1. "
            + (f"Men with an intact Y {pm(sx['men_intact_Y'])}; women at most {fmt(Y['F'].get('max'), 3)}. The {sx['n_mosaic_Y']} men below 0.85 copies have lost the Y in part "
               f"of their cell line." if women_ok else f"Men {pm(Y['M'])}; women {pm(Y['F'], 4)}."))
    P.chart("dj", dict(type="hist", col=kt["DJ_col"], xlabel="copies", ref=[dict(x=10, label="expected 10")], xfmt=2), "Known copy number: the distal junction (expected 10)",
            f"A 400-kb sequence present once on each of the ten acrocentric short arms, beside the rDNA arrays, measured by the same k-mer path as "
            f"the rDNA{' and calibrated the same way' if kt['DJ_col'] == 'DJ.cn' else ''}, in each of {DJ.get('n', 0):,} {cohort_}. Bars count genomes; "
            f"the line marks 10. Mean ± SD {pm(DJ)}; genomes a whole copy away carry a structural variant of a short arm (section 3.2).")
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
        tot_ = dj.get("transmitted", 0) + dj.get("not_transmitted", 0)
        P.chart("djstep", dict(type="hist", col="DJ.step", xlabel="distal-junction copies relative to the cohort's level", ref=[dict(x=k, label=str(k)) for k in (-2, -1, 0, 1)], xfmt=1, bins=40),
                "The distal junction changes in whole copies",
                f"Distal-junction copy number in each of {DJ.get('n', 0):,} genomes of the 1000 Genomes 30× cohort, relative to the cohort's level "
                f"({fmt(dj['median'], 2)} copies); lines at whole copies. A rearranged acrocentric short arm leaves one copy fewer, a Robertsonian "
                f"translocation two. Within ±0.3 of a step: {nr(-2)} genomes at −2, {nr(-1)} at −1, {nr(1)} at +1"
                + (f"; where a carrier parent and a child were both counted, the step was passed on in {dj['transmitted']} of {tot_}." if tot_ else "."))
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
        P.chart("fetch45", dict(type="hist", col="fetch_ratio.rDNA45S", xlabel="fetch / scan, 45S copies", ref=[dict(x=1, label="1")], xfmt=4),
                "One minute of the file gives the whole-file answer: 45S, targeted fetch against whole-file scan",
                f"For each of {md['n_both']:,} genomes of the 1000 Genomes 30× cohort counted both ways from the same CRAM: the 45S estimate from the "
                f"targeted fetch (the control regions and the intervals where the aligner places rDNA reads, about 0.5 GB of a 15-GB file) divided by "
                f"the estimate from reading the whole file. Bars count genomes; the line marks 1. Median {fmt(modes45.get('median'), 4)}, range "
                f"{fmt(modes45.get('min'), 4)}–{fmt(modes45.get('max'), 4)}.")
    else:
        P.h("<p>No genome has been counted in both modes yet; this section fills in when the per-sample jobs, which do both, land.</p>")
    if fc and fc.get("agreement"):
        S_ = {t["column"]: t for t in tr["table"]}
        cls_rows = [(col, d) for col, d in fc["agreement"].items() if not d.get("identical")]
        P.h(f"""<h3>Does the fetch carry the same information?</h3>
<p>The cohort layer (window calibration, control-region PCs) and the trio test were run a second time on the fetch counts of the
{fc["n"]:,} genomes alone, without reference to the scans. Columns made from the same reads in both modes (the known truths, the
dosage regions, the library's properties) come out identical and are not listed; the classes differ by what the sinks miss and by a
calibration learned twice.</p>""")
        P.table([[d["label"], d["n"], f'{fmt(d["ratio_median"], 4)} ({fmt(d["q10"], 4)}–{fmt(d["q90"], 4)})', fmt(d.get("r"), 4),
                  (fmt(min(S_[col]["R"], 1.0), 2) + (f' ({fmt(S_[col]["R_lo"], 2)}–{fmt(S_[col]["R_hi"], 2)})' if "R_lo" in S_[col] else "")) if col in S_ else "–",
                  (fmt(min(fR[col]["R"], 1.0), 2) + (f' ({fmt(fR[col]["R_lo"], 2)}–{fmt(fR[col]["R_hi"], 2)})' if "R_lo" in fR[col] else "")) if col in fR else "–"]
                 for col, d in cls_rows],
                ["metric", "genomes", "fetch / scan, median (10–90%)", "r across genomes", "reliability from scan (95% CI)", "reliability from fetch (95% CI)"], numeric={1, 2, 3, 4, 5})
        P.h('<p class="small">The full comparison, every column: <code>data/fetch_check.tsv</code>.</p>')
        # every genome's estimate both ways, one figure per fetchable class
        what = {"rDNA45S.cn": ("45S rDNA copy number", "copies", "45S rDNA"), "rDNA5S.cn": ("5S rDNA copy number", "copies", "5S rDNA"),
                "DJ.cn": ("Distal-junction copy number", "copies", "distal-junction"), "TEL.mass_Mb": ("Telomeric-repeat mass", "Mb", "telomeric")}
        fpts = {col: p for col, p in (fc.get("points") or {}).items() if col in what and col in fc["agreement"]}
        if fpts:
            P.h('<div class="grid2">')
            for col, pts in fpts.items():
                d, (name, unit, reads) = fc["agreement"][col], what[col]
                rel = (f"; reliability in {S_[col]['n_trios']} trios {fmt(min(S_[col]['R'], 1.0), 2)} from the scan and {fmt(min(fR[col]['R'], 1.0), 2)} from the fetch"
                       + (" (capped at 1)" if max(S_[col]["R"], fR[col]["R"]) > 1 else "") if col in S_ and col in fR else "")
                r_txt = "r > 0.99999" if d.get("r", 0) >= 0.999995 else f"r = {fmt(d.get('r'), 5)}"
                P.chart(f"fs_{col.split('.')[0]}", dict(type="scatter", points=[dict(x=x, y=y, label=s) for s, x, y in pts], xlabel=f"whole-file scan, {unit}",
                                                        ylabel=f"targeted fetch, {unit}", identity=True),
                        f"{name}, NGS-DOSE: targeted fetch against whole-file scan, per genome",
                        f"Each dot is one of {d['n']:,} genomes of the 1000 Genomes 30× cohort, counted twice from the same CRAM: by reading the whole "
                        f"file (x) and by the targeted fetch (y), which reads only the control regions and the intervals where the aligner places "
                        f"{reads} reads, about 0.5 GB of a 15-GB file. The whole NGS-DOSE pipeline was run separately on each set of counts. Diagonal: "
                        f"agreement. {r_txt}; fetch / scan median {fmt(d['ratio_median'], 4)} (10–90% {fmt(d['q10'], 4)}–{fmt(d['q90'], 4)}){rel}.")
            P.h("</div>")
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
depth, alignment pipeline and a GC response the reverse of NovaSeq's. Agreement between the two libraries of a person, per estimator: NGS-DOSE's
(the calibrated estimate with anchor windows chosen with the person's family held out, so the cross-technology level is out of sample) and,
for comparison, the 18S depth ratio of published studies computed from the same reads:</p>''')
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
        spct = lambda v: ("+" if v >= 0 else "−") + fmt(abs(v), 0, pct=True)
        # which of the agreement measures favour NGS-DOSE over the published ratio: the claim is made only as far as they do
        wins = {"offset": abs(rc["offset"]) < abs(rf["offset"]), "pair SD": rc["sd_log_ratio"] < rf["sd_log_ratio"], "Pearson r": rc["r"] > rf["r"],
                "intraclass correlation": rc["icc"] > rf["icc"], "within-person CV": rc["within_cv"] < rf["within_cv"]}
        P.chart("replicates", dict(type="scatter", points=[dict(x=p["x"], y=p["y"], label=p["sample"], si=p["si"]) for p in rep["points"]],
                                   legend=["NGS-DOSE (this method)", "18S depth ratio (published method, for comparison)"], xlabel="45S copies, NovaSeq 2×150 (2019)",
                                   ylabel="45S copies, HiSeq 2×100 / 2×126 (2012–15)", identity=True),
                "The same people on two sequencing technologies: NGS-DOSE reproduces, the published 18S ratio does not" if all(wins.values())
                else "The same people on two sequencing technologies: NGS-DOSE and the published 18S ratio",
                f"Each of {rep['n']} people of the 1000 Genomes cohort (four parent–child trios) was sequenced twice from the same cell line: by the New York Genome Center on "
                f"NovaSeq 2×150 in 2019 (x) and on HiSeq 2000 or 2500 (2×100 or 2×126) in 2012–15 (y). Each person appears twice. Blue: NGS-DOSE's "
                f"calibrated 45S estimate, with its anchor windows chosen with the person's family held out, so the level is not fitted to them. Orange: "
                f"the 18S read-depth ratio of published studies (no GC model, no calibration), computed from the same reads for comparison; it is not "
                f"NGS-DOSE's estimate. On the diagonal both sequencing runs give the same answer. NGS-DOSE: offset {spct(rc['offset'])}, pair SD "
                f"{fmt(rc['sd_log_ratio'], 1, pct=True)}, intraclass correlation {fmt(rc['icc'], 2)}. 18S ratio: offset {spct(rf['offset'])}, pair SD "
                f"{fmt(rf['sd_log_ratio'], 1, pct=True)}, intraclass correlation {fmt(rf['icc'], 2)} ({fmt(rfc['icc'], 2)} once its offset is removed).")
        P.h(f'''<p>The intraclass correlation measures absolute agreement, so an offset between technologies counts against it: the published depth ratio's
{spct(rf["offset"])} offset is a property of the library, not of the person, and removing it as a batch effect leaves an ICC of
{fmt(rfc["icc"], 2)} against {fmt(rc["icc"], 2)} for NGS-DOSE, whose offset is {spct(rc["offset"])}. The two DNA batches were
drawn from different cultures of each line, so these figures bound measurement error from above. This is the comparison that separates
NGS-DOSE from the published depth ratio{", and NGS-DOSE does better on every measure in the table: " + ", ".join(wins) if all(wins.values())
else "; NGS-DOSE does better on " + (", ".join(k for k, v in wins.items() if v) or "none") + " of the measures in the table"}. Within one chemistry
(3.6) the difference cannot be seen.</p>''')
    else:
        P.h("<p>Appears when the pilot's replicate tables are given (<code>--pilot</code>).</p>")
    P.end()

    # ---------------------------------------------------------------- 3.5 GC model
    P.section("gcmodel", "3.5 What the fragment-GC model removes", "GC model")
    if gcb.get("flat_vs_gc", {}).get("n", 0) >= 10:
        fv, mv, g65 = gcb["flat_vs_gc"], gcb["modelled_vs_gc"], gcb["gc65"]
        P.h(f'''<p>Libraries differ in GC bias even within one chemistry: the rate at which 65%-GC fragments were sequenced, relative to each library's
mean, runs from {fmt(g65.get("q10"), 2)} to {fmt(g65.get("q90"), 2)} across the middle 80% of genomes, and the rDNA is GC-rich. The published 18S depth
ratio divided by NGS-DOSE's calibrated estimate of the same genome follows that bias with <strong>r = {ci(fv)}</strong> (n = {fv["n"]:,}); the same 18S
region under NGS-DOSE's fragment-GC model, r = {ci(mv)}. Within the cohort the effect is a few percent (SD of the log ratio {fmt(gcb["flat_sd_log"], 3)}),
small next to the {fmt(bio.get("cn45_cv"), 0, pct=True)} by which people differ, which is why it does not show in the trios; between technologies
(3.4) it is {fmt(abs(rf["offset"]), 0, pct=True) if rf else "large"}.</p>''')
        # both estimators of one region in one figure, each genome twice: NGS-DOSE's GC model blue, the published ratio orange
        gc_pts = [dict(x=r["gc_rel_65"], y=r[c], label=r["sample"], si=si) for si, c in enumerate(("rDNA45S.18S_over_cn", "rDNA45S.18S.flat_over_cn"))
                  for r in rows if isinstance(r.get("gc_rel_65"), float) and isinstance(r.get(c), float) and math.isfinite(r["gc_rel_65"]) and math.isfinite(r[c])]
        P.chart("gc", dict(type="scatter", points=gc_pts, legend=["18S under NGS-DOSE's fragment-GC model", "18S depth ratio (published method, for comparison)"],
                           xlabel="library GC bias: rate at 65% GC relative to the library's mean", ylabel="18S estimate / NGS-DOSE 45S estimate", fit=True,
                           fit_labels=[f"NGS-DOSE GC model: r = {fmt(mv.get('r'), 2)}", f"published 18S ratio: r = {fmt(fv.get('r'), 2)}"]),
                "Library GC bias: the published 18S ratio follows it, NGS-DOSE's GC model removes it" if abs(fv["r"]) > abs(mv["r"]) else "Library GC bias: the published 18S ratio and NGS-DOSE's GC model",
                f"One sequencing chemistry: each of {fv['n']:,} genomes of the 1000 Genomes 30× cohort (NovaSeq 2×150, one sequencing centre) appears twice. "
                f"x: the library's GC bias, the rate at which 65%-GC fragments were sequenced relative to the library's mean, measured on its 800 "
                f"single-copy control regions. y: an estimate from the 18S gene divided by the same genome's NGS-DOSE calibrated 45S estimate, which "
                f"removes the differences between people (CV {fmt(bio.get('cn45_cv'), 0, pct=True)}) and leaves what each estimator adds of its own. Orange: the "
                f"18S read-depth ratio of published studies, with no GC model: r = {fmt(fv.get('r'), 2)}. Blue: the same 18S reads under NGS-DOSE's "
                f"fragment-GC model: r = {fmt(mv.get('r'), 2)}. Lines: least-squares fits."
                + (f" Between technologies, whose GC responses differ more, the published ratio is offset by {fmt(abs(rf['offset']), 0, pct=True)} (section 3.4)." if rf else ""))
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
        P.h(f"""<p><strong>How to read the numbers.</strong> Correlations are Pearson's, on values centred within population. The child–midparent
correlation cannot reach 1 even for a perfectly measured heritable trait: half of a child's variance is segregation, which the midparent does
not predict, so r is bounded by about √((1 + ρ)/2) — {fmt((0.5 * (1 + t45["spousal_r"])) ** 0.5, 2) if t45 else "0.71"} at the 45S's spousal
correlation of {fmt(t45["spousal_r"], 2) if t45 else "0"}. The midparent <em>slope</em> b has no such ceiling (it is 1 when every unit of
parental variance reappears in the children), which is why the slope, corrected for the spousal correlation to R = b − ρ(1 − b), is the
number that answers "what share of the measured variance is real". Inference rests on the family-bootstrap interval of R; the permutation
p-value (children shuffled among families) says whether a slope of that size arises by chance.</p>""")
        ts = next((t for t in tr["table"] if t["column"] == tr["scatter_column"]), None)
        trio_stats = (f" Computed on values centred within population: midparent slope b = {fmt(ts['slope'], 2)} ± {fmt(ts['slope_se'], 2)}, spousal correlation "
                      f"ρ = {fmt(ts['spousal_r'], 2)}, reliability R = b − ρ(1 − b) = {fmt(min(ts['R'], 1.0), 2)}"
                      + (" (" + "; ".join((["capped at 1"] if ts["R"] > 1 else []) + ([f"95% CI {fmt(ts['R_lo'], 2)} to {fmt(ts['R_hi'], 2)}"] if "R_lo" in ts else [])) + ")"
                         if ts["R"] > 1 or "R_lo" in ts else "")
                      + ": the share of the measured differences between people that is inherited, 1 for a perfectly measured heritable trait and 0 for pure "
                        "measurement error.") if ts else ""
        P.chart("trio", dict(type="scatter", points=[dict(x=p["mid"], y=p["c"], label=p["child"], extra=[f"father {fmt(p['f'], 0)}, mother {fmt(p['m'], 0)}", p["pop"]]) for p in tr["scatter"]],
                             xlabel="midparent copies", ylabel="child copies", identity=True, fit=True),
                f"Inherited: child against midparent, {ts['label'] if ts else tr['scatter_column']}",
                f"Each dot is one of {len(tr['scatter'])} complete parent–child trios of the 1000 Genomes 30× cohort: x the mean of the two parents' "
                f"copy numbers, y the child's{' (NGS-DOSE)' if 'NGS-DOSE' in (ts or {}).get('label', '') else ''}. Grey diagonal: child equals midparent; "
                f"line: least-squares fit to the points shown.{trio_stats}")
        # the heatmap: every metric the trios were asked about, grouped by what the answer must be
        from .report import TRIO_GROUPS
        by_col = {t["column"]: t for t in tr["table"]}
        stats = [("r_father", "r father"), ("r_mother", "r mother"), ("r_mid", "r midparent"), ("spousal_r", "spousal r"), ("slope", "slope"), ("R", "R")]
        titles = ["child–father correlation", "child–mother correlation", "child–midparent correlation", "spousal correlation", "midparent slope", "reliability R = b − ρ(1 − b)"]
        hm_rows, hm_groups, hm_vals, hm_extra = [], [], [], []
        for key, gl, cols in TRIO_GROUPS:
            for col, label in cols:
                t = by_col.get(col)
                if not t:
                    continue
                hm_rows.append(label); hm_groups.append(gl)
                v = [t.get(k) for k, _ in stats]
                pp = t.get("perm_p")
                x = [f"n = {t['n_trios']}"] * 4 + [f"± {fmt(t['slope_se'], 2)} (SE)" + (f"; permutation p {'< 0.001' if pp < 0.001 else fmt(pp, 3)}" if pp is not None else ""),
                                                   (f"{fmt(t['R_lo'], 2)} to {fmt(t['R_hi'], 2)} (95%)" if "R_lo" in t else "")]
                if fR:
                    f_ = fR.get(col)
                    v.append(f_["R"] if f_ else None)
                    x.append((f"{fmt(f_['R_lo'], 2)} to {fmt(f_['R_hi'], 2)} (95%)" if f_ and "R_lo" in f_ else "") + (f"; n = {f_['n_trios']}" if f_ else ""))
                hm_vals.append(v); hm_extra.append(x)
        cols_hm = [l for _, l in stats] + (["R, fetch"] if fR else [])
        titles += ["reliability R from the fetch counts alone"] if fR else []
        P.h('''<p>Every metric the cohort measures was put to the same test, in four groups: the rDNA estimators, whose transmission is
the claim; the satellite arrays, whose mass is a property of the genome and must be inherited (positive controls); sequence of known
copy number, which has nothing to inherit except the distal junction's whole-copy steps (3.2); and the culture's and the library's properties,
which are not in the nuclear genome (negative controls). A method that measured library artefacts would light up the last group; one that measured nothing would light up
none. HSat1B lives mostly on Yq and passes from father to son only, so its midparent statistics are diluted by design; the split by the sex of parent and child (<a href="#bysex">below</a>) tests it as it is inherited.'''
            + (" The last column repeats the reliability with the cohort layer and the trio test run on the fetch counts alone (3.3).</p>" if fR else "</p>"))
        P.chart("heat", dict(type="heatmap", rows=hm_rows, groups=hm_groups, cols=cols_hm, col_titles=titles, values=hm_vals, extra=hm_extra, lo=0, hi=1, nd=2),
                f"Transmission of every metric through {tr['n_complete']:,} parent–child trios",
                "Rows: every metric measured in the 1000 Genomes 30× cohort, grouped by what its transmission has to be. NGS-DOSE's rDNA estimates are "
                "the claim, shown with the 18S depth ratio of published studies computed from the same reads for comparison; satellite-array mass is "
                "genomic and must be inherited (positive controls); sequence of known copy number has nothing to inherit; the cell line's and the "
                "library's properties are not in the nuclear genome (negative controls). Columns: Pearson correlations of the child with father, mother "
                "and midparent, the parents' (spousal) correlation, the midparent slope b and the reliability R = b − ρ(1 − b)"
                + (", and R with the whole analysis run on the targeted-fetch counts alone" if fR else "")
                + ". Values centred within population; colour from 0 (white) to 1 (blue), value printed in each cell. A real, heritable measurement "
                  "reads near 1 in slope and R; the negative controls near 0. On the page, hover a cell for n, the standard error, the permutation p "
                  "and the interval.")
        rows_t = []
        for t in tr["table"]:
            r_ci = f" ({fmt(t['R_lo'], 2)} to {fmt(t['R_hi'], 2)})" if "R_lo" in t else ""
            err = ("≤ " + fmt(t["error_cv_max"], 1, pct=True)) if "error_cv_max" in t else fmt(t["error_cv"], 1, pct=True)
            pp = t.get("perm_p")
            rows_t.append([t["label"], t.get("group", ""), t["n_trios"], fmt(t["r_mid"], 3), fmt(t["slope"], 3) + " ± " + fmt(t["slope_se"], 3), fmt(t["spousal_r"], 3), fmt(min(t["R"], 1.0), 3) + r_ci,
                           ("< 0.001" if pp is not None and pp < 0.001 else fmt(pp, 3)), fmt(min(t["R_single"], 1.0), 3), err])
        heads = ["metric", "trios", "child–midparent r (Pearson)", "midparent slope b", "spousal r", "reliability R (95% CI)", "permutation p", "single-parent R", "error CV the interval allows"]
        main_rows = [r for r in rows_t if r[1] in ("rDNA", "truth", "culture")]
        P.table([r[:1] + r[2:] for r in main_rows], heads, numeric={1, 2, 3, 4, 5, 6, 7, 8})
        if len(rows_t) > len(main_rows):
            P.h("<details><summary>The same for the satellite arrays (the full table is <code>data/transmission.tsv</code>)</summary>")
            P.table([r[:1] + r[2:] for r in rows_t if r[1] == "satellites"], heads, numeric={1, 2, 3, 4, 5, 6, 7, 8})
            P.h("</details>")
        P.h('<p class="small">Reliability is capped at 1; a slope above 1 is noise around 1, and the interval says how much. A spousal correlation far from zero means members of a family share something other than DNA (a batch), and every reliability in the table is inflated by about as much. The negative-control rows should read near zero.</p>')
        if t45 and have_ci:
            cv = bio.get("cn45_cv")
            P.h(f'''<p>The 45S reliability is <strong>{fmt(min(t45["R"], 1.0), 2)}</strong> ({fmt(t45["R_lo"], 2)} to {fmt(t45["R_hi"], 2)}); the interval allows a
measurement error of at most {fmt(t45.get("error_cv_max"), 1, pct=True)} of a person's value. The variation measured is therefore inherited. Because
people differ in 45S copy number by a CV of {fmt(cv, 0, pct=True) if cv else "about 20%"} and every estimator errs by a few percent, every estimator
has a reliability near 1 within one cohort and one pipeline, and the trios cannot rank them: the paired differences below are the test, and they
are small. Ranking rests on the comparison across technologies (3.4).</p>''')
        if tr["compare"]:
            P.h("<p>Paired family bootstrap of each estimator's reliability minus that of the published 18S depth ratio:</p>")
            P.table([[c["label"], fmt(c["delta"], 3), f"{fmt(c['lo'], 3)} to {fmt(c['hi'], 3)}", fmt(c["p_better"], 3)] for c in tr["compare"]],
                    ["estimator", "ΔR vs the published 18S ratio", "95% CI", "P(better)"], numeric={1, 2, 3})
        # ---- transmission by the sex of parent and child
        nb = tr.get("by_sex") or {}
        bs = nb.get("table") or {}
        if bs:
            from .report import TRIO_GROUPS
            pair_cols = [("father_son", "father → son"), ("father_daughter", "father → daughter"), ("mother_son", "mother → son"), ("mother_daughter", "mother → daughter")]
            hs_rows, hs_groups, hs_vals, hs_extra = [], [], [], []
            for _, gl, cols in TRIO_GROUPS:
                for col, label in cols:
                    d = bs.get(col)
                    if not d:
                        continue
                    hs_rows.append(label); hs_groups.append(gl)
                    hs_vals.append([d[k]["r"] if k in d else None for k, _ in pair_cols])
                    hs_extra.append([(f"{d[k]['n']} pairs; 95% CI {fmt(d[k]['r_lo'], 2)} to {fmt(d[k]['r_hi'], 2)}; slope {fmt(d[k]['slope'], 2)} ± {fmt(d[k]['slope_se'], 2)}"
                                      if k in d else "") for k, _ in pair_cols])
            half = sorted((d[k]["r_hi"] - d[k]["r_lo"]) / 2 for d in bs.values() for k, _ in pair_cols if k in d)
            half = half[len(half) // 2] if half else None
            contrast = lambda col: (bs[col].get("father_contrast") or {}).get("z", 0) or 0
            # sex linkage is a property of genomic sequence: only the rDNA and the satellite arrays are read that way
            genomic = lambda col: bs[col]["group"] in ("rDNA", "satellites")
            ylinked = [col for col in bs if genomic(col) and contrast(col) > 3]
            xlinked = [col for col in bs if genomic(col) and contrast(col) < -3]
            other_far = [col for col in bs if not genomic(col) and abs(contrast(col)) > 3]
            hint = [col for col in bs if 2 < abs(contrast(col)) <= 3]
            lab = lambda col: bs[col]["label"]
            P.h(f"""<h3 id="bysex">By the sex of parent and child</h3>
<p>The midparent test assumes that each parent passes on half of what they carry, as they do for sequence on the autosomes. Sequence on the
sex chromosomes is inherited differently: a father passes his Y whole to his sons and none of it to his daughters, and his X whole to his
daughters and none of it to his sons, while a mother passes one of her two X chromosomes to every child. Split by the sex of parent and child,
the correlation of child with parent is therefore about one half in all four pairings for a well-measured autosomal quantity; close to 1
from father to son and 0 from father to daughter for a Y-linked one; 0 from father to son and high from father to daughter for an X-linked
one. Correlations, not slopes, are compared here: a slope also carries the ratio of the child's spread to the parent's, which differs between
men and women for a sex-linked quantity. Values are centred within population and sex, so that the difference between men and women is not
read as transmission. {nb['n_sons']} of the trios have a son and {nb['n_daughters']} a daughter, so each pairing rests on about half the
trios{f" and a correlation is uncertain by about ±{fmt(half, 2)} (95%)" if half else ""}.</p>""")
            P.chart("heat_sex", dict(type="heatmap", rows=hs_rows, groups=hs_groups, cols=[l for _, l in pair_cols], col_titles=[f"slope of child on parent, {l}" for _, l in pair_cols],
                                     values=hs_vals, extra=hs_extra, lo=0, hi=1, nd=2),
                    f"Transmission by the sex of parent and child: {nb['n_sons']} sons, {nb['n_daughters']} daughters",
                    "Rows: every metric of the trio test above, grouped as there. Columns: the Pearson correlation of child with parent for each pairing "
                    "of the parent's and the child's sex, on values centred within population and sex. A well-measured autosomal quantity reads about 0.5 "
                    "in every column; a Y-linked one close to 1 from father to son and 0 from father to daughter; an X-linked one 0 from father to son. "
                    "Colour from 0 (white) to 1 (blue), value printed" + (f"; each value is uncertain by about ±{fmt(half, 2)} (95%)" if half else "")
                    + ". On the page, hover a cell for the pairs, the interval and the slope.")
            notes = []
            for group, word in ((ylinked, "sequence on the Y chromosome: fathers pass more of it to their sons than to their daughters"),
                                (xlinked, "sequence on the X chromosome: fathers pass more of it to their daughters than to their sons")):
                if group:
                    notes.append(" ".join(f"{lab(c)}: father–son r = {fmt(bs[c]['father_son']['r'], 2)}, father–daughter r = {fmt(bs[c]['father_daughter']['r'], 2)} "
                                          f"(the difference is {fmt(contrast(c), 1)} standard errors on Fisher's scale)." for c in group)
                                 + f" By that contrast {'this metric carries' if len(group) == 1 else 'these metrics carry'} {word}, and the midparent reliability "
                                   "above understates how faithfully it is inherited.")
            why_not = {"truth": "sequence of known copy number, which varies between people only by the rare whole-copy steps of a few families (3.2), so a "
                                "handful of carrier parents who happen to have sons or daughters decides it",
                       "culture": "a property of the cell line or the library, not of the nuclear genome"}
            for c in other_far:
                notes.append(f"{sentence_case(lab(c))}: father–son r = {fmt(bs[c]['father_son']['r'], 2)}, father–daughter r = {fmt(bs[c]['father_daughter']['r'], 2)} "
                             f"({fmt(contrast(c), 1)} standard errors), which is no sign of the sex chromosomes: it is {why_not.get(bs[c]['group'], 'not genomic sequence')}.")
            if hint:
                notes.append(f"Between two and three standard errors, as one or two of {len(bs)} metrics would be by chance (the held-out autosomal "
                             "sequence, which has nothing to inherit, is a guide): " + "; ".join(
                    f"{lab(c)}, r {fmt(bs[c]['father_son']['r'], 2)} against {fmt(bs[c]['father_daughter']['r'], 2)} (z = {fmt(contrast(c), 1)})" for c in hint) + ".")
            notes.append(("No other metric's" if ylinked or xlinked or other_far else "No metric's") + " father–son and father–daughter correlations differ by more than three standard errors.")
            P.h("<p>" + " ".join(notes) + "</p>")
            pts_all = tr.get("points") or {}
            flagged = [c for c in ylinked + xlinked if c in pts_all]
            if flagged:
                P.h('<div class="grid2">')
                for col in flagged:
                    d, unit = bs[col], unit_of(col)
                    pts = [dict(x=p[3], y=p[2], label=p[0], si=0 if p[5] == "M" else 1, extra=[f"midparent {fmt(p[1], 2)}"]) for p in pts_all[col] if p[5] in ("M", "F")]
                    P.chart(f"sexfather_{re.sub(r'[^A-Za-z0-9]', '_', col)}", dict(type="scatter", points=pts, legend=["son", "daughter"], xlabel=f"father, {unit}",
                                                                                   ylabel=f"child, {unit}", diagonal=True, fit=True, fit_labels=["sons", "daughters"]),
                            f"{sentence_case(d['label'])}: children against their fathers, sons and daughters apart",
                            f"Each dot is one complete parent–child trio of the 1000 Genomes 30× cohort: x the father's {d['label']}, y the child's, in {unit}; "
                            f"sons blue ({d['father_son']['n']}), daughters orange ({d['father_daughter']['n']}), each with its least-squares line; grey: "
                            f"equality, where it falls inside the plot. On values centred within population and sex, the slope of child on father is "
                            f"{fmt(d['father_son']['slope'], 2)} for sons (r = {fmt(d['father_son']['r'], 2)}) and {fmt(d['father_daughter']['slope'], 2)} for "
                            f"daughters (r = {fmt(d['father_daughter']['r'], 2)}).")
                P.h("</div>")
        # ---- class by class: every metric's transmission in the same terms, and what this design can and cannot tell apart
        if have_ci and tr["table"] and "sd_ratio_lo" in tr["table"][0]:
            from .report import TRIO_GROUPS
            from .trios import PAIRS
            byc = {t["column"]: t for t in tr["table"]}
            bsx = (tr.get("by_sex") or {}).get("table") or {}
            zf = lambda col: ((bsx.get(col) or {}).get("father_contrast") or {}).get("z", 0) or 0
            linked = {col: ("Y" if zf(col) > 0 else "X") for col, d in bsx.items() if d["group"] in ("rDNA", "satellites") and abs(zf(col)) > 3}
            ylinked = [c for c in linked if linked[c] == "Y" and c in byc]
            cvt = lambda t: t["parent_sd"] / t["parent_mean"] if t.get("parent_mean") else float("nan")
            sg = lambda x, nd=2: fmt(x, nd).replace("-", "−")
            iv = lambda t, k, nd=2: f"{sg(t[k + '_lo'], nd)} to {sg(t[k + '_hi'], nd)}" if t.get(k + "_lo") is not None else ""
            signed = lambda x: ("+" if x >= 0 else "−") + fmt(abs(x), 1, pct=True)
            level = lambda t: (f"{fmt(abs(t['mean_ratio'] - 1), 1, pct=True)} {'more' if t['mean_ratio'] > 1 else 'less'} than their parents "
                               f"({signed(t['mean_ratio_lo'] - 1)} to {signed(t['mean_ratio_hi'] - 1)})")
            # the part of a person's value not inherited, CV x sqrt(1 - R), at the lower of the two reliabilities (a batch's scale must not hide it)
            pc = lambda x: fmt(x, 1 if x < 0.1 else 0, pct=True)
            lost = lambda t, hi=False: cvt(t) * max(0.0, 1 - (min(t["R_lo"], t["R_rescaled_lo"]) if hi else min(t["R"], t["R_rescaled"]))) ** 0.5
            het = lambda col: (bsx.get(col) or {}).get("heterogeneity") or {}
            hp = lambda col: ("≤ 0.001" if het(col)["p"] <= 0.001 else fmt(het(col)["p"], 3)) if het(col) else "–"
            pr = lambda col, k: fmt(bsx[col][k]["r"], 2) if k in (bsx.get(col) or {}) else "–"
            name = {col: label for _, _, cols in TRIO_GROUPS for col, label in cols}
            nm = lambda col: esc(name.get(col, col))
            bt = tr.get("batches") or {}
            top = lambda d: max(d.items(), key=lambda kv: kv[1]) if d else (None, 0)
            (kb, kn), (pb, pn) = top(bt.get("child") or {}), top(bt.get("parent") or {})
            apart = bool(bt) and kb != pb and kn >= 0.9 * bt["n"] and pn >= 0.9 * 2 * bt["n"]
            blab = lambda b: f"{int(b):,}" if str(b).isdigit() else esc(str(b))
            # the children's spread against their parents', where the interval excludes 1; largest departure first
            moved = sorted((t for t in tr["table"] if t["sd_ratio_lo"] > 1 or t["sd_ratio_hi"] < 1), key=lambda t: -abs(math.log(t["sd_ratio"])))
            # the four sex pairings: the metrics that differ beyond chance, the sex-linked ones aside and the three 45S estimators (one set of reads) counted once
            tested = [c for c in bsx if c not in linked and het(c)]
            uneven = [c for c in tested if het(c)["p"] < 0.05]
            once = lambda cols: len([c for c in cols if not c.startswith("rDNA45S")]) + any(c.startswith("rDNA45S") for c in cols)
            u45 = [c for c in uneven if c.startswith("rDNA45S")]
            brief = lambda c: "the 45S" if c.startswith("rDNA45S") else "the 5S" if c.startswith("rDNA5S") else nm(c).replace(" array", "")
            uneven_list = ", ".join(([f"the 45S by {'all three' if len(u45) == 3 else len(u45)} of its estimators"] if len(u45) > 1 else [brief(c) for c in u45])
                                    + [brief(c) for c in uneven if not c.startswith("rDNA45S")])
            said = []
            spelt = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five", 6: "Six", 7: "Seven", 8: "Eight", 9: "Nine"}
            uneven_txt = (f"{spelt.get(once(uneven), once(uneven))} of the {once(tested)} metrics (the 45S's three estimators, which share their reads, counted once) reach p &lt; 0.05 "
                          f"by the same test ({uneven_list}), where chance alone would give about {spelt.get(max(1, round(0.05 * once(tested))), '').lower() or round(0.05 * once(tested))}"
                          + ("; the held-out autosomal sequence, which has nothing to inherit, is among them" if "truth.auto" in uneven else "")
                          + ". It is a lead for the full cohort, not a finding")

            def pairings(col):
                s = bsx.get(col) or {}
                if not all(k in s for k, _, _ in PAIRS):
                    return ""
                txt = (f"split by the sex of parent and child, the four correlations are {pr(col, 'father_son')} (father–son), {pr(col, 'father_daughter')} "
                       f"(father–daughter), {pr(col, 'mother_son')} (mother–son) and {pr(col, 'mother_daughter')} (mother–daughter)")
                if not het(col):
                    return txt + "."
                if het(col)["p"] >= 0.05:
                    return txt + f", no more different than chance makes them (p = {hp(col)})."
                lo_k = min((k for k, _, _ in PAIRS), key=lambda k: s[k]["r"])
                hi_k = max((k for k, _, _ in PAIRS), key=lambda k: s[k]["r"])
                first = not said
                said.append(col)
                return (txt + f", more different than chance makes them (p = {hp(col)}): {lo_k.replace('_', '–')} {pr(col, lo_k)} against "
                        f"{hi_k.replace('_', '–')} {pr(col, hi_k)}. Inheritance on the autosomes gives no reason for that, since a child's sex does not "
                        "depend on which of a parent's chromosomes it received" + (f". {uneven_txt}." if first else ", and the caution under the 45S applies."))

            P.h("""<h3 id="byclass">Class by class</h3>
<p>Every metric's transmission in the same terms: how much people differ (the parents' coefficient of variation, CV); the share of that
difference the children inherit (R, with its 95% interval); the part of a person's value that is not inherited, as a percentage of the value
(CV × √(1 − R): measurement error, or change in the cell line; the interval allows at most the figure in brackets); the children's spread and
level against their parents'; and whether the four pairings of the parent's and the child's sex differ, which for a quantity carried on the
autosomes they should not.</p>""")
            if apart:
                dep = byc.get("depth")
                wide_ok = [t for t in moved if t["sd_ratio"] > 1 and t["R_lo"] > 1]            # a slope above 1: new variation cannot raise it
                narrow = [t for t in moved if t["sd_ratio"] < 1 and t["group"] in ("rDNA", "satellites")]   # a narrower spread: new variation cannot narrow it
                P.h(f"""<p><strong>Parents and children were sequenced apart.</strong> {"All" if kn == bt["n"] else f"{kn:,} of the"} {bt["n"]:,} children are among the
{blab(kb)} related genomes that the 1000 Genomes 30× release added to its original {blab(pb)} (the release batch in NGS-PCA's table), and {pn:,} of
their {2 * bt["n"]:,} parents are among the {blab(pb)}; {bt["shared"]:,} {"family has" if bt["shared"] == 1 else "families have"} a parent in the child's batch. A batch shared
within a family, which could imitate inheritance, is therefore all but absent. But a difference between the batches is a difference between the
generations. Where the later batch reads a quantity on a slightly larger or smaller scale, the children's spread changes, and the midparent slope
and R change with it by the same factor, while R computed with the children's values first put on their parents' scale does not. New
variation arising in the children would also widen their spread, and lowers the rescaled R instead; the two values bracket the reliability.
They part where the children's spread departs furthest from their parents': """
                    + "; ".join(f"{nm(t['column'])}, {fmt(t['sd_ratio'], 2)}× ({iv(t, 'sd_ratio')})" for t in moved[:6]) + "."
                    + (" " + " ".join([f"The {nm(t['column'])}'s wider spread cannot be new variation arising in the children, because its slope rises with it "
                                        f"(R {iv(t, 'R')})." for t in wide_ok]
                                       + [f"{'Nor can' if wide_ok else 'Neither can'} the {nm(t['column'])}'s narrower one." for t in narrow])
                       + " New variation widens the children's spread and leaves the slope alone; a batch that reads the children on a slightly different "
                         "scale moves both together." if wide_ok or narrow else "")
                    + (" The Mendelian reliability shown in earlier versions of this page's table, 1.5 − Var(child − midparent) / Var(parent), equals "
                       "R + 1 + ρ/2 − s², with s the ratio of spreads, so on this design it measured the batch"
                       + (f": sequencing depth, which no child inherits (R = {sg(dep['R'])}), reads {sg(dep['R_mendel'])} by it" if dep else "")
                       + ". It is no longer shown; <code>data/transmission.tsv</code> keeps it.</p>"))
            # the forest: R and its interval, and R on the parents' scale, for every metric in its group
            f_rows, f_groups, est, lo, hi, est2, extra, rows_c = [], [], [], [], [], [], [], []
            for _, gl, cols in TRIO_GROUPS:
                for col, label in cols:
                    t = byc.get(col)
                    if not t or t.get("R_lo") is None:
                        continue
                    lk = linked.get(col)
                    f_rows.append(label + (f" ({lk}-linked)" if lk else "")); f_groups.append(gl)
                    est.append(t["R"]); lo.append(t["R_lo"]); hi.append(t["R_hi"]); est2.append(t["R_rescaled"])
                    extra.append(f"people differ by a CV of {fmt(cvt(t), 1, pct=True)}; the children spread {fmt(t['sd_ratio'], 2)}× as much as their "
                                 f"parents and read {fmt(abs(t['mean_ratio'] - 1), 1, pct=True)} {'more' if t['mean_ratio'] > 1 else 'less'}")
                    rows_c.append([label + (f" ({lk}-linked)" if lk else ""), "–" if lk else fmt(cvt(t), 1, pct=True), sg(t["R"]), sg(t["R_rescaled"]),
                                   "–" if lk else f"{pc(lost(t))} (≤ {pc(lost(t, True))})",
                                   f"{fmt(t['sd_ratio'], 2)} ({fmt(t['sd_ratio_lo'], 2)}–{fmt(t['sd_ratio_hi'], 2)})",
                                   f"{fmt(t['mean_ratio'], 2)} ({fmt(t['mean_ratio_lo'], 2)}–{fmt(t['mean_ratio_hi'], 2)})",
                                   (f"{lk}-linked, " if lk else "") + hp(col)])
            span = [v for v in lo + hi + est + est2 if isinstance(v, float) and math.isfinite(v)] + [0.0, 1.0]
            P.chart("forest", dict(type="forest", rows=f_rows, groups=f_groups, est=est, lo=lo, hi=hi, est2=est2, extra=extra, refs=[0, 1],
                                   xmin=math.floor(min(span) * 4) / 4, xmax=math.ceil(max(span) * 4) / 4,
                                   xlabel="reliability R: the share of the differences between people that is inherited",
                                   legend=["R from the midparent slope, with its 95% interval", "R with the children put on their parents' scale"]),
                    f"Reliability of every metric, class by class, through {tr['n_complete']:,} trios",
                    "Each row is one metric measured in the 1000 Genomes 30× cohort, grouped by what its transmission has to be: NGS-DOSE's rDNA "
                    "estimates (the claim, with the 18S depth ratio of published studies computed from the same reads), satellite-array mass (genomic, "
                    "so it must be inherited), sequence of known copy number, and the cell line's and the library's properties (not in the nuclear "
                    f"genome). Filled dot and bar: the reliability R = b − ρ(1 − b) from the slope b of child on midparent in {tr['n_complete']:,} complete "
                    "parent–child trios, corrected for the spousal correlation ρ, with its 95% family-bootstrap interval: 1 for a perfectly measured "
                    "heritable quantity, 0 for one the children do not inherit. Open ring: the same with each child's value first rescaled to the "
                    "parents' spread"
                    + (", because every child was sequenced in a later batch than its parents and a batch can read on a slightly different scale" if apart else "")
                    + ". Values above 1 are noise around 1, or a larger scale in the children's measurements. Lines at 0 and 1; values centred within "
                    "population."
                    + (" " + "; ".join(f"The {name.get(c, c)} lies mostly on the Y chromosome and passes from father to son (r = {pr(c, 'father_son')}), "
                                       "which the midparent does not measure" for c in ylinked) + "." if ylinked else "")
                    + " On the page, hover a row for the spread between people and the children's spread and level against their parents'.")
            P.table(rows_c, ["metric", "CV", "R", "R, rescaled", "not inherited", "children's spread", "children's level", "sex pairings: p"],
                    numeric={1, 2, 3, 4, 5, 6, 7}, wrap=True, cls="compact")
            nperm = next((het(c)["n_perm"] for c in bsx if het(c)), None)
            P.h('<p class="small">CV: how much people differ, the parents\' standard deviation over their mean. R: from the midparent slope, and rescaled: with the '
                "children's values first put on their parents' spread (both with their 95% intervals in the figure). Not inherited: CV × √(1 − R) at the "
                "lower of the two, as a percentage of a person's value, and in brackets the most the intervals allow. Children's spread: their standard "
                "deviation over their parents'; children's level: their mean over their parents'; both on values centred within population, with 95% "
                "family-bootstrap intervals. Sex pairings: how often the four child–parent correlations (father–son, father–daughter, mother–son, "
                f"mother–daughter) differ as much as these when the children's sexes are shuffled among the families and each family's parents swap roles{f' ({nperm:,} shuffles)' if nperm else ''}. The "
                "textbook test, which assumes normal values, is far too small for skewed ones"
                + (f" (the distal junction's four pairings: {'&lt; 0.001' if het('DJ.cn')['p_normal'] < 0.001 else fmt(het('DJ.cn')['p_normal'], 3)} by it, "
                   f"{hp('DJ.cn')} by the shuffles)" if het("DJ.cn") else "")
                + ".</p>")
            # ---- the rDNA
            rep = (data.get("replicates") or {}).get("table") or {}
            t, t5 = byc.get("rDNA45S.cn"), byc.get("rDNA5S.cn")
            if t:
                ceiling = (0.5 * (1 + t["spousal_r"])) ** 0.5
                rc = rep.get("calibrated") or {}
                alike = abs(t["r_father"] - t["r_mother"]) < 0.15
                big = max((x for x in moved if x["group"] in ("rDNA", "satellites")), key=lambda x: abs(x["sd_ratio"] ** 2 - 1), default=None)
                room = max(0.0, t["sd_ratio_hi"] ** 2 - 1)
                P.h(f"""<p><strong>45S rDNA.</strong> People differ by a CV of {fmt(cvt(t), 0, pct=True)}. The children inherit
{"all of it that the trios can resolve" if t["R_lo"] <= 1 <= t["R_hi"] else "most of it" if t["R_hi"] < 1 else "all of it"}: R = {sg(t["R"])} ({iv(t, "R")}) from the slope and {sg(t["R_rescaled"])}
({iv(t, "R_rescaled")}) with the children on their parents' scale, and at most {fmt(lost(t, True), 1, pct=True)} of a person's value is not inherited{f" (the pilot's replicates across technologies put the measurement's own error at {fmt(rc['sd_log_ratio'] / 2 ** 0.5, 1, pct=True)}, section 3.4)" if rc.get("sd_log_ratio") else ""}. The
child–midparent correlation, {fmt(t["r_mid"], 2)} ({iv(t, "r_mid")}), {"reaches" if t["r_mid_hi"] >= ceiling else "falls short of"} the ceiling of {fmt(ceiling, 2)} that segregation sets
for a perfectly measured trait at the spousal correlation of {sg(t["spousal_r"])} ({iv(t, "spousal")}). {"The father and the mother pass it on alike" if alike else "The two parents pass it on"}
(child–father r = {fmt(t["r_father"], 2)}, child–mother r = {fmt(t["r_mother"], 2)}){", as they should for arrays on the short arms of five autosomes" if alike else ""};
{pairings("rDNA45S.cn")} The children spread {fmt(t["sd_ratio"], 2)}× as much as their parents ({iv(t, "sd_ratio")}) and read {level(t)}{": a difference of generation or of batch, which this design cannot separate, and one the slopes and correlations, unchanged when every child is shifted by the same amount, do not see" if not t["mean_ratio_lo"] <= 1 <= t["mean_ratio_hi"] else ""}.
An array is rearranged in more than one meiosis in ten ({doi("Stults et al. 2008", "10.1101/gr.6858507")}). Rearrangement as likely to lengthen
an array as to shorten it scatters the children about the midparent without flattening the slope, and widens their spread: if the batch changed
nothing, by at most {fmt(room, 0, pct=True)} of the parents' variance, the top of the interval{f". The batch alone moves the variance of the {nm(big['column'])} by {fmt(abs(big['sd_ratio'] ** 2 - 1), 0, pct=True)}, so the trios bound that variation; they do not measure it" if big and abs(big["sd_ratio"] ** 2 - 1) >= 0.5 * room else ""}.</p>""")
            if t5:
                r5 = rep.get("rDNA5S") or {}
                similar = t is not None and abs(cvt(t5) - cvt(t)) < 0.05
                meas = r5["sd_log_ratio"] / 2 ** 0.5 if r5.get("sd_log_ratio") else None
                P.h(f"""<p><strong>5S rDNA.</strong> People differ {"as much as in the 45S" if similar else "by"} (CV {fmt(cvt(t5), 0, pct=True)}),
{"but the children inherit less of it, or less certainly" if t and t5["R"] < t["R"] - 0.1 else "and the children inherit"}: R = {sg(t5["R"])} ({iv(t5, "R")}), {sg(t5["R_rescaled"])}
({iv(t5, "R_rescaled")}) on the parents' scale, at a spousal correlation of {sg(t5["spousal_r"])}.{f" The interval runs from an array transmitted as faithfully as the 45S to one with {fmt(1 - t5['R_lo'], 0, pct=True)} of its differences not inherited." if t5["R_hi"] >= 1 and t5["R_lo"] < 0.8 else ""}
New variation arising in the children at random would scatter them about the midparent line without flattening it (they spread
{fmt(t5["sd_ratio"], 2)}× as much as their parents, {iv(t5, "sd_ratio")}). A slope below 1 means that part of the parents' measured differences is not passed
on: error, change in the cell line, or arrays that drift toward the average as they are transmitted; {fmt(lost(t5), 0, pct=True)} of a person's value at the
estimate, at most {fmt(lost(t5, True), 0, pct=True)}.{f" The pilot's replicates across technologies put the measurement's own error near {fmt(meas, 1, pct=True)} (section 3.4), so at the estimate the part not inherited would lie in the cell lines or in transmission rather than in the measurement." if meas and lost(t5) > 2 * meas else ""}
The array is on chromosome 1 (1q42); {pairings("rDNA5S.cn")}</p>""")
            # ---- the satellite arrays, the telomeric repeat aside
            sats = [byc[c] for _, _, cols in TRIO_GROUPS for c, _ in cols if c in byc and byc[c]["group"] == "satellites" and c != "TEL.mass_Mb"]
            if sats:
                plain = [t for t in sats if t["column"] not in linked and t not in moved]
                low = sorted((t for t in plain if cvt(t) < 0.08), key=cvt)
                wide = sorted((t for t in plain if cvt(t) >= 0.1), key=cvt)
                errs = [lost(t) for t in plain]
                scaled = [t for t in sats if t["column"] not in linked and t in moved]
                pmax = max((t["perm_p"] for t in sats if t["column"] not in linked and t.get("perm_p") is not None), default=None)
                short = lambda col: nm(col).replace(" array", "")
                names = lambda ts: ", ".join(short(t["column"]) for t in ts[:-1]) + (" and " if len(ts) > 1 else "") + short(ts[-1]["column"])
                span = lambda ts, f, nd=0: f"{fmt(min(map(f, ts)), nd, pct=True)} to {fmt(max(map(f, ts)), nd, pct=True)}"
                sat_uneven = [c for c in uneven if byc.get(c, {}).get("group") == "satellites"]
                P.h("<p><strong>Satellite arrays.</strong> "
                    + (f"Every class is inherited (permutation p ≤ {fmt(pmax, 3)}). " if pmax is not None and pmax < 0.05 else "")
                    + (f"In every class whose children spread as their parents do, the part of a person's value that is not inherited lies between "
                       f"{fmt(min(errs), 0, pct=True)} and {fmt(max(errs), 0, pct=True)}, and how much of R that costs depends on how much people differ. " if errs else "")
                    + (f"The classes that vary least between people ({names(low)}, CV {span(low, cvt, 1)}) read the lowest R "
                       f"({', '.join(sg(t['R']) for t in low[:-1])} and {sg(low[-1]['R'])}): the same few percent is a larger share of a smaller spread. " if len(low) > 1 else "")
                    + (f"Those that vary most ({names(wide)}, CV {span(wide, cvt)}) read {sg(min(t['R'] for t in wide))} to {sg(max(t['R'] for t in wide))}. " if len(wide) > 1 else "")
                    + (f"{'Two classes are' if len(scaled) == 2 else 'One class is' if len(scaled) == 1 else f'{len(scaled)} classes are'} read on a different scale in "
                       "the children (above): " + "; ".join(f"{short(t['column'])} spreads {fmt(t['sd_ratio'], 2)}× as much in the children as in their parents "
                                                           f"({iv(t, 'sd_ratio')}), and its R reads {sg(t['R'])} from the slope but {sg(t['R_rescaled'])} on the "
                                                           "parents' scale" for t in scaled) + "." if scaled else "")
                    + (" " + " ".join(f"{short(c)} lies mostly on the Y chromosome: it passes from father to son at r = {pr(c, 'father_son')} and to daughters at "
                                      f"{pr(c, 'father_daughter')}, which the midparent statistics do not measure (<a href=\"#bysex\">by sex</a>, above)." for c in ylinked) if ylinked else "")
                    + (" By sex, " + "; ".join(f"{short(c)}'s four pairings differ at p = {hp(c)}" for c in sat_uneven)
                       + (", and the caution under the 45S applies." if said else f". {uneven_txt}.") if sat_uneven else "")
                    + "</p>")
            t = byc.get("TEL.mass_Mb")
            if t:
                genomic = [x for x in tr["table"] if x["group"] in ("rDNA", "satellites") and x["column"] not in linked]
                P.h(f"""<p><strong>Telomeric repeat.</strong> {"The least inherited genomic quantity" if t["R"] <= min(x["R"] for x in genomic) else "Weakly inherited"}:
R = {sg(t["R"])} ({iv(t, "R")}), at {"the largest" if cvt(t) >= max(cvt(x) for x in genomic) else "a large"} spread between people (CV {fmt(cvt(t), 0, pct=True)}),
so {fmt(lost(t), 0, pct=True)} of a person's value is not inherited. Telomeres shorten with age and change in culture; the children, younger than
their parents{" and sequenced in the later batch" if apart else ""}, read {level(t)}.</p>""")
            ta, tdj = byc.get("truth.auto"), byc.get("DJ.cn")
            if ta or tdj:
                P.h("<p><strong>Known copy number.</strong> "
                    + (f"Held-out autosomal sequence has nothing to inherit, and the children inherit nothing of it: R = {sg(ta['R'])} ({iv(ta, 'R')}), with "
                       f"people differing by {fmt(cvt(ta), 1, pct=True)}, the measurement's floor. " if ta else "")
                    + ("The distal junction is 10 copies in most genomes; the whole-copy steps some carry (3.2) run in families"
                       + (" and are inherited" if tdj["R_lo"] > 0 else ", and the trios cannot yet tell whether they are inherited")
                       + f": R = {sg(tdj['R'])} ({iv(tdj, 'R')}), people differing by {fmt(cvt(tdj), 1, pct=True)}." if tdj else "")
                    + "</p>")
            cul = [byc[c] for _, _, cols in TRIO_GROUPS for c, _ in cols if c in byc and byc[c]["group"] == "culture"]
            if cul:
                above = [t for t in cul if t["R_lo"] > 0]
                below = [t for t in cul if t["R_hi"] < 0]
                P.h("<p><strong>Culture and library.</strong> "
                    + ("No property of the culture or the library is inherited" if not above else "Only " + " and ".join(nm(t["column"]) for t in above) + " read above zero")
                    + f": R runs from {sg(min(t['R'] for t in cul))} to {sg(max(t['R'] for t in cul))}"
                    + (", and every interval includes zero." if not above and not below else
                       ", and every interval but " + " and ".join(f"that of {nm(t['column'])} ({iv(t, 'R')})" for t in above + below) + " includes zero"
                       + ("; an interval wholly below zero is not inheritance" if below else "")
                       + f", and one such interval among {len(tr['table'])} metrics is what chance gives at 95%.")
                    + " These are also the metrics whose spread differs most between the generations"
                    + (f" (the children's sequencing depth varies {fmt(byc['depth']['sd_ratio'], 2)}× as much as their parents')" if "depth" in byc else "")
                    + ": properties of the batch, not of the family.</p>")
        P.h('<p class="small">A printable assessment built from these tables, with a child-against-midparent scatter for every metric, the fetch check and the assembly comparison: <a href="trio_report.pdf">trio_report.pdf</a>. Every trio\'s values: <code>data/trios.tsv</code>; the split by the sex of parent and child: <code>data/transmission_by_sex.tsv</code>'
            + ('; every metric\'s child against midparent: <a href="#supp">supplementary figures</a>' if (tr.get("points") or {}) and tr["n_complete"] >= 10 else "") + ".</p>")
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
        P.h(f'''<p>On the {hall["n"]:,} genomes shared so far with Hall, Turner &amp; Queitsch (2021), their 18S value against the same 18S depth ratio computed
from NGS-DOSE's counts (no GC model), halved to their per-haploid scale: r = <strong>{fmt(hall["flat"].get("r"), 3)}</strong>, their values {fmt(hall["flat_ratio"], 3)}× ours. Re-applying
their exclusion of duplicate-flagged reads to our counts brings the ratio to {fmt(hall["dup_corrected_ratio"], 3)} (SD {fmt(hall["dup_corrected_ratio_sd"], 3)}):
the offset is the duplicate flag, which marks fewer reads inside the collapsed rDNA than outside it (section 4). Against the calibrated
estimate ({esc(hall["calibrated_column"])}/2): r = {fmt(hall["calibrated"].get("r"), 3)}, ratio {fmt(hall["calibrated_ratio"], 3)}.</p>''')
        g65m = (gcb.get("gc65") or {}).get("median")
        P.chart("hall", dict(type="scatter", points=[dict(x=p["cal"], y=p["theirs"], label=p["sample"]) for p in hall["points"]], xlabel="NGS-DOSE, calibrated 45S / 2", ylabel="Hall et al. 2021, 18S", identity=True, fit=True),
                "The same files, measured independently: NGS-DOSE against Hall et al. (2021)",
                f"Each dot is one of {hall['n']:,} genomes whose rDNA copy number Hall, Turner &amp; Queitsch (Sci Rep 2021) published from the same "
                f"1000 Genomes CRAMs. x: NGS-DOSE's calibrated 45S estimate, halved to their per-haploid scale; y: their 18S copy number, read depth "
                f"relative to chromosome 1 with duplicate-flagged reads excluded. Diagonal: equality; line: least-squares fit. r = "
                f"{fmt(hall['calibrated'].get('r'), 3)}; their values run {fmt(hall['calibrated_ratio'], 2)}× NGS-DOSE's, because their depth ratio has "
                f"no GC model or calibration{f' (these libraries sequence 65%-GC fragments at {fmt(g65m, 2)}× their mean rate, and the rDNA is GC-rich)' if g65m else ''} and "
                f"drops duplicate-flagged reads. The same 18S ratio computed from NGS-DOSE's counts matches theirs at r = {fmt(hall['flat'].get('r'), 3)}, "
                f"{fmt(hall['flat_ratio'], 2)}×, and {fmt(hall['dup_corrected_ratio'], 2)}× once their duplicate exclusion is applied.")
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
                "Mitochondrial genomes per cell: NGS-DOSE against NGS-PCA, same files",
                f"Each dot is one of {nq['n']:,} genomes of the 1000 Genomes 30× cohort measured two ways from the same CRAM. x: NGS-DOSE, fragment "
                f"ends in a region of the mitochondrial genome under the same library model as the rDNA; y: NGS-PCA, twice the mean chrM coverage over the median "
                f"autosomal coverage (mosdepth, 1-kb bins, duplicate-flagged reads excluded). Diagonal: equality; line: least-squares fit. r = "
                f"{fmt(mt['r'], 3)}; NGS-PCA reads {fmt(mt['ratio']['median'], 2)}× NGS-DOSE, the duplicate flag, which saturates on a 16.6-kb genome "
                f"sequenced to thousands-fold depth.")
        P.chart("qc_chrX", dict(type="scatter", x="truth.chrX", y="ngspca.chrX", group=dict(col="sex_inferred", levels=sex_levels), xlabel="NGS-DOSE, chrX copies (60 regions)", ylabel="NGS-PCA, 2 × chrX coverage ratio", identity=True),
                "chrX copies: NGS-DOSE against NGS-PCA, same files",
                f"Each dot is one of {nq['n']:,} genomes of the 1000 Genomes 30× cohort, men blue and women orange. x: NGS-DOSE, 60 chrX regions; "
                f"y: NGS-PCA, twice the chrX-to-autosome coverage ratio. Diagonal: equality. r = {fmt(xx['r'], 4)}. Women below the cluster have lost "
                f"an X in part of their cell line, and read so by both routes.")
        P.h("</div>")
    P.end()

    # ---------------------------------------------------------------- 3.8 assemblies
    P.section("assemblies", "3.8 Satellite arrays against long-read assemblies", "Assemblies")
    hp = sat.get("hprc")
    if hp:
        st = hp["stats"]
        good = [cls for cls, st_ in st.items() if st_.get("n", 0) >= 4 and st_.get("pearson", 0) >= 0.95]
        and_ = lambda xs: ", ".join(xs[:-1]) + (" and " if len(xs) > 1 else "") + xs[-1]
        relative = [cls for cls in st if SAT_FAMILY.get(cls, ("", 1.0))[1] < 0.8]
        P.h(f'''<p>The satellite arrays are measured by the same k-mer machinery as the rDNA, and unlike the rDNA they have a truth:
{hp["n_samples"]} of the genomes counted so far belong to people with an HPRC release-2 assembly, built from long reads for both haplotypes,
whose CenSat annotation gives the length of every satellite array. Summed over the two haplotypes, that is the person's array mass for each
family, in megabases. NGS-DOSE estimates the same mass from the short reads, from the reads that carry the family's k-mers. The rDNA itself
cannot be compared this way: the assemblies do not close its arrays (<a href="#why">section 1</a>).</p>
<p>If both measurements are right, a person's two values are equal. Where the k-mer panel sees only part of a family, NGS-DOSE reads the
family low by about the same share in everyone, and the people lie on a line below equality. The panel's <em>recall</em> is the share of a
family's reads it can assign in CHM13, the genome it was built from{(": " + and_(relative) + " have a recall well below 1 and are relative measures, comparable between people but not in absolute megabases") if relative else ""}.
What tests the method is how closely each person sits on their family's line.</p>''')
        P.h(f'''<dl class="methods">
<dt>genomes</dt><dd>People compared in the family: counted here, with an HPRC release-2 assembly, and not left out.</dd>
<dt>left out (gaps)</dt><dd>An assembly does not always finish an array. Where it stops, the annotation marks a gap of unknown length
labelled with the family (for example GAP,HSat2), and the assembly's mass for that family is then only a lower bound. A person is left out of
a family's comparison when such gaps amount to more than {fmt(MAX_GAPPED, 0, pct=True)} of the family's annotated mass in their assembly. An
array the assembler collapsed or lost without leaving a gap cannot be seen this way.</dd>
<dt>median estimate / assembly</dt><dd>The typical ratio of NGS-DOSE's mass to the assembly's: near 1 for a family the panel sees fully,
lower for a relative one.</dd>
<dt>SD of log ratio, robust SD</dt><dd>How far individual people scatter around that ratio, roughly the per-person disagreement as a
fraction (0.05 is about 5%). The robust version, 1.4826 × the median absolute deviation, is not moved by a few outliers.</dd>
<dt>between-person CV, assembly</dt><dd>How much people differ in the family, by the assembly's measure: the variation the comparison has
to reproduce.</dd>
<dt>Pearson r, Spearman</dt><dd>How well the order and spacing of people is reproduced (Spearman: the order only). r can be high only
where people differ by much more than the two measurements disagree. <strong>r ≥ 0.95</strong> is the mark the summary uses to say that a
family <em>tracks</em> the assemblies: a round threshold set when six people could be compared, not a statistical test.{(" By it, " + and_(good) + " track the assemblies.") if good else " No family reaches it yet."}</dd>
</dl>''')
        P.table([[cls, s_["n"], s_.get("n_gapped", 0), fmt(s_.get("ratio_median"), 2), fmt(s_.get("sd_log"), 3), fmt(s_.get("sd_log_robust"), 3), fmt(s_.get("cv_assembly"), 1, pct=True),
                  fmt(s_.get("pearson"), 2), fmt(s_.get("spearman"), 2)] for cls, s_ in st.items() if s_.get("n")],
                ["class", "genomes", "left out (gaps)", "median estimate / assembly", "SD of log ratio", "robust SD of log ratio", "between-person CV, assembly",
                 "Pearson r", "Spearman"], numeric={1, 2, 3, 4, 5, 6, 7, 8})
        # what r can say depends on how much people differ against how far the two measurements disagree per genome; and the outliers' direction
        judged = {cls: s_ for cls, s_ in st.items() if s_.get("n", 0) >= 10 and s_.get("sd_log_robust") and s_.get("cv_assembly") is not None}
        narrow = [cls for cls, s_ in judged.items() if s_["sd_log_robust"] <= 0.08 and s_["sd_log_robust"] <= s_["cv_assembly"] < 2 * s_["sd_log_robust"]
                  and (s_.get("pearson") or 0) < 0.9]
        untestable = [cls for cls, s_ in judged.items() if s_["cv_assembly"] < s_["sd_log_robust"]]
        n_far = sum(s_.get("n_far", 0) for s_ in st.values())
        n_short = sum(s_.get("n_far_assembly_short", 0) for s_ in st.values())
        far_in = [cls for cls, s_ in st.items() if s_.get("n_far")]
        notes = []
        if narrow:
            notes.append("A modest r does not mean poor agreement where people barely differ: " + "; ".join(
                f"{cls} arrays differ between people by {fmt(st[cls]['cv_assembly'], 1, pct=True)} and the two measurements agree per genome to "
                f"{fmt(st[cls]['sd_log_robust'], 1, pct=True)} (robust SD), so with so little to separate people r stays at {fmt(st[cls].get('pearson'), 2)}"
                for cls in narrow) + ".")
        if untestable:
            notes.append(f"For {and_(untestable)} the two measurements disagree per genome by more than people differ, so this comparison cannot test "
                         f"{'that panel' if len(untestable) == 1 else 'those panels'}.")
        if n_far:
            notes.append(f"Of the {n_far} comparisons more than three robust SDs from their family's median (in {and_(far_in)}), {n_short} are genomes "
                         "whose assembly holds less than the reads show" + (", the direction expected where part of an array is missing from an assembly "
                                                                             "without a marked gap." if n_short > n_far / 2 else "."))
        if notes:
            P.h('<p class="small">' + " ".join(notes) + "</p>")
        # one clean scatter per family: equality, the family's median ratio, and the people far from it
        P.h('<div class="grid2">')
        for cls, s_ in st.items():
            k, rsd = s_.get("ratio_median"), s_.get("sd_log_robust") or 0.0
            if s_.get("n", 0) < 3 or not k:
                continue
            desc, recall = SAT_FAMILY.get(cls, (cls, 1.0))
            pts_c = [r for r in hp["rows"] if r["cls"] == cls and r["assembly_Mb"] > 0 and r["ngsdose_Mb"] > 0
                     and r["assembly_gapped_Mb"] <= MAX_GAPPED * (r["assembly_Mb"] + r["assembly_gapped_Mb"])]
            far = [rsd > 0 and abs(math.log(r["ngsdose_Mb"] / r["assembly_Mb"] / k)) > 3 * rsd for r in pts_c]
            n_out = sum(far)
            n_out_short = sum(1 for r, f in zip(pts_c, far) if f and r["ngsdose_Mb"] / r["assembly_Mb"] > k)
            spec = dict(type="scatter", points=[dict(x=r["assembly_Mb"], y=r["ngsdose_Mb"], label=r["sample"], si=int(f)) for r, f in zip(pts_c, far)],
                        xlabel="assembly, Mb (both haplotypes)", ylabel="NGS-DOSE, Mb", diagonal=True, slope_ref=dict(k=k, label=f"median ratio {fmt(k, 2)}"))
            if n_out:
                spec["legend"] = ["person", "more than 3 robust SDs from the median ratio"]
            # whether the equality line falls inside the plot, by the renderer's own padding (6% of the x range, 8% of the y range)
            xv, yv = [r["assembly_Mb"] for r in pts_c], [r["ngsdose_Mb"] for r in pts_c]
            px, py = (max(xv) - min(xv) or 1) * 0.06, (max(yv) - min(yv) or 1) * 0.08
            eq_shown = max(min(xv) - px, min(yv) - py) < min(max(xv) + px, max(yv) + py)
            cap = (f"Each dot is one of {s_['n']} people of the 1000 Genomes 30× cohort with an HPRC release-2 assembly"
                   + (f" ({s_['n_gapped']} more left out: their assembly leaves part of the array as a gap)" if s_.get("n_gapped") else "")
                   + f": x the {cls} array mass in the assembly (both haplotypes), y NGS-DOSE's estimate from the short reads. "
                   + ("Solid grey line: equality. " if eq_shown else "")
                   + f"Dashed line: the median ratio, {fmt(k, 2)}"
                   + (f" (the panel's recall for this family in CHM13 is {fmt(recall, 0, pct=True)})" if recall < 0.95 else "")
                   + ("" if eq_shown else "; equality lies off the plot")
                   + f". People differ by {fmt(s_.get('cv_assembly'), 1, pct=True)} (CV of the assembly mass); per person the two measurements agree to "
                     f"{fmt(rsd, 1, pct=True)} (robust SD of the log ratio); Pearson r = {fmt(s_.get('pearson'), 2)}."
                   + ((f" Orange: {n_out} {'person' if n_out == 1 else 'people'} more than three robust SDs from the median ratio"
                       + (", all with less in the assembly than in the reads." if n_out_short == n_out else ".")) if n_out else "")
                   + SAT_NOTE.get(cls, ""))
            P.chart(f"sat_{cls}", spec, f"{cls}, {desc}: NGS-DOSE against the assemblies", cap)
        P.h("</div>")
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
        adj = pcs.get("adjusted")
        if var:
            P.chart("scree", dict(type="lines", series=[dict(name="variance explained", x=list(range(1, len(var) + 1)), y=var)], xlabel="component", ylabel="fraction of variance"),
                    "Technical structure: variance explained by the control regions' principal components",
                    f"Principal components of the residual depth of 800 single-copy control regions after each library's GC model, across {n:,} "
                    f"genomes of the 1000 Genomes 30× cohort: the share of that variance each component explains."
                    + (f" The first {adj['k']} clear the Marchenko–Pastur edge of the noise and are the technical covariates tested in this section." if adj else ""))
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
            names = {"truth.auto": "held-out autosomal (2 copies)", "truth.chrX": "chrX", "truth.chrY": "chrY", kt["DJ_col"]: "distal junction"}
            for col in ("truth.auto", "truth.chrX", "truth.chrY", kt["DJ_col"]):
                rr = [r for r in sw["rows"] if r["column"] == col and "sd_log_robust" in r]
                if rr:
                    base = rr[0]["sd_log_robust"] or 1
                    series.append(dict(name=names[col], x=[r["n_pc"] for r in rr], y=[r["sd_log_robust"] / base for r in rr]))
            if series:
                P.chart("sweep_truth", dict(type="lines", series=series[:3], xlabel="control PCs regressed out", ylabel="error relative to none", ref=1),
                        "Known truths: measurement error against the number of control PCs regressed out",
                        f"For the sequence of known copy number, in {n:,} genomes of the 1000 Genomes 30× cohort: the cross-validated robust SD of "
                        f"log(estimate / expected copies) after regressing out 0 to {sw['max_pc']} of the control regions' principal components, relative to "
                        f"regressing out none. Below the line at 1 the components remove measurement error."
                        + (" Cross-validated one-standard-error picks: " + ", ".join(f"{names[c]} {rec[c]['pick']}" for c in names if c in rec) + " components." if any(c in rec for c in names) else ""))
            series = []
            names = {"rDNA45S.cn": "45S, NGS-DOSE calibrated", "rDNA45S.18S.flat": "45S, 18S depth ratio (published)", "rDNA45S.cn_single": "45S, NGS-DOSE single-sample"}
            for col in ("rDNA45S.cn", "rDNA45S.18S.flat", "rDNA45S.cn_single"):         # colours as everywhere: NGS-DOSE blue, the published ratio orange
                rr = [r for r in sw["rows"] if r["column"] == col and "R_midparent" in r]
                if rr:
                    series.append(dict(name=names[col], x=[r["n_pc"] for r in rr], y=[r["R_midparent"] for r in rr], lo=[r.get("R_lo", r["R_midparent"]) for r in rr], hi=[r.get("R_hi", r["R_midparent"]) for r in rr]))
            if series:
                near1 = min(min(s["y"]) for s in series) >= 0.85
                P.chart("sweep_R", dict(type="lines", series=series, xlabel="control PCs regressed out", ylabel="transmission reliability"),
                        "Transmission reliability against the number of control PCs regressed out",
                        f"Reliability R in {sw['n_trios']} trios of the 1000 Genomes 30× cohort for NGS-DOSE's two 45S estimates and the 18S depth ratio "
                        f"of published studies computed from the same reads, after regressing out 0 to {sw['max_pc']} of the control regions' principal "
                        f"components. Bands: family-bootstrap 95% intervals."
                        + (" People differ so much in 45S copy number that every estimator reads near 1 here: within one chemistry the trios do not rank "
                           "the estimators (3.6); the two technologies do (3.4)." if near1 else ""))
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
    c5 = rd[col5]
    dist = lambda d, nd=0: f"Median {fmt(d.get('median'), nd)}; 10–90% {fmt(d.get('q10'), nd)}–{fmt(d.get('q90'), nd)}; range {fmt(d.get('min'), nd)}–{fmt(d.get('max'), nd)}."
    P.chart("cn45", dict(type="hist", col=col45, xlabel="45S copies per diploid genome", xfmt=0), "45S rDNA copy number, NGS-DOSE",
            f"Copies of the 45S rDNA unit per diploid genome (NGS-DOSE, {'calibrated' if col45 == 'rDNA45S.cn' else 'single-sample'} estimate) in each of "
            f"{c45.get('n', 0):,} genomes of the 1000 Genomes 30× cohort, from lymphoblastoid-cell-line DNA. Bars count genomes. {dist(c45)}")
    P.chart("cn5", dict(type="hist", col=col5, xlabel="5S copies per diploid genome", xfmt=0), "5S rDNA copy number, NGS-DOSE",
            f"Copies of the 5S rDNA unit per diploid genome (NGS-DOSE, {'calibrated' if col5 == 'rDNA5S.cn' else 'single-sample'} estimate) in each of "
            f"{c5.get('n', 0):,} genomes of the 1000 Genomes 30× cohort. Bars count genomes. {dist(c5)}")
    P.h("</div>")
    if superpop_order:
        P.chart("pop", dict(type="strip", col=col45, by="superpop", order=superpop_order, labels=SUPERPOP_NAMES, ylabel="45S copies"), "45S rDNA copy number by super-population",
                "Each dot is one genome of the 1000 Genomes 30× cohort (NGS-DOSE, 45S copies per diploid genome); the bar is the median. Counted so far: "
                + "; ".join(f"{esc(g['group'])} {g['n']:,}" + (f" of {m['by_superpop'][g['group']]['total']:,}" if g["group"] in (m.get("by_superpop") or {}) else "")
                            for g in rd["by_superpop"])
                + (" (" + " and ".join(c for c in SUPERPOPS if c in (m.get("by_superpop") or {}) and c not in {g["group"] for g in rd["by_superpop"]}) + " not yet)"
                   if any(c in (m.get("by_superpop") or {}) and c not in {g["group"] for g in rd["by_superpop"]} for c in SUPERPOPS) else "")
                + ". Whether differences between populations survive adjustment for technical structure is a question for the complete cohort.")
        if len(rd["by_pop"]) > 1:
            P.h("<details><summary>By population</summary>")
            P.table([[g["group"], g["n"], fmt(g["median"], 0), fmt(g.get("q10"), 0) + "–" + fmt(g.get("q90"), 0)] for g in rd["by_pop"]], ["population", "n", "median 45S", "10–90%"], numeric={1, 2, 3})
            P.h("</details>")
    cx = bio["cn45_vs_5S"]
    if cx.get("n", 0) >= 3:
        P.chart("c45v5", dict(type="scatter", x=col45, y=col5, xlabel="45S copies", ylabel="5S copies", fit=True),
                "45S against 5S rDNA copy number",
                f"Each dot is one of {cx['n']:,} genomes of the 1000 Genomes 30× cohort: 45S (x) and 5S (y) copies per diploid genome, both NGS-DOSE; line: "
                f"least-squares fit. Pearson r = {ci(cx)}, Spearman {fmt(cx.get('spearman'), 2)}. Gibbons et al. (2015) reported the two arrays' copy numbers "
                f"to be correlated; Hall et al. (2021) did not see it in these genomes. Here both are measured the same way, with no shared denominator.")
    P.h("<h3>The cell line</h3>")
    P.h(f'''<p>Every genome was sequenced from a lymphoblastoid cell line, whose state leaves marks on coverage. Mitochondrial genomes per cell: median
<strong>{fmt(bio["chrM"].get("median"), 0)}</strong> (10–90%: {fmt(bio["chrM"].get("q10"), 0)}–{fmt(bio["chrM"].get("q90"), 0)}); EBV episomes per cell: median
<strong>{fmt(bio["chrEBV"].get("median"), 1)}</strong> ({fmt(bio["chrEBV"].get("q10"), 1)}–{fmt(bio["chrEBV"].get("q90"), 1)}). Both vary far more between people
than the rDNA does and neither is inherited through the nuclear genome, which is why they serve as negative controls in 3.6.</p>''')
    P.h('<div class="grid2">')
    P.chart("chrM", dict(type="hist", col="chrM.copies", xlabel="mitochondrial genomes per cell", xfmt=0), "Mitochondrial genomes per cell",
            f"Mitochondrial genomes per cell in each of {bio['chrM'].get('n', 0):,} genomes of the 1000 Genomes 30× cohort (lymphoblastoid cell lines), "
            f"NGS-DOSE. Bars count genomes. {dist(bio['chrM'])} A property of the culture, not inherited through the nuclear genome: a negative "
            f"control for the trio test (3.6).")
    P.chart("ebv", dict(type="hist", col="chrEBV.copies", xlabel="EBV episomes per cell", xfmt=0), "Epstein–Barr virus episomes per cell",
            f"EBV genomes per cell in each of {bio['chrEBV'].get('n', 0):,} genomes of the 1000 Genomes 30× cohort, NGS-DOSE; the cell lines were "
            f"immortalised with EBV. Bars count genomes. {dist(bio['chrEBV'], 1)} A property of the culture, not inherited: a negative control for "
            f"the trio test (3.6).")
    P.h("</div>")
    cm, sph = bio["log_cn45_vs_log_chrM"], bio["DJ_vs_chrX_female"]
    P.h('<div class="grid2">')
    if cm.get("n", 0) >= 3:
        P.chart("c45vM", dict(type="scatter", x="chrM.copies", y=col45, xlabel="mitochondrial genomes per cell", ylabel="45S copies", log="xy", fit=True),
                "45S copy number against mitochondrial content",
                f"Each dot is one of {cm['n']:,} genomes of the 1000 Genomes 30× cohort: mitochondrial genomes per cell (x) and 45S copies per diploid "
                f"genome (y), both NGS-DOSE; log scales; line: least-squares fit. r = {ci(cm)}, on the log scale. Gibbons et al. (2014) reported rDNA copy "
                f"number to be coupled with mitochondrial DNA abundance in lymphoblastoid lines.")
    if sph.get("n", 0) >= 3:
        P.chart("sphase", dict(type="scatter", x="truth.chrX", y=kt["DJ_col"], where=dict(sex_inferred="F"), xlabel="chrX copies (women)", ylabel="distal junction copies", fit=False, xref=2, yref=10),
                "Two late-replicating controls, in women",
                f"Each dot is one of {sph['n']:,} women of the 1000 Genomes 30× cohort whose cell line has kept both X chromosomes (chrX 1.85–2.15): chrX "
                f"copies (x; expected 2) and distal-junction copies (y; expected 10), both NGS-DOSE; lines at the expected values. r = {ci(sph)}. The "
                f"inactive X and the acrocentric short arms replicate late; a culture with more cells in S phase should under-represent both together. "
                f"An r near zero says the two deficits are not one thing.")
    P.h("</div>")
    du = bio["dup"]
    P.chart("dup", dict(type="scatter", x="ctrl_dup_frac", y="rDNA45S.dup_flag_frac", xlabel="duplicate-flagged, control reads", ylabel="duplicate-flagged, 45S reads", identity=True, xfmt=3),
            "The duplicate flag inside and outside the rDNA",
            f"Each dot is one of {du['ratio'].get('n', 0):,} genomes of the 1000 Genomes 30× cohort: the share of reads the sequencing centre's pipeline "
            f"flagged as duplicates among single-copy control reads (x) and among 45S rDNA reads (y); diagonal: equal rates. Median "
            f"{fmt(du['control'].get('median'), 1, pct=True)} of control reads against {fmt(du['rDNA'].get('median'), 1, pct=True)} of 45S reads (ratio "
            f"{fmt(du['ratio'].get('median'), 2)}, range {fmt(du['ratio'].get('min'), 2)}–{fmt(du['ratio'].get('max'), 2)})"
            + (": inside a collapsed array duplicates escape the marker. A pipeline that drops flagged reads therefore reads the rDNA high, by a "
               "different amount in every genome; NGS-DOSE counts all primary reads." if (du["ratio"].get("median") or 1) < 1 else
               ". A pipeline that drops flagged reads reads the rDNA by a different amount in every genome; NGS-DOSE counts all primary reads."))
    if sat["classes"]:
        P.h("<h3>Satellite arrays and the telomeric repeat</h3>")
        P.h('''<p>Satellite arrays are dispersed over the alignment, so only a whole-file scan measures them. Diploid mass per family; what each
panel can see was measured on the genome it was built from (four families are relative measures, under-read by their k-mer recall).</p>''')
        rows_s = [[cls, d["n"], fmt(d["median"], 1), fmt(d.get("q10"), 1) + "–" + fmt(d.get("q90"), 1), fmt(d["by_sex"]["M"].get("median"), 1), fmt(d["by_sex"]["F"].get("median"), 1)] for cls, d in sat["classes"].items()]
        P.table(rows_s, ["class", "n", "median Mb (diploid)", "10–90%", "men", "women"], numeric={1, 2, 3, 4, 5})
        P.h('<div class="grid2">')
        sc = sat["classes"]
        sat_cap = lambda cls, what, tail="", nd=1: (f"{what} per diploid genome, NGS-DOSE from whole-file scans, in each of {(sc.get(cls) or {}).get('n', 0):,} genomes of "
                                                    f"the 1000 Genomes 30× cohort. Bars count genomes. {dist(sc.get(cls) or {}, nd)}{tail}")
        P.chart("hsat3", dict(type="hist", col="HSat3.mass_Mb", xlabel="Mb per diploid genome", xfmt=0), "HSat3 satellite: array mass",
                sat_cap("HSat3", "Megabases of HSat3 satellite arrays", " Compared with long-read assemblies in 3.8."))
        m1b = lambda s: fmt(((sc.get("HSat1B") or {}).get("by_sex", {}).get(s) or {}).get("median"), 1)
        P.chart("hsat1b", dict(type="hist", col="HSat1B.mass_Mb", group=dict(col="sex_inferred", levels=sex_levels), xlabel="Mb per diploid genome", xfmt=1),
                "HSat1B satellite, by sex: mostly on the Y chromosome",
                sat_cap("HSat1B", "Megabases of HSat1B satellite arrays", f" Men blue, women orange: median {m1b('M')} Mb in men and {m1b('F')} Mb in women, "
                                                                          "since most HSat1B lies on the long arm of chrY."))
        P.chart("ahor", dict(type="hist", col="aSatHOR.mass_Mb", xlabel="Mb per diploid genome", xfmt=0), "α-satellite higher-order repeats: array mass",
                sat_cap("aSatHOR", "Megabases of centromeric α-satellite higher-order-repeat arrays", " Compared with long-read assemblies in 3.8."))
        P.chart("tel", dict(type="hist", col="TEL.mass_Mb", xlabel="Mb of (TTAGGG)n-bearing reads, diploid", xfmt=2), "Telomeric repeat: a relative measure",
                sat_cap("TEL", "Megabases of reads carrying the telomeric repeat (TTAGGG)n",
                        " A relative measure of telomeric content, not a telomere length: exact 31-mers under-count error-bearing and variant repeats.", 2))
        P.h("</div>")
    P.end()

    # ---------------------------------------------------------------- 5. limitations
    h2 = hpst.get("HSat2") or {}
    hsat2_limit = (f"HSat2, compared in the {h2['n']} assemblies that close its arrays, gives r = {fmt(h2.get('pearson'), 2)} (SD of the log ratio "
                   f"{fmt(h2.get('sd_log'), 2)}){' and is not yet a usable measure' if (h2.get('pearson') or 0) < 0.8 else ''}." if h2.get("n", 0) >= 10
                   else "HSat2 is unjudged until assemblies without gaps in it have been compared.")
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
<li><strong>The satellite panels</strong> were built from one genome (CHM13). Four of the ten families are relative measures; {hsat2_limit}
The telomere class is a relative measure of (TTAGGG)n content, not a telomere length.</li>
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
<code>data/transmission.tsv</code>, <code>data/transmission_by_sex.tsv</code> and <code>data/trios.tsv</code> (every trio's values), <code>data/fetch_check.tsv</code>, <code>data/pcsweep.tsv</code>,
<code>data/satellites_hprc.tsv</code>,{" <code>data/rdna_hprc.tsv</code>," if rd.get("assemblies") else ""} <code>data/flags.tsv</code>; the numbers
in the prose, <code>report.json</code>; the trio assessment as a document, <code>trio_report.pdf</code>. The counts were made by <code>ngs-dose count</code> ({eng}) from the 1000 Genomes 30× CRAMs
(Byrska-Bishop et al., <em>Cell</em> 2022; AWS Open Data) with the {esc(m.get("bundle"))} resource bundle. Method, design document and
audit: <a href="https://github.com/jlanej/NGS-DOSE">github.com/jlanej/NGS-DOSE</a>.</p>''')
    P.end()

    # ---------------------------------------------------------------- 7. every sample
    P.section("samples", "7. Every genome", "Samples")
    cols = [("sample", "sample"), ("pop", "pop"), ("sex", "sex (ped)"), ("sex_inferred", "sex (reads)"), ("depth", "depth"), ("truth.auto", "auto"), ("truth.chrX", "chrX"),
            ("truth.chrY", "chrY"), (kt["DJ_col"], "DJ"), (col45, "45S"), ("rDNA45S.cn_single", "45S single"), ("rDNA45S.18S.flat", "18S ratio (published)"),
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

    # ---------------------------------------------------------------- supplementary: every trio metric, child against midparent
    pts_all = tr.get("points") or {}
    if pts_all and tr["n_complete"] >= 10:
        from .report import TRIO_GROUPS
        by_col = {t["column"]: t for t in tr["table"]}
        bs = (tr.get("by_sex") or {}).get("table") or {}
        P.section("supp", "Supplementary figures: every trio metric, child against midparent", "Supplementary")
        P.h(f"""<p>One figure for each metric of the trio test (<a href="#trios">3.6</a>), in its groups: x the mean of a child's two parents, y the child,
on the natural scale; sons blue, daughters orange, each with its least-squares line, and the grey diagonal where child equals midparent. The
statistics in each caption are those of 3.6, on values centred within population (and, for the split by sex, within sex too). {tr['n_complete']:,}
complete trios.</p>""")
        note = {"rDNA": "", "satellites": " Genomic, so it must be inherited: a positive control.",
                "truth": " Known copy number: nothing to inherit but the distal junction's whole-copy steps.",
                "culture": " Not in the nuclear genome, so a negative control: expected near 0."}
        for key, gl, cols in TRIO_GROUPS:
            present = [(c, l) for c, l in cols if c in pts_all and c in by_col]
            if not present:
                continue
            P.h(f"<h3>{esc(gl[0].upper() + gl[1:])}</h3>")
            P.h('<div class="grid2">')
            for col, label in present:
                t, unit, d = by_col[col], unit_of(col), bs.get(col) or {}
                pts = [dict(x=p[1], y=p[2], label=p[0], si=0 if p[5] == "M" else 1, extra=[f"father {fmt(p[3], 2)}, mother {fmt(p[4], 2)}"])
                       for p in pts_all[col] if p[5] in ("M", "F")]
                split = ", ".join(f"{name} {fmt(d[k]['r'], 2)}" for k, name in (("father_son", "father → son"), ("father_daughter", "father → daughter"),
                                                                                      ("mother_son", "mother → son"), ("mother_daughter", "mother → daughter")) if k in d)
                cap = (f"Each dot is one of {len(pts)} complete parent–child trios of the 1000 Genomes 30× cohort: x the mean of the two parents' {label}, "
                       f"y the child's, in {unit}; sons blue, daughters orange, each with its least-squares line; grey diagonal: child equals midparent. "
                       f"Midparent slope {fmt(t['slope'], 2)} ± {fmt(t['slope_se'], 2)}, spousal r = {fmt(t['spousal_r'], 2)}, reliability "
                       f"R = {fmt(min(t['R'], 1.0), 2)}"
                       + (" (" + "; ".join((["capped at 1"] if t["R"] > 1 else []) + ([f"95% CI {fmt(t['R_lo'], 2)} to {fmt(t['R_hi'], 2)}"] if "R_lo" in t else [])) + ")"
                          if t["R"] > 1 or "R_lo" in t else "") + "."
                       + (f" Correlation of child with parent by sex (values centred within population and sex): {split}." if split else "") + note.get(key, ""))
                P.chart(f"trio_{re.sub(r'[^A-Za-z0-9]', '_', col)}", dict(type="scatter", points=pts, legend=["son", "daughter"], xlabel=f"midparent, {unit}",
                                                                          ylabel=f"child, {unit}", identity=True, fit=True, fit_labels=["sons", "daughters"]),
                        f"{sentence_case(label)}: child against midparent", cap, supp=True)
            P.h("</div>")
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
