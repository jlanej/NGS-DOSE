"""`ngsdose report`: one static, self-contained page from whatever counts files exist.

The page states the rationale and the method, then shows - computed afresh from the counts
files it was given - the validation: sequence of known copy number in every sample, the
distal junction's whole-copy steps and their transmission, the targeted fetch against the
whole-file scan, the same people on two sequencing technologies (the pilot), what the GC model
removes, transmission through trios, an independent pipeline on the same files, and satellite
arrays against assemblies; then the descriptive results and the limitations. It is meant to be regenerated
as a cohort run proceeds (a sample's estimate is cached and reused if its counts file has not
changed) and published as it stands, partial or complete.

    ngsdose report --scan counts_scan/ --fetch counts_fetch/ -p pedigree.txt -o docs/

Everything the page shows is also written as TSV under `<out>/data/`, and the numbers behind the
prose as `<out>/report.json`.
"""
from __future__ import annotations

import datetime as dt
import glob
import html
import json
import math
import os
import re
import sys
from pathlib import Path

import numpy as np

from . import __version__, cohort, estimate, hprc, io, pcselect, resources, sinks
from . import trios as T
from .tables import dump, load_result, num, summary_row, write_table

REPORT_VERSION = 4
POSITIONAL = ("rDNA45S", "rDNA5S", "DJ")
# the classes fetch mode retrieves through the sinks, and the column their two modes are compared on
FETCHABLE = [("rDNA45S", "rDNA45S.cn_single"), ("rDNA5S", "rDNA5S.cn_single"), ("DJ", "DJ.cn_single"), ("TEL", "TEL.mass_Mb")]
# the reported estimate of each fetchable class, plotted per genome against the scan once the whole pipeline has run on each mode
FETCH_SCATTER = ("rDNA45S.cn", "rDNA5S.cn", "DJ.cn", "TEL.mass_Mb")
SATELLITES = hprc.CLASSES + ("TEL",)
# NGS-DOSE's estimators, and the 18S depth ratio of published studies computed from the same reads as the comparator (not NGS-DOSE's)
ESTIMATORS = [("rDNA45S.cn", "45S, NGS-DOSE calibrated"), ("rDNA45S.cn_single", "45S, NGS-DOSE single-sample"),
              ("rDNA45S.18S.flat", "45S, 18S depth ratio (published)"), ("rDNA5S.cn", "5S, NGS-DOSE calibrated"), ("DJ.cn", "distal junction, NGS-DOSE calibrated")]
NEGATIVE_CONTROLS = [("truth.auto", "held-out autosomal (no variance but error)"), ("chrM.copies", "mitochondrial genomes per cell (culture)"),
                     ("chrEBV.copies", "EBV episomes per cell (culture)")]
# every metric the trios are asked about, in the groups a sceptic would want to see side by side: what is claimed to be
# inherited, what must be inherited (the satellite arrays are genomic), what has no variation to inherit, and what is
# not in the nuclear genome at all
TRIO_GROUPS = [
    ("rDNA", "rDNA copy number: NGS-DOSE (the claim), and the published 18S ratio for comparison", [e for e in ESTIMATORS if e[0] != "DJ.cn"]),
    ("satellites", "satellite arrays, mass (genomic: positive controls)", [(f"{c}.mass_Mb", f"{c} array") for c in hprc.CLASSES] + [("TEL.mass_Mb", "telomeric repeat")]),
    ("truth", "known copy number (nothing to inherit but the distal junction's whole-copy steps)", [("truth.auto", "held-out autosomal (2)"), ("DJ.cn", "distal junction (10)")]),
    ("culture", "culture and library (not in the nuclear genome)", [("chrM.copies", "mitochondrial genomes per cell"), ("chrEBV.copies", "EBV episomes per cell"),
                                                                    ("depth", "sequencing depth"), ("ctrl_dup_frac", "duplicate-flagged fraction"),
                                                                    ("gc_rel_65", "library GC bias"), ("insert_median", "insert size")]),
]
TRIO_COLUMNS = [(col, label) for _, _, cols in TRIO_GROUPS for col, label in cols]
TRIO_GROUP_OF = {col: key for key, _, cols in TRIO_GROUPS for col, _ in cols}
ASSETS = Path(__file__).parent / "report_assets"


# ------------------------------------------------------------------------------------------------
# per-sample estimates, cached
# ------------------------------------------------------------------------------------------------
_E: dict = {}


def _init(bundle_dir):
    res = resources.Bundle(bundle_dir)
    bed = [(p[0], int(p[1]), int(p[2]), p[3].strip()) for p in (line.split("\t") for line in open(res.sinks) if line.strip())]
    _E.update(res=res, panel=io.load_panel(res.panel), units=res.units(), feats=res.features(), anchors=res.anchors(),
              lengths=res.contig_lengths(), regions=res.regions(), sinks=bed, tables={})


def _flat_summary(r: dict) -> dict:
    row = summary_row(r)
    x = r.get("report_extras", {})
    for k in ("elapsed_sec", "primary", "records", "placement_bin", "input", "panel_sha256", "controls_sha256", "sinks_sha256", "n_classes"):
        row[k] = x.get(k)
    for cls, v in x.get("capture", {}).items():
        row[f"capture.{cls}"] = v
    for c, v in x.get("flat_copies", {}).items():
        row[f"flat.{c}"] = v
    return row


def _one(job):
    path, out = job
    if out.exists() and os.path.getmtime(out) >= os.path.getmtime(path):
        r = load_result(out)
        if r.get("report_extras", {}).get("report_version") == REPORT_VERSION:
            return _flat_summary(r)
    counts = io.load_counts(path)
    L = estimate.nearest_table(counts, None)["l"]
    if L not in _E["tables"]:
        _E["tables"][L] = _E["res"].region_tables(L)
    r = estimate.estimate_sample(counts, _E["panel"], _E["units"], _E["feats"], region_tables=_E["tables"][L], anchors=_E["anchors"],
                                 contig_lengths=_E["lengths"], regions=_E["regions"])
    x = dict(report_version=REPORT_VERSION, elapsed_sec=counts.get("elapsed_sec"), primary=counts.get("primary"), records=counts.get("records"),
             placement_bin=counts.get("placement_bin"), input=counts.get("input"), n_classes=len(counts["classes"]),
             panel_sha256=",".join(h[:8] for h in counts.get("panel_sha256", [])), controls_sha256=(counts.get("controls_sha256") or "")[:8],
             sinks_sha256=(counts.get("sinks_sha256") or "")[:8])
    if counts["mode"] == "scan":
        sunk = {b[3] for b in _E["sinks"]}
        x["capture"] = {cls: round(inside / max(tot, 1), 6) for cls, (tot, inside) in sinks.capture(counts, _E["sinks"]).items() if cls in sunk}
        rate = r["ctrl_rate"]
        x["flat_copies"] = {c["name"]: round(c["reads"] / (c["len"] * rate), 3) for c in counts.get("contigs", [])
                            if c.get("reads") is not None and c["name"] in ("chrM", "chrEBV", "chrX", "chrY", "chr1") and rate > 0}
    r["report_extras"] = x
    out.parent.mkdir(parents=True, exist_ok=True)
    dump(r, out)
    return _flat_summary(r)


def estimate_all(paths: list[Path], cache: Path, bundle_dir, jobs: int, log) -> list[dict]:
    jobs_list = [(p, cache / (p.name.split(".")[0] + ".estimate.json.gz")) for p in paths]
    if not jobs_list:
        return []
    log(f"[report] {len(jobs_list)} counts files in {paths[0].parent} ({sum(1 for p, o in jobs_list if o.exists())} cached)")
    if jobs > 1 and len(jobs_list) > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(jobs, initializer=_init, initargs=(bundle_dir,)) as pool:
            return pool.map(_one, jobs_list, chunksize=8)
    _init(bundle_dir)
    return [_one(j) for j in jobs_list]


# ------------------------------------------------------------------------------------------------
# small statistics
# ------------------------------------------------------------------------------------------------
def fin(v) -> np.ndarray:
    a = np.asarray(v, float)
    return a[np.isfinite(a)]


def describe(v) -> dict:
    a = fin(v)
    if len(a) == 0:
        return dict(n=0)
    d = dict(n=int(len(a)), mean=float(a.mean()), median=float(np.median(a)), min=float(a.min()), max=float(a.max()))
    if len(a) > 1:
        d["sd"] = float(a.std(ddof=1))
        d["mad_sd"] = float(1.4826 * np.median(np.abs(a - np.median(a))))
        d["q10"], d["q90"] = float(np.percentile(a, 10)), float(np.percentile(a, 90))
    return d


def corr(x, y) -> dict:
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3 or x.std() == 0 or y.std() == 0:
        return dict(n=int(len(x)))
    rank = lambda v: np.argsort(np.argsort(v))
    r = float(np.corrcoef(x, y)[0, 1])
    # Fisher interval on r
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    se = 1 / math.sqrt(max(len(x) - 3, 1))
    return dict(n=int(len(x)), r=r, r_lo=float(np.tanh(z - 1.96 * se)), r_hi=float(np.tanh(z + 1.96 * se)),
                spearman=float(np.corrcoef(rank(x), rank(y))[0, 1]), slope=float(np.polyfit(x, y, 1)[0]))


def rnd(o, nd=4):
    """Round floats throughout a nested structure (a deterministic page is a diffable page)."""
    if isinstance(o, float):
        return None if not math.isfinite(o) else round(o, nd)
    if isinstance(o, dict):
        return {k: rnd(v, nd) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [rnd(v, nd) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return rnd(o.item(), nd)
    if isinstance(o, np.ndarray):
        return rnd(o.tolist(), nd)
    return o


# ------------------------------------------------------------------------------------------------
# pedigree, sex, flags
# ------------------------------------------------------------------------------------------------
def load_pedigree(path):
    """1000 Genomes pedigree: sample -> dict(sex, pop, superpop, father, mother) and the complete trios."""
    info = {}
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 5 or p[1] in ("SampleID", "IID", "sampleID"):
                continue
            info[p[1]] = dict(family=p[0], father=p[2], mother=p[3], sex={"1": "M", "2": "F"}.get(p[4], ""), pop=p[5] if len(p) > 5 else "",
                              superpop=p[6] if len(p) > 6 else "")
    trio_list, population = T.load_pedigree(path)
    return info, trio_list, population


def inferred_sex(row) -> str:
    """A Y is a Y even when most of a culture has lost it (men read as low as 0.3 here); a woman
    has none (the female maximum in this cohort is 0.004)."""
    y, x = num(row, "truth.chrY"), num(row, "truth.chrX")
    if np.isfinite(y):
        return "M" if y > 0.1 else "F"
    if np.isfinite(x):
        return "M" if x < 1.5 else "F"
    return ""


def flags_for(row: dict, majority_engine: str | None) -> list[str]:
    f = []
    if row.get("eof_marker") not in (None, "present"):
        f.append("truncated input")
    if row.get("flagged_chromosomes"):
        f.append("aneuploid: " + row["flagged_chromosomes"])
    d = num(row, "depth")
    if np.isfinite(d) and d < 15:
        f.append(f"depth {d:.1f}")
    if abs(num(row, "truth.auto") - 2) > 0.06:
        f.append(f"held-out autosomal {num(row, 'truth.auto'):.3f}")
    sx, si = row.get("sex", ""), row.get("sex_inferred", "")
    if sx and si and sx != si:
        f.append(f"sex: pedigree {sx}, reads {si}")
    x, y = num(row, "truth.chrX"), num(row, "truth.chrY")
    if si == "F" and np.isfinite(x) and x < 1.85:
        f.append(f"chrX {x:.2f} (mosaic X loss?)")
    if si == "M" and np.isfinite(y) and y < 0.85:
        f.append(f"chrY {y:.2f} (mosaic Y loss?)")
    if si == "M" and np.isfinite(x) and x > 1.5:
        f.append(f"chrX {x:.2f} with a Y")
    step = row.get("DJ.step")
    if step is not None and abs(step) >= 0.75:
        f.append(f"DJ {step:+.2f} copies")
    if majority_engine and row.get("engine") and row["engine"] != majority_engine:
        f.append(f"engine {row['engine']}")
    if num(row, "gc_curve_max_se") > 0.5:
        f.append("GC curve poorly determined")
    return f


# ------------------------------------------------------------------------------------------------
# the analyses
# ------------------------------------------------------------------------------------------------
def known_truth(rows) -> dict:
    out = {}
    out["auto"] = describe([num(r, "truth.auto") for r in rows])
    for lab, col in (("chrX", "truth.chrX"), ("chrY", "truth.chrY")):
        out[lab] = {sex: describe([num(r, col) for r in rows if r.get("sex_inferred") == sex]) for sex in ("M", "F")}
    dj_col = "DJ.cn" if any(np.isfinite(num(r, "DJ.cn")) for r in rows) else "DJ.cn_single"
    out["DJ"] = describe([num(r, dj_col) for r in rows])
    out["DJ_col"] = dj_col
    out["sex"] = dict(n_pedigree=sum(1 for r in rows if r.get("sex")), n_inferred=sum(1 for r in rows if r.get("sex_inferred")),
                      mismatch=[r["sample"] for r in rows if r.get("sex") and r.get("sex_inferred") and r["sex"] != r["sex_inferred"]])
    xm = fin([num(r, "truth.chrX") for r in rows if r.get("sex_inferred") == "M"])
    xf = fin([num(r, "truth.chrX") for r in rows if r.get("sex_inferred") == "F"])
    if len(xm) and len(xf):
        out["sex"]["chrX_men_max"], out["sex"]["chrX_women_min"] = float(xm.max()), float(xf.min())
        out["sex"]["women_intact"] = describe(xf[xf >= 1.85])
        ym = fin([num(r, "truth.chrY") for r in rows if r.get("sex_inferred") == "M"])
        out["sex"]["men_intact_Y"] = describe(ym[ym >= 0.85])
        out["sex"]["n_mosaic_X"], out["sex"]["n_mosaic_Y"] = int((xf < 1.85).sum()), int((ym < 0.85).sum())
    return out


def dj_steps(rows, ped) -> dict:
    """The distal junction has ten copies, one per acrocentric short arm, in nearly everyone; a
    person with a rearranged short arm has nine or eight (a Robertsonian translocation loses two).
    Copy number relative to the cohort's median is therefore near an integer, and a step, being a
    structural variant, should be transmitted to half of a carrier's children and appear de novo
    in none - which the trios can check."""
    col = "DJ.cn" if any(np.isfinite(num(r, "DJ.cn")) for r in rows) else "DJ.cn_single"
    v = np.array([num(r, col) for r in rows])
    if not np.isfinite(v).any():
        return {}
    med = float(np.nanmedian(v))
    by = {}
    for r, x in zip(rows, v):
        r["DJ.step"] = None if not np.isfinite(x) else round(float(x - med), 3)
        by[r["sample"]] = r
    steps = v - med
    near = {k: int(np.sum(np.abs(steps - k) < 0.3)) for k in (-2, -1, 0, 1, 2)}
    between = int(np.sum(np.isfinite(steps) & (np.abs(steps - np.round(steps)) > 0.35)))
    main = steps[np.abs(steps) < 0.5]
    spread = float(1.4826 * np.median(np.abs(main - np.median(main)))) if len(main) > 2 else float("nan")
    carriers = []
    for r in rows:
        st = r.get("DJ.step")
        if st is None or abs(st) < 0.75:
            continue
        p = ped.get(r["sample"], {})
        rel = []
        for who in ("father", "mother"):
            o = p.get(who, "0")
            if o in by and by[o].get("DJ.step") is not None:
                rel.append(dict(who=who, sample=o, step=by[o]["DJ.step"]))
        for o, q in ped.items():
            if o in by and (q.get("father") == r["sample"] or q.get("mother") == r["sample"]) and by[o].get("DJ.step") is not None:
                rel.append(dict(who="child", sample=o, step=by[o]["DJ.step"]))
        carriers.append(dict(sample=r["sample"], pop=r.get("pop"), sex=r.get("sex_inferred"), step=st, relatives=rel))
    # transmission: a carrier parent and a counted child (the other parent need not be counted); a child's
    # step within half a copy of the parent's is the step transmitted, a child near zero is the step not
    # transmitted. De novo: a child with a step when both counted parents have none.
    transmitted = not_transmitted = 0
    for c in carriers:
        for rel in c["relatives"]:
            if rel["who"] != "child":
                continue
            if abs(rel["step"] - c["step"]) < 0.5:
                transmitted += 1
            elif abs(rel["step"]) < 0.5:
                not_transmitted += 1
    de_novo = []
    for r in rows:
        st = r.get("DJ.step")
        p = ped.get(r["sample"], {})
        f, m = by.get(p.get("father", "0")), by.get(p.get("mother", "0"))
        if st is None or f is None or m is None or f.get("DJ.step") is None or m.get("DJ.step") is None:
            continue
        if abs(st) >= 0.75 and abs(f["DJ.step"]) < 0.5 and abs(m["DJ.step"]) < 0.5:
            de_novo.append(r["sample"])
    # the satellite families of the acrocentric short arms, in the carriers, relative to the cohort: a
    # lost arm takes its satellites with it, and the pan-centromeric alpha satellite stays
    arm = ("ACRO", "SST1", "bSat", "HSat3", "CER", "HSat1A", "aSatHOR")
    ok = [r for r in rows if r.get("DJ.step") is not None and abs(r["DJ.step"]) < 0.3]
    arm_ref = {}
    for cls in arm:
        v = np.array([num(r, f"{cls}.mass_Mb") for r in ok])
        if np.isfinite(v).sum() >= 10:
            arm_ref[cls] = dict(median=float(np.nanmedian(v)), sd_rel=float(np.nanstd(v / np.nanmedian(v), ddof=1)))
    for c in carriers:
        r = by[c["sample"]]
        c["arm_content"] = {cls: round(float(num(r, f"{cls}.mass_Mb") / arm_ref[cls]["median"]), 3) for cls in arm_ref if np.isfinite(num(r, f"{cls}.mass_Mb"))}
    return dict(column=col, median=med, near=near, between=between, spread=spread, carriers=sorted(carriers, key=lambda c: c["step"]),
                transmitted=transmitted, not_transmitted=not_transmitted, de_novo=de_novo, arm_ref=arm_ref)


def mode_agreement(S: dict) -> dict:
    """fetch / scan for the samples counted both ways, and the sink capture of every scan."""
    cols = [("rDNA45S.cn_single", "45S"), ("rDNA5S.cn_single", "5S"), ("DJ.cn_single", "distal junction"), ("truth.auto", "held-out autosomal"),
            ("chrM.copies", "chrM"), ("TEL.mass_Mb", "telomeric repeat, mass")]
    both = [s for s, d in S.items() if "fetch" in d and "scan" in d]
    out = dict(n_both=len(both), columns={})
    for col, label in cols:
        ratios = [num(S[s]["fetch"], col) / num(S[s]["scan"], col) for s in both if num(S[s]["scan"], col) > 0]
        d = describe(ratios)
        if d["n"]:
            lr = np.log(fin(ratios))
            d["sd_log"] = float(lr.std(ddof=1)) if len(lr) > 1 else 0.0
            d["label"] = label
        out["columns"][col] = d
    out["capture"] = {}
    for cls, _ in FETCHABLE:
        v = [num(S[s]["scan"], f"capture.{cls}") for s in S if "scan" in S[s]]
        d = describe(v)
        d["below_99"] = [s for s in S if "scan" in S[s] and num(S[s]["scan"], f"capture.{cls}") < 0.99]
        out["capture"][cls] = d
    return out


def trio_analysis(rows, trio_list, population, columns) -> dict:
    vals = {col: {r["sample"]: num(r, col) for r in rows if np.isfinite(num(r, col))} for col, _ in columns}
    have = set(r["sample"] for r in rows)
    complete = [t for t in trio_list if {t.child, t.father, t.mother} <= have]
    out = dict(n_complete=len(complete), n_total=len(trio_list), table=[], compare=[], scatter=[])
    if len(complete) < 3:
        return out
    n = len(complete)
    base = "rDNA45S.18S.flat"
    for col, label in columns:
        try:
            t = T.transmission(vals[col], complete, population, n_perm=1000 if n >= 10 else 0, n_boot=1000 if n >= 20 else 0)
        except ValueError:
            continue
        row = dict(column=col, label=label, group=TRIO_GROUP_OF.get(col.split(".adj")[0], ""), n_trios=t["n_trios"], R=t["reliability_midparent"],
                   R_single=t["reliability_single_parent"], R_mendel=t["reliability_mendel"], spousal_r=t["spousal_r"], slope=t["midparent_slope"],
                   slope_se=t["midparent_slope_se"], r_mid=t["r_midparent"], r_father=t["r_father"], r_mother=t["r_mother"],
                   error_cv=t["error_cv"], perm_null_sd=t.get("perm_null_sd"), perm_p=t.get("perm_p"))
        if "reliability_midparent_ci95" in t:
            row["R_lo"], row["R_hi"] = t["reliability_midparent_ci95"]
            row["spousal_lo"], row["spousal_hi"] = t["spousal_r_ci95"]
            row["r_mid_lo"], row["r_mid_hi"] = t["r_midparent_ci95"]
            # the error the interval allows: from its lower bound (a slope above 1 is noise around 1)
            row["error_cv_max"] = float(np.sqrt(max(0.0, 1 - row["R_lo"]) * t["parent_sd"] ** 2) / t["parent_mean"])
        row["parent_sd"], row["child_minus_midparent_sd"], row["parent_mean"] = t["parent_sd"], t["child_minus_midparent_sd"], t["parent_mean"]
        out["table"].append(row)
        if n >= 20 and col != base and base in vals and col.split(".adj")[0] in {c for c, _ in ESTIMATORS}:
            try:
                c = T.compare(vals[col], vals[base], complete, population, n_boot=2000)
                out["compare"].append(dict(column=col, label=label, delta=c["delta"], lo=c["ci95"][0], hi=c["ci95"][1], p_better=c["p_a_better"]))
            except ValueError:
                pass
    # the scatter shows the calibrated 45S estimate (or its adjusted form) when it is among the columns
    pref = [c for c, _ in columns if c.split(".adj")[0] == "rDNA45S.cn" and c in vals] or [c for c, _ in columns if c in vals]
    col = pref[0]
    v = vals[col]
    for t in complete:
        if all(s in v for s in (t.child, t.father, t.mother)):
            out["scatter"].append(dict(child=t.child, mid=(v[t.father] + v[t.mother]) / 2, c=v[t.child], f=v[t.father], m=v[t.mother], pop=t.population))
    out["scatter_column"] = col
    out["values"] = []                                    # one row per complete trio: every column's child, father and mother values
    for t in complete:
        row = dict(child=t.child, father=t.father, mother=t.mother, population=t.population)
        for c, _ in columns:
            v = vals.get(c, {})
            row[f"{c}.child"], row[f"{c}.father"], row[f"{c}.mother"] = v.get(t.child), v.get(t.father), v.get(t.mother)
        out["values"].append(row)
    return out


def transmission_by_sex(rows, trio_list, population, columns) -> dict:
    """Every trio metric's transmission split by the sex of parent and child (`trios.by_sex`): the four
    pairings' slopes and correlations on values centred within population and sex, and the father-to-son
    minus father-to-daughter contrast, which a Y-linked quantity makes large and positive. The child's sex is
    the pedigree's, or the one the reads show where the pedigree has none."""
    sex = {r["sample"]: r.get("sex") or r.get("sex_inferred") or "" for r in rows}
    have = {r["sample"] for r in rows}
    complete = [t for t in trio_list if {t.child, t.father, t.mother} <= have and sex.get(t.child) in ("M", "F")]
    table = {}
    for col, label in columns:
        vals = {r["sample"]: num(r, col) for r in rows if np.isfinite(num(r, col))}
        res = T.by_sex(vals, complete, population, sex)
        if res:
            table[col] = dict(label=label, group=TRIO_GROUP_OF.get(col, ""), **res)
    return dict(n_sons=sum(sex[t.child] == "M" for t in complete), n_daughters=sum(sex[t.child] == "F" for t in complete), table=table)


def trio_points(values: list[dict], rows, columns) -> dict:
    """For the page's child-against-midparent figures: per metric, one [child, midparent, child value, father, mother,
    child's sex] per complete trio, on the natural scale."""
    sex = {r["sample"]: r.get("sex") or r.get("sex_inferred") or "" for r in rows}
    fin_ = lambda x: isinstance(x, (int, float)) and np.isfinite(x)
    out = {}
    for col, _ in columns:
        pts = [[t["child"], round((t[f"{col}.father"] + t[f"{col}.mother"]) / 2, 4), round(t[f"{col}.child"], 4), round(t[f"{col}.father"], 4),
                round(t[f"{col}.mother"], 4), sex.get(t["child"], "")]
               for t in values if all(fin_(t.get(f"{col}.{w}")) for w in ("child", "father", "mother"))]
        if len(pts) >= 10:
            out[col] = pts
    return out


def fetch_check(rows, S, cache, res, trio_list, population) -> dict | None:
    """The same cohort layer and the same trio test on the fetch-mode counts alone, independently of
    the scan: does the targeted fetch give the same estimate for every genome, and does it carry the
    same inherited variation (reliability from fetch beside reliability from scan)? Columns made
    from the same reads in both modes (controls, dosage regions, library properties) come out
    identical; the classes differ by what the sinks miss and by a calibration learned twice."""
    f_samples = sorted(s for s in S if "fetch" in S[s])
    if len(f_samples) < 10:
        return None
    rows_f, _, _ = cohort.cohort_table((load_result(cache / "fetch" / f"{s}.estimate.json.gz") for s in f_samples), res.anchors(), log=lambda m: None)
    for r in rows_f:
        r.update({k: v for k, v in S[r["sample"]]["fetch"].items() if k not in r})
    by_scan = {r["sample"]: r for r in rows}
    agree, points = {}, {}
    for col, label in TRIO_COLUMNS + [("truth.chrX", "chrX"), ("truth.chrY", "chrY")]:
        x = np.array([num(by_scan[r["sample"]], col) if r["sample"] in by_scan else np.nan for r in rows_f])
        y = np.array([num(r, col) for r in rows_f])
        ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        if ok.sum() >= 10:
            q = y[ok] / x[ok]
            d = corr(x[ok], y[ok])
            d.update(label=label, n=int(ok.sum()), ratio_median=float(np.median(q)), q10=float(np.percentile(q, 10)), q90=float(np.percentile(q, 90)),
                     sd_log=float(np.log(q).std(ddof=1)), identical=bool(np.all(np.abs(q - 1) < 1e-9)))
            agree[col] = d
            if col in FETCH_SCATTER and not d["identical"]:          # every genome's pair, for the page's fetch-against-scan scatters
                points[col] = [[r["sample"], round(float(a), 4), round(float(b), 4)] for r, a, b, k in zip(rows_f, x, y, ok) if k]
    trios_f = trio_analysis(rows_f, trio_list, population, TRIO_COLUMNS) if trio_list else dict(n_complete=0, n_total=0, table=[], compare=[], scatter=[], values=[])
    trios_f.pop("values", None)
    trios_f.pop("scatter", None)
    return dict(n=len(rows_f), agreement=agree, points=points, trios=trios_f)


def pc_analysis(rows, info, trio_list, population, pcs_file, log) -> dict:
    out = dict(control=info or {})
    cols = [c for c, _ in ESTIMATORS + NEGATIVE_CONTROLS if any(np.isfinite(num(r, c)) for r in rows)]
    ctrl_cols = sorted([c for c in rows[0] if re.fullmatch(r"ctrlPC\d+", c)], key=lambda c: int(c[6:])) if rows else []
    k = int(info.get("mp", 0)) if info else 0
    if ctrl_cols and k > 0:
        keep = [r for r in rows if all(np.isfinite(num(r, c)) for c in ctrl_cols[:k])]
        X = np.array([[num(r, c) for c in ctrl_cols[:k]] for r in keep])
        adj_info = {}
        for col in cols:
            y = np.array([num(r, col) for r in keep])
            adj, r2 = cohort.adjust_for_covariates(y, X, log=True)
            for r, v in zip(keep, adj):
                r[f"{col}.adj"] = None if not np.isfinite(v) else round(float(v), 4)
            n = int(np.isfinite(adj).sum())
            adj_info[col] = dict(r2=r2, chance=k / max(n - 1, 1), r2_adj=1 - (1 - r2) * (n - 1) / max(n - k - 1, 1))
        out["adjusted"] = dict(k=k, source="control-region PCs", columns=adj_info)
    # NGS-PCA's coverage PCs, when their SVD is at hand
    if pcs_file and Path(pcs_file).exists():
        from .cli import load_pcs
        pcs = load_pcs(pcs_file)
        d = Path(pcs_file).parent
        sv_path, bins_path = d / "svd.singularvalues.txt", d / "svd.bins.txt"
        if sv_path.exists() and bins_path.exists():
            sv = [float(line.split()[-1]) for line in open(sv_path) if line.split() and line.split()[-1].replace(".", "", 1).replace("e", "", 1).replace("-", "", 1).isdigit()]
            n_bins = sum(1 for line in open(bins_path) if line.strip())
            sel = pcselect.mp_select(sv, len(pcs), n_bins)
            have = [r for r in rows if r["sample"] in pcs]
            out["ngspca"] = dict(n_with_pcs=len(have), n_pc=sel.n_pc, describe=sel.describe())
            if len(have) >= sel.n_pc + 10 and sel.n_pc > 0:
                X = np.array([pcs[r["sample"]][:sel.n_pc] for r in have])
                for col in cols:
                    y = np.array([num(r, col) for r in have])
                    adj, r2 = cohort.adjust_for_covariates(y, X, log=True)
                    for r, v in zip(have, adj):
                        r[f"{col}.adj_ngspca"] = None if not np.isfinite(v) else round(float(v), 4)
                    out["ngspca"].setdefault("columns", {})[col] = dict(r2=r2, chance=sel.n_pc / max(len(have) - 1, 1))
    # the sweep: what each further control PC does to the known truths and to transmission
    if ctrl_cols and len(rows) >= 60:
        samples = [r["sample"] for r in rows]
        max_pc = min(len(ctrl_cols), max(2 * k, 10), max(len(rows) // 5, 1))
        P = np.array([[num(r, c) for c in ctrl_cols[:max_pc]] for r in rows])
        sex = {r["sample"]: r.get("sex_inferred") for r in rows}
        male = np.array([sex[s] == "M" for s in samples])
        female = np.array([sex[s] == "F" for s in samples])
        truths = {"truth.auto": np.full(len(rows), 2.0), "truth.chrX": np.where(male, 1.0, np.where(female, 2.0, np.nan)),
                  "truth.chrY": np.where(male, 1.0, np.nan)}
        dj = "DJ.cn" if "DJ.cn" in cols else "DJ.cn_single"
        truths[dj] = np.full(len(rows), 10.0)
        classes = [c for c in ("rDNA45S.cn", "rDNA45S.cn_single", "rDNA45S.18S.flat", "rDNA5S.cn") if c in cols]
        table = {c: np.array([num(r, c) for r in rows]) for c in list(truths) + classes}
        have = set(samples)
        complete = [t for t in trio_list if {t.child, t.father, t.mother} <= have]
        res = pcselect.sweep(table, P, max_pc, truths, classes, samples, complete if len(complete) >= 20 else None, population, n_boot=200)
        rec = pcselect.recommend(res)
        out["sweep"] = dict(max_pc=max_pc, rows=res, recommend=rec, n_trios=len(complete))
        log(f"[report] PC sweep: 0..{max_pc} control PCs, {len(complete)} trios; picks: " + ", ".join(f"{c}={v['pick']}" for c, v in rec.items()))
    return out


def satellite_analysis(rows, censat_dir) -> dict:
    out = dict(classes={}, hprc=None)
    for cls in SATELLITES:
        col = f"{cls}.mass_Mb"
        v = [num(r, col) for r in rows]
        if not any(np.isfinite(x) for x in v):
            continue
        d = describe(v)
        d["by_sex"] = {sex: describe([num(r, col) for r in rows if r.get("sex_inferred") == sex]) for sex in ("M", "F")}
        out["classes"][cls] = d
    if censat_dir and out["classes"]:
        cmp_rows, stats = hprc.compare(rows, censat_dir)
        if cmp_rows:
            out["hprc"] = dict(rows=cmp_rows, stats=stats, n_samples=len({r["sample"] for r in cmp_rows}))
    return out


def rdna_assemblies(rows, censat_dir, col: str, unit_bp: int) -> dict | None:
    """What long-read assemblies hold of the rDNA: sequence annotated as rDNA in the HPRC assemblies of
    the counted genomes that have one, against the 45S copy number measured here times the unit length.
    The ten arrays of a diploid genome (five acrocentrics, two homologues) give the mean array size."""
    if not censat_dir or not unit_bp:
        return None
    out = []
    for r in rows:
        cn = num(r, col)
        if not np.isfinite(cn) or cn <= 0:
            continue
        bp, longest, n_files = hprc.rdna_in_assembly(r["sample"], censat_dir)
        if n_files != 2:
            continue
        implied = cn * unit_bp
        out.append(dict(sample=r["sample"], assembly_rDNA_Mb=round(bp / 1e6, 3), longest_stretch_Mb=round(longest / 1e6, 3), cn45=round(cn, 1),
                        implied_Mb=round(implied / 1e6, 2), fraction=round(bp / implied, 3), mean_array_Mb=round(implied / 10 / 1e6, 2)))
    if not out:
        return None
    return dict(rows=out, n=len(out), column=col, unit_bp=unit_bp, fraction=describe([x["fraction"] for x in out]),
                longest_stretch_Mb=describe([x["longest_stretch_Mb"] for x in out]), mean_array_Mb=describe([x["mean_array_Mb"] for x in out]))


def hall_comparison(rows, hall_path) -> dict | None:
    if not hall_path or not Path(hall_path).exists():
        return None
    import csv
    H = {r["Sample"]: r for r in csv.DictReader(open(hall_path), delimiter="\t")}
    shared = [r for r in rows if r["sample"] in H and np.isfinite(num(r, "rDNA45S.18S.flat"))]
    if len(shared) < 3:
        return dict(n=len(shared))
    theirs = np.array([float(H[r["sample"]]["HC.18S.CN"]) for r in shared])
    ours_flat = np.array([num(r, "rDNA45S.18S.flat") / 2 for r in shared])
    dup_c = np.array([num(r, "ctrl_dup_frac") for r in shared])
    dup_r = np.array([num(r, "rDNA45S.dup_flag_frac") for r in shared])
    pred = ours_flat * (1 - dup_r) / (1 - dup_c)
    cal_col = "rDNA45S.cn" if all(np.isfinite(num(r, "rDNA45S.cn")) for r in shared) else "rDNA45S.cn_single"
    ours_cal = np.array([num(r, cal_col) / 2 for r in shared])
    return dict(n=len(shared), flat=corr(theirs, ours_flat), flat_ratio=float(np.mean(theirs / ours_flat)),
                dup_corrected_ratio=float(np.mean(theirs / pred)), dup_corrected_ratio_sd=float(np.std(theirs / pred, ddof=1)),
                calibrated=corr(theirs, ours_cal), calibrated_ratio=float(np.mean(theirs / ours_cal)), calibrated_column=cal_col,
                points=[dict(sample=r["sample"], theirs=float(t), flat=float(o), cal=float(c)) for r, t, o, c in zip(shared, theirs, ours_flat, ours_cal)])


def replicate_analysis(pilot_dir) -> dict | None:
    """The pilot's twelve people, each sequenced twice: the NYGC library (NovaSeq 2x150, 2019) and
    an older library of the same cell line (HGSVC HiSeq 2500 2x126, 2015, or Illumina Platinum
    HiSeq 2000 2x100, 2012-13). Test-retest agreement of each estimator across the two
    technologies, as a one-way intraclass correlation (absolute agreement, so an offset between
    technologies counts against it), the pair SD of the log ratio and the within-person CV. The
    calibrated 45S values are those with the anchor windows chosen with the person's family held
    out (`evaluate_pilot.py`), so the cross-technology level is out of sample."""
    if not pilot_dir:
        return None
    import csv
    d = Path(pilot_dir)
    held, rep = d / "pilot_heldout.tsv", d / "pilot_replicates.tsv"
    if not held.exists():
        return None
    H = list(csv.DictReader(open(held), delimiter="\t"))
    R = {r["sample"]: r for r in csv.DictReader(open(rep), delimiter="\t")} if rep.exists() else {}
    pairs = {"calibrated": ("45S, NGS-DOSE calibrated (anchors chosen with the family held out)", [(float(r["nygc"]), float(r["replicate"]), r["sample"]) for r in H]),
             "flat": ("45S, 18S depth ratio (published method; no GC model)", [(float(r["nygc_18S_flat"]), float(r["replicate_18S_flat"]), r["sample"]) for r in H])}
    for key, label, a, b in (("rDNA5S", "5S, NGS-DOSE (fragment-GC model)", "rDNA5S, fragment-GC model [NYGC]", "rDNA5S, fragment-GC model [replicate]"),
                             ("DJ", "distal junction, NGS-DOSE (fragment-GC model)", "DJ, fragment-GC model [NYGC]", "DJ, fragment-GC model [replicate]")):
        if R and all(a in R[s] and b in R[s] for s in R):
            pairs[key] = (label, [(float(R[s][a]), float(R[s][b]), s) for s in sorted(R)])

    def stats(xy):
        x, y = np.array([q[0] for q in xy]), np.array([q[1] for q in xy])
        lr = np.log(y / x)
        n, mean = len(x), (x + y) / 2
        msb, msw = 2 * float(mean.var(ddof=1)), float(((x - y) ** 2).sum() / (2 * n))
        return dict(n=n, mean_log_ratio=float(lr.mean()), offset=float(np.exp(lr.mean()) - 1), sd_log_ratio=float(lr.std(ddof=1)),
                    r=float(np.corrcoef(x, y)[0, 1]), icc=float((msb - msw) / (msb + msw)), within_cv=float(np.sqrt(msw) / mean.mean()))

    out = dict(n=len(H), table={}, points=[])
    for key, (label, xy) in pairs.items():
        out["table"][key] = dict(label=label, **stats(xy))
    # the depth ratio with its offset between technologies removed: what a batch correction would leave
    off = float(np.exp(out["table"]["flat"]["mean_log_ratio"]))
    out["table"]["flat_centred"] = dict(label="45S, 18S depth ratio (published), its offset between technologies removed", **stats([(x, y / off, q) for x, y, q in pairs["flat"][1]]))
    for si, key in enumerate(("calibrated", "flat")):
        out["points"] += [dict(sample=q, x=x, y=y, si=si) for x, y, q in pairs[key][1]]
    return out


def ngspca_qc_comparison(rows, qc_path) -> dict | None:
    """NGS-PCA's per-sample QC for the same CRAMs (mosdepth, 1-kb bins, duplicate-flagged reads
    excluded by mosdepth's default filter): mitochondrial copies per cell as 2 x chrM mean
    coverage / high-quality median autosomal coverage, the chrX and chrY coverage ratios, the
    autosomal depth and the sex inferred from them - a coverage route to four quantities measured
    here by fragment-end counting under the GC model."""
    if not qc_path or not Path(qc_path).exists():
        return None
    import csv
    Q = {r["SAMPLE_ID"]: r for r in csv.DictReader(open(qc_path), delimiter="\t")}
    shared = [r for r in rows if r["sample"] in Q]
    if len(shared) < 3:
        return dict(n=len(shared), n_qc=len(Q))
    g = lambda r, c: float(Q[r["sample"]][c]) if Q[r["sample"]].get(c) not in (None, "", "NA") else float("nan")
    for r in shared:
        r["ngspca.MTDNA_CN"], r["ngspca.chrX"], r["ngspca.chrY"] = g(r, "MTDNA_CN"), 2 * g(r, "X_COV_RATIO"), 2 * g(r, "Y_COV_RATIO")
        r["ngspca.depth"], r["ngspca.sex"], r["ngspca.batch"] = g(r, "MEAN_AUTOSOMAL_COV"), Q[r["sample"]].get("INFERRED_SEX", ""), Q[r["sample"]].get("RELEASE_BATCH", "")

    def against(ours, theirs):
        x, y = np.asarray(ours, float), np.asarray(theirs, float)
        ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        d = corr(x, y)
        if ok.sum() >= 2:
            q = y[ok] / x[ok]
            d["ratio"] = dict(median=float(np.median(q)), q10=float(np.percentile(q, 10)), q90=float(np.percentile(q, 90)), sd_log=float(np.log(q).std(ddof=1)))
        return d

    out = dict(n=len(shared), n_qc=len(Q))
    ours_m, theirs_m = [num(r, "chrM.copies") for r in shared], [r["ngspca.MTDNA_CN"] for r in shared]
    out["mtdna"] = against(ours_m, theirs_m)
    if any(np.isfinite(num(r, "flat.chrM")) for r in shared):
        out["mtdna_flat"] = against([num(r, "flat.chrM") for r in shared], theirs_m)
    with np.errstate(invalid="ignore", divide="ignore"):
        out["mtdna_ratio_vs_depth"] = corr([num(r, "depth") for r in shared], np.log(np.asarray(theirs_m, float) / np.asarray(ours_m, float)))
    out["chrX"] = against([num(r, "truth.chrX") for r in shared], [r["ngspca.chrX"] for r in shared])
    men = [r for r in shared if r.get("sex_inferred") == "M"]
    out["chrY_men"] = against([num(r, "truth.chrY") for r in men], [r["ngspca.chrY"] for r in men])
    out["chrY_intact_men"] = describe([r["ngspca.chrY"] for r in men if num(r, "truth.chrY") >= 0.85])
    out["depth"] = against([num(r, "depth") for r in shared], [r["ngspca.depth"] for r in shared])
    both = [r for r in shared if r["ngspca.sex"] and r.get("sex_inferred")]
    out["sex_n"], out["sex_agree"] = len(both), sum(1 for r in both if r["ngspca.sex"] == r["sex_inferred"])
    out["mosaic_X"] = [dict(sample=r["sample"], ours=num(r, "truth.chrX"), theirs=r["ngspca.chrX"]) for r in shared if r.get("sex_inferred") == "F" and num(r, "truth.chrX") < 1.85]
    out["mosaic_Y"] = [dict(sample=r["sample"], ours=num(r, "truth.chrY"), theirs=r["ngspca.chrY"]) for r in men if num(r, "truth.chrY") < 0.85]
    return out


def biology(rows) -> dict:
    """The relations the cohort can test: 45S with 5S, 45S with the culture's mitochondrial content and
    EBV load, the late-replicating controls with each other (the S-phase hypothesis), and the
    duplicate-flag rates inside and outside the rDNA."""
    g = lambda c: np.array([num(r, c) for r in rows])
    cn45 = g("rDNA45S.cn") if any(np.isfinite(g("rDNA45S.cn"))) else g("rDNA45S.cn_single")
    out = dict(cn45_column="rDNA45S.cn" if any(np.isfinite(g("rDNA45S.cn"))) else "rDNA45S.cn_single")
    out["cn45_vs_5S"] = corr(cn45, g("rDNA5S.cn") if any(np.isfinite(g("rDNA5S.cn"))) else g("rDNA5S.cn_single"))
    out["cn45_vs_chrM"] = corr(cn45, g("chrM.copies"))
    out["cn45_vs_EBV"] = corr(cn45, g("chrEBV.copies"))
    out["log_cn45_vs_log_chrM"] = corr(np.log(cn45), np.log(g("chrM.copies")))
    # women whose culture has not lost an X: a mosaic loss is not an S-phase effect
    female = np.array([r.get("sex_inferred") == "F" for r in rows]) & (g("truth.chrX") >= 1.85) & (g("truth.chrX") <= 2.15)
    dj = g("DJ.cn") if any(np.isfinite(g("DJ.cn"))) else g("DJ.cn_single")
    out["DJ_vs_chrX_female"] = corr(np.where(female, dj, np.nan), np.where(female, g("truth.chrX"), np.nan))
    out["cn45_vs_chrX_female"] = corr(np.where(female, cn45, np.nan), np.where(female, g("truth.chrX"), np.nan))
    out["cn45_cv"] = float(np.nanstd(cn45, ddof=1) / np.nanmean(cn45))
    out["DJ_vs_chrM"] = corr(dj, g("chrM.copies"))
    # what the fragment-GC model removes: the uncorrected 18S ratio follows each library's GC bias
    flat, s18, gc65 = g("rDNA45S.18S.flat"), g("rDNA45S.18S"), g("gc_rel_65")
    with np.errstate(invalid="ignore", divide="ignore"):
        for r, a, b in zip(rows, flat / cn45, s18 / cn45):
            r["rDNA45S.18S.flat_over_cn"] = None if not np.isfinite(a) else round(float(a), 4)
            r["rDNA45S.18S_over_cn"] = None if not np.isfinite(b) else round(float(b), 4)
        out["gc_bias"] = dict(flat_vs_gc=corr(np.log(flat / cn45), gc65), modelled_vs_gc=corr(np.log(s18 / cn45), gc65),
                              flat_sd_log=float(np.nanstd(np.log(flat / cn45), ddof=1)), modelled_sd_log=float(np.nanstd(np.log(s18 / cn45), ddof=1)),
                              gc65=describe(gc65), flat_vs_cn=corr(cn45, flat))
    out["dup"] = dict(control=describe(g("ctrl_dup_frac")), rDNA=describe(g("rDNA45S.dup_flag_frac")),
                      ratio=describe(g("rDNA45S.dup_flag_frac") / g("ctrl_dup_frac")))
    out["chrM"] = describe(g("chrM.copies"))
    out["chrEBV"] = describe(g("chrEBV.copies"))
    return out


def _ratio(a, b) -> float:
    return a / b if isinstance(a, (int, float)) and isinstance(b, (int, float)) and b > 0 else float("nan")


def by_group(rows, col, key) -> list[dict]:
    groups = {}
    for r in rows:
        v = num(r, col)
        if np.isfinite(v) and r.get(key):
            groups.setdefault(r[key], []).append(v)
    return [dict(group=k, **describe(v)) for k, v in sorted(groups.items())]


# ------------------------------------------------------------------------------------------------
# the run
# ------------------------------------------------------------------------------------------------
def build(a, log=lambda m: print(m, file=sys.stderr)) -> dict:
    out = Path(a.out)
    (out / "data").mkdir(parents=True, exist_ok=True)
    cache = Path(a.cache) if a.cache else out / "cache"
    bundle_dir = a.resources or resources.default_bundle()
    res = resources.Bundle(bundle_dir)
    counts = {}
    for mode, d in (("scan", a.scan), ("fetch", a.fetch)):
        if d:
            counts[mode] = sorted(Path(p) for p in glob.glob(os.path.join(d, "*.json.gz")))
    if not counts:
        raise SystemExit("nothing to report: give --scan and/or --fetch directories of counts files")
    summaries = {mode: estimate_all(paths, cache / mode, bundle_dir, a.jobs, log) for mode, paths in counts.items()}
    S: dict[str, dict] = {}
    for mode, rows in summaries.items():
        for r in rows:
            S.setdefault(r["sample"], {})[mode] = r
    primary = "scan" if "scan" in counts and counts["scan"] else "fetch"
    ped, trio_list, population = load_pedigree(a.pedigree) if a.pedigree else ({}, [], {})

    # the cohort layer on the primary mode's estimates (streamed from the cache)
    prim_samples = sorted(s for s in S if primary in S[s])
    est_paths = [cache / primary / f"{s}.estimate.json.gz" for s in prim_samples]
    rows, eff, info = cohort.cohort_table((load_result(p) for p in est_paths), res.anchors(), log=log)
    for r in rows:
        r.update({k: v for k, v in S[r["sample"]][primary].items() if k not in r})
        r["mode_used"] = primary
        p = ped.get(r["sample"], {})
        r["sex"], r["pop"], r["superpop"] = p.get("sex", ""), p.get("pop", ""), p.get("superpop", "")
        r["sex_inferred"] = inferred_sex(r)
        f = S[r["sample"]].get("fetch")
        if f is not None and primary == "scan":
            for cls, col in FETCHABLE:
                sv, fv = num(r, col), num(f, col)
                r[f"fetch_ratio.{cls}"] = round(fv / sv, 6) if sv > 0 and np.isfinite(fv) else None
    dj = dj_steps(rows, ped)
    engines = {}
    for r in rows:
        engines[r.get("engine")] = engines.get(r.get("engine"), 0) + 1
    majority_engine = max(engines, key=engines.get) if engines else None

    data = dict(meta=dict(title=a.title, as_of=a.as_of or dt.date.today().isoformat(), total=a.total, generator=f"ngsdose {__version__}",
                          primary_mode=primary, n_scan=len(counts.get("scan", [])), n_fetch=len(counts.get("fetch", [])),
                          n_both=sum(1 for s in S if "scan" in S[s] and "fetch" in S[s]), n=len(rows), bundle=res.meta.get("name"),
                          engines=engines, panel_sha=sorted({r.get("panel_sha256") for r in rows if r.get("panel_sha256")}),
                          controls_sha=sorted({r.get("controls_sha256") for r in rows if r.get("controls_sha256")}),
                          sinks_sha=sorted({S[s]["fetch"].get("sinks_sha256") for s in S if "fetch" in S[s] and S[s]["fetch"].get("sinks_sha256")}),
                          eof_absent=[r["sample"] for r in rows if r.get("eof_marker") not in (None, "present")]))
    # progress by superpopulation, against the pedigree's totals
    if ped:
        tot, done = {}, {}
        for s, p in ped.items():
            tot[p["superpop"]] = tot.get(p["superpop"], 0) + 1
            if s in S and primary in S[s]:
                done[p["superpop"]] = done.get(p["superpop"], 0) + 1
        data["meta"]["by_superpop"] = {k: dict(total=tot[k], done=done.get(k, 0)) for k in sorted(tot)}
        data["meta"]["pedigree_total"] = len(ped)
    data["qc"] = dict(depth=describe([num(r, "depth") for r in rows]), insert=describe([num(r, "insert_median") for r in rows]),
                      dup=describe([num(r, "ctrl_dup_frac") for r in rows]), read_length=sorted({r.get("read_length") for r in rows}),
                      gc_max_se=describe([num(r, "gc_curve_max_se") for r in rows]),
                      elapsed=describe([num(S[s]["scan"], "elapsed_sec") for s in S if "scan" in S[s]]),
                      elapsed_fetch=describe([num(S[s]["fetch"], "elapsed_sec") for s in S if "fetch" in S[s]]),
                      flagged_chromosomes=[(r["sample"], r["flagged_chromosomes"]) for r in rows if r.get("flagged_chromosomes")],
                      placement_bins=sorted({str(r.get("placement_bin")) for r in rows}))
    data["known_truth"] = known_truth(rows)
    data["known_truth"]["DJ_steps"] = dj
    data["modes"] = mode_agreement(S)
    data["rdna"] = {c: describe([num(r, c) for r in rows]) for c in ("rDNA45S.cn", "rDNA45S.cn_single", "rDNA45S.18S.flat", "rDNA5S.cn", "DJ.cn")}
    data["rdna"]["by_superpop"] = by_group(rows, "rDNA45S.cn" if data["rdna"]["rDNA45S.cn"]["n"] else "rDNA45S.cn_single", "superpop")
    data["rdna"]["by_pop"] = by_group(rows, "rDNA45S.cn" if data["rdna"]["rDNA45S.cn"]["n"] else "rDNA45S.cn_single", "pop")
    data["biology"] = biology(rows)
    data["trios"] = trio_analysis(rows, trio_list, population, TRIO_COLUMNS) if trio_list else dict(n_complete=0, n_total=0, table=[], compare=[], scatter=[], values=[])
    if trio_list:
        data["trios"]["by_sex"] = transmission_by_sex(rows, trio_list, population, TRIO_COLUMNS)
        data["trios"]["points"] = trio_points(data["trios"].get("values") or [], rows, TRIO_COLUMNS)
    data["fetch_check"] = fetch_check(rows, S, cache, res, trio_list, population) if primary == "scan" and "fetch" in counts else None
    data["pcs"] = pc_analysis(rows, info, trio_list, population, a.pcs, log)
    if data["pcs"].get("adjusted") and data["trios"]["n_complete"] >= 3:
        adj_cols = [(f"{c}.adj", f"{l}, adjusted") for c, l in ESTIMATORS + NEGATIVE_CONTROLS if any(np.isfinite(num(r, f"{c}.adj")) for r in rows)]
        data["trios_adjusted"] = trio_analysis(rows, trio_list, population, adj_cols)
    data["satellites"] = satellite_analysis(rows, a.censat)
    unit45 = len(res.units().get("rDNA45S", "")) if a.censat else 0
    data["rdna"]["assemblies"] = rdna_assemblies(rows, a.censat, "rDNA45S.cn" if data["rdna"]["rDNA45S.cn"]["n"] else "rDNA45S.cn_single", unit45)
    data["hall"] = hall_comparison(rows, a.hall)
    data["replicates"] = replicate_analysis(a.pilot)
    data["ngspca_qc"] = ngspca_qc_comparison(rows, a.qc)
    for r in rows:
        r["flags"] = "; ".join(flags_for(r, majority_engine))
    data["flags"] = [(r["sample"], r["flags"]) for r in rows if r["flags"]]

    # tables on disk
    write_table(rows, out / "data" / "cohort.tsv")
    if "fetch" in summaries:
        write_table(summaries["fetch"], out / "data" / "fetch.tsv")
    if data["modes"]["n_both"]:
        write_table([dict(sample=s, **{f"fetch_over_scan.{cls}": _ratio(S[s]["fetch"].get(col), S[s]["scan"].get(col)) for cls, col in FETCHABLE},
                          **{f"capture.{cls}": S[s]["scan"].get(f"capture.{cls}") for cls, _ in FETCHABLE})
                     for s in sorted(S) if "scan" in S[s] and "fetch" in S[s]], out / "data" / "modes.tsv")
    if data["trios"]["table"]:
        write_table(data["trios"]["table"] + data.get("trios_adjusted", {}).get("table", []), out / "data" / "transmission.tsv")
    if data["trios"].get("values"):
        write_table(data["trios"]["values"], out / "data" / "trios.tsv")
    bs = (data["trios"].get("by_sex") or {}).get("table") or {}
    if bs:
        write_table([dict(metric=col, label=d["label"], group=d["group"],
                          **{f"{k}.{f}": d[k][f] for k, _, _ in T.PAIRS if k in d for f in ("n", "slope", "slope_se", "r", "r_lo", "r_hi")},
                          **{f"father_contrast.{f}": v for f, v in (d.get("father_contrast") or {}).items()},
                          **{f"father_slope_contrast.{f}": v for f, v in (d.get("father_slope_contrast") or {}).items()}) for col, d in bs.items()],
                    out / "data" / "transmission_by_sex.tsv")
    if data.get("fetch_check"):
        fc, R_s, R_f = data["fetch_check"], {t["column"]: t for t in data["trios"]["table"]}, {t["column"]: t for t in data["fetch_check"]["trios"]["table"]}
        write_table([dict(metric=col, **{k: v for k, v in d.items() if k != "label"}, label=d["label"],
                          R_scan=R_s.get(col, {}).get("R"), R_scan_lo=R_s.get(col, {}).get("R_lo"), R_scan_hi=R_s.get(col, {}).get("R_hi"),
                          R_fetch=R_f.get(col, {}).get("R"), R_fetch_lo=R_f.get(col, {}).get("R_lo"), R_fetch_hi=R_f.get(col, {}).get("R_hi"),
                          r_mid_scan=R_s.get(col, {}).get("r_mid"), r_mid_fetch=R_f.get(col, {}).get("r_mid")) for col, d in fc["agreement"].items()],
                    out / "data" / "fetch_check.tsv")
    data["trios"].pop("values", None)                     # on disk, not in the page
    if data["pcs"].get("sweep"):
        write_table(data["pcs"]["sweep"]["rows"], out / "data" / "pcsweep.tsv")
    if data["satellites"].get("hprc"):
        write_table(data["satellites"]["hprc"]["rows"], out / "data" / "satellites_hprc.tsv")
    if data["rdna"]["assemblies"]:
        write_table(data["rdna"]["assemblies"]["rows"], out / "data" / "rdna_hprc.tsv")
    write_table([dict(sample=s, flags=f) for s, f in data["flags"]], out / "data" / "flags.tsv")
    Path(out / "data" / "efficiencies.json").write_text(json.dumps(eff))
    data["efficiencies"] = {cls: dict(start=e["start"], a=e["a"], anchor=e["anchor"], gc=e["gc"]) for cls, e in eff.items()}
    data["samples"] = [sample_record(r) for r in rows]
    data = rnd(data, 5)
    dump(data, out / "report.json")
    (out / "index.html").write_text(render(data, rows))
    log(f"[report] {len(rows)} samples -> {out / 'index.html'}")
    return data


SAMPLE_COLUMNS = ["sample", "sex", "sex_inferred", "pop", "superpop", "mode_used", "engine", "depth", "insert_median", "ctrl_dup_frac", "rDNA45S.dup_flag_frac",
                  "gc_curve_max_se", "truth.auto", "truth.chrX", "truth.chrY", "chrM.copies", "chrEBV.copies", "rDNA45S.cn", "rDNA45S.cn_single",
                  "rDNA45S.18S.flat", "rDNA5S.cn", "rDNA5S.cn_single", "DJ.cn", "DJ.cn_single", "rDNA45S.cn.adj", "rDNA45S.cn.adj_ngspca",
                  "elapsed_sec", "flags", "DJ.step", "gc_rel_65", "rDNA45S.18S.flat_over_cn", "rDNA45S.18S_over_cn", "flat.chrM",
                  "ngspca.MTDNA_CN", "ngspca.chrX", "ngspca.chrY", "ngspca.depth", "ngspca.batch"] + [f"{c}.mass_Mb" for c in SATELLITES] + [f"fetch_ratio.{c}" for c, _ in FETCHABLE] + [f"capture.{c}" for c, _ in FETCHABLE]


def sample_record(r: dict) -> dict:
    o = {}
    for c in SAMPLE_COLUMNS:
        v = r.get(c)
        if isinstance(v, float) and not math.isfinite(v):
            v = None
        o[c] = v
    return o


# ------------------------------------------------------------------------------------------------
# the page
# ------------------------------------------------------------------------------------------------
def esc(x) -> str:
    return html.escape(str(x))


def fmt(x, nd=2, pct=False) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "–"
    if pct:
        return f"{100 * x:.{nd}f}%"
    if isinstance(x, int):
        return f"{x:,}"
    return f"{x:,.{nd}f}"


def render(data: dict, rows: list[dict]) -> str:
    from .report_page import page
    return page(data, rows)


def add_arguments(ap):
    ap.add_argument("--scan", help="directory of scan-mode counts files")
    ap.add_argument("--fetch", help="directory of fetch-mode counts files")
    ap.add_argument("-p", "--pedigree", help="1000 Genomes pedigree (sex, population, trios)")
    ap.add_argument("--pcs", help="NGS-PCA svd.pcs.txt (with svd.singularvalues.txt and svd.bins.txt beside it)")
    ap.add_argument("--censat", help="directory of HPRC CenSat annotations (<sample>_<hap>_...cenSat.bed)")
    ap.add_argument("--hall", help="Hall, Turner & Queitsch 2021 Supplementary Data 1 (per-sample table for the same CRAMs)")
    ap.add_argument("--qc", help="NGS-PCA sample_qc.tsv for the same cohort (mosdepth-based mitochondrial copy number, X/Y coverage ratios, depth, inferred sex)")
    ap.add_argument("--pilot", help="the pilot directory (pilot_heldout.tsv, pilot_replicates.tsv): the same twelve people on two sequencing technologies")
    ap.add_argument("--total", type=int, default=3202, help="samples the run will have when complete")
    ap.add_argument("--title", default="NGS-DOSE · 1000 Genomes 30× cohort")
    ap.add_argument("--as-of", help="date stamp (default: today)")
    ap.add_argument("-r", "--resources", default=None)
    ap.add_argument("-o", "--out", required=True, help="output directory (index.html, report.json, data/)")
    ap.add_argument("--cache", help="where per-sample estimates are cached (default: <out>/cache)")
    ap.add_argument("-j", "--jobs", type=int, default=4)
    return ap


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(prog="ngsdose report", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    build(add_arguments(ap).parse_args(argv))


if __name__ == "__main__":
    main()
