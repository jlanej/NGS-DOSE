"""Per-sample estimation: counts JSON -> copy number / array mass per class, with diagnostics.

Estimator for a positional class (one with a unit consensus, e.g. the 45S rDNA unit):

    C_w = 2 * obs_w / exp_w,     exp_w = sum over position-strands in w of lambda(g(p))

where obs_w counts fragment 5' ends assigned to window w of the unit by k-mer placement and
lambda is the fragment-GC rate fitted on single-copy controls (gcmodel.py). A position-strand
enters obs and exp only if (i) a read starting there would hold at least `min_kmers` panel
k-mers, so that it is reliably recovered, and (ii) its fragment GC lies inside the range where
the control curve is supported. Everything else is masked from numerator and denominator alike.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from . import gcmodel
from .io import Panel, PanelClass

ANCHOR_GC = (0.40, 0.60)


def nearest_table(counts: dict, L: int | None) -> dict:
    tabs = counts["gc_tables"]
    if not tabs:
        raise ValueError("counts file holds no GC tables")
    target = L if L else counts["insert_median"]
    if not target:
        raise ValueError("no insert-size estimate in counts file; pass L explicitly")
    return min(tabs, key=lambda t: abs(t["l"] - target))


def callable_masks(pc: PanelClass, k: int, R: int, min_kmers: int) -> tuple[np.ndarray, np.ndarray]:
    """Position-strands at which a read of length R holds >= min_kmers panel k-mers."""
    U = pc.length
    has = np.zeros(U, np.int64)
    has[pc.kmer_pos] = 1
    nk = R - k + 1
    if nk < 1:
        raise ValueError(f"read length {R} shorter than k={k}")
    if pc.circular:
        ext = np.tile(has, 3)
        cs = np.r_[0, np.cumsum(ext)]
        p = np.arange(U) + U
        fwd = cs[p + nk] - cs[p]                       # k-mer starts p .. p+R-k
        rev = cs[p - k + 2] - cs[p - R + 1]            # k-mer starts x-R+1 .. x-k+1
    else:
        cs = np.r_[0, np.cumsum(has)]
        p = np.arange(U)
        fwd = np.where(p + nk <= U, cs[np.minimum(p + nk, U)] - cs[p], 0)
        rev = np.where(p - R + 1 >= 0, cs[np.clip(p - k + 2, 0, U)] - cs[np.clip(p - R + 1, 0, U)], 0)
    return fwd >= min_kmers, rev >= min_kmers


@dataclass
class ControlQC:
    rescale: float
    n_regions: int
    n_flagged_regions: int
    flagged_chromosomes: list[str]
    region_log_sd: float
    ratios: np.ndarray


def _undefined_window_factor(regions: list[dict], tables: np.ndarray) -> np.ndarray:
    """A position whose fragment window runs over an N has no GC and is in no table, but a read
    that starts there is in the region's observed count. Expected counts are therefore scaled
    from the tabulated position-strands to all 2 x length of them (the bundle's regions have
    none - its build refuses them - so this is exactly 1 there; GRCh38's chrM has one at 3,107)."""
    total = 2.0 * np.array([r["len"] for r in regions], float)
    return total / np.maximum(tables.sum(1), 1.0)


def _chrom_key(name: str) -> str:
    return name.rsplit(":", 1)[0]


def control_qc(counts: dict, curve: gcmodel.GCCurve, region_tables: np.ndarray | None,
               max_region_dev: float = 0.30, max_chrom_dev: float = 0.04, z_min: float = 5.0) -> ControlQC | None:
    """Robust denominator: drop control regions (CNV) and whole chromosomes (aneuploidy, common in
    cell lines) whose depth departs from the GC-model expectation, and rescale the curve."""
    if region_tables is None or not counts["regions"]:
        return None
    if len(counts["regions"]) != len(region_tables):
        raise ValueError(f"{counts.get('sample')}: counts hold {len(counts['regions'])} control regions, the bundle {len(region_tables)} - "
                         "these counts were made with a different controls file")
    is_ctrl = np.array([r.get("role", "control") == "control" for r in counts["regions"]])
    regs = [r for r, k in zip(counts["regions"], is_ctrl) if k]
    region_tables = region_tables[is_ctrl]
    obs = np.array([r["obs"] for r in regs], float)
    rate = np.where(np.isnan(curve.rate), 0.0, curve.rate)
    supported = ~np.isnan(curve.rate)
    exp = region_tables[:, supported] @ rate[supported]
    frac = region_tables[:, supported].sum(1) / np.maximum(region_tables.sum(1), 1)
    exp = exp / np.maximum(frac, 1e-9)                 # regions partly outside GC support
    exp = exp * _undefined_window_factor(regs, region_tables)
    ratio = obs / np.maximum(exp, 1e-9)
    lr = np.log(np.maximum(ratio, 1e-9))
    chroms = np.array([_chrom_key(r["name"]) for r in regs])
    # A region or chromosome is dropped only if its departure is both large AND significant: at
    # low depth sampling noise alone exceeds any fixed threshold, and trimming noise is biased
    # (log ratios of counts are left-skewed), which would push the denominator up.
    se = 1.0 / np.sqrt(np.maximum(exp, 1.0))
    centre = np.median(lr)
    phi = max(1.0, (1.4826 * np.median(np.abs(lr - centre)) / np.median(se)) ** 2)      # overdispersion
    keep = ~((np.abs(lr - centre) > max_region_dev) & (np.abs(lr - centre) > z_min * np.sqrt(phi) * se))
    flagged_chroms = []
    tot = np.log(obs[keep].sum() / exp[keep].sum())
    for c in np.unique(chroms):
        m = (chroms == c) & keep
        if m.sum() < 5:
            continue
        dev = abs(np.log(obs[m].sum() / exp[m].sum()) - tot)
        if dev > max_chrom_dev and dev > z_min * np.sqrt(phi / exp[m].sum()):
            flagged_chroms.append(str(c))
    keep &= ~np.isin(chroms, flagged_chroms)
    # The fitted curve already reproduces the pooled control counts, so the correction is the
    # *relative* change from dropping flagged regions: exactly 1 when nothing is flagged. (Using
    # obs/exp of the kept regions directly is biased whenever the curve's supported GC range is
    # narrow - low depth - because reads at unsupported, GC-rich control positions then count
    # in obs but are only averaged into exp.)
    rescale = float((obs[keep].sum() / exp[keep].sum()) / (obs.sum() / exp.sum()))
    mad = 1.4826 * np.median(np.abs(lr[keep] - np.median(lr[keep])))
    return ControlQC(rescale, len(regs), int((~keep).sum()), sorted(flagged_chroms, key=_natural), float(mad), ratio)


def test_regions(counts: dict, curve: gcmodel.GCCurve, region_tables: np.ndarray | None) -> dict:
    """Copy number of every labelled set of non-control regions, estimated exactly as a class is:
    2 * observed / expected under the control GC curve. Role `test` regions have a known answer
    (held-out autosomal sequence: 2; chrX: 1 or 2; chrY: 1 or 0); role `dosage` regions do not
    (mitochondrial genome, EBV episome: copies per cell). A set whose contig the alignment
    header lacks is left out rather than reported as zero."""
    if region_tables is None:
        return {}
    rate = np.where(np.isnan(curve.rate), 0.0, curve.rate) * curve.rescale
    out: dict[str, dict] = {}
    other = [(i, r) for i, r in enumerate(counts["regions"]) if r.get("role", "control") != "control" and not r.get("absent")]
    for label in sorted({r.get("label", "") for _, r in other}):
        idx = [i for i, r in other if r.get("label", "") == label]
        obs = np.array([counts["regions"][i]["obs"] for i in idx], float)
        tab = region_tables[idx]
        # reads at positions outside the curve's supported GC range are in obs; give them the
        # region's average supported rate rather than none
        exp = (tab @ rate) / np.maximum(tab[:, rate > 0].sum(1) / np.maximum(tab.sum(1), 1), 1e-9)
        exp = exp * _undefined_window_factor([counts["regions"][i] for i in idx], tab)
        per = 2 * obs / np.maximum(exp, 1e-9)
        out[label] = dict(role=counts["regions"][idx[0]].get("role"), n_regions=len(idx), cn=float(2 * obs.sum() / exp.sum()), cn_median=float(np.median(per)),
                          region_sd=float(1.4826 * np.median(np.abs(per - np.median(per)))))
    return out


def _natural(s: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def estimate_positional(cls: dict, pc: PanelClass, seq: str, k: int, curve: gcmodel.GCCurve, R: int,
                        window: int = 250, min_kmers: int = 20,
                        features: list[tuple[str, int, int]] | None = None,
                        anchors: list[tuple[int, int]] | None = None) -> dict:
    """`anchors`: unit intervals that set the absolute level (a bundle's empirically clean windows);
    without them the rule is fragment GC within ANCHOR_GC."""
    U, b = pc.length, cls["bin"]
    if len(seq) != U:
        raise ValueError(f"unit FASTA for {pc.name} has {len(seq)} bp, panel says {U}")
    ef, er, gf, gr = gcmodel.expected_per_position(seq, curve, pc.circular)
    cf, cr = callable_masks(pc, k, R, min_kmers)
    nb = len(cls["fwd"])
    obs = np.stack([np.array(cls["fwd"], float), np.array(cls["rev"], float)])          # strand x bin
    pad = nb * b - U
    def binned(x, fn):
        x = np.r_[x, np.full(pad, np.nan)] if pad else x
        return fn(x.reshape(nb, b), axis=1)
    ok_pos = np.stack([cf & ~np.isnan(ef), cr & ~np.isnan(er)])
    usable = np.stack([binned(ok_pos[s].astype(float), np.nanmin) == 1 for s in (0, 1)])  # every position ok
    exp = np.stack([binned(np.where(ok_pos[0], ef, 0.0), np.nansum), binned(np.where(ok_pos[1], er, 0.0), np.nansum)])
    gcb = np.stack([binned(gf, np.nanmean), binned(gr, np.nanmean)])
    obs_u, exp_u = np.where(usable, obs, 0.0), np.where(usable, exp, 0.0)

    def ratio(mask):
        o, e = obs_u[mask].sum(), exp_u[mask].sum()
        return (2.0 * o / e if e > 0 else float("nan")), int(o)

    npos = np.stack([binned(ok_pos[s].astype(float), np.nansum) for s in (0, 1)])
    flat_rate = curve.scale * curve.rescale

    def ratio_flat(mask):
        """The same windows with no GC model: what a plain depth ratio reports."""
        o, n = obs_u[mask].sum(), np.where(usable, npos, 0.0)[mask].sum()
        return 2.0 * o / (n * flat_rate) if n > 0 else float("nan")

    all_mask = usable
    anchor_mask = usable & (gcb >= ANCHOR_GC[0]) & (gcb <= ANCHOR_GC[1])
    if anchors:
        inside = np.zeros(nb, bool)
        for s0, e0 in anchors:
            inside[(s0 + b - 1) // b: e0 // b] = True
        anchor_mask &= inside[None, :]
    cn_all, n_all = ratio(all_mask)
    cn_anchor, n_anchor = ratio(anchor_mask)
    # windows
    per = max(1, window // b)
    wins = []
    for w0 in range(0, nb, per):
        sl = slice(w0, min(nb, w0 + per))
        u = usable[:, sl]
        o, e = obs_u[:, sl].sum(), exp_u[:, sl].sum()
        frac = float(u.mean())
        g = float(np.nanmean(np.where(u, gcb[:, sl], np.nan))) if u.any() else float(np.nanmean(gcb[:, sl]))
        wins.append(dict(start=w0 * b, end=min(U, (w0 + per) * b), gc=round(g, 4), obs=int(o), exp=round(float(e), 3),
                         usable=round(frac, 3), cn=(round(2.0 * o / e, 3) if e > 0 and frac >= 0.5 else None)))
    cns = np.array([w["cn"] for w in wins if w["cn"] is not None], float)
    # headline: anchor windows when the unit has them, otherwise every usable window
    has_anchor = n_anchor >= 1000
    out = dict(kind="positional", length=U, reads=cls["reads"], usable_fraction=round(float(usable.mean()), 4),
               cn=cn_anchor if has_anchor else cn_all, cn_basis="anchor" if has_anchor else "all",
               cn_all=cn_all, n_all=n_all, cn_anchor=cn_anchor, n_anchor=n_anchor, cn_all_flat=ratio_flat(all_mask),
               cn_median=float(np.median(cns)) if len(cns) else float("nan"),
               window_log_sd=float(np.std(np.log(cns[cns > 0]))) if (cns > 0).sum() > 2 else float("nan"),
               windows=wins, features={})
    for name, s, e in features or []:
        m = np.zeros_like(usable)
        m[:, s // b:(e + b - 1) // b] = True
        # bins straddling a feature edge are kept only when fully inside it
        m[:, :(s + b - 1) // b] = False
        m[:, e // b:] = False
        v, n = ratio(usable & m)
        out["features"][name] = dict(cn=v, cn_flat=ratio_flat(usable & m), n=n, gc=float(np.nanmean(np.where(m, gcb, np.nan))) if m.any() else float("nan"))
    return out


def estimate_compositional(cls: dict, curve_read: gcmodel.GCCurve) -> dict:
    """Diploid array mass: M = sum_g T[g] / lambda_R(g), reads binned by their own GC."""
    T = np.array(cls["gc_read"], float)
    lam = curve_read.rate * curve_read.rescale
    ok = ~np.isnan(lam) & (T > 0)
    covered = T[ok].sum() / max(T.sum(), 1.0)
    # no reads at all is a mass of zero; reads, none of them at a GC the curve supports, is no estimate
    mass = float(np.sum(T[ok] / lam[ok]) / max(covered, 1e-9)) if ok.any() else (0.0 if T.sum() == 0 else float("nan"))
    return dict(kind="compositional", reads=cls["reads"], mass_bp=mass, mass_Mb=mass / 1e6,
                gc_supported_fraction=round(float(covered), 4),
                mean_read_gc=float(np.sum(T * np.arange(101)) / max(T.sum(), 1.0) / 100.0))


def check_build(counts: dict, contig_lengths: dict[str, int] | None):
    """Refuse counts made from a file aligned to another reference build. Contig *names* do not
    tell builds apart (hg19 and GRCh38 both have a chr1); the control regions and sinks are
    coordinates, and on the wrong build they are coordinates of something else."""
    have = {c["name"]: c["len"] for c in counts.get("contigs", [])}
    bad = [(n, have[n], ln) for n, ln in (contig_lengths or {}).items() if n in have and have[n] != ln]
    if bad:
        n, got, want = bad[0]
        raise ValueError(f"{counts.get('sample')}: the alignment header gives {n} a length of {got:,}, the resource bundle was built "
                         f"on {want:,} ({len(bad)} contigs differ): this file is aligned to a different reference build")


def align_region_tables(counts: dict, region_tables: np.ndarray, regions: list[tuple[str, str]]) -> np.ndarray:
    """Rows of the bundle's region tables in the order of the counts file's regions, matched by
    name. The counts must hold exactly the bundle's `control` regions - they are the denominator
    and the GC curve, and a cohort is only one cohort if they are the same everywhere - and no
    region the bundle does not know. Other bundle regions may be missing: a known-truth or dosage
    set added to a bundle later does not invalidate counts made before it existed (the one layer
    that cannot be redone without re-reading the alignments); it is simply not reported."""
    row = {name: i for i, (name, _) in enumerate(regions)}
    have = [r["name"] for r in counts["regions"]]
    unknown = [n for n in have if n not in row]
    ctrl_bundle = {n for n, role in regions if role == "control"}
    ctrl_counts = {r["name"] for r in counts["regions"] if r.get("role", "control") == "control"}
    if unknown or ctrl_bundle != ctrl_counts:
        raise ValueError(f"{counts.get('sample')}: these counts were made with a different controls file ({len(unknown)} regions unknown to the "
                         f"bundle, {len(ctrl_bundle ^ ctrl_counts)} control regions not shared)")
    return region_tables[[row[n] for n in have]]


def estimate_sample(counts: dict, panel: Panel, units: dict[str, str], features: dict | None = None,
                    region_tables: np.ndarray | None = None, L: int | None = None, window: int = 250,
                    min_kmers: int = 20, anchors: dict | None = None, contig_lengths: dict | None = None,
                    regions: list[tuple[str, str]] | None = None) -> dict:
    """`regions`: the bundle's (name, role) list in the row order of `region_tables`; with it the
    tables are matched to the counts by name, without it they must correspond row for row."""
    check_build(counts, contig_lengths)
    if region_tables is not None and regions is not None:
        region_tables = align_region_tables(counts, region_tables, regions)
    tab = nearest_table(counts, L)
    curve = gcmodel.fit_gc_curve(tab["n"], tab["o"], tab["l"])
    qc = control_qc(counts, curve, region_tables)
    if qc is not None:
        curve.rescale = qc.rescale
    R = counts["read_length_mode"]
    tab_r = nearest_table(counts, R)
    curve_r = gcmodel.fit_gc_curve(tab_r["n"], tab_r["o"], tab_r["l"])
    if qc is not None:
        curve_r.rescale = qc.rescale
    ctrl = max(counts["ctrl_reads"], 1)
    res = dict(
        sample=counts["sample"], mode=counts["mode"], engine_version=counts["engine_version"], engine_build=counts.get("engine_build"),
        read_length=R, insert_median=counts["insert_median"], gc_L=curve.L,
        ctrl_reads=counts["ctrl_reads"], ctrl_rate=curve.scale * curve.rescale,
        # single-copy depth implied by the 5'-end rate: ends per position-strand x 2 strands x read length
        depth_equiv=curve.scale * curve.rescale * 2 * R,
        ctrl_dup_frac=counts["ctrl_dup_flagged"] / ctrl, ctrl_mapq0_frac=counts["ctrl_mapq0"] / ctrl,
        gc_support=(curve.lo / 100.0, curve.hi / 100.0), gc_dispersion=curve.dispersion,
        gc_curve_max_se=float(np.nanmax(curve.se_log)),
        gc_rel={f"{g}": (round(float(curve.relative()[g]), 4) if not np.isnan(curve.rate[g]) else None) for g in range(20, 85, 5)},
        truth_regions=test_regions(counts, curve, region_tables),
        eof_marker=counts.get("eof_marker"),
        control_qc=None if qc is None else dict(rescale=qc.rescale, n_regions=qc.n_regions, n_flagged_regions=qc.n_flagged_regions,
                                                flagged_chromosomes=qc.flagged_chromosomes, region_log_mad_sd=qc.region_log_sd,
                                                # log(observed / expected) per control region: what is left after the GC
                                                # model. Across a cohort its leading components are technical and
                                                # biological covariates (library, replication timing in cycling cells)
                                                # measured on sequence disjoint from every class.
                                                region_log_ratio=[round(float(x), 4) for x in np.log(np.maximum(qc.ratios, 1e-6))]),
        classes={},
    )
    for cls in counts["classes"]:
        pc = panel.classes.get(cls["name"])
        if cls["kind"] == "positional":
            # needs the class's k-mer positions (what is callable) and its unit sequence (expected GC)
            if pc is None or cls["name"] not in units:
                continue
            r = estimate_positional(cls, pc, units[cls["name"]], panel.k, curve, R, window, min_kmers,
                                    (features or {}).get(cls["name"]), (anchors or {}).get(cls["name"]))
        else:
            r = estimate_compositional(cls, curve_r)
        r["dup_flag_frac"] = cls["dup_flagged"] / max(cls["reads"], 1)
        res["classes"][cls["name"]] = r
    return res


def control_region_tables(controls_fasta, L: int) -> np.ndarray:
    """N_r[g] for every control region (rows in file order): position-strands per 1% GC bin."""
    from .io import _open
    rows, name, buf, flank = [], None, [], 0

    def flush():
        if name is None:
            return
        seq = "".join(buf)
        a = np.frombuffer(seq.encode(), dtype=np.uint8)
        isgc = (a == 71) | (a == 67)
        isn = ~((a == 65) | (a == 67) | (a == 71) | (a == 84))
        cg, cn = np.r_[0, np.cumsum(isgc)], np.r_[0, np.cumsum(isn)]
        m = len(seq) - 2 * flank
        p = np.arange(m) + flank
        h = np.zeros(101, np.int64)
        for lo_, hi_ in ((p, p + L), (p + 1 - L, p + 1)):
            g = cg[hi_] - cg[lo_]
            ok = (cn[hi_] - cn[lo_]) == 0
            h += np.bincount((g[ok] * 100 + L // 2) // L, minlength=101)
        rows.append(h)

    with _open(controls_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                name, buf = line[1:].split()[0], []
                flank = int(re.search(r"flank=(\d+)", line).group(1))
            else:
                buf.append(line.strip().upper())
    flush()
    return np.array(rows)
