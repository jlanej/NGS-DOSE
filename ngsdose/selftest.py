"""Simulation checks with known truth; no sequencing data required.

Each check simulates the generative model the estimator assumes *plus* the nuisance it is meant
to survive, and asserts recovery of the truth. They are the executable form of the claims in
docs/DESIGN.md.
"""
from __future__ import annotations

import numpy as np

from . import cohort, estimate, gcmodel, trios
from .io import PanelClass


def _random_unit(rng, n=12000, gc_blocks=((0.45, 3000), (0.75, 2500), (0.55, 3000), (0.35, 1500), (0.62, 2000))):
    parts = []
    for gc, ln in gc_blocks:
        p = np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])
        parts.append(rng.choice(list("ACGT"), size=ln, p=p))
    seq = "".join(np.concatenate(parts))
    return seq[:n]


def _bias(g, strength=1.0):
    """A unimodal library GC bias (relative rate), as in Benjamini & Speed 2012."""
    return np.exp(-strength * ((g - 0.48) / 0.22) ** 2)


def _control_tables(rng, L, depth_rate, strength, n_pos=6_000_000):
    g = np.arange(101) / 100.0
    w = np.exp(-0.5 * ((g - 0.41) / 0.07) ** 2) + 0.002 * ((g > 0.15) & (g < 0.85))
    N = np.round(n_pos * w / w.sum())
    O = rng.poisson(N * depth_rate * _bias(g, strength))
    return N, O


def _simulate_sample(rng, seq, pc, k, R, L, cn, strength, depth_rate=0.125, window_eff=None, bin_=50):
    N, O = _control_tables(rng, L, depth_rate, strength)
    gf, gr = gcmodel.window_gc_counts(seq, L, pc.circular)
    lam_f = depth_rate * _bias(gf / L, strength)
    lam_r = depth_rate * _bias(gr / L, strength)
    eff = np.ones(len(seq)) if window_eff is None else np.repeat(window_eff, 250)[: len(seq)]
    nb = -(-len(seq) // bin_)
    pad = nb * bin_ - len(seq)
    def binsum(x):
        return np.r_[x, np.zeros(pad)].reshape(nb, bin_).sum(1)
    fwd = rng.poisson(binsum(cn / 2 * lam_f * eff))
    rev = rng.poisson(binsum(cn / 2 * lam_r * eff))
    counts = dict(sample="sim", mode="sim", engine_version="sim", insert_median=L, read_length_mode=R,
                  ctrl_reads=int(O.sum()), ctrl_dup_flagged=0, ctrl_mapq0=0, regions=[],
                  gc_tables=[dict(l=L, n=N.tolist(), o=O.tolist()), dict(l=R, n=N.tolist(), o=O.tolist())],
                  classes=[dict(name=pc.name, kind="positional", bin=bin_, fwd=fwd.tolist(), rev=rev.tolist(),
                                reads=int(fwd.sum() + rev.sum()), dup_flagged=0)])
    return counts


def check_estimator(rng, verbose):
    """Copy number is recovered under a strong GC bias that a depth ratio gets badly wrong."""
    seq = _random_unit(rng)
    pc = PanelClass("unit", "positional", len(seq), True, np.arange(len(seq)))
    from .io import Panel
    panel = Panel(31, {"unit": pc})
    truth, L, R = 400.0, 450, 150
    c = _simulate_sample(rng, seq, pc, 31, R, L, truth, strength=1.5)
    r = estimate.estimate_sample(c, panel, {"unit": seq})["classes"]["unit"]
    naive = 2 * (np.sum(c["classes"][0]["fwd"]) + np.sum(c["classes"][0]["rev"])) / (2 * len(seq) * 0.125 * np.mean(_bias(np.array([0.41]), 1.5)))
    ok = abs(r["cn_all"] / truth - 1) < 0.02 and abs(r["cn_anchor"] / truth - 1) < 0.02
    if verbose:
        print(f"[1] GC-aware estimator     : truth {truth:.0f}, cn_all {r['cn_all']:.1f}, cn_anchor {r['cn_anchor']:.1f}, "
              f"uncorrected depth ratio {naive:.1f}  {'OK' if ok else 'FAIL'}")
    return ok


def check_calibration(rng, verbose):
    """Window efficiencies shared by all samples are separated from per-sample copy number."""
    seq = _random_unit(rng)
    pc = PanelClass("unit", "positional", len(seq), True, np.arange(len(seq)))
    from .io import Panel
    panel = Panel(31, {"unit": pc})
    nwin = -(-len(seq) // 250)
    gf, gr = gcmodel.window_gc_counts(seq, 450, True)
    wgc = np.array([np.mean(np.r_[gf[i * 250:(i + 1) * 250], gr[i * 250:(i + 1) * 250]]) / 450 for i in range(nwin)])
    eff = np.where(wgc > 0.65, rng.uniform(0.45, 0.9, nwin), 1.0)        # dropout confined to GC-rich windows
    truths = rng.uniform(200, 700, 40)
    results = []
    for t in truths:
        c = _simulate_sample(rng, seq, pc, 31, 150, 450, t, strength=rng.uniform(0.3, 1.5), window_eff=eff)
        r = estimate.estimate_sample(c, panel, {"unit": seq})
        r["sample"] = f"s{len(results)}"
        results.append(r)
    cal = cohort.calibrate(results, "unit")
    est = np.exp(cal.c)
    naive = np.array([r["classes"]["unit"]["cn_all"] for r in results])
    err, err_naive = np.abs(est / truths - 1).max(), np.abs(naive / truths - 1).mean()
    a_err = np.nanmax(np.abs(np.exp(cal.a) - eff[: len(cal.a)]))
    ok = err < 0.02 and a_err < 0.03
    if verbose:
        print(f"[2] window calibration     : max |error| {100 * err:.2f}% (all-window ratio is off by {100 * err_naive:.1f}% on average); "
              f"efficiencies recovered to {a_err:.3f}  {'OK' if ok else 'FAIL'}")
    return ok


def _trio_cohort(rng, n, s_hap, s_err, shared, rho_pop=0.0):
    vals, ped, pop = {}, [], {}
    for i in range(n):
        shift = rng.normal(0, rho_pop)                                     # population/assortment term shared by spouses
        hap = rng.normal(100 + shift / 2, s_hap, 4)
        T = dict(f=hap[0] + hap[1], m=hap[2] + hap[3], c=hap[rng.integers(0, 2)] + hap[2 + rng.integers(0, 2)])
        fam = rng.normal(0, s_err * np.sqrt(shared))
        for who in "fmc":
            sid = f"t{i}{who}"
            vals[sid] = T[who] + fam + rng.normal(0, s_err * np.sqrt(1 - shared))
            pop[sid] = "P"
        ped.append(trios.Trio(f"t{i}c", f"t{i}f", f"t{i}m", "P"))
    return vals, ped, pop


def _trio_cohort_lognormal(rng, n, cv_hap, cv_err):
    """Skewed array sizes and multiplicative error - what real dosage looks like."""
    vals, ped, pop = {}, [], {}
    sig = np.sqrt(np.log(1 + cv_hap ** 2))
    for i in range(n):
        hap = 200 * rng.lognormal(-sig ** 2 / 2, sig, 4)
        T = dict(f=hap[0] + hap[1], m=hap[2] + hap[3], c=hap[rng.integers(0, 2)] + hap[2 + rng.integers(0, 2)])
        for who in "fmc":
            sid = f"t{i}{who}"
            vals[sid] = T[who] * (1 + rng.normal(0, cv_err))
            pop[sid] = "P"
        ped.append(trios.Trio(f"t{i}c", f"t{i}f", f"t{i}m", "P"))
    return vals, ped, pop


def check_trios(rng, verbose):
    """Reliability is recovered; spousal correlation is corrected for; shared error is exposed."""
    s_hap, s_err = 18.0, 12.0
    R_true = 2 * s_hap ** 2 / (2 * s_hap ** 2 + s_err ** 2)
    v, ped, pop = _trio_cohort(rng, 4000, s_hap, s_err, shared=0.0)
    t0 = trios.transmission(v, ped, pop, n_perm=50, n_boot=0)
    v, ped, pop = _trio_cohort(rng, 4000, s_hap, s_err, shared=0.0, rho_pop=20.0)
    t1 = trios.transmission(v, ped, None, n_perm=0, n_boot=0)
    R1 = 1 - s_err ** 2 / np.var([x for k, x in v.items() if not k.endswith("c")])
    v, ped, pop = _trio_cohort(rng, 4000, s_hap, s_err, shared=0.8)
    t2 = trios.transmission(v, ped, pop, n_perm=0, n_boot=0)
    # skewed (log-normal) array sizes with multiplicative error: R is then the ratio of true to
    # observed variance with the *average* error variance, and the natural-scale estimators hold
    cv_hap, cv_err = 0.35, 0.10
    v, ped, pop = _trio_cohort_lognormal(rng, 6000, cv_hap, cv_err)
    t3 = trios.transmission(v, ped, pop, n_perm=0, n_boot=0)
    var_T = 2 * (200 * cv_hap) ** 2
    mean_T2 = var_T + 400.0 ** 2
    R3 = var_T / (var_T + cv_err ** 2 * mean_T2)
    # a worse estimator of the same quantity is told apart by the paired bootstrap
    noisy = {k: x * (1 + rng.normal(0, 0.08)) for k, x in v.items()}
    cmp_ = trios.compare(v, noisy, ped[:602], pop, n_boot=300)
    ok = (abs(t0["reliability_midparent"] - R_true) < 0.03 and abs(t0["reliability_mendel"] - R_true) < 0.04
          and abs(t3["reliability_midparent"] - R3) < 0.03 and abs(t3["reliability_single_parent"] - R3) < 0.04
          and cmp_["ci95"][0] > 0
          and abs(t0["reliability_single_parent"] - R_true) < 0.04 and abs(t0["perm_null_mean"]) < 0.02
          and abs(t1["reliability_midparent"] - R1) < 0.03 and t1["midparent_slope"] > t1["reliability_midparent"] + 0.02
          and t2["spousal_r"] > 0.15)
    if verbose:
        print(f"[3] transmission           : true R {R_true:.3f}; midparent {t0['reliability_midparent']:.3f}, single-parent "
              f"{t0['reliability_single_parent']:.3f}, Mendelian {t0['reliability_mendel']:.3f}; permuted null {t0['perm_null_mean']:+.3f}")
        print(f"    assortment/structure   : raw slope {t1['midparent_slope']:.3f} -> corrected {t1['reliability_midparent']:.3f} (true {R1:.3f})")
        print(f"    skewed sizes, mult. err: midparent {t3['reliability_midparent']:.3f}, single-parent {t3['reliability_single_parent']:.3f} (true {R3:.3f}); "
              f"paired bootstrap, 602 trios: better estimator wins by {cmp_['delta']:+.3f} CI95 ({cmp_['ci95'][0]:+.3f}, {cmp_['ci95'][1]:+.3f})")
        print(f"    family-shared error    : spousal r {t2['spousal_r']:+.3f} flags it (apparent R {t2['reliability_midparent']:.3f} "
              f"vs true {R_true:.3f})  {'OK' if ok else 'FAIL'}")
    return ok


def check_callable_mask(rng, verbose):
    """Positions whose reads cannot be recovered are removed from numerator and denominator alike."""
    seq = _random_unit(rng)
    U = len(seq)
    blind = np.zeros(U, bool)
    blind[4000:4600] = True                                                # e.g. a pseudogene copy elsewhere in the genome
    pc = PanelClass("unit", "positional", U, True, np.where(~blind)[0])
    from .io import Panel
    panel = Panel(31, {"unit": pc})
    c = _simulate_sample(rng, seq, pc, 31, 150, 450, 300.0, strength=0.5)
    cf, cr = estimate.callable_masks(pc, 31, 150, 20)
    # reads that start where too few panel k-mers remain are never assigned by the engine
    for strand, mask in (("fwd", cf), ("rev", cr)):
        x = np.array(c["classes"][0][strand], float)
        keep = np.r_[mask, np.ones(len(x) * 50 - U, bool)].reshape(len(x), 50).mean(1)
        c["classes"][0][strand] = rng.binomial(x.astype(int), keep).tolist()
    r = estimate.estimate_sample(c, panel, {"unit": seq})["classes"]["unit"]
    ok = abs(r["cn_all"] / 300.0 - 1) < 0.02 and r["usable_fraction"] < 0.99
    if verbose:
        print(f"[4] recoverability mask    : truth 300, estimate {r['cn_all']:.1f} with {100 * (1 - r['usable_fraction']):.1f}% of the unit masked  "
              f"{'OK' if ok else 'FAIL'}")
    return ok


def check_pc_selection(rng, verbose):
    """The number of PCs to regress out: components planted in noise of unequal variance are counted
    exactly at the fitted edge of the noise bulk, none are found in noise alone, and a sweep against a
    known truth and against transmission stops where the PCs stop carrying technical error."""
    from . import pcselect
    from .trios import Trio
    n, p, m = 500, 800, 6
    noise = lambda: rng.normal(size=(n, p)) * rng.uniform(0.6, 1.4, (n, 1)) * rng.uniform(0.7, 1.3, (1, p))
    X = noise()
    for j in range(m):
        X += (2.5 - 0.1 * j) * (np.sqrt(n) + np.sqrt(p)) * np.outer(rng.normal(size=n) / np.sqrt(n), rng.normal(size=p) / np.sqrt(p))
    sv = np.linalg.svd(X - X.mean(0), compute_uv=False)
    found, found_top, in_noise = (pcselect.mp_select(sv, n, p).n_pc, pcselect.mp_select(sv[:60], n, p).n_pc,
                                  pcselect.mp_select(np.linalg.svd(noise(), compute_uv=False), n, p).n_pc)
    # 150 trios; three of fifteen PCs carry technical error into a class and into a known truth
    nt = 150
    P = rng.normal(size=(3 * nt, 15))
    tech = P[:, :3] @ np.array([0.12, 0.10, 0.08])
    f, mo = rng.normal(400, 80, nt), rng.normal(400, 80, nt)
    true = np.empty(3 * nt)
    true[0::3], true[1::3], true[2::3] = (f + mo) / 2 + rng.normal(0, 80 / np.sqrt(2), nt), f, mo
    names = [f"S{i}" for i in range(3 * nt)]
    table = {"cls": true * np.exp(tech + rng.normal(0, 0.01, 3 * nt)), "truth": 2.0 * np.exp(0.3 * tech + rng.normal(0, 0.004, 3 * nt))}
    rows = pcselect.sweep(table, P, 15, {"truth": np.full(3 * nt, 2.0)}, ["cls"], names,
                          [Trio(names[3 * t], names[3 * t + 1], names[3 * t + 2], "P") for t in range(nt)], {s: "P" for s in names}, n_boot=100)
    rec = pcselect.recommend(rows)
    ok = found == m and found_top == m and in_noise <= 1 and rec["truth"]["pick"] == 3 and rec["cls"]["pick"] in (2, 3)
    if verbose:
        print(f"[5] number of PCs          : {m} planted in unequal noise -> {found} (whole spectrum), {found_top} (top 60), {in_noise} in noise alone; "
              f"sweep picks {rec['truth']['pick']} by known truth, {rec['cls']['pick']} by transmission (3 PCs carry the error)  {'OK' if ok else 'FAIL'}")
    return ok


def run(verbose=True, seed=11) -> bool:
    rng = np.random.default_rng(seed)
    checks = (check_estimator, check_calibration, check_trios, check_callable_mask, check_pc_selection)
    ok = all([f(rng, verbose) for f in checks])
    if verbose:
        print("\nall checks passed" if ok else "\nFAILURES - see above")
    return ok
