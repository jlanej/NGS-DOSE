"""How many PCs to regress out: the Marchenko-Pastur default, and the sweep that can overrule it."""
from pathlib import Path

import numpy as np

from ngsdose import pcselect
from ngsdose.trios import Trio

ROOT = Path(__file__).resolve().parents[1]


def planted(rng, n, p, m, strength=2.5, unequal=True):
    """m components, the weakest still well above the noise edge, in noise whose SD differs between
    rows by +-40% (samples of different depth) and between columns by +-30% (bins of different mappability)."""
    X = rng.normal(size=(n, p))
    if unequal:
        X *= rng.uniform(0.6, 1.4, (n, 1)) * rng.uniform(0.7, 1.3, (1, p))
    edge = np.sqrt(n) + np.sqrt(p)
    for j in range(m):
        X += (strength - 0.1 * j) * edge * np.outer(rng.normal(size=n) / np.sqrt(n), rng.normal(size=p) / np.sqrt(p))
    X -= X.mean(0)
    return np.linalg.svd(X, compute_uv=False)


def test_the_edge_finds_planted_components_in_unequal_noise():
    rng = np.random.default_rng(5)
    for n, p, m in ((600, 800, 7), (2000, 800, 12), (400, 20000, 5)):
        sv = planted(rng, n, p, m)
        whole, leading = pcselect.mp_select(sv, n, p), pcselect.mp_select(sv[:max(60, 5 * m)], n, p)
        assert whole.n_pc == m and leading.n_pc == m, (n, p, m, whole.describe(), leading.describe())
        assert whole.slope_ratio > 1.0                       # and it says the noise is not identically distributed


def test_pure_noise_has_no_components():
    rng = np.random.default_rng(6)
    counts = [pcselect.mp_select(planted(rng, n, 800, 0), n, 800).n_pc for n in (60, 60, 200, 200, 600, 600, 1500)]
    assert sum(c == 0 for c in counts) >= len(counts) - 1 and max(counts) <= 2, counts
    few = pcselect.mp_select(planted(rng, 12, 800, 0), 12, 800)
    assert few.n_pc == 0 and "too few" in few.describe()     # a dozen samples have no noise bulk to speak of


def test_the_count_on_the_1000_genomes_spectrum_does_not_depend_on_how_much_of_it_was_kept():
    """NGS-PCA keeps the top 200 singular values of 3,200. The textbook fit of the law to that tail
    gives 59, 66 and 80 components for the top 100, 150 and 200; the edge fit must not wander so."""
    sv = np.array([float(line.split()[1]) for line in open(ROOT / "example/1000G/ngspca/svd.singularvalues.txt").read().splitlines()[1:]])
    n_bins = sum(1 for _ in open(ROOT / "example/1000G/ngspca/svd.bins.txt"))
    counts = [pcselect.mp_select(sv[:k], 3200, n_bins).n_pc for k in (100, 120, 150, 200)]
    assert max(counts) - min(counts) <= 6 and 40 <= min(counts) and max(counts) <= 55, counts


def cohort_with_technical_factors(rng, n_trios=200, n_factors=3, n_pcs=20):
    """Trios whose estimates carry real, transmitted dosage plus error driven by the first
    `n_factors` of `n_pcs` coverage PCs; a known-truth column that carries the same error."""
    n = 3 * n_trios
    P = rng.normal(size=(n, n_pcs))
    tech = P[:, :n_factors] @ np.array([0.12, 0.10, 0.08][:n_factors])          # each one costs the class several points of reliability
    f, m = rng.normal(400, 80, n_trios), rng.normal(400, 80, n_trios)
    c = (f + m) / 2 + rng.normal(0, 80 / np.sqrt(2), n_trios)
    true = np.empty(n)
    true[0::3], true[1::3], true[2::3] = c, f, m
    names = [f"S{i}" for i in range(n)]
    table = {"cls": true * np.exp(tech + rng.normal(0, 0.01, n)), "truth.auto": 2.0 * np.exp(0.3 * tech + rng.normal(0, 0.004, n))}
    trios = [Trio(names[3 * t], names[3 * t + 1], names[3 * t + 2], "POP") for t in range(n_trios)]
    return table, P, names, trios


def test_the_sweep_shows_where_adjustment_stops_helping():
    rng = np.random.default_rng(8)
    table, P, names, trios = cohort_with_technical_factors(rng)
    rows = pcselect.sweep(table, P, 20, {"truth.auto": np.full(len(names), 2.0)}, ["cls"], names, trios, {s: "POP" for s in names}, n_boot=100)
    at = lambda col, k, key: next(r[key] for r in rows if r["column"] == col and r["n_pc"] == k)
    # known truth: the error falls while the PCs carry its cause, then stops - cross-validated, it does not keep falling
    assert at("truth.auto", 3, "sd_log_robust") < 0.5 * at("truth.auto", 0, "sd_log_robust")
    assert at("truth.auto", 20, "sd_log_robust") > 0.97 * at("truth.auto", 3, "sd_log_robust")
    # the class: more of what is left is transmitted, and the paired interval against no adjustment excludes zero
    assert at("cls", 3, "R_midparent") > at("cls", 0, "R_midparent") + 0.05 and at("cls", 3, "dR_lo") > 0
    rec = pcselect.recommend(rows)
    assert rec["truth.auto"]["pick"] == 3 and rec["cls"]["pick"] == 3, rec


def test_cross_validation_is_what_makes_the_sweep_honest():
    """Residual variance falls with every regressor, meaningful or not; out-of-fold residuals do not."""
    rng = np.random.default_rng(9)
    n, k = 120, 40
    y, P = rng.normal(size=n), rng.normal(size=(n, k))                       # nothing to find
    cv = pcselect.cv_adjusted(y, P, k)
    insample = [np.std(y - np.column_stack([np.ones(n), P[:, :j]]) @ np.linalg.lstsq(np.column_stack([np.ones(n), P[:, :j]]), y, rcond=None)[0]) for j in (0, k)]
    assert insample[1] < 0.85 * insample[0]                                  # looks like a third of the variance explained
    assert np.std(cv[k]) > 1.1 * np.std(cv[0])                               # is in fact worse than doing nothing
    assert np.allclose(cv[0], y) and np.isfinite(cv).all()
