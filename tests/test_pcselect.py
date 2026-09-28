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
    sv = np.array([float(line.split()[1]) for line in open(ROOT / "tests/data/ngspca_1000G.singularvalues.txt").read().splitlines()[1:]])
    n_bins = 142_070                       # the bins of NGS-PCA's 1000 Genomes run (svd.bins.txt, kept in NGS-DOSE-1000G/meta/ngspca)
    counts = [pcselect.mp_select(sv[:k], 3200, n_bins).n_pc for k in (100, 120, 150, 200)]
    assert counts == [38, 40, 40, 40], counts
    bare = [pcselect.mp_select(sv[:k], 3200, n_bins, margin=0).n_pc for k in (100, 120, 150, 200)]
    assert max(bare) - min(bare) <= 6 and 44 <= min(bare) and max(bare) <= 50, bare          # 45-49 without the margin
    assert len({pcselect.mp_select(sv, 3200, n_bins, z=z).n_pc for z in (3, 4, 5, 6)}) == 1   # the residual-SD multiple is immaterial


def _noise_spectrum(rng, n, p, column_log_sd):
    """Singular values of noise whose SD differs between rows by +-40% and between columns log-normally."""
    X = rng.standard_normal((n, p), dtype=np.float32) * rng.uniform(0.6, 1.4, (n, 1)).astype(np.float32)
    if column_log_sd:
        X *= np.exp(rng.normal(0, column_log_sd, (1, p))).astype(np.float32)
    X -= X.mean(0)
    return np.sqrt(np.clip(np.linalg.eigvalsh((X @ X.T).astype(np.float64))[::-1], 0, None))


def test_bins_with_heavy_tailed_variance_make_the_count_lean_high():
    """A bin whose noise variance is several times the typical one is a component of its own: structure
    in the matrix, of no interest. With log-normal bin SDs (sigma 0.3) noise alone yields a few
    components (about three per run without the margin, about one with it); with equal bins, none.
    This test documents the size of that lean - the reason the count is a default for `pcsweep` to
    confirm or overrule, not a measurement."""
    rng = np.random.default_rng(99)
    n, p = 1500, 12000
    equal = [pcselect.mp_select(_noise_spectrum(rng, n, p, 0.0)[:200], n, p).n_pc for _ in range(2)]
    assert equal == [0, 0], equal
    with_margin, bare = [], []
    for _ in range(4):
        sv = _noise_spectrum(rng, n, p, 0.3)
        with_margin.append(pcselect.mp_select(sv[:200], n, p).n_pc)
        bare.append(pcselect.mp_select(sv[:200], n, p, margin=0).n_pc)
    assert all(a <= b for a, b in zip(with_margin, bare)), (with_margin, bare)
    assert 1 <= np.mean(bare) <= 6 and np.mean(with_margin) <= 3 and max(with_margin) <= 5, (with_margin, bare)


def test_the_margin_costs_no_component_that_was_planted_just_above_the_noise():
    rng = np.random.default_rng(12)
    n, p = 1500, 12000
    X = rng.standard_normal((n, p), dtype=np.float32) * rng.uniform(0.6, 1.4, (n, 1)).astype(np.float32)
    top = np.sqrt(np.linalg.eigvalsh(((X - X.mean(0)) @ (X - X.mean(0)).T).astype(np.float64))[-1])
    for strength in (1.15, 1.10, 1.06, 1.03):                 # 3-15% above the largest noise value, before the noise inflates them
        X += (strength * top * (rng.standard_normal(n).astype(np.float32) / np.sqrt(n))[:, None]) * (rng.standard_normal(p).astype(np.float32) / np.sqrt(p))[None, :]
    X -= X.mean(0)
    sv = np.sqrt(np.clip(np.linalg.eigvalsh((X @ X.T).astype(np.float64))[::-1], 0, None))
    assert [pcselect.mp_select(sv[:200], n, p, margin=m).n_pc for m in (0.0, 0.01, 0.02)] == [4, 4, 4]


def test_the_one_standard_error_band_uses_the_standard_error_of_a_mad():
    """The MAD is 37% efficient: a MAD-based SD has 1.65 times the standard error of a sample SD."""
    rng = np.random.default_rng(3)
    n = 400
    x = rng.normal(0, 0.02, (3000, n))
    empirical = (1.4826 * np.median(np.abs(x - np.median(x, axis=1, keepdims=True)), axis=1)).std()
    rows = [dict(column="truth.auto", kind="truth", n_pc=k, n=n, sd_log_robust=0.02, rmse_log=0.02) for k in range(3)]
    assert abs(pcselect.recommend(rows)["truth.auto"]["se"] / empirical - 1) < 0.08


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


def test_a_truth_that_differs_between_samples_is_not_mistaken_for_error():
    """chrX is 1 in males and 2 in females. The sweep regresses the error, log(estimate / truth), on the
    PCs - not the estimate: the estimate carries SD(log truth) = 0.35 of variance that is no error, and
    regressing it on k irrelevant PCs adds 0.35 * sqrt(k / n) of estimation noise out of fold, which
    read as adjustment making chrX several times worse with every PC (0.010 -> 0.040 at 46 PCs here)."""
    rng = np.random.default_rng(1)
    n, K = 3200, 60
    P = rng.normal(size=(n, K))
    truth = np.where(rng.random(n) < 0.5, 1.0, 2.0)
    table = {"truth.chrX": truth * np.exp(rng.normal(0, 0.01, n))}
    rows = pcselect.sweep(table, P, K, {"truth.chrX": truth}, [])
    e = np.array([r["sd_log_robust"] for r in rows])
    assert abs(e[0] - 0.01) < 0.001 and e.max() < 1.05 * e[0], e[[0, 10, 20, 46, 60]]
    assert pcselect.recommend(rows)["truth.chrX"]["pick"] == 0
    # ... and when a PC does follow sex, the truth is not what it explains
    P[:, 0] = np.log(truth) + rng.normal(0, 0.05, n)
    e_sex = np.array([r["sd_log_robust"] for r in pcselect.sweep(table, P, 5, {"truth.chrX": truth}, [])])
    assert np.all(np.abs(e_sex / e_sex[0] - 1) < 0.05), e_sex


def test_recommend_survives_undefined_intervals_and_columns():
    """A bootstrap interval is undefined when one resample has parents all alike (a mostly constant
    column): the band is then 0 and the pick the best. A column undefined at every k is left out."""
    rows = [dict(column="x", kind="class", n_pc=k, n=66, R_midparent=r, R_lo=np.nan, R_hi=np.nan) for k, r in enumerate((0.50, 0.53, 0.51))]
    rows += [dict(column="y", kind="class", n_pc=k, n=66, R_midparent=np.nan) for k in range(3)]
    rows += [dict(column="t", kind="truth", n_pc=k, n=66, sd_log_robust=np.nan) for k in range(3)]
    rec = pcselect.recommend(rows)
    assert set(rec) == {"x"} and rec["x"]["best"] == rec["x"]["pick"] == 1 and rec["x"]["se"] == 0.0


def test_the_sweep_of_a_mostly_constant_column_gives_a_recommendation():
    """22 trios, a copy number that is 2 in all but a few samples: the case that crashed `pcsweep`."""
    trios = [Trio(f"c{i}", f"f{i}", f"m{i}") for i in range(22)]
    samples = [s for t in trios for s in (t.child, t.father, t.mother)]
    rng = np.random.default_rng(1)
    P = rng.normal(size=(len(samples), 2))
    p = rng.uniform(0.85, 0.99)
    y = np.where(rng.random(len(samples)) < p, 2.0, 3.0)
    rows = pcselect.sweep({"x": y}, P, 1, {}, ["x"], samples, trios, None, folds=5, n_boot=100, seed=1)
    rec = pcselect.recommend(rows)
    assert rec["x"]["pick"] in (0, 1) and np.isfinite(rec["x"]["se"])
