"""Unit tests of the modelling layer (no engine, no data)."""
import numpy as np
import pytest

from ngsdose import cohort, estimate, gcmodel, selftest, trios
from ngsdose.io import PanelClass


def test_selftest_passes():
    assert selftest.run(verbose=False)


def test_window_gc_counts_are_strand_symmetric_and_circular():
    seq = "GGGGGCCCCC" + "AAAAATTTTT" * 2
    fwd, rev = gcmodel.window_gc_counts(seq, 10, circular=True)
    assert fwd[0] == 10 and rev[9] == 10          # the same window seen from either fragment end
    assert fwd[25] == 5                           # [25, 35) wraps over the unit junction: TTTTT + GGGGG
    assert rev[4] == 5                            # [4-9, 4] wraps backwards over the junction
    lin_f, lin_r = gcmodel.window_gc_counts(seq, 10, circular=False)
    assert lin_f[-1] == -1 and lin_r[0] == -1     # undefined at the ends of a linear sequence
    assert (lin_f[:21] == fwd[:21]).all()


def test_gc_curve_recovers_rate_and_bounds_support():
    rng = np.random.default_rng(0)
    g = np.arange(101) / 100
    N = np.round(5e6 * np.exp(-0.5 * ((g - 0.41) / 0.06) ** 2))
    true = 0.12 * np.exp(-((g - 0.5) / 0.25) ** 2)
    O = rng.poisson(N * true)
    c = gcmodel.fit_gc_curve(N, O, 450)
    inside = ~np.isnan(c.rate)
    assert np.abs(c.rate[inside] / true[inside] - 1).max() < 0.08
    assert 0.15 <= c.lo / 100 < 0.35 and 0.5 < c.hi / 100 < 0.7    # no extrapolation beyond the data
    assert np.isnan(c.rate[95])
    with pytest.raises(ValueError):
        gcmodel.fit_gc_curve(N, np.zeros(101), 450)


def test_callable_mask_matches_bruteforce():
    rng = np.random.default_rng(1)
    U, k, R = 400, 31, 100
    pos = np.sort(rng.choice(U, 250, replace=False))
    pc = PanelClass("u", "positional", U, True, pos)
    f, r = estimate.callable_masks(pc, k, R, 40)
    has = np.zeros(U, int)
    has[pos] = 1
    for p in (0, 17, 199, 399):
        nf = sum(has[(p + j) % U] for j in range(R - k + 1))                 # k-mer starts p .. p+R-k
        nr = sum(has[(p - R + 1 + j) % U] for j in range(R - k + 1))         # k-mer starts x-R+1 .. x-k+1
        assert f[p] == (nf >= 40) and r[p] == (nr >= 40)


def test_reliability_algebra():
    rng = np.random.default_rng(2)
    vals, ped, pop = selftest._trio_cohort(rng, 6000, 18.0, 12.0, shared=0.0)
    t = trios.transmission(vals, ped, pop, n_perm=0)
    R = 2 * 18.0 ** 2 / (2 * 18.0 ** 2 + 12.0 ** 2)
    assert abs(t["reliability_midparent"] - R) < 0.03
    assert abs(t["reliability_mendel"] - R) < 0.03
    assert abs(t["spousal_r"]) < 0.05


def test_pedigree_parser(tmp_path):
    p = tmp_path / "ped.txt"
    p.write_text("FamilyID SampleID FatherID MotherID Sex Population Superpopulation\n"
                 "F1 kid dad mum 1 CEU EUR\nF1 dad 0 0 1 CEU EUR\nF1 mum 0 0 2 CEU EUR\n")
    ped, pop = trios.load_pedigree(p)
    assert [(t.child, t.father, t.mother) for t in ped] == [("kid", "dad", "mum")]
    assert pop["mum"] == "CEU"


def test_calibration_with_fixed_efficiencies_needs_no_cohort():
    rng = np.random.default_rng(3)
    a = np.where(np.arange(40) % 4 == 0, np.log(0.6), 0.0)
    gc = np.where(np.arange(40) % 4 == 0, 0.75, 0.5)

    def fake(cn, name):
        wins = [dict(start=i * 250, end=(i + 1) * 250, gc=float(gc[i]), obs=1, exp=1.0, usable=1.0,
                     cn=float(cn * np.exp(a[i] + rng.normal(0, 0.01)))) for i in range(40)]
        return dict(sample=name, classes={"u": dict(windows=wins)})

    cal = cohort.calibrate([fake(300, "one")], "u", a_fixed=a)
    assert abs(np.exp(cal.c[0]) / 300 - 1) < 0.01


def test_bundle_anchor_intervals_override_the_gc_rule():
    """Anchors pin the absolute level: efficiencies are expressed relative to them."""
    rng = np.random.default_rng(4)
    nwin = 40
    a = np.zeros(nwin)
    a[10:20] = np.log(0.8)                         # a dropout zone that the GC rule cannot see
    gc = np.full(nwin, 0.5)

    def fake(cn, name):
        wins = [dict(start=i * 250, end=(i + 1) * 250, gc=float(gc[i]), obs=1, exp=1.0, usable=1.0,
                     cn=float(cn * np.exp(a[i] + rng.normal(0, 0.005)))) for i in range(nwin)]
        return dict(sample=name, classes={"u": dict(windows=wins)})

    res = [fake(cn, f"s{i}") for i, cn in enumerate((200, 300, 400, 500, 600))]
    clean = cohort.calibrate(res, "u", anchors=[(0, 2500), (5000, 10000)])
    assert clean.anchor.sum() == 30 and not clean.anchor[10:20].any()
    assert np.allclose(np.exp(clean.c), (200, 300, 400, 500, 600), rtol=0.01)
    # anchored inside the dropout zone instead, every estimate is scaled by the zone's efficiency
    biased = cohort.calibrate(res, "u", anchors=[(2500, 5000)])
    assert np.allclose(np.exp(biased.c), 0.8 * np.array((200, 300, 400, 500, 600)), rtol=0.01)


def test_usable_windows_do_not_depend_on_depth():
    """One sample at 30x and at 0.5x must be averaged over the same windows. With dropout confined
    to GC-rich windows, a GC support that narrows at low depth would raise the all-window estimate."""
    from ngsdose.io import Panel
    rng = np.random.default_rng(9)
    seq = selftest._random_unit(rng)
    pc = PanelClass("unit", "positional", len(seq), True, np.arange(len(seq)))
    panel = Panel(31, {"unit": pc})
    nwin = -(-len(seq) // 250)
    gf, gr = gcmodel.window_gc_counts(seq, 450, True)
    wgc = np.array([np.mean(np.r_[gf[i * 250:(i + 1) * 250], gr[i * 250:(i + 1) * 250]]) / 450 for i in range(nwin)])
    eff = np.where(wgc > 0.6, 0.6, 1.0)
    out = {}
    for depth_rate in (0.125, 0.002):
        v, usable = [], set()
        for _ in range(12):
            c = selftest._simulate_sample(rng, seq, pc, 31, 150, 450, 400.0, strength=0.6, depth_rate=depth_rate, window_eff=eff)
            r = estimate.estimate_sample(c, panel, {"unit": seq})["classes"]["unit"]
            v.append(r["cn_all"])
            usable.add(r["usable_fraction"])
        out[depth_rate] = (np.mean(v), np.std(v, ddof=1) / np.sqrt(len(v)), usable)
    hi, lo = out[0.125], out[0.002]
    assert hi[2] == lo[2] and len(hi[2]) == 1, (hi[2], lo[2])          # identical window sets
    assert abs(lo[0] / hi[0] - 1) < 3 * np.hypot(lo[1], hi[1]) / hi[0] + 0.005, (hi, lo)


def test_an_N_inside_a_region_does_not_bias_its_copy_number():
    """A position whose fragment window runs over an N is in no GC table, but a read starting there
    is in the region's observed count (GRCh38's chrM has such an N at 3,107: the first chrM
    estimate was 18% high). Expected counts cover all position-strands of the region."""
    rng = np.random.default_rng(4)
    L, n_reg, length, rate, copies = 300, 6, 4000, 0.1, 40.0
    g = np.arange(101) / 100
    lam = rate * np.exp(1.5 * (g - 0.45))
    N = np.zeros(101)
    tables = np.zeros((n_reg + 1, 101))
    regions = []
    for i in range(n_reg):                                      # controls: every window defined
        h = np.bincount(np.clip(rng.normal(45, 6, 2 * 50000).astype(int), 20, 70), minlength=101)
        tables[i] = h
        N += h
        regions.append(dict(name=f"chr1:{i * 100000}-{i * 100000 + 50000}", role="control", label="", len=50000, obs=int(rng.poisson(h @ lam))))
    h = np.bincount(np.clip(rng.normal(44, 3, 2 * length).astype(int), 20, 70), minlength=101)
    undefined = 2 * L                                           # one N: L windows lost on each strand
    keep = h.astype(float) * (2 * length - undefined) / (2 * length)
    tables[n_reg] = keep
    regions.append(dict(name="chrM:1000-5000", role="dosage", label="chrM", len=length, obs=int(rng.poisson(copies / 2 * (h @ lam)))))
    O = np.array([rng.poisson(N[b] * lam[b]) for b in range(101)], float)
    counts = dict(sample="s", regions=regions)
    curve = gcmodel.fit_gc_curve(N, O, L, min_positions=500)
    t = estimate.test_regions(counts, curve, tables)["chrM"]
    assert t["role"] == "dosage" and abs(t["cn"] / copies - 1) < 0.03, t
    # and the control QC's expectation is untouched when nothing is undefined
    assert np.allclose(estimate._undefined_window_factor(regions[:n_reg], tables[:n_reg]), 1.0)
