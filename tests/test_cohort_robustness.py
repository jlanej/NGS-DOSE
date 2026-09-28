"""The cohort layer on inputs a real cohort produces sooner or later: an estimate without one of
the positional classes, a sample with no reads of a class, the same sample twice, a sample
without control residuals, an efficiency table learned on other windows. Each must either give
the other samples the numbers they would have had on their own, or stop with a message that
names what is wrong - never drop a column for everyone, or write a wrong number, without a word."""
import copy

import numpy as np
import pytest

from ngsdose import cohort, estimate, selftest
from ngsdose.io import Panel, PanelClass


def _results(n, seed=1, classes=("DJ", "rDNA45S"), n_regions=800):
    """n simulated estimates, each class a random unit at its own copy number, with control residuals."""
    rng = np.random.default_rng(seed)
    units = {c: selftest._random_unit(np.random.default_rng(100 + i), n=12000 if i == 0 else 6000) for i, c in enumerate(classes)}
    pcs = {c: PanelClass(c, "positional", len(s), True, np.arange(len(s))) for c, s in units.items()}
    panel = Panel(31, pcs)
    out = []
    for i in range(n):
        strength = rng.uniform(0.3, 1.2)
        counts = None
        for j, c in enumerate(classes):
            one = selftest._simulate_sample(rng, units[c], pcs[c], 31, 150, 450, (300.0, 400.0)[j % 2] * rng.uniform(0.8, 1.25),
                                            strength=strength)
            if counts is None:
                counts = one
            else:
                counts["classes"] += one["classes"]
        r = estimate.estimate_sample(counts, panel, units)
        r["sample"] = f"s{i:02d}"
        r["control_qc"] = dict(region_log_ratio=rng.normal(0, 0.03, n_regions).tolist(), region_log_mad_sd=0.03, flagged_chromosomes=[])
        out.append(r)
    return out


def _table(results, **kw):
    log = []
    rows, eff, info = cohort.cohort_table(copy.deepcopy(results), {}, log=log.append, **kw)
    return {r["sample"]: r for r in rows}, eff, info, "\n".join(log)


@pytest.fixture(scope="module")
def twelve():
    return _results(12)


def test_a_class_missing_from_one_estimate_is_calibrated_on_the_others(twelve):
    full, _, _, _ = _table(twelve)
    part = copy.deepcopy(twelve)
    del part[3]["classes"]["rDNA45S"]
    rows, eff, _, log = _table(part)
    alone, _, _, _ = _table([r for i, r in enumerate(twelve) if i != 3])
    assert "rDNA45S.cn" not in rows["s03"] and eff["rDNA45S"]["n_samples"] == 11
    for s, r in rows.items():
        assert r["DJ.cn"] == full[s]["DJ.cn"]                                        # the other class is untouched
        if s != "s03":
            assert r["rDNA45S.cn"] == alone[s]["rDNA45S.cn"]                         # as if the sample were not there
            assert r["rDNA45S.profilePC1"] == alone[s]["rDNA45S.profilePC1"]
    assert "WARNING: rDNA45S" in log and "11 of 12" in log and "s03 (class absent from the estimate)" in log
    part[3]["skipped_classes"] = [dict(name="rDNA45S", kind="positional", reads=0, reason="not in the bundle panel")]
    _, _, _, log = _table(part)
    assert "s03 (skipped: not in the bundle panel)" in log


def test_a_class_with_no_usable_estimate_says_so_and_skips_only_that_class(twelve):
    part = copy.deepcopy(twelve)
    for r in part:
        r["classes"]["rDNA45S"]["status"] = "no_sinks_in_fetch"
    rows, eff, _, log = _table(part)
    assert "rDNA45S" not in eff and "DJ" in eff and all("rDNA45S.cn" not in r and "DJ.cn" in r for r in rows.values())
    assert "rDNA45S: no sample has usable windows" in log and "status no_sinks_in_fetch" in log


def test_a_sample_with_no_reads_of_a_class_gets_NA_even_in_a_small_cohort():
    res = _results(7, seed=4, classes=("DJ",))
    for w in res[0]["classes"]["DJ"]["windows"]:
        w["cn"] = 0.0                                                             # what zero fwd/rev counts give
    rows, _, _, log = _table(res)
    assert "DJ.cn" not in rows["s00"] and "s00 (no usable window)" in log
    truth = np.array([r["classes"]["DJ"]["cn_all"] for r in res[1:]])
    got = np.array([rows[r["sample"]]["DJ.cn"] for r in res[1:]])
    assert np.all(np.abs(got / truth - 1) < 0.03), (got, truth)


def test_calibration_errors_name_the_class():
    with pytest.raises(ValueError, match="no usable windows for DJ"):
        cohort.calibrate_matrix(["a", "b"], np.full((2, 5), np.nan), np.arange(5) * 250, np.arange(1, 6) * 250, np.full(5, 0.5), cls="DJ")
    with pytest.raises(ValueError, match="DJ: the efficiency table has 4 windows"):
        cohort.calibrate_matrix(["a"], np.zeros((1, 5)), np.arange(5) * 250, np.arange(1, 6) * 250, np.full(5, 0.5), a_fixed=np.zeros(4), cls="DJ")


def test_the_same_sample_twice_is_refused(twelve):
    dup = copy.deepcopy(twelve[5])
    dup["mode"] = "fetch"
    with pytest.raises(ValueError, match=r"s05: given twice \(modes sim and fetch\)"):
        _table(twelve + [dup])


@pytest.fixture(scope="module")
def thirty():
    return _results(30, seed=7, classes=("DJ",))


def test_control_pcs_leave_out_only_the_samples_without_residuals(thirty):
    part = copy.deepcopy(thirty)
    part[4]["control_qc"] = None
    part[9]["control_qc"]["region_log_ratio"] = part[9]["control_qc"]["region_log_ratio"][:-1]
    rows, _, info, log = _table(part, n_control_pcs=5)
    alone, _, info_alone, _ = _table([r for i, r in enumerate(thirty) if i not in (4, 9)], n_control_pcs=5)
    assert "ctrlPC1" not in rows["s04"] and "ctrlPC_mp" not in rows["s09"]
    for s, r in alone.items():
        assert [rows[s].get(f"ctrlPC{k}") for k in range(1, 6)] == [r[f"ctrlPC{k}"] for k in range(1, 6)]
    assert info["excluded"] == ["s04", "s09"] and info["mp"] == info_alone["mp"]
    assert "1 sample has no control residuals" in log and "other than 800" in log and "s04" in log and "s09" in log


def test_control_pcs_asked_for_with_too_few_samples_stop_the_run(thirty):
    few = copy.deepcopy(thirty[:12])
    for r in few[:4]:
        r["control_qc"] = None
    with pytest.raises(ValueError, match="--control-pcs 5: control-region PCs need at least 10 samples"):
        _table(few, n_control_pcs=5)
    rows, _, info, log = _table(few)                                               # the default only says so
    assert info == {} and not any("ctrlPC1" in r for r in rows.values()) and "8 of 12 have them" in log
    rows, _, info, log = _table(few, n_control_pcs=0)
    assert info == {} and "control" not in log


def test_efficiencies_from_other_windows_are_refused(thirty):
    _, eff, _, _ = _table(thirty)
    rows, _, _, log = _table(thirty[:3], efficiencies=eff)                           # the intended use: a few new samples
    assert all(np.isfinite(r["DJ.cn"]) for r in rows.values()) and "WARNING" not in log
    rotated = copy.deepcopy(eff)
    for key in ("start", "gc", "anchor", "a", "window_sd"):
        rotated["DJ"][key] = rotated["DJ"][key][10:] + rotated["DJ"][key][:10]
    with pytest.raises(ValueError, match="--efficiencies: the DJ table .* window 1 starts at 2500 there and at 0 here"):
        _table(thirty[:3], efficiencies=rotated)
    halved = copy.deepcopy(eff)
    for key in ("start", "gc", "anchor", "a", "window_sd"):
        halved["DJ"][key] = halved["DJ"][key][::2]
    with pytest.raises(ValueError, match="learned on 24 windows, this cohort's estimates have 48"):
        _table(thirty[:3], efficiencies=halved)
    _, _, _, log = _table(thirty[:3], efficiencies={"rDNA45S": eff["DJ"]})
    assert "--efficiencies has no table for DJ" in log


def test_adjustment_refuses_a_fit_with_no_residual_degrees_of_freedom():
    rng = np.random.default_rng(2)
    X = rng.normal(size=(40, 3))
    y = np.exp(rng.normal(size=40))
    assert cohort.adjust_mask(y, X).sum() == 40
    y[4:] = np.nan
    assert cohort.adjust_mask(y, X).sum() == 4
    with pytest.raises(ValueError, match="4 usable values for 3 covariates"):
        cohort.adjust_for_covariates(y, X)
    y[4] = 1.0
    adj, _ = cohort.adjust_for_covariates(y, X)
    assert np.isfinite(adj).sum() == 5

