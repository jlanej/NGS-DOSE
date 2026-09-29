"""Class rules of the cohort layer (the core a level is set on, the pin of the scale, the polymorphic
intervals) and integer copy states along a unit (segments.py), on profiles whose truth is known."""
import json
import subprocess
import sys

import numpy as np
import pytest

from ngsdose import cohort, segments

W, U = 250, 400000
STARTS = np.arange(0, U, W)


def profile(rng, level=10.0, pieces=(), noise=0.7, scale=1.0, usable=0.75):
    """One genome's windows: `level` copies, `pieces` = (start, end, delta) on top, a scale and window noise."""
    x = np.full(len(STARTS), float(level))
    for a, b, d in pieces:
        x[(STARTS >= a) & (STARTS < b)] += d
    x = x * scale + rng.normal(0, noise, len(x))
    x[rng.random(len(x)) > usable] = np.nan
    return x


def call(x, poly=((0, 30000), (190000, 232000)), **kw):
    c = segments.segment(STARTS, x, unit_length=U, **kw)
    return segments.describe(c, list(poly), expected=10)


def test_a_flat_genome_is_one_segment_and_its_scale_is_found():
    rng = np.random.default_rng(1)
    for scale in (1.0, 0.985, 1.02):
        c = call(profile(rng, scale=scale))
        assert [s.state for s in c.segments] == [10] and c.copies == 10 and c.events == []
        assert abs(c.scale - scale) <= 0.006 and not c.uncertain and c.gap > 5 and not c.segments[0].off_integer
    c = call(profile(rng, level=9, scale=1.0))
    assert c.copies == 9 and c.events == [] and not c.uncertain
    between = call(profile(rng, scale=0.95))                       # a genome at 9.5 throughout: a whole number is called, and the other reading is close behind
    assert between.copies in (9, 10) and between.uncertain and between.alternative["bulk"] in (9, 10) and between.alternative["bulk"] != between.copies
    assert between.gap < segments.MIN_GAP


def test_a_partial_copy_is_found_with_its_breakpoint():
    rng = np.random.default_rng(2)
    c = call(profile(rng, pieces=[(0, 316000, 1)]))
    assert [s.state for s in c.segments] == [11, 10] and abs(c.segments[0].end - 316000) <= 1500
    assert c.copies == 10 and [(e.delta, e.kind) for e in c.events] == [(1, "partial copy")]
    assert segments.partial_copies(c) == [(1, 0, c.segments[0].end)]
    # the same genome at a level of 10.5 on a median (three quarters of it at 11): its integers are 10 and 11 all the same
    assert abs(c.scale - 1) <= 0.01
    # a copy that begins at 122 kb: ten copies, one of which lacks the distal 122 kb
    d = call(profile(rng, level=9, pieces=[(122000, 400000, 1)]))
    lost = segments.partial_copies(d, kind="partial loss")
    assert d.copies == 10 and segments.partial_copies(d) == [] and len(lost) == 1 and lost[0][:2] == (-1, 0) and abs(lost[0][2] - 122000) <= 1500
    # a genome that holds ten nowhere is described against what it holds most: eight and a copy of 122-400 kb
    e = call(profile(rng, level=8, pieces=[(122000, 400000, 1)]))
    assert e.copies == 9 and [(x.delta, x.kind) for x in e.events] == [(-1, "partial loss")]


def test_an_internal_deletion_is_a_loss_not_two_partial_copies():
    rng = np.random.default_rng(3)
    c = call(profile(rng, level=9, pieces=[(246000, 300000, -1)]))
    assert c.copies == 9 and [(e.delta, e.kind) for e in c.events] == [(-1, "loss")]
    assert abs(c.events[0].start - 246000) <= 1500 and abs(c.events[0].end - 300000) <= 1500
    assert segments.partial_copies(c) == []


def test_a_gain_interrupted_by_a_polymorphic_deletion_is_one_partial_copy():
    rng = np.random.default_rng(4)
    c = call(profile(rng, pieces=[(0, 316000, 1), (197000, 217000, -2)]))
    assert [s.state for s in c.segments] == [11, 9, 11, 10]
    assert c.copies == 10 and segments.events_string(c).count(":") == 3
    pc = segments.partial_copies(c, [(0, 30000), (190000, 232000)])
    assert len(pc) == 1 and pc[0][0] == 1 and pc[0][1] == 0 and abs(pc[0][2] - 316000) <= 1500


def test_wild_windows_and_short_blips_make_no_segment():
    rng = np.random.default_rng(5)
    x = profile(rng)
    ok = np.flatnonzero(np.isfinite(x))
    x[ok[rng.choice(len(ok), 25, replace=False)]] += rng.choice([-6, 6, 9], 25)          # wild single windows
    x[(STARTS >= 100000) & (STARTS < 101500)] += 1.5                                      # six windows a copy and a half up
    c = call(x)
    assert [s.state for s in c.segments] == [10]
    # what can be found: a change of state costs tau, so one copy over about 10 kb (thirty usable windows) or two over 4 kb
    assert [s.state for s in call(profile(rng, pieces=[(150000, 165000, -1)])).segments] == [10, 9, 10]
    assert [s.state for s in call(profile(rng, pieces=[(150000, 158000, -2)])).segments] == [10, 8, 10]
    assert [s.state for s in call(profile(rng, pieces=[(150000, 155000, -1)])).segments] == [10]


def test_a_lean_is_a_lean_and_a_step_a_step():
    """A profile that rises smoothly by a copy across the unit is one state with a lean; a step stays a step, with or
    without a lean under it; without the lean in the model the rise is broken into a step."""
    rng = np.random.default_rng(7)
    t = (STARTS + W / 2) / U - 0.5
    for g in (0.09, -0.06, 0.03):
        c = call(profile(rng) * np.exp(g * t))
        assert [s.state for s in c.segments] == [10] and abs(c.tilt - g) <= 0.015 and not c.uncertain and abs(c.segments[0].mean - 10) < 0.1
    c = call(profile(rng, pieces=[(0, 316000, 1)]))
    assert [s.state for s in c.segments] == [11, 10] and abs(c.tilt) <= 0.02 and abs(c.segments[0].end - 316000) <= 1500
    c = call(profile(rng, pieces=[(0, 316000, 1)]) * np.exp(0.06 * t))
    assert [s.state for s in c.segments] == [11, 10] and abs(c.tilt - 0.06) <= 0.02 and abs(c.segments[0].end - 316000) <= 1500
    c = call(profile(rng, pieces=[(0, 316000, 1), (197000, 217000, -2)]) * np.exp(-0.05 * t))
    assert [s.state for s in c.segments] == [11, 9, 11, 10] and abs(c.tilt + 0.05) <= 0.02
    steep = profile(rng) * np.exp(0.14 * t)                              # a copy and a half across the unit
    broken = call(steep, tilt=False)
    assert len(broken.segments) > 1 and broken.tilt == 0
    c = call(steep)
    assert [s.state for s in c.segments] == [10] and abs(c.tilt - 0.14) <= 0.02


def test_too_few_windows_give_no_call():
    x = np.full(len(STARTS), np.nan)
    x[:20] = 10.0
    assert segments.segment(STARTS, x) is None


def test_a_profile_without_noise_is_still_called():
    """A synthetic profile (every window the same, a step at 200 kb) has no window noise to measure; the noise has a floor."""
    x = np.full(len(STARTS), 10.0)
    c = segments.describe(segments.segment(STARTS, x, unit_length=U), [], expected=10)
    assert c.copies == 10 and not c.events and len(c.segments) == 1 and c.sigma > 0
    x[STARTS >= 200000] = 11.0
    c = segments.describe(segments.segment(STARTS, x, unit_length=U), [], expected=10)
    assert [s.state for s in c.segments] == [10, 11] and c.segments[0].end == 200000 and [str(e) for e in c.events] == ["+1:200-400kb"]


# ------------------------------------------------------------------ the rules in the median polish
def cohort_matrix(rng, n=300, m=400, anchors_read=0.975, del_frac=0.4, noise=0.05):
    """log C_iw of a cohort: copies 10 in most genomes, 9 and 11 in some; window efficiencies; the anchors read
    `anchors_read` of the truth; windows 40-60 lack one copy in `del_frac` of the genomes."""
    cn = rng.choice([9, 10, 11], n, p=[0.06, 0.88, 0.06]).astype(float)
    a = rng.normal(0, 0.08, m)
    state = np.tile(cn[:, None], (1, m))
    carrier = rng.random(n) < del_frac
    state[np.ix_(carrier, np.arange(40, 60))] -= 1
    Y = np.log(state) + a[None, :] + np.log(anchors_read) + rng.normal(0, noise, (n, m))
    starts = np.arange(m) * 250
    return Y, starts, starts + 250, np.full(m, 0.5), cn, carrier, a


RULES = dict(expected_copies=10, scale=dict(rule="mode", min_samples=50), level_exclude=[(40 * 250, 60 * 250)],
             polymorphic=[dict(name="the deletion", interval=(40 * 250, 60 * 250), offset=0.04)])


def test_the_level_is_set_on_the_core_and_pinned_to_the_mode():
    rng = np.random.default_rng(11)
    Y, starts, ends, gc, cn, carrier, a = cohort_matrix(rng)
    names = [f"s{i}" for i in range(len(cn))]
    log = []
    cal = cohort.calibrate_matrix(names, Y, starts, ends, gc, rules=RULES, log=log.append)
    lv = np.exp(cal.c)
    assert np.allclose(lv, cn, atol=0.12)                                           # pinned: ten copies read ten
    assert cal.scale["rule"] == "mode" and abs(cal.scale["factor"] - 1 / 0.975) < 0.004 and "pinned" in "\n".join(log)
    assert not cal.level[40:60].any() and cal.level[:40].all() and cal.level[60:].all()
    # a genome with the deletion has the level of one without: the interval is not in the level
    ten = cn == 10
    assert abs(np.median(lv[ten & carrier]) - np.median(lv[ten & ~carrier])) < 0.02
    # without rules the scale is the anchors' and the deletion pulls the carriers' level down
    plain = cohort.calibrate_matrix(names, Y, starts, ends, gc)
    assert abs(np.median(np.exp(plain.c)) - 9.75) < 0.03 and plain.scale["rule"] == "anchors" and plain.level.all()
    # a cohort too small for a pin keeps the anchors' scale, and says so
    log = []
    small = cohort.calibrate_matrix(names[:30], Y[:30], starts, ends, gc, rules=RULES, log=log.append)
    assert small.scale["rule"] == "anchors" and abs(np.median(np.exp(small.c)) - 9.75) < 0.08 and "fewer than the 50" in "\n".join(log)


def test_a_polymorphic_interval_is_put_on_whole_numbers():
    rng = np.random.default_rng(12)
    Y, starts, ends, gc, cn, carrier, a = cohort_matrix(rng)
    names = [f"s{i}" for i in range(len(cn))]
    cal = cohort.calibrate_matrix(names, Y, starts, ends, gc, rules=RULES)
    prof = np.exp(cal.c[:, None] + cal.resid)
    inside = np.nanmedian(prof[:, 40:60], axis=1)
    truth = np.where(carrier, cn - 1, cn)
    assert abs(np.median(inside - truth)) < 0.06 and np.mean(np.round(inside) == truth) > 0.97
    assert len(cal.offsets) == 1 and cal.offsets[0]["how"].startswith("this cohort") and 0.01 <= cal.offsets[0]["offset"] <= 0.07
    assert abs(cal.offsets[0]["states"][0] - 0.6) < 0.08 and abs(cal.offsets[0]["states"][-1] - 0.4) < 0.08
    # a cohort in which most genomes lack a copy of the interval: the median is a whole copy off, and the rule still finds the reference
    Y2, _, _, _, cn2, carrier2, _ = cohort_matrix(np.random.default_rng(15), del_frac=0.7)
    cal2 = cohort.calibrate_matrix(names, Y2, starts, ends, gc, rules=RULES)
    inside2 = np.nanmedian(np.exp(cal2.c[:, None] + cal2.resid)[:, 40:60], axis=1)
    assert np.mean(np.round(inside2) == np.where(carrier2, cn2 - 1, cn2)) > 0.97 and 0.08 <= cal2.offsets[0]["offset"] <= 0.12
    # and one in which nobody does: no offset
    Y3, _, _, _, cn3, _, _ = cohort_matrix(np.random.default_rng(16), del_frac=0.0)
    cal3 = cohort.calibrate_matrix(names, Y3, starts, ends, gc, rules=RULES)
    assert abs(cal3.offsets[0]["offset"]) <= 0.01
    # without the rule the cohort's median sits between the two groups and both read off a whole number
    no_poly = dict(RULES, polymorphic=[])
    plain = cohort.calibrate_matrix(names, Y, starts, ends, gc, rules=no_poly)
    off = np.nanmedian(np.exp(plain.c[:, None] + plain.resid)[:, 40:60], axis=1)
    assert np.median(off[~carrier] - cn[~carrier]) > 0.2 and np.mean(np.round(off) == truth) < 0.9


def test_a_saved_table_carries_the_pin_to_a_new_genome():
    rng = np.random.default_rng(13)
    Y, starts, ends, gc, cn, carrier, a = cohort_matrix(rng)
    cal = cohort.calibrate_matrix([f"s{i}" for i in range(len(cn))], Y, starts, ends, gc, rules=RULES)
    new = cohort.calibrate_matrix(["one"], Y[:1], starts, ends, gc, a_fixed=cal.a, rules=RULES)
    assert abs(np.exp(new.c[0]) - cn[0]) < 0.12 and new.scale["rule"] == "table" and new.offsets == []


def test_rules_for_another_unit_are_not_applied():
    from test_cohort_robustness import _results
    res = _results(12, classes=("DJ",))
    log = []
    rows, eff, _ = cohort.cohort_table(res, {}, log=log.append, n_control_pcs=0, rules={"DJ": dict(RULES, unit_length=400000, segments={})})
    assert "DJ.copies" not in rows[0] and "level" not in eff["DJ"] and "the rules are not applied" in "\n".join(log)


def test_cohort_table_calls_segments_and_fills_profiles():
    """Estimates of a 40-kb unit whose windows are made by hand: a cohort at ten copies, one genome with a partial copy."""
    rng = np.random.default_rng(14)
    m, n = 160, 80
    eff_true = rng.normal(0, 0.08, m)
    import copy
    from test_cohort_robustness import _results
    template = _results(1, classes=("U",))[0]                                          # a whole estimate; its class's windows are replaced
    res = []
    for i in range(n):
        cn = np.full(m, 10.0)
        if i == 0:
            cn[:120] += 1                                                               # a partial copy of the first 30 kb
        vals = cn * np.exp(eff_true) * 0.97 * np.exp(rng.normal(0, 0.05, m))
        r = copy.deepcopy(template)
        r["sample"] = f"g{i:02d}"
        r["classes"]["U"]["windows"] = [dict(start=j * 250, end=(j + 1) * 250, gc=0.5, obs=300, exp=30.0, usable=1.0, cn=round(float(v), 3)) for j, v in enumerate(vals)]
        res.append(r)
    rules = dict(U=dict(unit_length=40000, expected_copies=10, scale=dict(rule="mode", min_samples=50), level_exclude=[],
                        segments=dict(tau=16, min_windows=16, min_core=5000)))
    prof = {}
    rows, eff, _ = cohort.cohort_table(res, {}, n_control_pcs=0, rules=rules, profiles=prof)
    by = {r["sample"]: r for r in rows}
    assert by["g00"]["U.copies"] == 10 and by["g00"]["U.partial"].startswith("+1:0-3") and by["g05"]["U.partial"] == "none" and by["g05"]["U.copies"] == 10
    assert abs(by["g05"]["U.cn"] - 10) < 0.15 and eff["U"]["scale"]["rule"] == "mode" and abs(eff["U"]["scale"]["factor"] - 1 / 0.97) < 0.01
    assert prof["U"]["cn"].shape == (n, m) and set(prof["U"]["calls"]) == set(by) and prof["U"]["calls"]["g00"].copies == 10
    # the saved table gives one new genome the same call
    rows1, eff1, _ = cohort.cohort_table(res[:1], {}, n_control_pcs=0, rules=rules, efficiencies=eff)
    assert rows1[0]["U.partial"] == by["g00"]["U.partial"] and eff1["U"]["scale"]["rule"] == "table" and eff1["U"]["scale"]["table_rule"] == "mode"


def test_the_bundle_rules_load_and_lie_inside_the_unit():
    from ngsdose import resources
    B = resources.Bundle()
    rules = B.calibration()
    assert set(rules) == {"DJ"} and rules["DJ"]["expected_copies"] == B.meta["expected_copies"]["DJ"]
    U_dj = len(B.units()["DJ"])
    assert rules["DJ"]["unit_length"] == U_dj
    for a, b in rules["DJ"]["level_exclude"] + [p["interval"] for p in rules["DJ"]["polymorphic"]]:
        assert 0 <= a < b <= U_dj
    for p in rules["DJ"]["polymorphic"]:                                               # a polymorphic interval is left out of the level
        assert any(a <= p["interval"][0] and p["interval"][1] <= b for a, b in rules["DJ"]["level_exclude"])
        assert 0 <= p["offset"] < 0.12 and p["kind"] == "deletion" and 0 < p["max_gain"] < 0.5


def test_the_cohort_command_writes_the_segments(tmp_path):
    """Sixty estimates of a 400-kb unit named DJ, one genome with the partial copy of the first 316 kb: the bundle's
    rules apply (the unit is theirs), the table has the calls and --segments the states."""
    import copy
    from ngsdose.tables import dump
    from test_cohort_robustness import _results
    rng = np.random.default_rng(21)
    template = _results(1, classes=("DJ",))[0]
    eff_true = rng.normal(0, 0.08, len(STARTS))
    usable = rng.random(len(STARTS)) < 0.75
    files = []
    for i in range(60):
        cn = np.full(len(STARTS), 10.0)
        if i == 0:
            cn[STARTS < 316000] += 1
        if i % 3 == 0:
            cn[(STARTS >= 197000) & (STARTS < 217000)] -= 1                              # the common deletion, in a third of the genomes
        vals = cn * np.exp(eff_true) * 0.975 * np.exp(rng.normal(0, 0.07, len(STARTS)))
        r = copy.deepcopy(template)
        r["sample"] = f"g{i:02d}"
        r["classes"]["DJ"]["windows"] = [dict(start=int(a), end=int(a) + W, gc=0.5, obs=300, exp=30.0, usable=1.0, cn=(round(float(v), 3) if u else None))
                                         for a, v, u in zip(STARTS, vals, usable)]
        f = tmp_path / f"g{i:02d}.estimate.json.gz"
        dump(r, f)
        files.append(str(f))
    table, segs, eff = tmp_path / "cohort.tsv", tmp_path / "segments.tsv", tmp_path / "eff.json"
    run = subprocess.run([sys.executable, "-m", "ngsdose", "cohort", *files, "-t", str(table), "--segments", str(segs), "--save-efficiencies", str(eff), "--control-pcs", "0"],
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert "scale pinned to the cohort's mode" in run.stderr and "polymorphic intervals" in run.stderr
    import csv
    T = {r["sample"]: r for r in csv.DictReader(open(table), delimiter="\t")}
    assert T["g00"]["DJ.copies"] == "10" and T["g00"]["DJ.partial"].startswith("+1:0-31") and T["g01"]["DJ.partial"] == "none"
    assert abs(float(T["g01"]["DJ.cn"]) - 10) < 0.15 and "-1:19" in T["g03"]["DJ.variants"] and T["g01"]["DJ.variants"] == "none"
    S = [r for r in csv.DictReader(open(segs), delimiter="\t")]
    assert {r["class"] for r in S} == {"DJ"} and len({r["sample"] for r in S}) == 60
    g0 = [r for r in S if r["sample"] == "g00"]
    assert [r["state"] for r in g0] == ["11", "10", "11", "10"] and g0[0]["complete_copies"] == "10"
    E = json.loads(eff.read_text())["DJ"]
    assert E["scale"]["rule"] == "mode" and abs(E["scale"]["factor"] - 1 / 0.975) < 0.01 and len(E["polymorphic"]) == 3 and not all(E["level"])
    by_name = {x["name"]: x for x in E["polymorphic"]}
    assert abs(by_name["distal 5-15 kb"]["offset"]) <= 0.01 and 0.02 <= by_name["197-217 kb"]["offset"] <= 0.06      # found on this cohort, not taken from the bundle
    # without the rules: the old estimate, no calls
    run = subprocess.run([sys.executable, "-m", "ngsdose", "cohort", *files, "-t", str(table), "--no-class-rules", "--control-pcs", "0"], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    T = {r["sample"]: r for r in csv.DictReader(open(table), delimiter="\t")}
    assert "DJ.copies" not in T["g00"] and abs(float(T["g01"]["DJ.cn"]) - 9.75) < 0.15
