"""Class rules of the cohort layer (the core a level is set on, the pin of the scale, the polymorphic
intervals) and integer copy states along a unit (segments.py), on profiles whose truth is known."""
import json
import subprocess
import sys

import numpy as np

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


# ------------------------------------------------------------------ off the whole numbers
POLY = [(0, 30000), (190000, 232000)]
USABLE = np.random.default_rng(99).random(len(STARTS)) < 0.75            # the unit's windows that hold panel k-mers: the same in every genome
GC = np.where(STARTS < 105000, 0.43, 0.39) + np.random.default_rng(98).normal(0, 0.03, len(STARTS))    # richer below 105 kb, as the junction is


def cohort_of(rng, n, special=None, noise=0.65, scale_sd=0.013, lean_sd=0.015, gc_sd=0.0):
    """Calibrated profiles of `n` genomes at ten copies, each with its own scale, lean and (with `gc_sd`) slope on the
    windows' GC; `special`: genome -> (level, pieces) in its place, at a scale of one and without a lean. Returns
    (values, calls)."""
    t = (STARTS + W / 2) / U - 0.5
    V, C = [], []
    for i in range(n):
        level, pieces = (special or {}).get(i, (10.0, ()))
        x = np.full(len(STARTS), float(level))
        for a, b, d in pieces:
            x[(STARTS >= a) & (STARTS < b)] += d
        own = rng.normal(0, scale_sd) + rng.normal(0, lean_sd) * t + rng.normal(0, gc_sd) * (GC - GC.mean())
        x = x * np.exp(0.0 if i in (special or {}) else own) + rng.normal(0, noise, len(x))
        x = np.where(USABLE, x, np.nan)
        V.append(x)
        C.append(segments.describe(segments.segment(STARTS, x, unit_length=U, scale_sd=0.015), POLY, expected=10))
    return np.array(V), C


def test_a_level_off_its_whole_number_is_measured_and_judged_against_the_cohort():
    """The nearest whole number is called; how far the level lies from it is kept, and judged against the cohort's scales."""
    rng = np.random.default_rng(5)
    V, C = cohort_of(rng, 120, special={0: (9.6, ()), 1: (10.36, ()), 2: (9.9, ()), 3: (9.0, ())}, scale_sd=0.012)
    assert C[0].copies == 10 and not C[0].uncertain and C[0].status == "settled" and C[0].off is None       # before the fractions are looked for
    used = segments.find_fractions(C, STARTS, V, unit_length=U, leave_out=POLY, scale_sd=0.015, level_z=2.5)
    assert used["steps"] and used["spread_from"] == "this cohort" and 0.009 < used["spread"] < 0.017 and used["n"] == 120 and used["level_z"] == 2.5
    low, high, near, nine = C[:4]
    assert low.copies == 10 and abs(low.off + 0.4) < 0.08 and low.off_z < -2.5 and low.status == "fractional" and low.fractional and not low.fractions
    assert high.copies == 10 and abs(high.off - 0.36) < 0.08 and high.off_z > 2.5 and high.status == "fractional"
    assert near.copies == 10 and abs(near.off + 0.1) < 0.08 and abs(near.off_z) < 2.5 and near.status == "settled"   # 9.9: a scale like many another
    assert nine.copies == 9 and abs(nine.off) < 0.3 and nine.status == "settled"                             # nine copies, on its whole number
    rest = C[4:]
    assert sum(c.status == "fractional" for c in rest) <= 5 and all(c.ratio is not None and abs(c.ratio - 1) < 0.06 for c in rest)
    assert abs(np.median([c.off for c in rest])) < 0.05
    # a genome half way between two whole numbers is uncertain, whatever its offset
    V2, C2 = cohort_of(rng, 60, special={0: (9.5, ())})
    segments.find_fractions(C2, STARTS, V2, unit_length=U, leave_out=POLY)
    assert C2[0].uncertain and C2[0].status == "uncertain" and not C2[0].fractional and abs(abs(C2[0].off) - 0.5) < 0.15


def test_a_step_of_fractional_height_is_found_where_the_cohort_has_none():
    """A partial copy of the first 316 kb in half of the cells: the chain of whole numbers calls nothing or a whole copy;
    the step is found, placed and sized, and nothing is found in the genomes that have none."""
    rng = np.random.default_rng(8)
    V, C = cohort_of(rng, 150, special={0: (10.0, ((0, 316000, 0.5),)), 1: (10.0, ((100000, 180000, -0.45),))})
    used = segments.find_fractions(C, STARTS, V, unit_length=U, leave_out=POLY, gc=GC)
    assert used["steps"] and "gc_slope_sd" in used
    f = C[0].fractions
    assert len(f) == 1 and C[0].status == "fractional" and abs(f[0].z) >= 4
    end = f[0].end if f[0].start == 0 else f[0].start                                                          # the step, whichever side is described
    assert abs(end - 316000) <= 12000 and abs(abs(f[0].offset) - 0.5) < 0.15
    assert {f[0].start, f[0].end} & {0, U} and segments.fractions_string(C[0]).endswith("kb") and segments.fractions_string(C[0])[0] in "+-"
    g = C[1].fractions                                                                                         # a loss inside the unit, in 45% of the cells
    assert len(g) == 1 and abs(g[0].start - 100000) <= 12000 and abs(g[0].end - 180000) <= 12000 and abs(g[0].offset + 0.45) < 0.15 and g[0].z <= -4
    assert g[0].label(10).startswith("-0.4") and g[0].height(10) == g[0].offset                                # the whole number there is the ten it is described against
    assert sum(bool(c.fractions) for c in C[2:]) <= 1                                                          # nothing there: nothing found
    # a partial copy in 0.58 of the cells (in quieter libraries): whether the chain calls a whole copy there or none, the fraction gives its height
    V3, C3 = cohort_of(rng, 150, special={0: (10.0, ((0, 316000, 0.58),))}, noise=0.45)
    segments.find_fractions(C3, STARTS, V3, unit_length=U, leave_out=POLY, gc=GC)
    (h,) = C3[0].fractions
    assert C3[0].status in ("fractional", "uncertain") and h.start == 0 and abs(h.end - 316000) <= 12000
    assert h.state in (10, 11) and abs(h.height(C3[0].copies) - 0.58) < 0.12 and (h.offset < 0) == (h.state == 11)
    # in two thirds of the cells it is a third of a copy from a whole number, which is within what the cohort's profiles wander by: nothing is said
    V4, C4 = cohort_of(rng, 150, special={0: (10.0, ((0, 316000, 0.68),))})
    segments.find_fractions(C4, STARTS, V4, unit_length=U, leave_out=POLY, gc=GC)
    assert any(e.delta == 1 and e.start == 0 and abs(e.end - 316000) <= 8000 and e.kind == "partial copy" for e in C4[0].events)
    assert all(abs(f.offset) < 0.45 for f in C4[0].fractions)


def test_a_profile_that_follows_gc_is_no_step():
    """A library whose GC response the model left in: its profile follows the windows' GC, which here is richer below
    105 kb. With the GC beside the lean in the fit, the library has no step; without it, it has."""
    rng = np.random.default_rng(12)
    V, C = cohort_of(rng, 150, gc_sd=0.05)
    x = 10.0 * np.exp(1.2 * (GC - GC.mean())) + rng.normal(0, 0.65, len(STARTS))                            # a slope of 1.2, 24 of the cohort's SDs: half a copy across 105 kb
    V[0] = np.where(USABLE, x, np.nan)
    C[0] = segments.describe(segments.segment(STARTS, V[0], unit_length=U, scale_sd=0.015), POLY, expected=10)
    segments.find_fractions(C, STARTS, V, unit_length=U, leave_out=POLY)
    without = list(C[0].fractions)
    segments.find_fractions(C, STARTS, V, unit_length=U, leave_out=POLY, gc=GC)
    assert without and not C[0].fractions and C[0].status in ("settled", "fractional")
    assert without[0].start == 0 or without[0].end == U


def test_few_genomes_are_judged_by_the_rules_and_no_step_is_looked_for():
    rng = np.random.default_rng(3)
    V, C = cohort_of(rng, 20, special={0: (9.6, ()), 1: (10.0, ((0, 316000, 0.5),))})
    used = segments.find_fractions(C, STARTS, V, unit_length=U, leave_out=POLY, gc=GC, scale_sd=0.015, level_z=2.5)
    assert not used["steps"] and used["spread_from"] == "the class's rules" and used["spread"] == 0.015 and "gc_slope_sd" not in used
    assert C[0].status == "fractional" and abs(C[0].off_z - np.log(C[0].ratio) / 0.015) < 1e-6 and not any(c.fractions for c in C)
    used = segments.find_fractions(C, STARTS, V, unit_length=U, leave_out=POLY, spread=0.03, level_z=2.5)
    assert used["spread_from"] == "the efficiency table's cohort" and used["spread"] == 0.03 and C[0].status == "settled" and abs(C[0].off_z) < 2.5
    none = segments.find_fractions([None, None], STARTS, V[:2], unit_length=U)
    assert none["n"] == 0 and not none["steps"]


def test_the_levels_of_a_fit_are_the_stretches_means():
    """fit_levels and propose_steps on a profile whose truth is known: two levels and a lean."""
    rng = np.random.default_rng(4)
    n = 800
    t = np.linspace(-0.5, 0.5, n)
    y = np.where(np.arange(n) < 500, 0.0, 0.05) + 0.03 * t + rng.normal(0, 0.06, n)
    w = np.full(n, 1 / 0.06 ** 2)
    S = segments._Sums(y, w, t[:, None])
    lv, co, ll = segments.fit_levels(S, [500], [0.03])
    assert abs(lv[1] - lv[0] - 0.05) < 0.02 and abs(co[0] - 0.03) < 0.03
    _, _, l0 = segments.fit_levels(S, [], [0.03])
    assert ll - l0 > 16
    cuts = segments.propose_steps(y, w, t[:, None], 10.0, [0.03])
    assert len(cuts) == 1 and abs(cuts[0] - 500) <= 16
    assert segments.propose_steps(y, w, t[:, None], 10.0, [0.03], also=[500]) == [500]                          # where the whole numbers change, a step is tried
    flat = 0.03 * t + rng.normal(0, 0.06, n)
    assert segments.propose_steps(flat, w, t[:, None], 10.0, [0.03]) == []
    assert segments.propose_steps(y[:150], w[:150], t[:150, None], 10.0, [0.03]) == []                          # too few windows for a stretch on either side


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
    fr = rules["DJ"]["fractions"]
    assert fr["level_z"] >= 3 and fr["event_z"] >= fr["level_z"] and 0.2 <= fr["min_height"] <= 0.5 and fr["min_samples"] >= 50 and fr["min_windows"] >= 50
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
        if i == 7:
            cn *= 0.958                                                                  # a junction lost in four cells of ten
        if i == 8:
            cn[STARTS < 110000] -= 0.5                                                   # a copy that lacks the first 110 kb, in half of the cells
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
    assert [r["state"] for r in g0] == ["11", "10", "11", "10"] and g0[0]["complete_copies"] == "10" and {r["kind"] for r in g0} == {"segment"}
    assert abs(float(g0[0]["raw"]) / float(g0[0]["mean"]) - float(g0[0]["scale_f"])) < 0.03              # what the reads give, before the genome's scale
    # off the whole numbers: a level, and a stretch
    assert "off the whole numbers" in run.stderr and "carry a step of fractional height" in run.stderr
    assert T["g07"]["DJ.copies"] == "10" and T["g07"]["DJ.call"] == "fractional" and abs(float(T["g07"]["DJ.off"]) + 0.42) < 0.1 and float(T["g07"]["DJ.off_z"]) < -3
    assert T["g07"]["DJ.fractional"] == "none" and T["g07"]["DJ.fractional_z"] == "none"
    assert T["g01"]["DJ.call"] == "settled" and abs(float(T["g01"]["DJ.off"])) < 0.3 and T["g00"]["DJ.call"] == "settled" and T["g00"]["DJ.fractional"] == "none"
    assert T["g08"]["DJ.call"] == "fractional" and T["g08"]["DJ.partial"] == "none"
    h, (a, b) = T["g08"]["DJ.fractional"].split(":")[0], T["g08"]["DJ.fractional"].split(":")[1].removesuffix("kb").split("-")
    assert abs(abs(float(h)) - 0.5) < 0.15 and abs(float(T["g08"]["DJ.fractional_z"])) >= 4
    assert (h[0] == "-" and a == "0" and abs(int(b) - 110) <= 12) or (h[0] == "+" and abs(int(a) - 110) <= 12 and b == "400")   # the stretch that lacks it, or the rest against it
    f8 = [r for r in S if r["sample"] == "g08" and r["kind"] == "fraction"]
    assert len(f8) == 1 and f8[0]["state"] == "10" and f8[0]["off_integer"] == "True" and abs(float(f8[0]["z"])) >= 4 and f8[0]["call"] == "fractional"
    assert abs(abs(float(f8[0]["mean"]) - 10) - 0.5) < 0.15 and sum(r["kind"] == "fraction" for r in S) <= 3
    E = json.loads(eff.read_text())["DJ"]
    assert E["fractions"]["spread_from"] == "this cohort" and E["fractions"]["steps"] and 0.002 < E["fractions"]["spread"] < 0.02
    # one genome counted later, with the saved table: its level is judged against the spread of the table's cohort
    one = tmp_path / "one.tsv"
    again = subprocess.run([sys.executable, "-m", "ngsdose", "cohort", files[7], "-t", str(one), "--efficiencies", str(eff), "--control-pcs", "0"], capture_output=True, text=True)
    assert again.returncode == 0, again.stderr
    (r7,) = list(csv.DictReader(open(one), delimiter="\t"))
    assert "the efficiency table's cohort" in again.stderr and "steps of fractional height are not looked for" in again.stderr
    assert r7["DJ.call"] == "fractional" and abs(float(r7["DJ.off"]) - float(T["g07"]["DJ.off"])) < 0.05 and r7["DJ.fractional"] == "none"
    assert E["scale"]["rule"] == "mode" and abs(E["scale"]["factor"] - 1 / 0.975) < 0.01 and len(E["polymorphic"]) == 3 and not all(E["level"])
    by_name = {x["name"]: x for x in E["polymorphic"]}
    assert abs(by_name["distal 5-15 kb"]["offset"]) <= 0.01 and 0.02 <= by_name["197-217 kb"]["offset"] <= 0.06      # found on this cohort, not taken from the bundle
    # without the rules: the old estimate, no calls
    run = subprocess.run([sys.executable, "-m", "ngsdose", "cohort", *files, "-t", str(table), "--no-class-rules", "--control-pcs", "0"], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    T = {r["sample"]: r for r in csv.DictReader(open(table), delimiter="\t")}
    assert "DJ.copies" not in T["g00"] and abs(float(T["g01"]["DJ.cn"]) - 9.75) < 0.15
