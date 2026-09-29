"""Cohort-level calibration of positional classes.

Within one unit every window is present at the same copy number, so after the fragment-GC model

    log C_iw = c_i + a_w + e_iw

where c_i is the sample's log copy number, a_w a window efficiency shared by all samples
(sequence-specific dropout that no genome-wide GC curve captures, residual mappability, ...)
and e_iw the sample-specific distortion. The scale is pinned by the *anchor* windows - those
present in >= 90% of samples, with fragment GC 40-60% (where the control curve is best supported
and the correction smallest) and, where the bundle ships them (anchors.json; ignored with
--gc-rule-anchors), lying wholly inside its empirically chosen anchor intervals; when fewer than
five windows qualify (e.g. the 68%-GC 5S unit), every retained window is used - by requiring
median(a_w) = 0 over them. c_i then uses every window (precision) while the absolute level is
set by the anchors (accuracy).

A class may carry rules of its own (the bundle's calibration.json), for what one model of every
window at one copy number does not hold:

- `level_exclude`: intervals of the unit where copies differ between people (the distal 22 kb of the
  distal junction, which one junction copy in sixteen lacks). c_i and the anchors use the rest, the
  *core*; the excluded windows keep their efficiencies and stay in the profile.
- `scale: mode`: a class with a known copy number in nearly everyone (the junction: ten) has its
  cohort's main mode pinned to it, c_i - delta and a_w + delta, when the cohort has enough genomes.
  The fragment-GC model in the anchor windows reads such a class a few percent off (its residual
  in repeat context); the pin is recorded, and a saved efficiency table carries it to new samples.
- `polymorphic`: intervals where the cohort's median genome is not at the level (a common deletion),
  so that the median gives the wrong efficiency. A cohort finds each interval's offset on its own
  comb of integers, with the reference state the highest that many genomes share (`comb_offset`);
  the bundle records what its reference cohort had, checked against assemblies, and a cohort too
  small for a comb takes that.
- `segments`: integer copy states along the unit are called for every genome (segments.py).
"""
from __future__ import annotations

import warnings
from collections import Counter
from dataclasses import dataclass

import numpy as np

from .estimate import ANCHOR_GC


@dataclass
class Calibration:
    samples: list[str]
    window_start: np.ndarray
    window_gc: np.ndarray
    anchor: np.ndarray            # bool per window
    a: np.ndarray                 # window efficiencies (log), NaN for windows not retained
    c: np.ndarray                 # per-sample log copy number
    c_se: np.ndarray              # robust SE of c_i from the window residuals
    resid: np.ndarray             # samples x windows
    resid_sd: np.ndarray          # per-sample residual MAD-SD ("profile roughness")
    window_sd: np.ndarray         # per-window residual MAD-SD across samples
    level: np.ndarray | None = None       # bool per window: the windows c_i is the median of (the core)
    scale: dict | None = None             # how the scale was set: rule, and for a pin its delta and the mode it moved
    offsets: list | None = None           # the polymorphic intervals' offsets as applied


def _nanmedian(x, axis=None, keepdims=False):
    with warnings.catch_warnings():                       # windows masked in every sample are all-NaN by design
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(x, axis=axis, keepdims=keepdims)


def _madsd(x, axis=None):
    med = _nanmedian(x, axis=axis, keepdims=True)
    return 1.4826 * _nanmedian(np.abs(x - med), axis=axis)


def sample_windows(result: dict, cls: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """One sample's windows of one class: log C_w (NaN where masked), GC, starts, ends."""
    wins = result["classes"][cls]["windows"]
    y = np.array([np.log(w["cn"]) if w["cn"] is not None and w["cn"] > 0 else np.nan for w in wins], float)
    return (y, np.array([w["gc"] for w in wins], float), np.array([w["start"] for w in wins]),
            np.array([w["end"] for w in wins]))


def window_matrix(results: list[dict], cls: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """log C_iw (NaN where masked), window starts, ends, and GC."""
    rows = [sample_windows(r, cls) for r in results]
    if len({len(r[0]) for r in rows}) != 1:
        raise ValueError("samples were estimated with different window layouts")
    return np.array([r[0] for r in rows]), rows[0][2], rows[0][3], _nanmedian(np.array([r[1] for r in rows]), axis=0)


def calibrate(results: list[dict], cls: str, **kw) -> Calibration:
    """Convenience wrapper over `calibrate_matrix` for results held in memory."""
    Y, starts, ends, gc = window_matrix(results, cls)
    return calibrate_matrix([r["sample"] for r in results], Y, starts, ends, gc, cls=cls, **kw)


def calibrate_matrix(samples: list[str], Y: np.ndarray, starts: np.ndarray, ends: np.ndarray, gc: np.ndarray,
                     min_present: float = 0.9, a_fixed: np.ndarray | None = None, max_window_sd: float | None = None,
                     n_iter: int = 50, min_anchor: int = 5, anchors: list[tuple[int, int]] | None = None,
                     cls: str = "class", rules: dict | None = None, log=None) -> Calibration:
    """Robust additive fit by median polish. With `a_fixed` (a shipped efficiency table) only the
    sample effects are estimated, which is what a single new sample needs. A sample with no usable
    window (no reads of the class) gets NaN and plays no part in which windows are retained."""
    Y = np.asarray(Y, float)
    n, m = Y.shape
    m_win = m
    present = ~np.all(np.isnan(Y), axis=1)
    if not present.any():
        raise ValueError(f"no usable windows for {cls}: none of the {n} samples has one")
    keep = np.mean(~np.isnan(Y[present]), axis=0) >= min_present
    say = log or (lambda *a, **k: None)
    rules = rules or {}
    level = keep.copy()
    for s0, e0 in rules.get("level_exclude", []):
        level &= ~((np.asarray(ends) > s0) & (np.asarray(starts) < e0))
    if level.sum() < min_anchor:
        raise ValueError(f"{cls}: the class's level_exclude leaves {int(level.sum())} windows for the level")
    anchor = level & (gc >= ANCHOR_GC[0]) & (gc <= ANCHOR_GC[1])
    if anchors:
        anchor &= np.array([any(s0 <= a and b <= e0 for s0, e0 in anchors) for a, b in zip(starts, ends)])
    if a_fixed is not None:
        if len(a_fixed) != m:
            raise ValueError(f"{cls}: the efficiency table has {len(a_fixed)} windows, the samples {m}")
        a = np.where(keep, a_fixed, np.nan)
        keep &= ~np.isnan(a)
        anchor &= keep
        level &= keep
    else:
        a = np.zeros(m)
    if anchor.sum() < min_anchor:
        # a unit with no moderate-GC sequence (the 5S unit is 68% GC throughout): the scale then
        # rests on the fragment-GC model alone, over every retained window
        anchor = level.copy()
    if anchor.sum() == 0:
        raise ValueError(f"no usable windows for {cls}: none is present in {100 * min_present:.0f}% of the {int(present.sum())} samples "
                         f"that have any" + (" and has an efficiency in the table" if a_fixed is not None else "")
                         + (f" ({n - int(present.sum())} samples have none)" if not present.all() else ""))
    Yk = np.where(keep[None, :], Y, np.nan)
    c = _nanmedian(Yk[:, anchor], axis=1)
    for _ in range(n_iter):
        if a_fixed is None:
            a_new = _nanmedian(Yk - c[:, None], axis=0)
            a_new = a_new - _nanmedian(a_new[anchor])
        else:
            a_new = a
        c_new = _nanmedian((Yk - a_new[None, :])[:, level], axis=1)
        done = np.nanmax(np.abs(c_new - c)) < 1e-9 and np.nanmax(np.abs(np.nan_to_num(a_new - a))) < 1e-9
        a, c = a_new, c_new
        if done:
            break
    resid = Yk - c[:, None] - a[None, :]
    wsd = _madsd(resid, axis=0) if n > 2 else np.full(m, np.nan)
    if max_window_sd is not None and a_fixed is None and n > 2:
        # unstable windows carry sample-specific structure (unit polymorphism); drop and refit once
        bad = keep & (wsd > max_window_sd) & ~anchor
        if bad.any():
            keep2 = keep & ~bad
            level = level & keep2
            Yk = np.where(keep2[None, :], Y, np.nan)
            for _ in range(n_iter):
                a = _nanmedian(Yk - c[:, None], axis=0)
                a = a - _nanmedian(a[anchor])
                c = _nanmedian((Yk - a[None, :])[:, level], axis=1)
            resid = Yk - c[:, None] - a[None, :]
            wsd = _madsd(resid, axis=0)
            keep = keep2
    a = np.where(keep, a, np.nan)
    n_fit = int(present.sum())
    scale, offsets = dict(rule="table" if a_fixed is not None else "anchors"), []
    if a_fixed is None and rules:
        sc = rules.get("scale") or {}
        expected = rules.get("expected_copies")
        if sc.get("rule") == "mode" and expected:
            if n_fit >= int(sc.get("min_samples", 50)):
                delta, mode, n_main = mode_pin(c, float(expected))
                c, a = c - delta, a + delta
                scale = dict(rule="mode", expected=float(expected), delta=round(delta, 5), factor=round(float(np.exp(-delta)), 5),
                             mode_on_anchors=round(mode, 4), n_main=n_main, n=n_fit)
                say(f"[cohort] {cls}: scale pinned to the cohort's mode: the core reads {mode:.3f} on the anchors' scale in the {n_main} genomes of "
                    f"the main mode, {expected:g} expected; every estimate x {np.exp(-delta):.4f}")
            else:
                say(f"[cohort] {cls}: {n_fit} genomes, fewer than the {int(sc.get('min_samples', 50))} a pin of the scale to the cohort's mode needs; "
                    "the scale is the anchors'")
        pinned = scale["rule"] == "mode"
        for pr in rules.get("polymorphic", []):
            s0, e0 = pr["interval"]
            m = keep & (np.asarray(starts) >= s0) & (np.asarray(starts) < e0)
            if not m.any():
                continue
            o, how, found = float(pr.get("offset", 0.0)), "the bundle's, measured on its reference cohort", {}
            if pinned:
                ratio = np.exp(_nanmedian((Yk - a[None, :])[:, m], axis=1) - c)               # the interval against the genome's level
                o, found = comb_offset(ratio, np.exp(c), max_gain=float(pr.get("max_gain", 0.15)))
                how = "this cohort's integer comb"
            a[m] = a[m] + o
            offsets.append(dict(name=pr.get("name", f"{s0}-{e0}"), interval=[int(s0), int(e0)], offset=round(o, 5), reference=pr.get("offset"), how=how,
                                windows=int(m.sum()), **found))
        if offsets:
            say(f"[cohort] {cls}: polymorphic intervals (log offset of their efficiencies): " + "; ".join(
                f"{x['name']} {x['offset']:+.4f} ({x['how']}" + (f"; the bundle's reference cohort had {x['reference']:+.3f}" if x.get("reference") is not None and x["how"].startswith("this") else "") + ")"
                for x in offsets))
    resid = Yk - c[:, None] - a[None, :]
    wsd = _madsd(resid, axis=0) if n > 2 else np.full(m_win, np.nan)
    rsd = _madsd(resid, axis=1)
    nw = np.sum(~np.isnan(resid[:, level]), axis=1)
    c_se = 1.2533 * _madsd(resid[:, level], axis=1) / np.sqrt(np.maximum(nw, 1))
    return Calibration(list(samples), starts, gc, anchor, a, c, c_se, resid, rsd, wsd, level=level, scale=scale, offsets=offsets)


def mode_pin(c: np.ndarray, expected: float) -> tuple[float, float, int]:
    """delta = log(mode / expected) for per-sample log levels `c`: the mode is the median of the genomes within half a
    copy of the cohort's median (the main mode of a class that most people carry at `expected` copies)."""
    lv = np.exp(c[np.isfinite(c)])
    m0 = float(np.median(lv))
    main = lv[np.abs(lv - m0) < 0.5 * m0 / expected]
    mode = float(np.median(main))
    return float(np.log(mode / expected)), mode, int(len(main))


def comb_offset(ratio: np.ndarray, level: np.ndarray, max_gain: float = 0.15, lo: float = -0.03, hi: float = 0.2, step: float = 0.0025,
                near: float = 0.25) -> tuple[float, dict]:
    """The log offset of a polymorphic interval's efficiencies, from the cohort itself.

    `ratio`: each genome's estimate of the interval over its level; `level`: its level in copies. A genome at N copies
    whose interval holds k of them has ratio k / N, so ratio * N is a whole number once the interval's efficiency is
    right. Where a deletion of the interval is common the cohort's median genome lacks part of a copy, the median gives
    the windows an efficiency too low, and every genome reads the interval too high by the same factor exp(offset).

    The offsets that put the cohort on whole numbers (the minima of its distance to them) lie one copy apart, and a
    comb cannot tell them from each other. What can: the interval's reference state, every copy holding it, is the
    highest state that many genomes share, because a copy can lack the interval but seldom holds it twice. The offset
    is the smallest that leaves at most `max_gain` of the genomes above their own level. Genomes within `near` of a
    whole number take part. Returns the offset and what was found (the share of genomes per state)."""
    n_int = np.round(level)
    ok = np.isfinite(ratio) & np.isfinite(level) & (np.abs(level - n_int) < near) & (n_int > 0)
    if ok.sum() < 20:
        return 0.0, dict(n=int(ok.sum()), note="too few genomes near a whole number")
    v, N = ratio[ok] * n_int[ok], n_int[ok]
    grid = np.arange(lo, hi + step / 2, step)
    cost = np.array([float(np.mean(np.minimum(np.abs(v * np.exp(-o) - np.round(v * np.exp(-o))), 0.5) ** 2)) for o in grid])
    minima = [k for k in range(len(grid)) if (k == 0 or cost[k] <= cost[k - 1]) and (k == len(grid) - 1 or cost[k] < cost[k + 1])]
    gain = lambda o: float(np.mean(np.round(v * np.exp(-o)) > N))
    chosen = next((k for k in minima if gain(grid[k]) <= max_gain), int(np.argmin(cost)))
    o = float(grid[chosen])
    d = np.round(v * np.exp(-o)) - N
    states = {int(k): round(float(np.mean(d == k)), 4) for k in np.unique(d)}
    return o, dict(n=int(ok.sum()), states=states, rms=round(float(np.sqrt(cost[chosen])), 4), rms_at_zero=round(float(np.sqrt(cost[np.argmin(np.abs(grid))])), 4))


def profile_pcs(cal: Calibration, n_pc: int = 5) -> tuple[np.ndarray, np.ndarray]:
    """Principal components of the window residuals: sample-specific profile distortions.

    Scores are candidate technical covariates; they are computed from e_iw, which is orthogonal
    to c_i by construction, so adjusting for them cannot absorb copy number itself."""
    E = cal.resid[:, ~np.all(np.isnan(cal.resid), axis=0)]
    E = np.where(np.isnan(E), 0.0, E)
    E = E - E.mean(0, keepdims=True)
    U, S, _ = np.linalg.svd(E, full_matrices=False)
    k = min(n_pc, len(S))
    return U[:, :k] * S[:k], (S ** 2 / max((S ** 2).sum(), 1e-30))[:k]


def adjust_mask(y: np.ndarray, X: np.ndarray, log: bool = True) -> np.ndarray:
    """The samples `adjust_for_covariates` fits on: the value finite (and positive on the log
    scale) and every covariate finite."""
    y, X = np.asarray(y, float), np.asarray(X, float)
    X = X[:, None] if X.ndim == 1 else X
    return np.isfinite(y) & np.all(np.isfinite(X), axis=1) & ((y > 0) if log else True)


def adjust_for_covariates(y: np.ndarray, X: np.ndarray, log: bool = True) -> tuple[np.ndarray, float]:
    """Residualise a per-sample estimate on covariates (coverage PCs); returns adjusted values
    on the original scale (cohort mean restored) and the fraction of variance removed. Raises
    ValueError when fewer than p + 2 samples are usable (`adjust_mask`) for p covariates."""
    y, X = np.asarray(y, float), np.asarray(X, float)
    X = X[:, None] if X.ndim == 1 else X
    ok = adjust_mask(y, X, log)
    if ok.sum() <= X.shape[1] + 1:
        raise ValueError(f"{int(ok.sum())} usable values for {X.shape[1]} covariates: the fit needs at least {X.shape[1] + 2}")
    z = np.log(y[ok]) if log else y[ok]
    A = np.column_stack([np.ones(ok.sum()), X[ok] - X[ok].mean(0)])
    beta, *_ = np.linalg.lstsq(A, z, rcond=None)
    res = z - A @ beta
    sst = float(((z - z.mean()) ** 2).sum())
    r2 = 1.0 - float((res ** 2).sum()) / sst if sst > 0 else 0.0
    out = np.full_like(y, np.nan)
    out[ok] = np.exp(res + z.mean()) if log else res + z.mean()
    return out, r2


def control_pcs(results: list[dict], n_pc: int | None = 10) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[int, int]] | None:
    """Principal components of the control regions' log(observed/expected) across the cohort.

    The matrix is samples x control regions, each row centred (a sample's overall level is its
    depth, already divided out) and each region centred across samples. The regions are
    single-copy sequence disjoint from every class, so - like NGS-PCA's coverage PCs, of which
    this is a small internal version - the scores can absorb library and sample structure
    (GC residue, replication timing in DNA from cycling cells, aneuploid chromosomes) but not
    the dosage of a class. Needs far more samples than components to be meaningful.

    Returns (scores, fraction of variance, all singular values, matrix shape); the last two are what
    `pcselect.mp_select` chooses the number of components from. `n_pc=None` returns every component."""
    rows = [(r.get("control_qc") or {}).get("region_log_ratio") for r in results] if isinstance(results, list) else results
    if isinstance(rows, list):
        if any(x is None for x in rows) or len({len(x) for x in rows}) != 1:
            return None
    X = np.array(rows, float)
    X = X - np.median(X, axis=1, keepdims=True)
    X = np.clip(X - np.median(X, axis=0, keepdims=True), -0.5, 0.5)       # a deleted region must not make a PC
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    k = len(S) if n_pc is None else min(n_pc, len(S))
    return U[:, :k] * S[:k], (S ** 2 / max((S ** 2).sum(), 1e-30))[:k], S, X.shape


def _few(names: list[str], k: int = 5) -> str:
    return ", ".join(names[:k]) + (f" and {len(names) - k} more" if len(names) > k else "")


def _n(names: list[str]) -> str:
    return "1 sample has" if len(names) == 1 else f"{len(names)} samples have"


def fixed_efficiencies(table: dict, cls: str, starts: np.ndarray) -> np.ndarray:
    """A saved efficiency table's log efficiencies for `cls`, after checking that it was learned on
    the same windows as this cohort's (a table from another window size or another unit origin
    would otherwise fail to broadcast, or be applied to the wrong windows without a word)."""
    have = [int(v) for v in starts]
    start = [int(v) for v in table.get("start", [])]
    if len(table.get("a", [])) != len(start):
        raise ValueError(f"--efficiencies: the {cls} table has {len(table.get('a', []))} efficiencies for {len(start)} window starts")
    if start != have:
        i = next((i for i, (p, q) in enumerate(zip(start, have)) if p != q), min(len(start), len(have)))
        at = (f"window {i + 1} starts at {start[i]} there and at {have[i]} here" if i < min(len(start), len(have))
              else "the shorter layout ends there")
        raise ValueError(f"--efficiencies: the {cls} table was learned on {len(start)} windows, this cohort's estimates have {len(have)}; "
                         f"{at} (another window size or another unit): fit the efficiencies on this cohort, or re-estimate with the bundle "
                         "the table came from")
    return np.array([np.nan if v is None else v for v in table["a"]], float)


def cohort_table(results, anchors: dict, max_window_sd: float | None = None, n_profile_pcs: int = 3, n_control_pcs="mp",
                 mp_margin: float = 0.01, efficiencies: dict | None = None, log=None, rules: dict | None = None,
                 profiles: dict | None = None) -> tuple[list[dict], dict, dict]:
    """The cohort layer over per-sample estimates: window calibration of every positional class,
    profile PCs, and the control-region PCs with their Marchenko-Pastur count.

    `results` is an iterable of estimate results (dicts), read one at a time: only each sample's
    summary row, window vectors and control residuals are kept (about 50 kB a sample), so thousands
    of samples fit in a few hundred MB. Returns (rows in input order, efficiency tables per class,
    information about the control PCs). `n_control_pcs` is "mp" or a number of components to write.

    Each sample id may appear once. A positional class is calibrated on the samples that have
    usable windows for it; the others (class absent from the estimate, a status other than "ok",
    no reads) get NA and are named in the log. The control PCs are computed on the samples whose
    control residuals have the most common length; the others get NA, and are named too.

    `rules`: per class, the bundle's calibration rules (`resources.Bundle.calibration()`): the core the level is set
    on, the pin of the scale, the polymorphic intervals, and whether integer copy states are called along the unit
    (columns `<class>.copies`, the copies the unit is described against; `.partial`, the copies that hold or lack an
    end of the unit; `.variants`, every event; `.scale_f`; `.call`, settled or uncertain, with `.call_gap`, how far behind
    the best other reading is; `<class>.cn_unit` is the level over the whole unit, what `<class>.cn` was before the class
    had a core). `profiles`, a dict, is filled per class with every genome's
    calibrated window profile and its segments, for a caller that draws or tabulates them."""
    from . import pcselect
    from . import segments as seg
    from .tables import summary_row
    say = log or (lambda *a, **k: None)
    rows, order, classes, win, layout, lacking, ctrl = {}, [], [], {}, {}, {}, []
    ctrl_key = []                                          # which control regions, in which order, per sample
    for r in results:
        s = r["sample"]
        if s in rows:
            raise ValueError(f"{s}: given twice (modes {rows[s].get('mode')} and {r.get('mode')}); a cohort takes one estimate per sample - "
                             "put scan and fetch estimates of the same genomes in separate cohorts, and give every sample its own name")
        rows[s] = summary_row(r)
        order.append(s)
        for sk in r.get("skipped_classes") or []:
            if sk.get("kind") == "positional":
                lacking.setdefault(sk["name"], {})[s] = f"skipped: {sk.get('reason')}"
        for cls, v in r["classes"].items():
            if v["kind"] != "positional":
                continue
            if cls not in classes:
                classes.append(cls)
            status = v.get("status", "ok")
            if status != "ok":
                lacking.setdefault(cls, {})[s] = f"status {status}"
                continue
            y, gc, starts, ends = sample_windows(r, cls)
            if cls in layout and len(layout[cls][0]) != len(starts):
                raise ValueError(f"{s}: {cls} was estimated with a different window layout ({len(starts)} windows, "
                                 f"{len(layout[cls][0])} in the samples before it)")
            layout.setdefault(cls, (starts, ends))
            if np.isnan(y).all():
                lacking.setdefault(cls, {})[s] = "no usable window"
                continue
            win.setdefault(cls, []).append((s, y, gc))
        qc = r.get("control_qc") or {}
        x = qc.get("region_log_ratio")
        ctrl.append(None if x is None else np.asarray(x, float))
        # the hash `estimate` writes of the regions' names in order (estimates written before it: their count)
        ctrl_key.append(None if x is None else (qc.get("region_order_sha256") or f"{len(x)} regions"))
    eff, info = {}, {}
    for cls in classes:
        per = win.pop(cls, [])
        if len(per) < len(order):
            have = {p[0] for p in per}
            why = lacking.get(cls, {})
            miss = [f"{s} ({why.get(s, 'class absent from the estimate')})" for s in order if s not in have]
            if not per:
                say(f"[cohort] WARNING: {cls}: no sample has usable windows ({_few(miss)}); the class is not calibrated")
                continue
            absent = any(s not in have and s not in why for s in order)
            say(f"[cohort] WARNING: {cls}: calibrated on the {len(per)} of {len(order)} samples that have it; NA for {_few(miss)}"
                + (" (a class absent from an estimate usually means estimates made with bundles whose panel or units differ)" if absent else ""))
        names = [p[0] for p in per]
        a_fixed = None
        if efficiencies:
            if cls in efficiencies:
                a_fixed = fixed_efficiencies(efficiencies[cls], cls, layout[cls][0])
            else:
                say(f"[cohort] WARNING: --efficiencies has no table for {cls}: its efficiencies are fitted on this cohort")
        Y = np.array([p[1] for p in per])
        gc = _nanmedian(np.array([p[2] for p in per]), axis=0)
        del per
        rule = (rules or {}).get(cls) or {}
        if rule and rule.get("unit_length") and int(rule["unit_length"]) != int(layout[cls][1][-1]):
            say(f"[cohort] WARNING: {cls}: the bundle's calibration rules are for a unit of {int(rule['unit_length']):,} bp, these estimates' windows end at "
                f"{int(layout[cls][1][-1]):,}: the rules are not applied (estimates from another bundle's unit)")
            rule = {}
        cal = calibrate_matrix(names, Y, layout[cls][0], layout[cls][1], gc, a_fixed=a_fixed, max_window_sd=max_window_sd,
                               anchors=(anchors or {}).get(cls), cls=cls, rules=rule, log=say)
        del Y
        if a_fixed is not None and (efficiencies[cls].get("scale") or {}).get("rule") not in (None, "anchors", "table"):
            cal.scale = dict(efficiencies[cls]["scale"], rule="table", table_rule=efficiencies[cls]["scale"]["rule"])
            cal.offsets = efficiencies[cls].get("polymorphic") or []
        kept = ~np.isnan(cal.a)
        calibrated = np.exp(cal.c[:, None] + cal.resid)                        # samples x windows, copies; NaN where not retained
        unit_level = np.exp(cal.c + _nanmedian(np.where(kept[None, :], cal.resid, np.nan), axis=1))
        calls, fractions = {}, None
        if rule.get("segments") is not None and rule.get("segments") is not False:
            sp = rule["segments"] if isinstance(rule["segments"], dict) else {}
            ref = _nanmedian(cal.window_sd[cal.level]) if cal.level is not None and np.isfinite(cal.window_sd[cal.level]).any() else np.nan
            rel = np.clip(cal.window_sd / ref, 0.5, 1.5) if np.isfinite(ref) and ref > 0 else np.ones(len(kept))
            rel = np.where(cal.level, rel, 1.0)                                # a polymorphic window's cohort SD holds the polymorphism, not its noise
            poly = [tuple(iv) for iv in rule.get("level_exclude", [])]
            unit_len = int(layout[cls][1][-1])
            for i, smp in enumerate(cal.samples):
                call = seg.segment(layout[cls][0], np.where(kept, calibrated[i], np.nan), rel, window=int(layout[cls][1][0] - layout[cls][0][0]),
                                   unit_length=unit_len, tau=float(sp.get("tau", seg.TAU)), min_windows=int(sp.get("min_windows", seg.MIN_WINDOWS)),
                                   scale_sd=float(sp.get("scale_sd", seg.SCALE_SD)), tilt=bool(sp.get("tilt", True)), tilt_sd=float(sp.get("tilt_sd", seg.TILT_SD)))
                if call is not None:
                    calls[smp] = seg.describe(call, poly, int(sp.get("min_core", seg.MIN_CORE)), expected=(int(rule["expected_copies"]) if rule.get("expected_copies") else None))
                    calls[smp].readings = []                                   # the other scales' readings have served
            fr = rule.get("fractions")
            if calls and fr is not None and fr is not False:
                fp = fr if isinstance(fr, dict) else {}
                table = ((efficiencies or {}).get(cls) or {}).get("fractions") or {}
                fractions = seg.find_fractions([calls.get(smp) for smp in cal.samples], layout[cls][0], np.where(kept[None, :], calibrated, np.nan), rel,
                                               window=int(layout[cls][1][0] - layout[cls][0][0]), unit_length=unit_len, leave_out=poly, gc=cal.window_gc,
                                               scale_sd=float(sp.get("scale_sd", seg.SCALE_SD)), tilt_sd=float(sp.get("tilt_sd", seg.TILT_SD)), tau=float(sp.get("tau", seg.TAU)),
                                               level_z=float(fp.get("level_z", seg.LEVEL_Z)), event_z=float(fp.get("event_z", seg.EVENT_Z)),
                                               min_height=float(fp.get("min_height", seg.FRACTION_MIN_HEIGHT)), min_windows=int(fp.get("min_windows", seg.FRACTION_MIN_WINDOWS)),
                                               min_samples=int(fp.get("min_samples", seg.FRACTION_MIN_SAMPLES)), spread=table.get("spread") if a_fixed is not None else None)
                n_level = sum(1 for c in calls.values() if not c.uncertain and c.off_z is not None and abs(c.off_z) >= c.level_z)
                n_step = sum(1 for c in calls.values() if not c.uncertain and c.fractions)
                say(f"[cohort] {cls}: off the whole numbers: the scales' spread is {100 * fractions['spread']:.2f}% ({fractions['spread_from']}); "
                    f"{n_level} genomes' levels lie {fractions['level_z']:g} SDs or more from their whole number"
                    + (f", {n_step} carry a step of fractional height" if fractions["steps"] else
                       f"; steps of fractional height are not looked for in fewer than {int(fp.get('min_samples', seg.FRACTION_MIN_SAMPLES))} genomes")
                    + f"; {sum(1 for c in calls.values() if c.uncertain)} calls are uncertain")
        scores, var = profile_pcs(cal, n_profile_pcs) if len(names) > n_profile_pcs + 2 else (None, None)
        for i, s in enumerate(cal.samples):
            rows[s][f"{cls}.cn"] = round(float(np.exp(cal.c[i])), 2)
            rows[s][f"{cls}.cn_se_rel"] = round(float(cal.c_se[i]), 5)
            rows[s][f"{cls}.profile_sd"] = round(float(cal.resid_sd[i]), 4)
            if rule:
                rows[s][f"{cls}.cn_unit"] = round(float(unit_level[i]), 2)
            if s in calls:
                cl = calls[s]
                poly = [tuple(iv) for iv in rule.get("level_exclude", [])]
                rows[s][f"{cls}.copies"] = cl.copies
                rows[s][f"{cls}.partial"] = ";".join(f"{d:+d}:{a // 1000}-{b // 1000}kb" for d, a, b in
                                                     seg.partial_copies(cl, poly) + seg.partial_copies(cl, poly, kind="partial loss")) or "none"
                rows[s][f"{cls}.variants"] = seg.events_string(cl) or "none"
                rows[s][f"{cls}.scale_f"] = round(cl.scale, 3)
                rows[s][f"{cls}.tilt"] = round(cl.tilt, 3)
                rows[s][f"{cls}.call_gap"] = None if not np.isfinite(cl.gap) else round(cl.gap, 1)
                rows[s][f"{cls}.call"] = cl.status
                if fractions is not None:
                    rows[s][f"{cls}.off"] = None if cl.off is None else round(cl.off, 2)
                    rows[s][f"{cls}.off_z"] = None if cl.off_z is None else round(cl.off_z, 2)
                    rows[s][f"{cls}.fractional"] = seg.fractions_string(cl) or "none"
                    rows[s][f"{cls}.fractional_z"] = ";".join(f"{f.z:+.1f}" for f in cl.fractions) or "none"
            if scores is not None:
                for k in range(scores.shape[1]):
                    rows[s][f"{cls}.profilePC{k + 1}"] = round(float(scores[i, k]), 5)
        eff[cls] = dict(start=cal.window_start.tolist(), gc=[round(float(g), 4) for g in cal.window_gc],
                        anchor=cal.anchor.tolist(), a=[None if np.isnan(v) else round(float(v), 5) for v in cal.a],
                        window_sd=[None if np.isnan(v) else round(float(v), 5) for v in cal.window_sd],
                        n_samples=len(cal.samples))
        if rule or (cal.scale or {}).get("rule") == "table":
            eff[cls].update(level=cal.level.tolist(), scale=cal.scale, polymorphic=cal.offsets)
        if fractions is not None:
            eff[cls]["fractions"] = {k: (round(v, 5) if isinstance(v, float) else v) for k, v in fractions.items()}
        if profiles is not None:
            profiles[cls] = dict(samples=list(cal.samples), start=cal.window_start.tolist(), end=[int(e) for e in layout[cls][1]],
                                 cn=np.where(kept[None, :], calibrated, np.nan).astype(np.float32), level=cal.level.copy(), calls=calls,
                                 scale=cal.scale, offsets=cal.offsets)
    # control-region PCs: how many are structure is decided at the Marchenko-Pastur edge of the noise
    # bulk ("mp"); the table carries more than that, so that `pcsweep` can look beyond the choice
    want = None if str(n_control_pcs).lower() == "mp" else int(n_control_pcs)
    if want == 0:
        return [rows[s] for s in order], eff, info
    # one set of control regions in one order for all: the residual vectors are matched by position, so a sample
    # whose regions differ in set or order (another controls file, another bundle revision) must stay out
    keys = Counter(k for k in ctrl_key if k is not None)
    m = keys.most_common(1)[0][0] if keys else None
    use = [i for i, k in enumerate(ctrl_key) if k is not None and k == m]
    m_len = len(ctrl[use[0]]) if use else None
    no_qc = [order[i] for i, k in enumerate(ctrl_key) if k is None]
    other_n = [order[i] for i, k in enumerate(ctrl_key) if k is not None and k != m and len(ctrl[i]) != m_len]
    other_set = [order[i] for i, k in enumerate(ctrl_key) if k is not None and k != m and len(ctrl[i]) == m_len]
    if no_qc or other_n or other_set:
        say("[cohort] WARNING: control-region PCs: " + "; ".join(
            ([f"{_n(no_qc)} no control residuals (estimated with --no-control-qc, or counts without regions): {_few(no_qc)}"] if no_qc else [])
            + ([f"{_n(other_n)} a number of control regions other than {m_len} (another controls file): {_few(other_n)}"] if other_n else [])
            + ([f"{_n(other_set)} {m_len} control regions that are not the same regions in the same order as the other "
                f"{len(use)} samples' (another controls file or bundle revision): {_few(other_set)}"] if other_set else []))
            + f"; their ctrlPC columns are NA, the PCs are computed on the other {len(use)}")
    if len(use) < 10:
        msg = f"control-region PCs need at least 10 samples with control residuals of one length; {len(use)} of {len(order)} have them"
        if want is not None:
            raise ValueError(f"--control-pcs {want}: {msg}")
        say(f"[cohort] {msg}: no ctrlPC columns written")
        return [rows[s] for s in order], eff, info
    X = np.array([ctrl[i] for i in use], float)
    del ctrl
    scores, var, sv, shape = control_pcs(X, None)
    del X
    sel = pcselect.mp_select(sv, *shape, margin=mp_margin)
    most = max(len(use) // 5, 1)                           # at least five samples per component
    n_write = min(scores.shape[1], most, want if want is not None else max(2 * sel.n_pc, 20))
    if want is not None and n_write < want:
        say(f"[cohort] WARNING: --control-pcs {want}, but {len(use)} samples support {n_write} (five samples per component): {n_write} written")
    for j, i in enumerate(use):
        s = order[i]
        rows[s]["ctrlPC_mp"] = sel.n_pc
        for k in range(n_write):
            rows[s][f"ctrlPC{k + 1}"] = round(float(scores[j, k]), 5)
    info = dict(mp=sel.n_pc, describe=sel.describe(), n_written=n_write, variance=[round(float(v), 5) for v in var[:n_write]],
                singular_values=[round(float(v), 5) for v in sv], shape=list(shape))
    other = other_n + other_set
    if no_qc or other:
        info["excluded"] = no_qc + other
    say(f"[cohort] control-region PCs: {sel.describe()}; {n_write} written (ctrlPC1..), variance explained: "
        + " ".join(f"{v:.3f}" for v in var[:min(n_write, 12)]))
    return [rows[s] for s in order], eff, info
