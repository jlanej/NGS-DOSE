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
                     cls: str = "class") -> Calibration:
    """Robust additive fit by median polish. With `a_fixed` (a shipped efficiency table) only the
    sample effects are estimated, which is what a single new sample needs. A sample with no usable
    window (no reads of the class) gets NaN and plays no part in which windows are retained."""
    Y = np.asarray(Y, float)
    n, m = Y.shape
    present = ~np.all(np.isnan(Y), axis=1)
    if not present.any():
        raise ValueError(f"no usable windows for {cls}: none of the {n} samples has one")
    keep = np.mean(~np.isnan(Y[present]), axis=0) >= min_present
    anchor = keep & (gc >= ANCHOR_GC[0]) & (gc <= ANCHOR_GC[1])
    if anchors:
        anchor &= np.array([any(s0 <= a and b <= e0 for s0, e0 in anchors) for a, b in zip(starts, ends)])
    if a_fixed is not None:
        if len(a_fixed) != m:
            raise ValueError(f"{cls}: the efficiency table has {len(a_fixed)} windows, the samples {m}")
        a = np.where(keep, a_fixed, np.nan)
        keep &= ~np.isnan(a)
        anchor &= keep
    else:
        a = np.zeros(m)
    if anchor.sum() < min_anchor:
        # a unit with no moderate-GC sequence (the 5S unit is 68% GC throughout): the scale then
        # rests on the fragment-GC model alone, over every retained window
        anchor = keep.copy()
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
        c_new = _nanmedian(Yk - a_new[None, :], axis=1)
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
            Yk = np.where(keep2[None, :], Y, np.nan)
            for _ in range(n_iter):
                a = _nanmedian(Yk - c[:, None], axis=0)
                a = a - _nanmedian(a[anchor])
                c = _nanmedian(Yk - a[None, :], axis=1)
            resid = Yk - c[:, None] - a[None, :]
            wsd = _madsd(resid, axis=0)
            keep = keep2
    a = np.where(keep, a, np.nan)
    rsd = _madsd(resid, axis=1)
    nw = np.sum(~np.isnan(resid), axis=1)
    c_se = 1.2533 * rsd / np.sqrt(np.maximum(nw, 1))
    return Calibration(list(samples), starts, gc, anchor, a, c, c_se, resid, rsd, wsd)


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
                 mp_margin: float = 0.01, efficiencies: dict | None = None, log=None) -> tuple[list[dict], dict, dict]:
    """The cohort layer over per-sample estimates: window calibration of every positional class,
    profile PCs, and the control-region PCs with their Marchenko-Pastur count.

    `results` is an iterable of estimate results (dicts), read one at a time: only each sample's
    summary row, window vectors and control residuals are kept (about 50 kB a sample), so thousands
    of samples fit in a few hundred MB. Returns (rows in input order, efficiency tables per class,
    information about the control PCs). `n_control_pcs` is "mp" or a number of components to write.

    Each sample id may appear once. A positional class is calibrated on the samples that have
    usable windows for it; the others (class absent from the estimate, a status other than "ok",
    no reads) get NA and are named in the log. The control PCs are computed on the samples whose
    control residuals have the most common length; the others get NA, and are named too."""
    from . import pcselect
    from .tables import summary_row
    say = log or (lambda *a, **k: None)
    rows, order, classes, win, layout, lacking, ctrl = {}, [], [], {}, {}, {}, []
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
        x = (r.get("control_qc") or {}).get("region_log_ratio")
        ctrl.append(None if x is None else np.asarray(x, float))
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
        cal = calibrate_matrix(names, Y, layout[cls][0], layout[cls][1], gc, a_fixed=a_fixed, max_window_sd=max_window_sd,
                               anchors=(anchors or {}).get(cls), cls=cls)
        del Y
        scores, var = profile_pcs(cal, n_profile_pcs) if len(names) > n_profile_pcs + 2 else (None, None)
        for i, s in enumerate(cal.samples):
            rows[s][f"{cls}.cn"] = round(float(np.exp(cal.c[i])), 2)
            rows[s][f"{cls}.cn_se_rel"] = round(float(cal.c_se[i]), 5)
            rows[s][f"{cls}.profile_sd"] = round(float(cal.resid_sd[i]), 4)
            if scores is not None:
                for k in range(scores.shape[1]):
                    rows[s][f"{cls}.profilePC{k + 1}"] = round(float(scores[i, k]), 5)
        eff[cls] = dict(start=cal.window_start.tolist(), gc=[round(float(g), 4) for g in cal.window_gc],
                        anchor=cal.anchor.tolist(), a=[None if np.isnan(v) else round(float(v), 5) for v in cal.a],
                        window_sd=[None if np.isnan(v) else round(float(v), 5) for v in cal.window_sd],
                        n_samples=len(cal.samples))
    # control-region PCs: how many are structure is decided at the Marchenko-Pastur edge of the noise
    # bulk ("mp"); the table carries more than that, so that `pcsweep` can look beyond the choice
    want = None if str(n_control_pcs).lower() == "mp" else int(n_control_pcs)
    if want == 0:
        return [rows[s] for s in order], eff, info
    lengths = Counter(len(x) for x in ctrl if x is not None)
    m = lengths.most_common(1)[0][0] if lengths else None
    use = [i for i, x in enumerate(ctrl) if x is not None and len(x) == m]
    no_qc = [order[i] for i, x in enumerate(ctrl) if x is None]
    other = [order[i] for i, x in enumerate(ctrl) if x is not None and len(x) != m]
    if no_qc or other:
        say("[cohort] WARNING: control-region PCs: " + "; ".join(
            ([f"{_n(no_qc)} no control residuals (estimated with --no-control-qc, or counts without regions): {_few(no_qc)}"] if no_qc else [])
            + ([f"{_n(other)} a number of control regions other than {m} (another controls file): {_few(other)}"] if other else []))
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
    if no_qc or other:
        info["excluded"] = no_qc + other
    say(f"[cohort] control-region PCs: {sel.describe()}; {n_write} written (ctrlPC1..), variance explained: "
        + " ".join(f"{v:.3f}" for v in var[:min(n_write, 12)]))
    return [rows[s] for s in order], eff, info
