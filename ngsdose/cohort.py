"""Cohort-level calibration of positional classes.

Within one unit every window is present at the same copy number, so after the fragment-GC model

    log C_iw = c_i + a_w + e_iw

where c_i is the sample's log copy number, a_w a window efficiency shared by all samples
(sequence-specific dropout that no genome-wide GC curve captures, residual mappability, ...)
and e_iw the sample-specific distortion. The scale is pinned by the *anchor* windows - those
whose fragment GC lies where the control curve is best supported and the correction smallest -
by requiring median(a_w) = 0 over them. c_i then uses every window (precision) while the
absolute level is set by the anchors (accuracy).
"""
from __future__ import annotations

import warnings
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
    return calibrate_matrix([r["sample"] for r in results], Y, starts, ends, gc, **kw)


def calibrate_matrix(samples: list[str], Y: np.ndarray, starts: np.ndarray, ends: np.ndarray, gc: np.ndarray,
                     min_present: float = 0.9, a_fixed: np.ndarray | None = None, max_window_sd: float | None = None,
                     n_iter: int = 50, min_anchor: int = 5, anchors: list[tuple[int, int]] | None = None) -> Calibration:
    """Robust additive fit by median polish. With `a_fixed` (a shipped efficiency table) only the
    sample effects are estimated, which is what a single new sample needs."""
    Y = np.asarray(Y, float)
    n, m = Y.shape
    keep = np.mean(~np.isnan(Y), axis=0) >= min_present
    anchor = keep & (gc >= ANCHOR_GC[0]) & (gc <= ANCHOR_GC[1])
    if anchors:
        anchor &= np.array([any(s0 <= a and b <= e0 for s0, e0 in anchors) for a, b in zip(starts, ends)])
    cls = "class"
    if a_fixed is not None:
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
        raise ValueError(f"no usable windows for {cls}")
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


def adjust_for_covariates(y: np.ndarray, X: np.ndarray, log: bool = True) -> tuple[np.ndarray, float]:
    """Residualise a per-sample estimate on covariates (coverage PCs); returns adjusted values
    on the original scale (cohort mean restored) and the fraction of variance removed."""
    y = np.asarray(y, float)
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1) & ((y > 0) if log else True)
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


def cohort_table(results, anchors: dict, max_window_sd: float | None = None, n_profile_pcs: int = 3, n_control_pcs="mp",
                 mp_margin: float = 0.01, efficiencies: dict | None = None, log=None) -> tuple[list[dict], dict, dict]:
    """The cohort layer over per-sample estimates: window calibration of every positional class,
    profile PCs, and the control-region PCs with their Marchenko-Pastur count.

    `results` is an iterable of estimate results (dicts), read one at a time: only each sample's
    summary row, window vectors and control residuals are kept, so thousands of samples fit in
    a few hundred MB. Returns (rows in input order, efficiency tables per class, information
    about the control PCs). `n_control_pcs` is "mp" or a number of components to write."""
    from . import pcselect
    from .tables import summary_row
    say = log or (lambda *a, **k: None)
    rows, order, win, layout, ctrl = {}, [], {}, {}, []
    for r in results:
        rows[r["sample"]] = summary_row(r)
        order.append(r["sample"])
        for cls, v in r["classes"].items():
            if v["kind"] != "positional":
                continue
            y, gc, starts, ends = sample_windows(r, cls)
            if cls in layout and len(layout[cls][0]) != len(starts):
                raise ValueError(f"{r['sample']}: {cls} was estimated with a different window layout")
            layout.setdefault(cls, (starts, ends))
            win.setdefault(cls, []).append((y, gc))
        ctrl.append((r.get("control_qc") or {}).get("region_log_ratio"))
    eff, info = {}, {}
    for cls, per in win.items():
        if len(per) != len(order):
            continue
        a_fixed = None
        if efficiencies and cls in efficiencies:
            a_fixed = np.array([np.nan if v is None else v for v in efficiencies[cls]["a"]], float)
        Y = np.array([p[0] for p in per])
        gc = _nanmedian(np.array([p[1] for p in per]), axis=0)
        cal = calibrate_matrix(order, Y, layout[cls][0], layout[cls][1], gc, a_fixed=a_fixed, max_window_sd=max_window_sd,
                               anchors=(anchors or {}).get(cls))
        scores, var = profile_pcs(cal, n_profile_pcs) if len(order) > n_profile_pcs + 2 else (None, None)
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
    ok = len(order) >= 10 and all(x is not None for x in ctrl) and len({len(x) for x in ctrl}) == 1
    cp = control_pcs(np.array(ctrl, float), None) if ok and want != 0 else None
    if cp is not None:
        scores, var, sv, shape = cp
        sel = pcselect.mp_select(sv, *shape, margin=mp_margin)
        most = max(len(order) // 5, 1)                     # at least five samples per component
        n_write = min(scores.shape[1], most, want if want is not None else max(2 * sel.n_pc, 20))
        if want is not None and n_write < want:
            say(f"[cohort] WARNING: --control-pcs {want}, but {len(order)} samples support {n_write} (five samples per component): {n_write} written")
        for i, s in enumerate(order):
            rows[s]["ctrlPC_mp"] = sel.n_pc
            for k in range(n_write):
                rows[s][f"ctrlPC{k + 1}"] = round(float(scores[i, k]), 5)
        info = dict(mp=sel.n_pc, describe=sel.describe(), n_written=n_write, variance=[round(float(v), 5) for v in var[:n_write]],
                    singular_values=[round(float(v), 5) for v in sv], shape=list(shape))
        say(f"[cohort] control-region PCs: {sel.describe()}; {n_write} written (ctrlPC1..), variance explained: "
            + " ".join(f"{v:.3f}" for v in var[:min(n_write, 12)]))
    return [rows[s] for s in order], eff, info

