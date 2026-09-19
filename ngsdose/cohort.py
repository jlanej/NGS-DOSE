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
