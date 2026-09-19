"""How many coverage PCs to regress out: a default from random-matrix theory, and the evidence to overrule it.

**The default: the Marchenko-Pastur edge.** The singular values of an n x p matrix of pure noise
fill a bulk with a sharp upper edge (sigma * (sqrt(n) + sqrt(p)) for identically distributed
noise); a component standing above the edge is structure, one inside the bulk cannot be told
from noise. Counting the components above the edge needs the edge, and the textbook way to get
it - fit the law to the spectrum, e.g. match its median (Gavish & Donoho 2014) - assumes every
entry has the same noise variance. Coverage does not oblige: noise falls with a sample's depth
and varies with a bin's mappability. On simulated noise whose rows differ in SD by +-40%, the
textbook estimate called 34 components where 7 had been planted, and 148 for 5.

What survives unequal variances is the *shape* of the edge: the density of any such noise bulk
vanishes like a square root at its top, so the j-th largest noise value sits at E - a * j^(2/3)
(the universality that also gives the Tracy-Widom law). So the edge E is fitted, with a, to the
lower half of the leading singular values - the only part of the spectrum NGS-PCA's randomized
SVD keeps in any case - and the fit is iterated, because j counts from the first noise value
and that depends on how many components are signal. A component is selected if it clears E by
four residual SDs of the fit (plus the Tracy-Widom scale, which matters for dozens of samples
and not for thousands). On the simulations above this recovers 7 of 7 and 5 of 5, selects
nothing in pure noise in 95% of runs, and on the 1000 Genomes spectrum gives 46-51 components
whether 100, 150 or all 200 kept values are used (the textbook fit: 59, 66, 80).

**The evidence: a sweep against known truth** (`sweep`). Every sample carries sequence of known
copy number (held-out autosomal: 2; chrX and chrY by sex; the distal junction: 10), and trios say
how much of an estimate is transmitted. Adding PCs one at a time, the error of the known truths
- cross-validated, because residual variance falls with every regressor whether it means
anything or not - and the transmission reliability of the classes show directly where
adjustment stops removing noise and starts removing nothing, or signal.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def mp_singular_value_quantiles(n: int, p: int, probs, grid: int = 200_001) -> np.ndarray:
    """Quantiles of the singular values of an n x p matrix of iid unit-variance noise
    (Marchenko-Pastur law for the eigenvalues of X X^T / max(n, p), mapped to singular values)."""
    n, p = (n, p) if n <= p else (p, n)
    g = n / p
    lo, hi = (1 - np.sqrt(g)) ** 2, (1 + np.sqrt(g)) ** 2
    x = np.linspace(lo, hi, grid)
    dens = np.sqrt(np.maximum((hi - x) * (x - lo), 0.0)) / (2 * np.pi * g * x)
    cdf = np.cumsum(dens)
    cdf /= cdf[-1]
    return np.sqrt(np.interp(np.asarray(probs, float), cdf, x) * p)


@dataclass
class MPSelection:
    n_pc: int                      # components standing above the noise edge
    edge: float                    # fitted upper edge of the noise bulk, in singular-value units
    slope: float                   # a in E - a * j^(2/3)
    resid_sd: float                # residual SD of that fit: how closely the tail follows the edge law
    sigma: float                   # noise SD per entry that an identically distributed bulk with this edge would have
    slope_ratio: float             # fitted slope / the slope of that identically distributed bulk (1 = equal variances)
    n: int
    p: int
    n_values: int                  # singular values supplied
    n_fitted: int                  # leading values the fit looked at
    notes: list[str] = field(default_factory=list)

    def describe(self) -> str:
        if not np.isfinite(self.edge):
            return "Marchenko-Pastur edge: " + "; ".join(self.notes)
        return (f"Marchenko-Pastur edge: {self.n_pc} components above the noise bulk (n = {self.n:,} x p = {self.p:,}; edge {self.edge:.4g} "
                f"fitted on the leading {self.n_fitted} of {self.n_values} singular values; noise SD ~{self.sigma:.3g}; edge {self.slope_ratio:.2f}x as "
                f"broad as equal-variance noise would make it)" + "".join(f"; {x}" for x in self.notes))


def mp_select(singular_values, n_rows: int, n_cols: int, min_samples: int = 30, top_frac: float = 0.15, min_top: int = 40,
              z: float = 4.0) -> MPSelection:
    """Number of components above the edge of the noise bulk. `singular_values` are those of the
    (centred) n_rows x n_cols matrix, in any order: all of them, or only the largest."""
    s = np.sort(np.asarray(singular_values, float))[::-1]
    s = s[np.isfinite(s) & (s > 0)]
    n, p = min(n_rows, n_cols), max(n_rows, n_cols)
    nan = float("nan")
    if n < min_samples or len(s) < 12:
        return MPSelection(0, nan, nan, nan, nan, nan, n, p, len(s), 0,
                           [f"too few samples or singular values for a noise bulk (n = {n}, {len(s)} values): no components selected"])
    # the edge law holds near the top of the bulk: of a whole spectrum, look at the leading part only
    k = len(s) if len(s) < 0.9 * (n - 1) else min(len(s), max(min_top, int(round(top_frac * len(s)))))
    tw = 2.0 * (1 / np.sqrt(n) + 1 / np.sqrt(p)) ** (1 / 3) / (2 * (np.sqrt(n) + np.sqrt(p)))      # Tracy-Widom scale / edge
    m, E, a, sd, notes = 0, nan, nan, nan, []
    for _ in range(500):
        r = np.arange(max(m + 1, k // 2 + 1), k + 1)                   # ranks taken to be noise
        if len(r) < 6:
            notes.append("nearly every kept component is above the edge: the count is a lower bound - keep more components")
            break
        x = (r - m - 0.5) ** (2 / 3)
        A = np.column_stack([np.ones_like(x), -x])
        (E, a), *_ = np.linalg.lstsq(A, s[r - 1], rcond=None)
        sd = float((s[r - 1] - A @ [E, a]).std(ddof=2))
        m_new = int((s[:k] > E * (1 + tw) + z * sd).sum())
        if m_new <= m:                                                 # the count only grows as signal is set aside
            break
        m = m_new
    sigma = float(E / (np.sqrt(n) + np.sqrt(p)))
    # the slope an identically distributed bulk with this edge would have, fitted the same way
    j = np.arange(1, max(k - m, 8) + 1)
    q = sigma * mp_singular_value_quantiles(max(n - m, 2), p, 1 - (j - 0.5) / max(n - m, 2))
    use = j > len(j) // 2
    iid = np.linalg.lstsq(np.column_stack([np.ones(use.sum()), -(j[use] - 0.5) ** (2 / 3)]), q[use], rcond=None)[0][1]
    return MPSelection(m, float(E), float(a), sd, sigma, float(a / iid) if iid > 0 else nan, n, p, len(s), k, notes)


# ------------------------------------------------------------------------------------------------
# sweep: what each further PC does to known truths and to transmission
# ------------------------------------------------------------------------------------------------

def cv_adjusted(y: np.ndarray, pcs: np.ndarray, max_pc: int, folds: int = 10, seed: int = 1) -> np.ndarray:
    """Out-of-fold residuals of `y` on the first k PCs, for every k = 0..max_pc: array (max_pc + 1, n),
    each row with the mean of `y` restored. Rows of `y` that are not finite stay NaN.

    The models are nested, so one QR factorisation per fold serves every k."""
    y = np.asarray(y, float)
    ok = np.isfinite(y) & np.all(np.isfinite(pcs[:, :max_pc]), axis=1)
    idx = np.where(ok)[0]
    out = np.full((max_pc + 1, len(y)), np.nan)
    if len(idx) < max(3 * folds, max_pc + folds + 2):
        return out
    rng = np.random.default_rng(seed)
    fold = rng.permutation(len(idx)) % folds
    mean = y[idx].mean()
    for f in range(folds):
        tr, te = idx[fold != f], idx[fold == f]
        centre = pcs[tr, :max_pc].mean(0)
        A = np.column_stack([np.ones(len(tr)), pcs[tr, :max_pc] - centre])
        Q, R = np.linalg.qr(A)
        z = Q.T @ y[tr]
        B = pcs[te, :max_pc] - centre
        for k in range(max_pc + 1):
            beta = np.linalg.solve(R[:k + 1, :k + 1], z[:k + 1]) if k else z[:1] / R[0, 0]
            out[k, te] = y[te] - (B[:, :k] @ beta[1:] if k else 0.0) - beta[0] + mean
    out[0, idx] = y[idx]                                   # no PCs is no adjustment: the values as they came
    return out


def sweep(table: dict[str, np.ndarray], pcs: np.ndarray, max_pc: int, truths: dict[str, np.ndarray],
          classes: list[str], samples: list[str] | None = None, trios=None, population=None,
          folds: int = 10, n_boot: int = 300, seed: int = 1) -> list[dict]:
    """For k = 0..max_pc coverage PCs regressed out of log(estimate), cross-validated:

    * known-truth columns (`truths[col]` is the per-sample truth; NaN where there is none):
      the root-mean-square and the robust SD of log(adjusted / truth);
    * class columns: the spread that is left, and - with `trios` - the transmission reliability
      of the adjusted values, its bootstrap interval, and its paired difference from k = 0.
    """
    from . import trios as T
    rows = []
    mad_sd = lambda v: float(1.4826 * np.median(np.abs(v - np.median(v)))) if len(v) else float("nan")
    for col in list(truths) + [c for c in classes if c not in truths]:
        y = np.asarray(table[col], float)
        truth = truths.get(col)
        logy = np.where(y > 0, np.log(np.where(y > 0, y, 1.0)), np.nan)
        if truth is not None:
            logy = np.where(np.isfinite(truth) & (truth > 0), logy, np.nan)
        adj = cv_adjusted(logy, pcs, max_pc, folds, seed)
        base = None
        for k in range(max_pc + 1):
            v = adj[k]
            fin = np.isfinite(v)
            row = dict(column=col, kind="truth" if truth is not None else "class", n_pc=k, n=int(fin.sum()))
            if not fin.any():
                rows.append(row)
                continue
            if truth is not None:
                d = v[fin] - np.log(truth[fin])
                row.update(rmse_log=float(np.sqrt(np.mean(d ** 2))), sd_log_robust=mad_sd(d), bias_log=float(np.mean(d)))
            else:
                row.update(sd_log=float(np.std(v[fin], ddof=1)), sd_log_robust=mad_sd(v[fin]))
            if trios and samples is not None and truth is None:
                vals = {s: float(np.exp(x)) for s, x in zip(samples, v) if np.isfinite(x)}
                try:
                    t = T.transmission(vals, trios, population, n_perm=0, n_boot=n_boot, seed=seed)
                    row.update(n_trios=t["n_trios"], R_midparent=t["reliability_midparent"], spousal_r=t["spousal_r"])
                    if "reliability_midparent_ci95" in t:
                        row.update(R_lo=t["reliability_midparent_ci95"][0], R_hi=t["reliability_midparent_ci95"][1])
                    if k == 0:
                        base = vals
                    elif base is not None and t["n_trios"] >= 20:
                        c = T.compare(vals, base, trios, population, n_boot=n_boot, seed=seed)
                        row.update(dR_vs_0=c["delta"], dR_lo=c["ci95"][0], dR_hi=c["ci95"][1])
                except ValueError:
                    pass
            rows.append(row)
    return rows


def recommend(rows: list[dict]) -> dict[str, dict]:
    """Per column, the number of PCs the sweep supports, by the one-standard-error rule: the fewest
    PCs that do as well as the best number does, to within the sampling error of "best" - for a
    known truth, the cross-validated robust SD of log(estimate / truth); for a class with trios,
    the transmission reliability (standard error from its bootstrap interval). A minimum at 8 PCs
    that is 3% below the value at 0 PCs in 66 samples is noise, and this rule says 0."""
    out: dict[str, dict] = {}
    for col in dict.fromkeys(r["column"] for r in rows):
        rs = [r for r in rows if r["column"] == col]
        if rs[0]["kind"] == "truth" and all("sd_log_robust" in r for r in rs):
            e = np.array([r["sd_log_robust"] for r in rs])
            best = int(np.nanargmin(e))
            se = e[best] * 1.25 / np.sqrt(2 * max(rs[best]["n"], 2))          # SE of a MAD-based SD
            pick = int(np.where(e <= e[best] + se)[0][0])
            out[col] = dict(kind="truth", best=best, pick=pick, at_0=float(e[0]), at_pick=float(e[pick]), at_best=float(e[best]), se=float(se))
        elif all("R_midparent" in r for r in rs):
            R = np.array([r["R_midparent"] for r in rs])
            best = int(np.nanargmax(R))
            se = (rs[best]["R_hi"] - rs[best]["R_lo"]) / 3.92 if "R_hi" in rs[best] else 0.0
            pick = int(np.where(R >= R[best] - se)[0][0])
            out[col] = dict(kind="class", best=best, pick=pick, at_0=float(R[0]), at_pick=float(R[pick]), at_best=float(R[best]), se=float(se))
    return out
