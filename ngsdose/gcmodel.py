"""Fragment-GC rate model.

The engine tabulates, over single-copy control regions, the number of position-strands N[g] whose
downstream fragment-length window has GC content g, and the number of fragment 5' ends O[g]
observed at those position-strands. The rate

    lambda(g) = E[5' ends per position-strand | window GC = g]        (diploid, two copies)

is what a multi-copy class is compared against: a class present in C copies per diploid genome
yields C/2 * lambda(g(p)) ends at unit position p. Benjamini & Speed (NAR 2012) showed that the
GC of the whole fragment, not of the read, is what predicts Illumina coverage; this repository
re-derived that on 1000 Genomes 30x data (docs/DESIGN.md, section 2, finding 2; the model itself
is described in section 6).

The curve is a Poisson regression of O on a natural cubic spline in g with offset log N.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

GC_BINS = 101


def _ncs_basis(x: np.ndarray, knots: np.ndarray) -> np.ndarray:
    """Natural cubic spline basis (Hastie, Tibshirani & Friedman, eq. 5.4); linear beyond the end knots."""
    K = len(knots)

    def d(k):
        return (np.clip(x - knots[k], 0, None) ** 3 - np.clip(x - knots[K - 1], 0, None) ** 3) / (knots[K - 1] - knots[k])

    cols = [np.ones_like(x), x]
    for k in range(K - 2):
        cols.append(d(k) - d(K - 2))
    return np.column_stack(cols)


@dataclass
class GCCurve:
    L: int
    rate: np.ndarray            # lambda(g) per 1% GC bin; NaN outside the supported range
    se_log: np.ndarray          # standard error of log lambda(g), overdispersion included
    lo: int                     # supported bin range [lo, hi], inclusive
    hi: int
    scale: float                # pooled rate: sum(O) / sum(N)
    dispersion: float           # Pearson chi2 / df of the binned fit
    n_positions: int
    n_ends: int
    rescale: float = 1.0        # robust scale factor applied after control-region trimming
    notes: list[str] = field(default_factory=list)

    def relative(self) -> np.ndarray:
        """Rate relative to the pooled mean: the conventional 'GC bias curve'."""
        return self.rate / self.scale

    def lookup(self, gc_count: np.ndarray) -> np.ndarray:
        """lambda for windows holding `gc_count` G/C bases out of L; NaN outside support."""
        b = (gc_count.astype(np.int64) * 100 + self.L // 2) // self.L
        return self.rate[b] * self.rescale


def _dropout(N, O, use, k: int = 3, limit: float = 50.0):
    """The first run of fitted GC bins with no fragment end whose expected ends exceed `limit`,
    as (first bin, last bin, expected ends), or None. Expected ends are the run's positions times
    the pooled rate of the `k` nearest fitted bins with ends on each side. The neighbours overstate
    it where the rate falls steeply at the edge of the support, hence the wide margin: thinned
    1000 Genomes controls (0.1-1.5x) reach at most about 12 in bins that simply got no end."""
    idx = np.where(use)[0]
    held = idx[O[idx] > 0]
    i = 0
    while i < len(idx):
        if O[idx[i]] > 0:
            i += 1
            continue
        j = i
        while j + 1 < len(idx) and O[idx[j + 1]] == 0:
            j += 1
        a, b = idx[i], idx[j]
        near = np.r_[held[held < a][-k:], held[held > b][:k]]
        if len(near):
            expect = N[a:b + 1][use[a:b + 1]].sum() * O[near].sum() / N[near].sum()
            if expect > limit:
                return int(a), int(b), float(expect)
        i = j + 1
    return None


def fit_gc_curve(N, O, L: int, knot_step: float = 0.05, max_se: float = float("inf"), min_positions: int = 2000) -> GCCurve:
    """Fit lambda(g).

    The supported GC range is set by the control *positions*, not by how many reads fell on them:
    it runs from the lowest to the highest 1% bin with at least `min_positions` control
    position-strands (bins inside it with fewer are left out of the fit and take the spline's
    interpolated value), so that the set of usable class windows is the same at every depth. A
    support that narrows with depth silently changes which windows an estimate averages over:
    with the original rule (log lambda known to 5%) the all-window 45S estimate of one sample
    rose by 2.4% at 1x and 6.5% at 0.4x as its GC-rich, low-efficiency windows dropped out. How
    well the curve is determined is reported instead (`se_log`; the estimate carries its maximum
    over the support as `gc_curve_max_se`); `max_se` can still be set to trim the support by
    precision when that is wanted."""
    N = np.asarray(N, float)
    O = np.asarray(O, float)
    if N.shape != (GC_BINS,) or O.shape != (GC_BINS,):
        raise ValueError("N and O must have 101 entries (1% GC bins)")
    if O.sum() < 10000:
        raise ValueError(f"only {int(O.sum())} control fragment ends: too few to fit a GC curve")
    g = np.arange(GC_BINS) / 100.0
    use = N >= min_positions
    idx = np.where(use)[0]
    if not len(idx):
        raise ValueError(f"no 1% GC bin holds {min_positions} control position-strands (the largest holds {int(N.max())}): "
                         "the controls are too small to fit a GC curve")
    lo, hi = idx.min(), idx.max()
    use &= (np.arange(GC_BINS) >= lo) & (np.arange(GC_BINS) <= hi)
    knots = np.arange(np.ceil(g[lo] / knot_step) * knot_step, g[hi] + 1e-9, knot_step)
    knots = np.unique(np.r_[g[lo], knots, g[hi]])
    if len(knots) < 4:
        raise ValueError("control GC range too narrow for a spline fit")
    gap = _dropout(N, O, use)
    if gap is not None:
        a, b, expect = gap
        where = f"GC bin {a}% holds" if a == b else f"GC bins {a}-{b}% hold"
        raise ValueError(f"GC curve fit at L={L}: the {where} control positions but no fragment end, where "
                         f"their neighbours give about {expect:.0f} (complete GC dropout): the curve cannot be fitted there")
    B = _ncs_basis(g, knots)
    Bu, Nu, Ou = B[use], N[use], O[use]
    beta = np.zeros(B.shape[1])
    beta[0] = np.log(Ou.sum() / Nu.sum())
    ridge = 1e-8 * np.eye(B.shape[1])
    # the linear predictor is clipped so that a bin with no ends at all (complete GC dropout) cannot
    # drive the rate to zero and the Hessian singular
    mean = lambda b: Nu * np.exp(np.clip(Bu @ b, -30.0, 30.0))
    dev = lambda m: 2.0 * np.sum(np.where(Ou > 0, Ou * np.log(np.maximum(Ou, 1e-300) / m), 0.0) - (Ou - m))
    mu = mean(beta)
    d = dev(mu)
    converged = False
    for _ in range(100):
        H = Bu.T @ (Bu * mu[:, None]) + ridge
        try:
            step = np.linalg.solve(H, Bu.T @ (Ou - mu))
        except np.linalg.LinAlgError as e:
            raise ValueError(f"GC curve fit failed at L={L} ({e})") from e
        # Newton steps are halved while they increase the Poisson deviance (by more than rounding:
        # near the optimum the steps are noise of about 1e-9 and the deviance moves by 1e-10)
        for _ in range(30):
            mu_new = mean(beta + step)
            d_new = dev(mu_new)
            if d_new <= d + 1e-8 * max(d, 1.0):
                break
            step = step / 2
        beta = beta + step
        mu, d = mu_new, d_new
        if np.max(np.abs(step)) < 1e-10:
            converged = True
            break
    # Near the optimum the coefficients can wander along a flat direction (rounding noise at 30x)
    # while the fitted rates stay put: that is a fit. A bin with no ends lets its rate slide toward
    # zero for ever, so the test is on the bins that hold ends: rates still moving by 1% an
    # iteration after 100 are not a fit (bins with no ends at all were checked above).
    if not np.all(np.isfinite(beta)):
        raise ValueError(f"GC curve fit failed at L={L}: non-finite coefficients")
    held = Ou > 0
    moving = float(np.max(np.abs(Bu[held] @ step)))
    if not converged and moving > 0.01:
        raise ValueError(f"GC curve fit did not converge at L={L} (deviance {d:.4g}, fitted log rates still moving by {moving:.3g})")
    mu = mean(beta)
    df = max(1, use.sum() - B.shape[1])
    dispersion = float(np.sum((Ou - mu) ** 2 / np.maximum(mu, 1e-9)) / df)
    cov = np.linalg.inv(Bu.T @ (Bu * mu[:, None]) + ridge) * max(dispersion, 1.0)
    se = np.sqrt(np.einsum("ij,jk,ik->i", B, cov, B))
    rate = np.exp(B @ beta)
    ok = (se <= max_se) & (np.arange(GC_BINS) >= lo) & (np.arange(GC_BINS) <= hi)
    # supported range: the contiguous stretch around the mode of N that satisfies the SE bound
    mode = int(np.argmax(N))
    a = mode
    while a - 1 >= 0 and ok[a - 1]:
        a -= 1
    b = mode
    while b + 1 < GC_BINS and ok[b + 1]:
        b += 1
    rate_out = np.full(GC_BINS, np.nan)
    rate_out[a:b + 1] = rate[a:b + 1]
    se_out = np.full(GC_BINS, np.nan)
    se_out[a:b + 1] = se[a:b + 1]
    return GCCurve(L=L, rate=rate_out, se_log=se_out, lo=a, hi=b, scale=float(O.sum() / N.sum()),
                   dispersion=dispersion, n_positions=int(N.sum()), n_ends=int(O.sum()))


def window_gc_counts(seq: str, L: int, circular: bool) -> tuple[np.ndarray, np.ndarray]:
    """G/C count of the L-window anchored at every position, for forward and reverse 5' ends.

    forward 5' end at p  -> window [p, p+L)
    reverse 5' end at x  -> window [x-L+1, x]
    Positions whose window is undefined (linear sequence edge, or a non-ACGT base) get -1.
    """
    U = len(seq)
    a = np.frombuffer(seq.upper().encode(), dtype=np.uint8)
    if circular:
        c = -(-L // U)                  # copies on each side, enough for a window longer than the unit
        a = np.tile(a, 2 * c + 1)
        off = U * c
    else:
        off = 0
    isgc = (a == 71) | (a == 67)
    isn = ~((a == 65) | (a == 67) | (a == 71) | (a == 84))
    cg = np.r_[0, np.cumsum(isgc)]
    cn = np.r_[0, np.cumsum(isn)]
    p = np.arange(U) + off
    fwd = np.full(U, -1, np.int64)
    rev = np.full(U, -1, np.int64)
    okf = p + L <= len(a)
    fwd[okf] = cg[p[okf] + L] - cg[p[okf]]
    fwd[okf] = np.where(cn[p[okf] + L] - cn[p[okf]] > 0, -1, fwd[okf])
    okr = p + 1 - L >= 0
    rev[okr] = cg[p[okr] + 1] - cg[p[okr] + 1 - L]
    rev[okr] = np.where(cn[p[okr] + 1] - cn[p[okr] + 1 - L] > 0, -1, rev[okr])
    return fwd, rev


def expected_per_position(seq: str, curve: GCCurve, circular: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """lambda at every unit position for both strands, plus the window GC fractions."""
    gf, gr = window_gc_counts(seq, curve.L, circular)
    ef = np.where(gf >= 0, curve.lookup(np.maximum(gf, 0)), np.nan)
    er = np.where(gr >= 0, curve.lookup(np.maximum(gr, 0)), np.nan)
    return ef, er, gf / curve.L, gr / curve.L
