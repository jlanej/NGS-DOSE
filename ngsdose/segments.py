"""Integer copy states along a positional unit, from one genome's calibrated window profile.

A unit present in N copies reads N in every window, and the cohort layer makes the windows
comparable (each divided by its efficiency). Copies that hold only part of the unit, or lack part
of it, leave the profile piecewise constant at whole numbers: a partial copy of the first 316 kb
of the distal junction reads N + 1 up to 316 kb and N beyond. `segment` finds the runs and their
integers; `describe` turns them into complete copies, partial copies and local gains and losses.

The model is a hidden Markov chain over the windows in unit order. Its states are the integers
1..max_state; a window in state k is expected at k * f, where f is the genome's own scale (its
level is known to a percent or two, so an integer's worth of windows sits at 0.98 or 1.02 of the
integer, not at it); the window's deviation is judged against its noise, with a cap so that a
single wild window cannot buy a change of state; every change of state costs `tau`. The most
probable path is found for each f of a grid, a Gaussian prior of SD `scale_sd` about 1 keeps f from
explaining a copy away, and the best f wins. Runs shorter than `min_windows` are given to the
neighbour that fits them better.

Some genomes' profiles lean: they rise or fall smoothly along the unit by a few percent, some by ten
(a property of the library or of the culture; batch, depth and GC response do not account for it). A chain
of whole numbers would break a lean into a step, so the lean is part of the model: `tilt`, the log
change across the unit. The leans of a grid are scored like the scales are, the chain found for
each pair, with a Gaussian prior of SD `TILT_SD`, and the best pair is kept. A smooth rise of one
copy across the unit costs less as a lean than as a change of state; a step, whose windows a lean
fits badly on both sides, stays a step; a step on a leaning profile is both.

A genome whose level lies between two whole numbers throughout has two readings, N copies at a
scale a few percent high or N + 1 at one a few percent low, and the data cannot choose. The call
says so: `alternative` is the best reading with other complete copies and `gap` how much less
probable it is (log units); a call with a gap below `MIN_GAP` is `uncertain`.

Whole numbers are what a germ line holds. A change in part of the cells (a culture that is
losing a junction, a rearrangement in some of a donor's blood) leaves the profile a fraction of a
copy off them, and a chain of whole numbers would hide it: the nearest whole number is called
and the difference goes into the scale. So the difference is measured and kept (`find_fractions`).
A genome's scale is compared with the scales of the cohort: `off` is its level less its whole
number, in copies, and `off_z` the same in the cohort's robust SDs. And steps of fractional
height are looked for in what the whole numbers leave: a level per stretch and a lean are fitted
together (a step and a lean are partly alike), a step is proposed where it improves the fit as
much as a change of state must, and it is kept where it stands `EVENT_Z` robust SDs from what
the same fit finds in the cohort's other genomes, whose profiles wander more than counting
alone allows. A call that is not uncertain is then `settled`, or `fractional` where the level
or a stretch sits off its whole number. What a fraction is, the measurement cannot say: a
change in part of the cells, or a library unlike the cohort's.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

TAU = 16.0                     # cost of a change of state, in units of log-likelihood
MIN_WINDOWS = 16               # the shortest run kept (about 5 kb of a unit with three quarters of its windows usable)
Z_CAP = 3.5                    # a window further than this many SDs from a state counts as this far
SCALE_SD = 0.02                # prior SD of a genome's scale about 1
SCALE_GRID = tuple(np.round(np.arange(0.93, 1.0701, 0.005), 4))
MAX_STATE = 24
OFF_INTEGER = 0.35             # a segment whose mean is further than this from its state (and 3 SE) is flagged
MIN_CORE = 40000               # a segment at least this long, outside the polymorphic intervals, counts towards the complete copies
MIN_GAP = 3.0                  # a call whose best alternative is less than this many log units behind is uncertain
TILT_SD = 0.03                 # prior SD of a profile's lean (log change across the unit)
TILT_MAX = 0.15                # the leans tried lie within this of none
LEVEL_Z = 3.0                  # a level this many of the cohort's robust SDs from its whole number is fractional
EVENT_Z = 4.0                  # a step of fractional height is kept this many robust SDs from what the cohort's other genomes show there
FRACTION_MIN_HEIGHT = 0.25     # and at least this many copies high
FRACTION_MIN_WINDOWS = 100     # with at least this many windows of the core on either side (about 30 kb of the junction)
FRACTION_GRID = 16             # windows between the places where a step is tried
FRACTION_MAX_STEPS = 4
FRACTION_MIN_SAMPLES = 50      # genomes a cohort needs for its own spread to judge a fraction by


@dataclass
class Segment:
    start: int
    end: int
    windows: int
    state: int
    mean: float                # the windows' weighted mean, in copies on the genome's scale (divided by f)
    se: float
    raw: float = float("nan")  # the same mean as the reads give it: neither levelled nor divided by the scale

    @property
    def off_integer(self) -> bool:
        return abs(self.mean - self.state) > max(OFF_INTEGER, 3 * self.se)

    def as_dict(self) -> dict:
        return dict(start=self.start, end=self.end, windows=self.windows, state=self.state, mean=round(self.mean, 3), se=round(self.se, 3),
                    raw=None if not np.isfinite(self.raw) else round(self.raw, 3))


@dataclass
class Event:
    delta: int                 # copies relative to the complete copies
    start: int
    end: int
    kind: str                  # "partial copy" (a gain that reaches an end of the unit), "gain", "loss"

    def __str__(self) -> str:
        return f"{self.delta:+d}:{self.start // 1000}-{self.end // 1000}kb"


@dataclass
class Fraction:
    """A stretch of the unit that reads a fraction of a copy off the whole number called there."""
    start: int
    end: int
    state: int                 # the whole number called over most of it
    offset: float              # copies above (+) or below (-) that number
    z: float                   # the stretch against the rest, in robust SDs of the same contrast in the cohort's other genomes
    windows: int               # the windows of the core it rests on

    def height(self, copies: int) -> float:
        """Relative to the copies the genome is described against."""
        return self.state - copies + self.offset

    def label(self, copies: int) -> str:
        return f"{self.height(copies):+.2f}:{self.start // 1000}-{self.end // 1000}kb"


@dataclass
class Call:
    segments: list[Segment]
    scale: float               # f
    sigma: float               # the genome's window noise, copies
    copies: int | None = None  # the copies the unit is described against: `expected` where the genome holds it, else its commonest state
    events: list[Event] = field(default_factory=list)
    score: float = float("nan")
    tilt: float = 0.0          # the profile's lean: log change across the unit, taken out before the states were called
    readings: list = field(default_factory=list)      # (score, scale, segments) of every scale of the grid within reach of the best
    alternative: dict | None = None                   # the best reading with other copies: copies, scale, gap
    gap: float = float("inf")
    ratio: float | None = None                        # the level over the whole numbers called, measured (no prior): the genome's scale
    off: float | None = None                          # the level less its whole number, copies
    off_z: float | None = None                        # the same in robust SDs of the cohort's scales
    fractions: list[Fraction] = field(default_factory=list)
    level_z: float = LEVEL_Z

    @property
    def uncertain(self) -> bool:
        """Another whole number of copies explains the profile nearly as well (a level between two whole numbers)."""
        return self.gap < MIN_GAP

    @property
    def fractional(self) -> bool:
        """The nearest whole numbers are clear and the genome is not on them: its level, or a stretch of its unit."""
        return not self.uncertain and (bool(self.fractions) or (self.off_z is not None and abs(self.off_z) >= self.level_z))

    @property
    def status(self) -> str:
        return "uncertain" if self.uncertain else "fractional" if self.fractional else "settled"

    @property
    def scale_uncertain(self) -> bool:
        return self.uncertain

    def state_at(self, pos: int) -> int | None:
        for s in self.segments:
            if s.start <= pos < s.end:
                return s.state
        return None


def window_noise(values: np.ndarray) -> float:
    """A genome's window noise from the differences of successive windows (steps between segments do not enter a median)."""
    d = np.diff(values)
    if len(d) < 4:
        return float("nan")
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2))


def _runs(path: np.ndarray) -> list[list[int]]:
    cut = np.flatnonzero(np.diff(path)) + 1
    lo = np.r_[0, cut]
    hi = np.r_[cut, len(path)] - 1
    return [[int(a), int(b), int(path[a])] for a, b in zip(lo, hi)]


def _merge_short(path: np.ndarray, x: np.ndarray, sigma: np.ndarray, mu: np.ndarray, min_windows: int) -> np.ndarray:
    path = path.copy()
    while True:
        R = _runs(path)
        short = [r for r in R if r[1] - r[0] + 1 < min_windows]
        if not short or len(R) == 1:
            return path
        i, j, s = min(short, key=lambda r: r[1] - r[0])
        k = R.index([i, j, s])
        cand = [R[k + d][2] for d in (-1, 1) if 0 <= k + d < len(R)]
        cost = [float(np.sum(np.minimum(((x[i:j + 1] - mu[c]) / sigma[i:j + 1]) ** 2, Z_CAP ** 2))) for c in cand]
        path[i:j + 1] = cand[int(np.argmin(cost))]


def segment(starts: np.ndarray, values: np.ndarray, rel_noise: np.ndarray | None = None, window: int = 250, unit_length: int | None = None,
            tau: float = TAU, min_windows: int = MIN_WINDOWS, scale_grid=SCALE_GRID, scale_sd: float = SCALE_SD,
            max_state: int = MAX_STATE, tilt: bool = True, tilt_sd: float = TILT_SD) -> Call | None:
    """`starts`, `values`: a genome's windows (unit coordinate, calibrated copies; NaN where unusable); `rel_noise`: each
    window's noise relative to the typical window's (1 where not given). With `tilt` the profile's lean is a parameter
    of the chain, like the scale (the segments' means are then of the levelled profile). None when fewer than
    2 * min_windows windows have a value. The window noise has a floor of a thousandth of the profile's median."""
    starts, values = np.asarray(starts), np.asarray(values, float)
    ok = np.isfinite(values) & (values > 0)
    if ok.sum() < 2 * min_windows:
        return None
    s, v = starts[ok], values[ok]
    r = np.ones(len(v)) if rel_noise is None else np.where(np.isfinite(np.asarray(rel_noise, float)[ok]), np.asarray(rel_noise, float)[ok], 1.0)
    sig = window_noise(v)
    if not np.isfinite(sig):
        return None
    sig = max(sig, 1e-3 * float(np.median(v)))          # a profile without noise (a synthetic one) is still a profile
    sigma = sig * r
    span = float(unit_length or (s[-1] + window))
    t = (s + window / 2) / span - 0.5                                       # position along the unit, -0.5 .. 0.5
    F = np.asarray(scale_grid, float)
    K = np.arange(1, max(int(np.ceil(np.nanmax(v) / (F.min() * np.exp(-0.5 * TILT_MAX)))) + 2, 4))
    K = K[K <= max_state].astype(float)
    n, nk = len(v), len(K)
    ar = np.arange(nk)[None, :]

    def solve(G: np.ndarray, Fc: np.ndarray):
        """The best chain's score for every candidate (lean G[c], scale Fc[c]), with what is needed to trace it."""
        nc = len(G)
        back = np.zeros((n, nc, nk), np.int8)
        rows = np.arange(nc)

        def loglik(q):
            z = (v[q] - (Fc * np.exp(G * t[q]))[:, None] * K[None, :]) / sigma[q]
            return -0.5 * np.minimum(z * z, Z_CAP * Z_CAP)
        score = loglik(0)
        for q in range(1, n):
            best = np.argmax(score, axis=1)
            move = score[rows, best] - tau
            take = move[:, None] > score
            back[q] = np.where(take, best[:, None], ar)
            score = np.where(take, move[:, None], score) + loglik(q)
        total = score.max(axis=1) - 0.5 * ((Fc - 1) / scale_sd) ** 2 - 0.5 * (G / tilt_sd) ** 2
        return total, score, back

    def reading(c, G, Fc, score, back) -> list[Segment]:
        path = np.zeros(n, int)
        path[-1] = int(np.argmax(score[c]))
        for q in range(n - 1, 0, -1):
            path[q - 1] = back[q, c, path[q]]
        lv = v * np.exp(-G[c] * t)                                           # the profile levelled
        path = _merge_short(path, lv, sigma, Fc[c] * K, min_windows)
        segs = []
        for a, b, k in _runs(path):
            w = 1 / sigma[a:b + 1] ** 2
            segs.append(Segment(int(s[a]), int(s[b]) + window, b - a + 1, int(K[k]), float(np.sum(w * lv[a:b + 1]) / w.sum()) / Fc[c], float(np.sqrt(1 / w.sum())) / Fc[c],
                                raw=float(np.sum(w * v[a:b + 1]) / w.sum())))
        for a, b in zip(segs[:-1], segs[1:]):                               # a breakpoint lies midway between the windows either side
            a.end = b.start = (a.end + b.start) // 2
        segs[0].start = 0
        if unit_length:
            segs[-1].end = int(unit_length)
        return segs
    g_best = 0.0
    if tilt:
        # first the lean, coarsely, on every second scale; then the scales in full about the best lean
        g1, f1 = np.arange(-TILT_MAX, TILT_MAX + 1e-9, 0.03), F[::2]
        G, Fc = np.repeat(g1, len(f1)), np.tile(f1, len(g1))
        total, _, _ = solve(G, Fc)
        g_best = float(G[int(np.argmax(total))])
        g2 = g_best + np.array([-0.02, -0.01, 0.0, 0.01, 0.02])
    else:
        g2 = np.array([0.0])
    G, Fc = np.repeat(g2, len(F)), np.tile(F, len(g2))
    total, score, back = solve(G, Fc)
    c_best = int(np.argmax(total))
    same = np.flatnonzero((G == G[c_best]) & (total >= total[c_best] - 25.0))   # the other scales, under the same lean
    readings = sorted(((float(total[c]), float(Fc[c]), reading(c, G, Fc, score, back)) for c in same), key=lambda x: -x[0])
    return Call(readings[0][2], readings[0][1], sig, score=readings[0][0], tilt=float(G[c_best]), readings=readings)


def _short_polymorphic(g, poly, slack: int = 5000) -> bool:
    """A segment that lies within a polymorphic interval (to within `slack` either side)."""
    return any(g.start >= a - slack and g.end <= b + slack for a, b in poly)


def _copies(segs: list[Segment], poly, min_core: int, expected: int | None) -> int:
    """The copies a profile is described against: `expected` where a segment of the core holds it (ten junctions are the
    norm, and a genome with a partial copy or a partial loss still holds ten over part of the unit); otherwise the state
    that holds most of the core."""
    cs = [g for g in segs if g.end - g.start >= min_core and not _short_polymorphic(g, poly)] or [max(segs, key=lambda g: g.end - g.start)]
    if expected is not None and any(g.state == expected for g in cs):
        return int(expected)
    length: dict[int, int] = {}
    for g in cs:
        length[g.state] = length.get(g.state, 0) + g.end - g.start
    return int(max(length, key=lambda k: (length[k], -abs(k - (expected or k)))))


def _bulk(segs: list[Segment], poly) -> int:
    """The state that holds most of the unit outside the polymorphic intervals."""
    length: dict[int, int] = {}
    for g in segs:
        if not _short_polymorphic(g, poly):
            length[g.state] = length.get(g.state, 0) + g.end - g.start
    if not length:
        return segs[0].state
    return int(max(length, key=length.get))


def describe(call: Call, polymorphic: list[tuple[int, int]] | None = None, min_core: int = MIN_CORE, expected: int | None = None) -> Call:
    """Copies and events. Every segment at another state than the copies is an event: a gain that reaches an end of the
    unit (through polymorphic segments or other gains) is a partial copy, a loss that does a partial loss (a copy that
    lacks that end), anything else a local gain or loss. The best reading with other copies, and how far behind it
    is, say whether the whole numbers are settled."""
    poly = polymorphic or []
    segs = call.segments
    n = _copies(segs, poly, min_core, expected)
    events: list[Event] = []
    unit_end = segs[-1].end
    for g in segs:
        d = g.state - n
        if d == 0:
            continue
        if events and events[-1].delta == d and events[-1].end == g.start:
            events[-1].end = g.end
        else:
            events.append(Event(d, g.start, g.end, ""))
    in_poly = lambda g: _short_polymorphic(g, poly)
    for e in events:
        e.kind = "loss" if e.delta < 0 else "gain"
        if e.end - e.start < min_core:
            continue
        sign = 1 if e.delta > 0 else -1
        to_start = e.start == 0 or all(in_poly(g) or (g.state - n) * sign > 0 for g in segs if g.end <= e.start)
        to_end = e.end == unit_end or all(in_poly(g) or (g.state - n) * sign > 0 for g in segs if g.start >= e.end)
        if to_start or to_end:
            e.kind = "partial copy" if sign > 0 else "partial loss"
    call.copies, call.events = int(n), events
    call.alternative, call.gap = None, float("inf")
    bulk = _bulk(segs, poly)
    for sc, f, other in call.readings[1:]:
        if _bulk(other, poly) != bulk:                                       # other whole numbers over most of the core
            call.alternative = dict(copies=_copies(other, poly, min_core, expected), bulk=_bulk(other, poly), scale=round(f, 3), gap=round(call.score - sc, 2))
            call.gap = float(call.score - sc)
            break
    return call


def partial_copies(call: Call, polymorphic: list[tuple[int, int]] | None = None, min_core: int = MIN_CORE, kind: str = "partial copy") -> list[tuple[int, int, int]]:
    """The partial copies (or, with `kind`, the partial losses) as (count, start, end), joining the pieces of one that
    polymorphic segments interrupt."""
    poly = polymorphic or []
    out: list[list[int]] = []
    for e in call.events:
        if e.kind != kind:
            continue
        if out and out[-1][0] == e.delta and any(a - 5000 <= out[-1][2] <= b + 5000 and a - 5000 <= e.start <= b + 5000 for a, b in poly):
            out[-1][2] = e.end
        else:
            out.append([e.delta, e.start, e.end])
    return [tuple(x) for x in out]


def events_string(call: Call) -> str:
    return ";".join(str(e) for e in call.events)


def states_at(call: Call, size: int, unit_length: int) -> np.ndarray:
    """The call's state per block of `size`: that of the segment holding most of the block."""
    n = unit_length // size
    out = np.zeros(n, int)
    for k in range(n):
        a0, b0 = k * size, (k + 1) * size
        out[k] = max(call.segments, key=lambda g: max(0, min(g.end, b0) - max(g.start, a0))).state
    return out


# ------------------------------------------------------------------ off the whole numbers
class _Sums:
    """The running sums a fit of levels needs: of the weights, of weight x value and of weight x each smooth term (a
    leading zero each), and the totals among the smooth terms. `X`: windows x terms (the position along the unit, for
    the lean; the window's GC, for a library whose GC response the model left in)."""

    def __init__(self, y: np.ndarray, w: np.ndarray, X: np.ndarray):
        z = lambda a: np.concatenate([np.zeros((1,) + a.shape[1:]), np.cumsum(a, axis=0)])
        self.n, self.q = len(y), X.shape[1]
        self.w, self.wy, self.wx = z(w), z(w * y), z(w[:, None] * X)
        self.xx, self.xy, self.yy = X.T @ (w[:, None] * X), X.T @ (w * y), float(np.sum(w * y * y))


def fit_levels(S: _Sums, cuts, prior_sd) -> tuple[np.ndarray, np.ndarray, float]:
    """A level per stretch and one coefficient per smooth term, by weighted least squares. `cuts`: the windows at which
    a new stretch begins; `prior_sd`: the SD of each smooth term's Gaussian prior about none. Returns (levels,
    coefficients, log-likelihood)."""
    e = np.r_[0, np.asarray(cuts, int), S.n]
    n, q = len(e) - 1, S.q
    M = np.zeros((n + q, n + q))
    b = np.zeros(n + q)
    d = np.arange(n)
    M[d, d] = S.w[e[1:]] - S.w[e[:-1]]
    M[:n, n:] = S.wx[e[1:]] - S.wx[e[:-1]]
    M[n:, :n] = M[:n, n:].T
    M[n:, n:] = S.xx + np.diag(1 / np.asarray(prior_sd, float) ** 2)
    M[d, d] += 1e-9 * max(float(S.w[-1]), 1.0)                         # a stretch without a usable window has no level of its own: it takes none
    b[:n] = S.wy[e[1:]] - S.wy[e[:-1]]
    b[n:] = S.xy
    beta = np.linalg.solve(M, b)
    return beta[:n], beta[n:], -0.5 * float(S.yy - b @ beta)


def propose_steps(y: np.ndarray, w: np.ndarray, X: np.ndarray, copies: float, prior_sd, also=(), tau: float = TAU,
                  min_height: float = FRACTION_MIN_HEIGHT, min_windows: int = FRACTION_MIN_WINDOWS, grid: int = FRACTION_GRID,
                  max_steps: int = FRACTION_MAX_STEPS) -> list[int]:
    """Where a genome's profile steps by a fraction of a copy. `y`: the log of each window over its called whole number;
    `w`: its weight; `X`: the smooth terms (see `_Sums`), fitted beside the levels because a step and a lean are partly
    alike. Steps are added one at a time, the one that improves the fit most, for as long as it improves it by `tau`
    (what a change of state costs the chain of whole numbers) and is `min_height` copies high. Where no single step
    does, a pair is tried, for a stretch inside the unit that one step describes badly: it must improve the fit by
    twice `tau` and stand `min_height` from both its neighbours. A step is tried every `grid` windows (a pair every
    second one of them) and at the windows of `also` (where the whole numbers change), never within `min_windows` of
    an end or of another step. Returns the windows at which a new stretch begins."""
    n = len(y)
    if n < 2 * min_windows:
        return []
    S = _Sums(y, w, X)
    cand = sorted({int(p) for p in list(range(min_windows, n - min_windows + 1, grid)) + list(also) if min_windows <= p <= n - min_windows})
    pairs = sorted(set(cand[::2]) | {int(p) for p in also if min_windows <= p <= n - min_windows})
    cuts: list[int] = []
    _, _, ll = fit_levels(S, cuts, prior_sd)
    high = lambda lv, j: copies * abs(np.exp(lv[j + 1]) - np.exp(lv[j])) >= min_height
    free = lambda p, cs: all(abs(p - c) >= min_windows for c in cs)
    while len(cuts) < max_steps:
        best = None
        for p in cand:
            if not free(p, cuts):
                continue
            cs = sorted(cuts + [p])
            lv, _, l2 = fit_levels(S, cs, prior_sd)
            if best is None or l2 > best[0]:
                best = (l2, cs, lv, cs.index(p))
        if best is not None and best[0] - ll >= tau and high(best[2], best[3]):
            cuts, ll = best[1], best[0]
            continue
        best = None
        if len(cuts) + 2 <= max_steps:
            for a, p in enumerate(pairs):
                if not free(p, cuts):
                    continue
                for q in pairs[a + 1:]:
                    if q - p < min_windows or not free(q, cuts):
                        continue
                    cs = sorted(cuts + [p, q])
                    lv, _, l2 = fit_levels(S, cs, prior_sd)
                    if best is None or l2 > best[0]:
                        best = (l2, cs, lv, cs.index(p))
        if best is None or best[0] - ll < 2 * tau or not (high(best[2], best[3]) and high(best[2], best[3] + 1)):
            break
        cuts, ll = best[1], best[0]
    return cuts


def find_fractions(calls: list, starts: np.ndarray, values: np.ndarray, rel_noise: np.ndarray | None = None, window: int = 250,
                   unit_length: int | None = None, leave_out=None, gc: np.ndarray | None = None, scale_sd: float = SCALE_SD,
                   tilt_sd: float = TILT_SD, tau: float = TAU, level_z: float = LEVEL_Z, event_z: float = EVENT_Z,
                   min_height: float = FRACTION_MIN_HEIGHT, min_windows: int = FRACTION_MIN_WINDOWS, grid: int = FRACTION_GRID,
                   min_samples: int = FRACTION_MIN_SAMPLES, spread: float | None = None) -> dict:
    """What the whole numbers leave, for a cohort: each call's `ratio`, `off`, `off_z` and `fractions` are set.

    `calls`: one per row of `values` (genomes x windows, calibrated copies; None where a genome has no call). The windows
    of `leave_out` (the stretches where copies differ between people, which the level leaves out) take no part. `gc`:
    the windows' GC; with it a genome's profile may follow GC, as the profile of a library does whose GC response the
    model left in, and what follows GC is no step. The priors on the lean and on that slope are the cohort's own
    spreads of them (`tilt_sd` where there is no cohort to take them from). With
    fewer than `min_samples` called genomes there is no cohort to judge by: the level is judged against `spread` (the
    spread of the scales in the cohort a saved efficiency table was made on) or, without one, against `scale_sd` (that
    of the cohort the class's rules were made on), and no step is looked for. Returns what was used: the genomes, the
    windows, the spread of the scales and where it came from."""
    starts, V = np.asarray(starts), np.asarray(values, float)
    span = float(unit_length or (starts[-1] + window))
    use = np.array([not any(a <= x < b for a, b in (leave_out or [])) for x in starts]) & (np.isfinite(V) & (V > 0)).any(axis=0)
    if gc is not None:
        use &= np.isfinite(np.asarray(gc, float))
    have = [i for i, c in enumerate(calls) if c is not None]
    out = dict(n=len(have), windows=int(use.sum()), spread=float(spread or scale_sd), spread_from="the efficiency table's cohort" if spread else "the class's rules",
               steps=False, level_z=float(level_z), event_z=float(event_z), min_height=float(min_height))
    if not have or use.sum() < 2 * min_windows:
        return out
    s = starts[use]
    t = (s + window / 2) / span - 0.5
    r = np.ones(len(s)) if rel_noise is None else np.where(np.isfinite(np.asarray(rel_noise, float)[use]), np.asarray(rel_noise, float)[use], 1.0)
    Y = np.zeros((len(have), len(s)))
    W = np.zeros((len(have), len(s)))
    K = np.zeros((len(have), len(s)))
    for g, i in enumerate(have):
        c = calls[i]
        k = np.zeros(len(s))
        for sg in c.segments:
            k[(s >= sg.start) & (s < sg.end)] = sg.state
        v = V[i, use]
        ok = np.isfinite(v) & (v > 0) & (k > 0)
        rel = c.sigma * r / np.where(k > 0, k * c.scale * np.exp(c.tilt * t), 1.0)      # a window's noise, relative: by what is expected there, not by what was read
        y = np.zeros(len(s))
        y[ok] = np.log(v[ok] / k[ok])
        mid = float(np.median(y[ok])) if ok.any() else 0.0
        y = np.where(ok, np.clip(y, mid - Z_CAP * rel, mid + Z_CAP * rel), mid)        # a wild window is drawn in; a window without a value takes no part (weight 0)
        Y[g], W[g], K[g] = y, np.where(ok, 1 / rel ** 2, 0.0), np.where(k > 0, k, np.nan)
    judged = len(have) >= min_samples
    quiet = np.array([not calls[i].uncertain for i in have])
    wc = np.median(W, axis=0)                                                          # one weight per window, for the fits that are compared between genomes

    def mad_sd(x):
        m = float(np.median(x))
        return m, float(1.4826 * np.median(np.abs(x - m)))
    X, prior = t[:, None], [float(tilt_sd)]
    if judged and wc.sum() > 0:
        # the smooth terms: a lean, and with the windows' GC a slope on it. Their priors are the cohort's own spreads of
        # them, from a fit of every genome without a prior: a lean or a slope like the cohort's is cheap, one unlike
        # it is dear, and a step is then the better account
        x = None if gc is None else np.asarray(gc, float)[use] - float(np.sum(wc * np.asarray(gc, float)[use]) / wc.sum())
        if x is not None and float(np.sum(wc * x * x) / wc.sum()) <= 1e-8:                 # windows that do not differ in GC
            x = None
        A = np.c_[np.ones(len(s)), t] if x is None else np.c_[np.ones(len(s)), t, x]
        co = Y @ np.linalg.lstsq(A.T @ (A * wc[:, None]), (A * wc[:, None]).T, rcond=None)[0].T
        lean_sd = max(mad_sd(co[quiet, 1])[1], 0.005)
        prior = [lean_sd]
        out.update(lean_sd=lean_sd)
        if x is not None:
            sd = mad_sd(co[quiet, 2])[1]
            if sd > 0:
                X, prior = np.c_[t, x], [lean_sd, sd]
                out.update(gc_slope_sd=sd)

    def cohort_levels(cuts):
        """The same fit in every genome, with the cohort's weights: levels, genomes x stretches."""
        e = np.r_[0, np.asarray(cuts, int), len(s)]
        A = np.zeros((len(s), len(e) - 1 + X.shape[1]))
        for j in range(len(e) - 1):
            A[e[j]:e[j + 1], j] = 1.0
        A[:, len(e) - 1:] = X
        P = np.r_[np.full(len(e) - 1, 1e-9 * max(float(wc.sum()), 1.0)), 1 / np.asarray(prior) ** 2]
        B = np.linalg.solve(A.T @ (A * wc[:, None]) + np.diag(P), (A * wc[:, None]).T)
        return (Y @ B.T)[:, :len(e) - 1]
    logr = np.zeros(len(have))
    found: dict[int, list] = {}
    for g, i in enumerate(have):
        c = calls[i]
        S = _Sums(Y[g], W[g], X)
        cuts: list[int] = []
        if judged:
            edges = [int(np.searchsorted(s, sg.start)) for sg in c.segments[1:]]
            cuts = propose_steps(Y[g], W[g], X, float(c.copies or np.nanmedian(K[g])), prior, also=edges, tau=tau, min_height=min_height,
                                 min_windows=min_windows, grid=grid)
        events: list[Fraction] = []
        while cuts:
            lv, _, _ = fit_levels(S, cuts, prior)
            L = cohort_levels(cuts)
            ref = int(np.argmin(np.abs(lv)))                                           # the stretch whose level puts the scale nearest one
            e = [0] + cuts + [len(s)]
            events, keep = [], set()
            for j in range(len(e) - 1):
                if j == ref:
                    continue
                kk = int(np.nanmedian(K[g, e[j]:e[j + 1]]))
                offset = float(kk * (np.exp(lv[j] - lv[ref]) - 1))
                contrast = L[:, j] - L[:, ref]
                m, sd = mad_sd(contrast[quiet & (np.arange(len(have)) != g)])
                zz = float((contrast[g] - m) / sd) if sd > 0 else 0.0
                if abs(zz) >= event_z and abs(offset) >= min_height:
                    a = 0 if e[j] == 0 else int((s[e[j] - 1] + window + s[e[j]]) // 2)
                    b = int(span) if e[j + 1] == len(s) else int((s[e[j + 1] - 1] + window + s[e[j + 1]]) // 2)
                    events.append(Fraction(a, b, kk, round(offset, 3), round(zz, 2), int(e[j + 1] - e[j])))
                    keep |= {x for x in (e[j], e[j + 1]) if 0 < x < len(s)}
            if keep == set(cuts):
                break
            cuts = sorted(keep)                                                        # the steps that did not stand are taken out, and the rest fitted again
            if not cuts:
                events = []
        lv, _, _ = fit_levels(S, cuts, prior)
        logr[g] = float(lv[int(np.argmin(np.abs(lv)))])
        found[i] = events
    m, sd = mad_sd(logr[quiet]) if judged and quiet.sum() >= min_samples else (0.0, 0.0)
    if judged and quiet.sum() >= min_samples and sd > 0:
        out.update(spread=sd, spread_from="this cohort", steps=True, median=m)
    else:
        m, sd = 0.0, float(spread or scale_sd)
    for g, i in enumerate(have):
        c = calls[i]
        c.ratio = float(np.exp(logr[g]))
        c.off = float((c.copies or np.nanmedian(K[g])) * (c.ratio / np.exp(m) - 1))
        c.off_z = float((logr[g] - m) / sd)
        c.fractions, c.level_z = found[i], float(level_z)
    return out


def fractions_string(call: Call) -> str:
    return ";".join(f.label(call.copies) for f in call.fractions)
