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


@dataclass
class Segment:
    start: int
    end: int
    windows: int
    state: int
    mean: float                # the windows' weighted mean, in copies on the genome's scale (divided by f)
    se: float

    @property
    def off_integer(self) -> bool:
        return abs(self.mean - self.state) > max(OFF_INTEGER, 3 * self.se)

    def as_dict(self) -> dict:
        return dict(start=self.start, end=self.end, windows=self.windows, state=self.state, mean=round(self.mean, 3), se=round(self.se, 3))


@dataclass
class Event:
    delta: int                 # copies relative to the complete copies
    start: int
    end: int
    kind: str                  # "partial copy" (a gain that reaches an end of the unit), "gain", "loss"

    def __str__(self) -> str:
        return f"{self.delta:+d}:{self.start // 1000}-{self.end // 1000}kb"


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

    @property
    def uncertain(self) -> bool:
        """Another whole number of copies explains the profile nearly as well (a level between two whole numbers)."""
        return self.gap < MIN_GAP

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
            segs.append(Segment(int(s[a]), int(s[b]) + window, b - a + 1, int(K[k]), float(np.sum(w * lv[a:b + 1]) / w.sum()) / Fc[c], float(np.sqrt(1 / w.sum())) / Fc[c]))
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
