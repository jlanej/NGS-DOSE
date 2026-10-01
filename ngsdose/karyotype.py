"""Chromosomes in copies: whole chromosomes, arms and segments, from the single-copy regions.

Every counts file holds single-copy regions: the controls every class is a ratio against, the held-out
known-truth sets on the autosomes, chrX and chrY, and the karyotype windows spread along every arm. Each
region's log(observed / expected) under the genome's own fragment-GC curve says how many copies of that
sequence the genome holds, relative to its autosomes. A chromosome is its regions read together, so a gain
or a loss of a chromosome, of an arm or of a stretch is a level shared by the regions that lie on it.

Three things stand between the regions and a number of copies, and a `Model` learned on a reference cohort
carries all three to any other genome, one genome at a time:

- a region's *efficiency*: what it reads in a genome that holds two copies (its median over the cohort);
- the libraries' *shared modes*: a few components along which the regions of one library depart together
  (late-replicating, AT-rich sequence reads low in DNA from cycling cells; GC residue). They are learned,
  and a genome is scored on them, from what the regions of a chromosome do relative to each other: every
  region less its chromosome's mean. Components of the whole matrix would take a common trisomy out with
  them; contrasts within chromosomes cannot, because a chromosome gained or lost moves all its regions
  alike and leaves no contrast. What a genome's scores predict for each region, the chromosome's share
  of it included, is then taken out;
- a region's *spread* in the cohort, which is its weight: a region that a common copy-number variant
  moves is left with little say, without being named in advance.

One thing the contrasts cannot see, because it has no pattern within a chromosome: in some libraries the
chromosomes rich in GC (19, 22, 17, 16) read a little low or high together, in proportion to the
chromosome's GC content and whatever the GC of the region itself. Each chromosome is therefore set against
what the OTHER autosomes of the same genome say: their levels are fitted, robustly, by a common shift and
a slope on the chromosomes' GC, and the fit's value for this chromosome is taken out of it. A chromosome
gained or lost takes no part in its own correction, and one that departs from the others' line (a gain, a
loss) is given no say in theirs.

A level's error is its regions' noise in this genome (their spread in the cohort, times the genome's own
noise factor, with the scores' share), and the cohort then says how far that is from the whole of it: the
levels of a cohort's chromosomes spread a little more than their regions' noise predicts, by a factor (the
regions' tails) and by a floor in copies that more regions do not lower (what a chromosome's regions share
and no pattern within the chromosome accounts for). Both are measured on the cohort (`calibrate`), across
chromosomes with few regions and with many, and z is taken against the wider error.

Along each chromosome the regions are then read by a chain: its states are levels in copies on a grid, a
region's distance from a level is judged against its spread (capped, so one region that a variant moves
cannot buy a change of level), and every change of level costs `tau`. One level throughout is a whole
chromosome; a change at the centromere is an arm; anything else a segment. Each level is reported as it is
measured, with its distance from the nearest whole number of copies and from the cohort's spread: a level
between whole numbers is a change in part of the cells (or two people's DNA in one sample), and it is
kept in view rather than rounded away.

The sex chromosomes are read on scales of their own. A region of chrX is calibrated on the genomes that
hold one X; a second X reads `u` of the first (0.95 in lymphoblastoid lines, where the inactive X
replicates late; the cohort measures it, per region), and it follows the shared modes more than the first
does. chrX is reported in copies on the scale that puts the cohort's one-X and two-X genomes at one and
two; chrY on the scale of the genomes that hold one.

Regions flagged in a genome (a chromosome or a segment off its whole number) take no part in that genome's
level or scores: the read is repeated until the set is stable.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field, replace

import numpy as np

AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
K_MAX = 20                 # components at most
MIN_FIT = 50               # genomes a cohort needs to learn a model of its own
MIN_REPIN = 10             # two-X genomes a small cohort needs to measure a second X's factor itself
MIN_PER_COMPONENT = 20     # genomes per component
TAU = 14.0                 # cost of a change of level along a chromosome, in units of log-likelihood
Z_CAP = 4.0                # the chain: a place further than this many SDs from a level counts as this far
MIN_LOCI = 4               # the shortest stretch kept, in places along the chromosome
LOCUS_GAP = 50_000         # regions with less than this between them are one place: a window's pieces, which one variant can move together
LOCUS_SPAN = 150_000       # and a place is no longer than this
GRID_STEP = 0.05           # copies between the chain's levels
MAX_COPIES = 5.0
C_FLOOR = 0.1              # a level's noise does not fall below that of this share of the reference copies
LEVEL_Z = 5.0              # a level this many SDs from its whole number is fractional ...
MIN_OFF = 0.03             # ... when it is also this many copies from it
MAX_SE = 0.10              # a level known no better than this (copies) is uncertain
MIN_CHROM_REGIONS = 3      # fewer regions than this on a chromosome: not read
ARM_MIN_LOCI = 5           # an arm with fewer places than this is not read on its own
MAX_REGION_SD = 3.0        # a region whose spread is this many times the typical one is not used
MAX_ITER = 5
CLIP = 4.0                 # learning: a region's values beyond this many of its SDs are pulled in
MAX_CHROM_SHARE = 0.5      # a component with more than this share of its weight on one chromosome is no library mode
MIN_CALIBRATE = 30         # levels a chromosome's (or a sex-chromosome class's) calibration needs
MIN_TILT = 8               # autosomes a genome's chromosome-level fit needs
FAR_Z = 5.0                # a place (or a piece within its place) this many SDs from its level is not noise: it has no say in the level
MIN_FAR, FAR_SHARE = 2, 0.15   # this many such places, and this share of a chromosome's, and its level is not settled
TILT_HUBER, TILT_BIWEIGHT = 1.5, 4.685   # the fit's robust weights: Huber's to start, Tukey's biweight to finish (a level this many SEs off has no say)
MAX_LOCAL = 0.05           # the largest SD (log) accepted for what the pieces of one place share
NAME = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
ONE = (0.75, 1.25)         # a raw reading in this range is one copy of a sex chromosome (for learning the scales)
TWO = (1.70, 2.20)
NONE_Y = 0.10


def _mad(x, axis=None):
    med = np.nanmedian(x, axis=axis, keepdims=True)
    return 1.4826 * np.nanmedian(np.abs(x - med), axis=axis)


def _natural(s: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


@dataclass
class Table:
    """The regions of one layout, in the order given: where each lies."""
    names: list[str]
    chrom: np.ndarray
    start: np.ndarray          # -1 for a set read as one value (an estimate written before regions were kept singly)
    end: np.ndarray
    arm: np.ndarray            # "p", "q", or "" where not known
    kind: np.ndarray           # "A" autosomal, "X", "Y", "" for a region on no chromosome read here

    @property
    def mid(self) -> np.ndarray:
        return (self.start + self.end) / 2.0

    def chromosomes(self) -> list[str]:
        return sorted({c for c, k in zip(self.chrom, self.kind) if k}, key=_natural)


def aggregate_name(chrom: str) -> str:
    """The name of a chromosome's known-truth set read as one value."""
    return f"{chrom}:set"


def table(names: list[str], arms: dict[str, int] | None = None) -> Table:
    """`arms`: per chromosome, the position where the short arm ends and the long arm begins."""
    chrom, start, end, arm, kind = [], [], [], [], []
    for n in names:
        m = NAME.match(n)
        c = m["chrom"] if m else n.rsplit(":", 1)[0]
        s, e = (int(m["start"]), int(m["end"])) if m else (-1, -1)
        chrom.append(c)
        start.append(s)
        end.append(e)
        b = (arms or {}).get(c)
        arm.append("" if b is None or s < 0 else "p" if (s + e) / 2 < b else "q")
        kind.append("A" if c in AUTOSOMES else "X" if c == "chrX" else "Y" if c == "chrY" else "")
    return Table(list(names), np.array(chrom), np.array(start), np.array(end), np.array(arm), np.array(kind))


def gather(result: dict, control_names: list[str] | None = None) -> tuple[list[str], np.ndarray] | None:
    """One genome's single-copy regions: names and log(observed / expected) under its GC curve.

    An estimate written by this version keeps them singly (`single_copy`). One written before keeps the control
    regions (`control_qc.region_log_ratio`, in the bundle's order: `control_names`) and, of chrX and chrY, only
    each known-truth set's pooled value: that becomes one region per chromosome (`aggregate_name`), which gives the
    chromosome's level and no more (no arms, no segments)."""
    sc = result.get("single_copy")
    if sc:
        return list(sc["names"]), np.array([np.nan if v is None else v for v in sc["log_ratio"]], float)
    qc = result.get("control_qc") or {}
    lr = qc.get("region_log_ratio")
    if lr is None or control_names is None or len(control_names) != len(lr):
        return None
    names, vals = list(control_names), [float(v) for v in lr]
    for c in ("chrX", "chrY"):
        t = (result.get("truth_regions") or {}).get(c)
        if t and t.get("cn") is not None:
            names.append(aggregate_name(c))
            vals.append(float(np.log(max(t["cn"], 1e-4) / 2.0)))
    return names, np.array(vals, float)


@dataclass
class Model:
    """What a reference cohort teaches about the regions (module docstring). Arrays are per region, in `names` order."""
    names: list[str]
    a: np.ndarray              # log efficiency: per two copies on an autosome, per one copy on chrX and chrY
    sd: np.ndarray             # spread (log) at those copies in a typical genome of the cohort; NaN where not learned
    V: np.ndarray              # regions x components: loadings on the shared modes
    lam2: np.ndarray           # prior variance of each component's score (its variance in the cohort)
    u: np.ndarray              # chrX regions: what a second X reads of the first; NaN elsewhere
    x2: np.ndarray             # components: how much more the second X follows each mode than the first (log)
    phi: dict[str, float]      # per chromosome (the sex chromosomes by whole number, "chrX:2"): the factor on a level's error
    arms: dict[str, int]
    n: int = 0                 # genomes it was learned on
    info: dict = field(default_factory=dict)
    floor: dict[str, float] = field(default_factory=dict)   # keyed like phi: the error's floor, in copies at the reference level
    local: float = 0.0         # SD (log) of what the pieces of one place share beyond the modes: a place's error does not fall below it
    gc: dict[str, float] = field(default_factory=dict)      # per chromosome: its GC content (empty: no chromosome-level fit)

    @property
    def k(self) -> int:
        return self.V.shape[1]

    def to_json(self) -> dict:
        r = lambda v, d=5: [None if not np.isfinite(x) else round(float(x), d) for x in v]
        return dict(format="ngsdose-karyotype-model-1", names=self.names, a=r(self.a), sd=r(self.sd), V=[r(row, 4) for row in self.V.T], lam2=r(self.lam2, 6),
                    u=r(self.u, 4), x2=r(self.x2, 5), phi={c: round(float(v), 3) for c, v in self.phi.items()},
                    floor={c: round(float(v), 5) for c, v in self.floor.items()}, local=round(float(self.local), 5), gc=self.gc, arms=self.arms, n=self.n,
                    info=self.info)

    @staticmethod
    def from_json(d: dict) -> "Model":
        if d.get("format") != "ngsdose-karyotype-model-1":
            raise ValueError(f"not a karyotype model (format {d.get('format')!r})")
        f = lambda v: np.array([np.nan if x is None else x for x in v], float)
        m = len(d["names"])
        V = np.array([f(row) for row in d["V"]], float).T if d["V"] else np.zeros((m, 0))
        return Model(list(d["names"]), f(d["a"]), f(d["sd"]), np.nan_to_num(V), f(d["lam2"]), f(d["u"]), f(d["x2"]) if d.get("x2") else np.zeros(V.shape[1]),
                     dict(d.get("phi") or {}), {c: int(v) for c, v in (d.get("arms") or {}).items()}, int(d.get("n", 0)), dict(d.get("info") or {}),
                     dict(d.get("floor") or {}), float(d.get("local") or 0.0), {c: float(v) for c, v in (d.get("gc") or {}).items()})


@dataclass
class Seg:
    """A stretch of a chromosome at one level."""
    chrom: str
    first: int                 # its regions, as positions among the chromosome's regions in order
    last: int
    start: int                 # bp; a breakpoint lies midway between the regions either side
    end: int
    regions: int
    copies: float
    se: float
    span: str = "whole"        # "whole", "p", "q", "pter", "qter", "inner"
    loci: int = 0              # the places its regions lie at (a window is one)

    def as_dict(self) -> dict:
        return dict(chrom=self.chrom, start=self.start, end=self.end, regions=self.regions, loci=self.loci, copies=round(self.copies, 3), se=round(self.se, 4), span=self.span)


@dataclass
class Chrom:
    chrom: str
    regions: int
    copies: float              # the level of the stretch that holds most of the chromosome
    se: float                  # its error: the regions' noise in this genome, widened as the model's cohort found it (phi, floor)
    mean: float                # over all its regions, whatever the stretches
    segments: list[Seg] = field(default_factory=list)
    whole: int | None = None   # the nearest whole number of copies
    off: float | None = None   # copies less that
    z: float | None = None     # off in SDs (the genome's own, widened by the cohort's where that is wider)
    status: str = "not read"   # "settled", "fractional", "uncertain", "not read"
    arms: dict = field(default_factory=dict)   # "p", "q" -> (copies, SE, regions): each arm read on its own, where it has ARM_MIN_LOCI places
    se0: float | None = None   # the error before the cohort's widening: what `calibrate` sets the cohort's spread against
    places: int = 0            # the places its regions lie at
    far: int = 0               # how many of them lie FAR_Z SDs from their stretch's level
    tilt: float = 0.0          # what the fit across chromosomes took out of it (log): its level without the fit is copies x exp(tilt)


@dataclass
class Event:
    """A stretch that is not at the whole number its chromosome is expected at (two; the sex chromosomes' own)."""
    chrom: str
    span: str
    start: int
    end: int
    regions: int
    copies: float
    se: float
    delta: float               # copies less the expected whole number
    z: float
    whole: bool                # on a whole number (a gain or loss in every cell), else in part of the cells

    def label(self) -> str:
        c = self.chrom[3:]
        d = 0 if self.end - self.start >= 5e6 else 1
        where = c if self.span == "whole" else f"{c}{self.span}" if self.span in ("p", "q") else f"{c}({self.start / 1e6:.{d}f}-{self.end / 1e6:.{d}f}Mb)"
        sign = "+" if self.delta > 0 else "-"
        n = int(round(abs(self.delta)))
        if self.whole:
            return ",".join([f"{sign}{where}"] * max(n, 1))
        return f"{sign}{where}[{abs(self.delta):.2f}]"


@dataclass
class Reading:
    """One genome read against a model."""
    chromosomes: dict[str, Chrom]
    events: list[Event]
    level: float               # the genome's autosomal level less the model's (log): what was taken out
    noise: float               # its regions' spread in units of the cohort's (1: a typical genome)
    regions: int
    scores: np.ndarray
    x: float | None = None     # chrX copies on the pinned scale
    y: float | None = None
    complement: str = ""       # the whole numbers of X and Y, as letters
    note: str = ""
    off: np.ndarray | None = None       # per model region: in an event (left out of level and scores)
    iterations: int = 0
    tilt: tuple | None = None           # the chromosome-level fit over the autosomes: (shift, slope per unit of GC, the slope's SE), log

    @property
    def gc_tilt(self) -> float | None:
        """By how many copies a chromosome ten points richer in GC reads higher in this genome (at two copies)."""
        return None if self.tilt is None else float(2.0 * self.tilt[1] * 0.1)

    @property
    def gc_tilt_z(self) -> float | None:
        """The slope in its own standard errors: how far this library's GC-rich chromosomes depart from what the fit can measure."""
        return None if self.tilt is None or not self.tilt[2] > 0 else float(self.tilt[1] / self.tilt[2])

    @property
    def status(self) -> str:
        st = [c.status for c in self.chromosomes.values()]
        if "fractional" in st or any(not e.whole for e in self.events):
            return "fractional"
        return "uncertain" if "uncertain" in st else "settled" if "settled" in st else "not read"

    def karyotype(self) -> str:
        """Written like a karyotype, measured by depth: counts of chromosomes, not their structure. A change in part of
        the cells carries the share in brackets, a stretch its span in megabases."""
        if not self.complement:
            return ""
        whole = [e for e in self.events if e.whole and e.span == "whole" and e.chrom in AUTOSOMES]
        n = 44 + sum(int(round(e.delta)) for e in whole) + len(self.complement)
        parts = [f"{n}", self.complement]
        rest = sorted((e for e in self.events if not (e.chrom in ("chrX", "chrY") and e.whole and e.span == "whole")),
                      key=lambda e: (e.chrom not in ("chrX", "chrY"), _natural(e.chrom), e.start))
        parts += [e.label() for e in rest]
        return ",".join(parts)


def _round(v: float) -> float:
    """To the nearest whole number, halves up."""
    return float(np.floor(v + 0.5))


def _wmedian(x: np.ndarray, w: np.ndarray) -> float:
    o = np.argsort(x)
    cw = np.cumsum(w[o])
    return float(x[o][min(int(np.searchsorted(cw, 0.5 * cw[-1])), len(x) - 1)])


def _near_mean(x: np.ndarray, sd: np.ndarray, at: float, far_z: float = FAR_Z) -> tuple[float, np.ndarray]:
    """The weighted mean of the values within `far_z` SDs of it, started at `at` and repeated until the set is stable; and
    which values those are. A value far from the rest (a copy-number variant under one region) is left out, not pulled
    in: pulled in to four SDs it would still move a chromosome of a few regions by more than its error."""
    w = 1.0 / sd ** 2
    c, keep = float(at), np.abs(x - at) <= far_z * sd
    for _ in range(6):
        if not keep.any():
            # nothing near the start: the nearest half decides
            keep = np.abs(x - c) / sd <= np.median(np.abs(x - c) / sd)
        c = float(np.sum(w[keep] * x[keep]) / w[keep].sum())
        new = np.abs(x - c) <= far_z * sd
        if (new == keep).all():
            break
        keep = new
    if not keep.any():
        keep = np.ones(len(x), bool)
        c = float(np.sum(w * x) / w.sum())
    return c, keep


def _ridge(V: np.ndarray, lam2: np.ndarray, e: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """A genome's scores on the shared modes from regions with loadings V, residuals e and weights w; and (V'WV + 1/lam2)^-1."""
    if V.shape[1] == 0:
        return np.zeros(0), np.zeros((0, 0))
    Vw = V * w[:, None]
    Ainv = np.linalg.inv(V.T @ Vw + np.diag(1.0 / lam2))
    return Ainv @ (Vw.T @ e), Ainv


def _contrast(chrom: np.ndarray, base: np.ndarray, w: np.ndarray, e: np.ndarray, V: np.ndarray, min_regions: int = MIN_CHROM_REGIONS):
    """For scoring on contrasts within chromosomes: which regions take part (those of `base` on a chromosome with at
    least `min_regions` of them), and each region's chromosome mean of `e` and of the loadings `V` over those (weighted)."""
    sel = np.zeros(len(e), bool)
    cm, Vm = np.zeros(len(e)), np.zeros(V.shape)
    for c in np.unique(chrom[base]):
        on = base & (chrom == c)
        if on.sum() < min_regions:
            continue
        sel |= on
        ws = w[on].sum()
        cm[on] = float((w[on] * e[on]).sum() / ws)
        Vm[on] = (w[on, None] * V[on]).sum(0) / ws
    return sel, cm, Vm


def _tilt(mu: np.ndarray, se: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Chromosome levels `mu` (log, SEs `se`) fitted by a common shift and a slope on `x` (the chromosomes' GC, centred),
    robustly: Huber's weights first, then Tukey's biweight, under which a chromosome far off the others' line (a gain,
    a loss) has no say. Returns (coefficients, the normal matrix, its right-hand side, each chromosome's weight)
    for leaving one chromosome out in closed form."""
    X = np.c_[np.ones(len(mu)), x]
    b = np.array([float(np.median(mu)), 0.0])
    eps = 1e-12 * np.eye(2)
    W = 1.0 / se ** 2
    for i in range(40):
        z = np.abs(mu - X @ b) / se
        if i < 10:
            rw = np.where(z <= TILT_HUBER, 1.0, TILT_HUBER / np.maximum(z, 1e-12))
        else:
            rw = np.where(z < TILT_BIWEIGHT, (1.0 - (z / TILT_BIWEIGHT) ** 2) ** 2, 0.0)
            if (rw > 0).sum() < max(MIN_TILT, len(mu) // 2):      # the line has lost the bulk: keep Huber's
                rw = np.where(z <= TILT_HUBER, 1.0, TILT_HUBER / np.maximum(z, 1e-12))
        W = rw / se ** 2
        A = (X * W[:, None]).T @ X + eps
        new = np.linalg.solve(A, (X * W[:, None]).T @ mu)
        done = i >= 12 and np.abs(new - b).max() < 1e-7
        b = new
        if done:
            break
    A = (X * W[:, None]).T @ X + eps
    return b, A, (X * W[:, None]).T @ mu, W


def loci(start: np.ndarray, end: np.ndarray, x: np.ndarray, sd: np.ndarray, gap: int = LOCUS_GAP, span: int = LOCUS_SPAN, far_z: float = FAR_Z,
         local: float = 0.0):
    """Regions in order -> places: runs of regions with less than `gap` between one's end and the next one's start, no
    run longer than `span`, each read as one value. A window's pieces lie within 100 kb, and one copy-number variant
    can move several of them; read as one place, capped like any other, it cannot pass for a stretch of the chromosome.
    Within a place a piece further than `far_z` SDs from the others is left out of the place's value.

    `local`: the SD (in x's units) of what the pieces of one place share. Each piece's own SD holds it already, so the
    mean of n pieces is not n times better known: its variance is 1/sum(w) plus the shared part's, less the share of it
    that 1/sum(w) counts (one n-th, for equal pieces).
    Returns (first region of each place, value, SD, regions per place)."""
    n = len(start)
    if n and start[0] >= 0:
        first = [0]
        for i in range(1, n):
            if start[i] - end[i - 1] > gap or end[i] - start[first[-1]] > span:
                first.append(i)
        first = np.array(first)
    else:
        first = np.arange(n)
    last = np.r_[first[1:], n]
    xv, sv = np.empty(len(first)), np.empty(len(first))
    for j, (a, b) in enumerate(zip(first, last)):
        if b - a == 1:
            xv[j], sv[j] = x[a], sd[a]
            continue
        w = 1.0 / sd[a:b] ** 2
        xv[j], keep = _near_mean(x[a:b], sd[a:b], _wmedian(x[a:b], w), far_z)
        w = w[keep]
        sv[j] = float(np.sqrt(1.0 / w.sum() + local ** 2 * (1.0 - np.sum(w * w) / w.sum() ** 2)))
    return first, xv, sv, last - first


def chain(x: np.ndarray, sd_ref: np.ndarray, ref: float, tau: float = TAU, step: float = GRID_STEP, z_cap: float = Z_CAP,
          min_loci: int = MIN_LOCI, max_copies: float = MAX_COPIES) -> np.ndarray:
    """The most probable piecewise-constant level along a chromosome. `x`: the places in order, in copies; `sd_ref`: each
    one's SD in copies when the level is `ref` (noise grows as the square root of the level, with a floor). Returns
    each place's level on the grid, runs shorter than `min_loci` given to the neighbour that fits them better."""
    grid = np.arange(0.0, max_copies + step / 2, step)
    n, K = len(x), len(grid)
    rel = np.sqrt(np.maximum(grid, C_FLOOR * ref) / ref)                       # sigma(c) / sigma(ref)
    z = (x[:, None] - grid[None, :]) / (sd_ref[:, None] * rel[None, :])
    ll = -0.5 * np.minimum(z * z, z_cap * z_cap) - np.log(rel)[None, :]
    back = np.zeros((n, K), np.int32)
    ar = np.arange(K)
    score = ll[0].copy()
    for q in range(1, n):
        best = int(np.argmax(score))
        move = score[best] - tau
        take = move > score
        back[q] = np.where(take, best, ar)
        score = np.where(take, move, score) + ll[q]
    path = np.zeros(n, int)
    path[-1] = int(np.argmax(score))
    for q in range(n - 1, 0, -1):
        path[q - 1] = back[q, path[q]]
    while True:                                                                # short runs go to a neighbour
        cut = np.flatnonzero(np.diff(path)) + 1
        lo, hi = np.r_[0, cut], np.r_[cut, n]
        if len(lo) == 1:
            break
        ln = hi - lo
        j = int(np.argmin(ln))
        if ln[j] >= min_loci:
            break
        cand = [path[lo[j + d]] for d in (-1, 1) if 0 <= j + d < len(lo)]
        cost = [float(-ll[lo[j]:hi[j], c].sum()) for c in cand]
        path[lo[j]:hi[j]] = cand[int(np.argmin(cost))]
    return grid[path]


def _level(x: np.ndarray, sd_ref: np.ndarray, ref: float, at: float) -> tuple[float, float, int]:
    """The level of a stretch's places, its SE (before the genome's noise factor's and the scores' shares), and how many
    places lie too far from it to have a say (`_near_mean`)."""
    sd = sd_ref * np.sqrt(max(at, C_FLOOR * ref) / ref)
    c, keep = _near_mean(x, sd, at)
    return c, float(np.sqrt(1.0 / np.sum(1.0 / sd[keep] ** 2))), int((~keep).sum())


def _merge(bounds: list, xl: np.ndarray, sl: np.ndarray, ref: float, path: np.ndarray, tau: float) -> list:
    """The chain's runs of places, re-judged on their measured levels: the chain's states lie on a grid, and a level between
    two grid steps can be cut in two, each half on its nearer step. A change is kept only where the levels either side
    differ by what the chain's own price asks, (a - b)^2 / (se_a^2 + se_b^2) > 2 tau; the closest pair is joined first, and
    the levels are measured again. Returns [(first place, end place, level, SE, places left out)]."""
    segs = [(a, b, *_level(xl[a:b], sl[a:b], ref, float(path[a]))) for a, b in bounds]
    while len(segs) > 1:
        chi2 = [(u[2] - v[2]) ** 2 / (u[3] ** 2 + v[3] ** 2) for u, v in zip(segs[:-1], segs[1:])]
        i = int(np.argmin(chi2))
        if chi2[i] > 2.0 * tau:
            break
        a, b = segs[i][0], segs[i + 1][1]
        at = float(np.round(_wmedian(xl[a:b], 1.0 / sl[a:b] ** 2) / GRID_STEP) * GRID_STEP)
        segs[i:i + 2] = [(a, b, *_level(xl[a:b], sl[a:b], ref, at))]
    return segs


def _span(t_arm: np.ndarray, first: int, last: int) -> str:
    n = len(t_arm)
    if first == 0 and last == n - 1:
        return "whole"
    inside = np.zeros(n, bool)
    inside[first:last + 1] = True
    for a in ("p", "q"):
        on = t_arm == a
        na, other = int(on.sum()), int((~on).sum())
        if na and other and (inside & on).sum() >= max(na - 1, 0.9 * na) and (inside & ~on).sum() <= max(1, 0.1 * other):
            return a
    if first == 0 and t_arm[0] == "p":
        return "pter"
    if last == n - 1 and t_arm[-1] == "q":
        return "qter"
    return "inner"


def read(model: Model, names: list[str], lr: np.ndarray, k: int | None = None, tau: float = TAU, level_z: float = LEVEL_Z, min_off: float = MIN_OFF,
         tab: Table | None = None, index: np.ndarray | None = None) -> Reading | None:
    """One genome against the model. `names`, `lr`: its regions and their log(observed / expected) (`gather`). Regions the
    model does not know are left out. None when fewer than 20 autosomal regions can be read.

    `tab`, `index`: the model's table and each model region's position in `names` (-1: absent), for a caller that reads
    many genomes with the same regions."""
    tab = tab or table(model.names, model.arms)
    if index is None:
        pos = {n: i for i, n in enumerate(names)}
        index = np.array([pos.get(n, -1) for n in model.names])
    y = np.where(index >= 0, np.asarray(lr, float)[np.maximum(index, 0)], np.nan)
    k = model.k if k is None else min(k, model.k)
    V, lam2 = model.V[:, :k], model.lam2[:k]
    ok = np.isfinite(y) & np.isfinite(model.a) & np.isfinite(model.sd) & (model.sd > 0)
    auto = ok & (tab.kind == "A")
    if auto.sum() < 20:
        return None
    w = np.where(ok, 1.0 / np.where(ok, model.sd, 1.0) ** 2, 0.0)
    e = np.where(ok, y - model.a, 0.0)
    off = np.zeros(len(y), bool)
    # each chromosome's usable regions, in order along it (a chromosome without one is not read)
    order = {c: np.flatnonzero(ok & (tab.chrom == c))[np.argsort(tab.mid[ok & (tab.chrom == c)], kind="stable")] for c in tab.chromosomes()}
    order = {c: idx for c, idx in order.items() if len(idx)}
    out = None
    for it in range(MAX_ITER):
        base = auto & ~off
        if base.sum() < 20:
            base = auto
        level = _wmedian(e[base], w[base])
        # a region far from its chromosome's others in this genome (a copy-number variant) is pulled in before the scores are
        # taken: CLIP of its SDs, as widened by this genome's modes and noise
        d = e - level
        g0 = max(float(_mad((d / np.where(ok, model.sd, 1.0))[base])), 0.5)
        for idx in order.values():
            j = idx[base[idx]]
            if len(j) >= MIN_CHROM_REGIONS:
                med, lim = _wmedian(d[j], w[j]), CLIP * g0 * model.sd[j]
                d[j] = np.clip(d[j], med - lim, med + lim)
        # the shared modes are scored on what the regions of a chromosome do relative to each other: a whole chromosome
        # gained or lost moves all its regions alike and so leaves the scores where they are
        sel, cm, Vm = _contrast(tab.chrom, base, w, d, V)
        s, Ainv = _ridge(V[sel] - Vm[sel], lam2, (d - cm)[sel], w[sel])
        r = e - np.where(tab.kind == "A", level, level - np.log(2.0)) - V @ s
        infl = {}
        for c, idx in order.items():
            vbar = (w[idx, None] * V[idx]).sum(0) / w[idx].sum()
            infl[c] = (float(vbar @ Ainv @ vbar) if k else 0.0, s)
        s_all = s
        g = float(_mad((r / np.where(ok, model.sd, 1.0))[base]))
        g = max(g, 0.5)
        # across chromosomes: each one against what the other autosomes say (a common shift, a slope on the chromosomes' GC)
        tilt, tilt_var, coef, slope_se = {}, {}, None, 0.0
        if model.gc:
            fc = [c for c, idx in order.items() if tab.kind[idx[0]] == "A" and c in model.gc and (~off[idx]).sum() >= MIN_CHROM_REGIONS]
            if len(fc) >= MIN_TILT:
                centre = float(np.mean([model.gc[c] for c in AUTOSOMES if c in model.gc]))
                mu, se_c = np.empty(len(fc)), np.empty(len(fc))
                for q, c in enumerate(fc):
                    j = order[c][~off[order[c]]]
                    mu[q], near = _near_mean(r[j], model.sd[j] * g, _wmedian(r[j], w[j]))
                    # its error as the cohort found it, where the model knows (log: a floor in copies at two is half that)
                    se_c[q] = float(np.hypot(g / np.sqrt(w[j][near].sum()) * max(model.phi.get(c, 1.0), 1.0), model.floor.get(c, 0.0) / 2.0))
                xs = np.array([model.gc[c] - centre for c in fc])
                coef, A, rhs, W = _tilt(mu, se_c, xs)
                slope_se = float(np.sqrt(max(np.linalg.inv(A)[1, 1], 0.0)))
                at = {c: q for q, c in enumerate(fc)}
                for c in order:
                    if c not in model.gc:
                        continue
                    xc = np.array([1.0, model.gc[c] - centre])
                    q = at.get(c)
                    Ac = A if q is None else A - W[q] * np.outer(xc, xc)
                    bc = rhs if q is None else rhs - W[q] * xc * mu[q]
                    try:
                        Ai = np.linalg.inv(Ac)
                    except np.linalg.LinAlgError:
                        continue
                    tilt[c], tilt_var[c] = float(xc @ Ai @ bc), float(xc @ Ai @ xc)
        chroms, events, new_off, notes = {}, [], np.zeros(len(y), bool), []
        xy = {}
        for c, idx in order.items():
            kind = tab.kind[idx[0]]
            var_s, s = infl[c]
            rc, var_t = r[idx] - tilt.get(c, 0.0), tilt_var.get(c, 0.0)
            if kind == "A":
                ref, x, sd_ref, shared = 2.0, 2.0 * np.exp(rc), 2.0 * model.sd[idx] * g, 2.0 * model.local
            elif kind == "X":
                # a second X reads u of the first, and follows the modes more: copies = 1 + (reading - 1) / u
                u = np.where(np.isfinite(model.u[idx]), model.u[idx], 1.0) * np.exp(float(model.x2[:k] @ s) if k else 0.0)
                ref, x, sd_ref, shared = 1.0, 1.0 + (np.exp(rc) - 1.0) / u, model.sd[idx] * g / u, model.local / float(np.median(u))
            else:
                ref, x, sd_ref, shared = 1.0, np.exp(rc), model.sd[idx] * g, model.local
            first, xl, sl, nreg = loci(tab.start[idx], tab.end[idx], x, sd_ref, local=shared)
            if len(first) < MIN_CHROM_REGIONS and not (len(idx) >= 1 and tab.start[idx[0]] < 0):
                chroms[c] = Chrom(c, len(idx), float("nan"), float("nan"), float("nan"))
                continue
            if len(first) >= 2 * MIN_LOCI:
                path = chain(xl, sl, ref, tau=tau)
            else:
                path = np.full(len(first), float(np.round(_wmedian(xl, 1 / sl ** 2) / GRID_STEP) * GRID_STEP))
            cut = np.flatnonzero(np.diff(path)) + 1
            last = np.r_[first[1:], len(idx)] - 1                                      # each place's last region
            segs, far = [], 0
            for a_, b_, lv, se, out_ in _merge(list(zip(np.r_[0, cut], np.r_[cut, len(first)])), xl, sl, ref, path, tau):
                far += out_
                # the scores' share of the error (the regions' mean loading, through the scores' covariance), and the share of
                # the fit across chromosomes
                se = float(np.sqrt(se ** 2 + max(lv, C_FLOOR * ref) ** 2 * (g * g * var_s + var_t)))
                segs.append(Seg(c, int(first[a_]), int(last[b_ - 1]), int(tab.start[idx[first[a_]]]), int(tab.end[idx[last[b_ - 1]]]),
                                int(nreg[a_:b_].sum()), lv, se, loci=int(b_ - a_)))
            for s0, s1 in zip(segs[:-1], segs[1:]):
                s0.end = s1.start = (s0.end + s1.start) // 2
            arm = tab.arm[idx]
            for sg in segs:
                sg.span = _span(arm, sg.first, sg.last)
            main = max(segs, key=lambda sg: (sg.loci, sg.regions))
            # over the chromosome's length, whatever its stretches: what the sex chromosomes' whole numbers are taken from
            # (an isochromosome's arms differ, and neither is the number of chromosomes)
            span = np.array([max(sg.end - sg.start, 1) for sg in segs], float)
            mean = float(np.sum(span * np.array([sg.copies for sg in segs])) / span.sum())
            chroms[c] = Chrom(c, len(idx), main.copies, main.se, mean, segs, places=len(first), far=far, tilt=tilt.get(c, 0.0))
            # each arm on its own: a whole chromosome gained or lost shows in both alike
            la = arm[first]
            for a_name in ("p", "q"):
                on = la == a_name
                if on.sum() >= ARM_MIN_LOCI and (la != a_name).sum() >= ARM_MIN_LOCI:
                    lv, se, _ = _level(xl[on], sl[on], ref, float(np.round(_wmedian(xl[on], 1 / sl[on] ** 2) / GRID_STEP) * GRID_STEP))
                    chroms[c].arms[a_name] = (lv, float(np.sqrt(se ** 2 + max(lv, C_FLOOR * ref) ** 2 * (g * g * var_s + var_t))), int(nreg[on].sum()))
            xy[c] = (idx, segs, main, ref)
        # the whole numbers: two on an autosome; on the sex chromosomes the nearest to the mean over the chromosome
        for c, (idx, segs, main, ref) in xy.items():
            ch = chroms[c]
            sexchr = tab.kind[idx[0]] != "A"
            expected = int(np.clip(_round(ch.mean), 0, 4)) if sexchr else 2
            # a sex chromosome's class is its whole number; one the cohort held too few of takes the nearest that was measured
            key = next((q for q in ([f"{c}:{n_}" for n_ in (expected, min(expected, 2), 1)] if sexchr else [c]) if q in model.phi), c)
            ph, fl = max(model.phi.get(key, 1.0), 1.0), max(model.floor.get(key, 0.0), 0.0)
            for sg in segs:
                # the error as the model's cohort found it: wider by a factor, and no smaller than a floor that scales with the level
                se = float(np.hypot(sg.se * ph, fl * max(sg.copies, C_FLOOR * ref) / ref))
                wn = int(np.clip(_round(sg.copies), 0, 6))
                d_whole, d_exp = sg.copies - wn, sg.copies - expected
                frac = abs(d_whole) >= min_off and abs(d_whole) >= level_z * se
                if sg is main:
                    ch.whole, ch.off, ch.z, ch.se, ch.se0 = wn, d_whole, d_whole / se, se, sg.se
                    ch.status = "uncertain" if se > MAX_SE else "fractional" if frac else "settled"
                    # places that no stretch accounts for: too few to be a stretch of their own, too many and too far to be
                    # noise or a copy-number variant. One level does not describe this chromosome
                    if ch.far >= MIN_FAR and ch.far >= FAR_SHARE * ch.places:
                        notes.append(f"{c}: {ch.far} of its {ch.places} places lie far from its level, a structure these regions do not resolve")
                        if ch.status == "settled":
                            ch.status = "uncertain"
                if (wn != expected and not frac) or (frac and abs(d_exp) >= min_off and abs(d_exp) >= level_z * se):
                    events.append(Event(c, sg.span, sg.start, sg.end, sg.regions, sg.copies, se, d_exp, d_exp / se, whole=not frac))
                    new_off[idx[sg.first:sg.last + 1]] = True
                    if frac and ch.status == "settled":
                        ch.status = "fractional"
        out = Reading(chroms, events, level, g, int(ok.sum()), s_all, off=new_off, iterations=it + 1,
                      tilt=None if coef is None else (float(coef[0]), float(coef[1]), slope_se), note="; ".join(notes))
        if not (new_off & ~off).any():
            break
        off = off | new_off
    if "chrX" in out.chromosomes and np.isfinite(out.chromosomes["chrX"].copies):
        out.x = out.chromosomes["chrX"].copies
    if "chrY" in out.chromosomes and np.isfinite(out.chromosomes["chrY"].copies):
        out.y = out.chromosomes["chrY"].copies
    if out.x is not None:
        mean = lambda c: out.chromosomes[c].mean
        nx, ny = int(np.clip(_round(mean("chrX")), 0, 5)), int(np.clip(_round(mean("chrY")), 0, 4)) if out.y is not None else 0
        out.complement = "X" * nx + "Y" * ny
        if out.y is not None:
            dx, dy = out.x - nx, out.y - ny
            sx, sy = out.chromosomes["chrX"].se, out.chromosomes["chrY"].se
            if abs(dx) >= min_off and abs(dy) >= min_off and dx * dy < 0 and abs(dx + dy) < 3 * np.hypot(sx, sy) + min_off:
                out.note = "; ".join(filter(None, [out.note, "X and Y are off their whole numbers by the same amount in opposite directions: cells of two kinds, or two "
                                                              "people's DNA"]))
    return out


def _lstsq_rows(T: np.ndarray, e: np.ndarray) -> np.ndarray:
    return np.linalg.lstsq(T, e, rcond=None)[0] if len(e) > T.shape[1] + 2 else np.zeros(T.shape[1])


def raw_sex(tab: Table, L: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Each genome's raw reading of chrX and chrY in copies (two times the median region's ratio to the autosomal level):
    good enough to say which genomes hold one X, two, one Y, none - the ones a model learns the sex chromosomes' scales on."""
    import warnings
    out = []
    for kind in ("X", "Y"):
        m = tab.kind == kind
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            out.append(2.0 * np.exp(np.nanmedian(L[:, m], axis=1)) if m.any() else np.full(len(L), np.nan))
    return out[0], out[1]


def fit(names: list[str], Y: np.ndarray, arms: dict[str, int], k_max: int = K_MAX, rounds: int = 2, tau: float = TAU, log=None,
        gc: dict[str, float] | None = None) -> Model:
    """Learn a model on a cohort. `Y`: genomes x regions, log(observed / expected), NaN where a genome lacks a region.

    Efficiencies are medians and spreads robust SDs, so the cohort need not be free of aneuploid genomes; and what a first
    read of every genome finds off its whole number (a chromosome, a stretch) is left out of the next fit, so that a
    trisomy common in the cohort (chromosome 12 in lymphoblastoid lines) neither shifts its regions' efficiencies nor
    becomes a component."""
    import warnings
    from . import pcselect
    say = log or (lambda *a, **k: None)
    Y = np.asarray(Y, float)
    n, m = Y.shape
    tab = table(names, arms)
    A, X, Yk = tab.kind == "A", tab.kind == "X", tab.kind == "Y"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lev = np.nanmedian(Y[:, A], axis=1)
    L = Y - lev[:, None]
    xr, yr = raw_sex(tab, L)
    one_x, two_x, one_y = (xr > ONE[0]) & (xr < ONE[1]), (xr > TWO[0]) & (xr < TWO[1]), (yr > ONE[0]) & (yr < ONE[1])
    info = dict(one_x=int(one_x.sum()), two_x=int(two_x.sum()), one_y=int(one_y.sum()))
    # what each region is learned on: every genome for an autosomal region, the one-X genomes for chrX, the one-Y for chrY
    rows = np.ones((n, m), bool)
    x_from = "one X"
    if one_x.sum() >= 10:
        rows[:, X] = one_x[:, None]
    elif two_x.sum() >= 10:
        rows[:, X] = two_x[:, None]                    # no men: the two-X genomes, a second X taken to read as the first
        x_from = "two X (no genomes with one)"
    else:
        rows[:, X] = False
    rows[:, Yk] = one_y[:, None] if one_y.sum() >= 10 else False
    shift = np.where(A, 0.0, np.log(2.0) if x_from == "one X" else 0.0)
    shift = np.where(Yk, np.log(2.0), shift)
    Z = L + shift[None, :]                             # one copy of a sex chromosome at 0, like two of an autosome
    K = int(max(0, min(k_max, n // MIN_PER_COMPONENT)))
    mask = np.zeros((n, m), bool)                      # genome x region: in an event, left out
    model = None
    for rnd in range(rounds + 1):
        use = rows & ~mask & np.isfinite(Z)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            Zu = np.where(use, Z, np.nan)
            a = np.nanmedian(Zu, axis=0)
            E = Zu - a[None, :]
            sd = _mad(E, axis=0)
        few = use.sum(0) < 10
        a[few], sd[few] = np.nan, np.nan
        V, lam2, T, dropped, above = np.zeros((m, 0)), np.zeros(0), np.zeros((n, 0)), [], 0
        resid = sd.copy()                              # each region's spread within the cohort (sd: as a genome outside it will show)
        if rnd > 0 and K > 0:
            typical = np.nanmedian(sd[A])
            core = A & (np.mean(np.isfinite(Y), axis=0) >= 0.95) & np.isfinite(sd) & (sd <= MAX_REGION_SD * typical)
            # a region's values beyond CLIP of its own SDs are pulled in: a copy-number variant that a few percent of the
            # cohort carry in one region would otherwise be a component of its own
            Ec = np.clip(np.nan_to_num(E[:, core]), -CLIP * sd[core][None, :], CLIP * sd[core][None, :])
            held = use[:, core]
            # contrasts within chromosomes: each genome's regions less their chromosome's mean (a chromosome gained or
            # lost in some of a cohort's genomes is then no component, however common it is)
            for c in np.unique(tab.chrom[core]):
                on = tab.chrom[core] == c
                cnt = held[:, on].sum(1, keepdims=True)
                mu = np.where(cnt > 0, (Ec[:, on] * held[:, on]).sum(1, keepdims=True) / np.maximum(cnt, 1), 0.0)
                Ec[:, on] = np.where(held[:, on] & (cnt >= MIN_CHROM_REGIONS), Ec[:, on] - mu, 0.0)
            Ec = Ec - Ec.mean(0)
            U, S, Wt = np.linalg.svd(Ec, full_matrices=False)
            Wt = Wt.T
            # a mode of the libraries is spread over the genome; a component that lives on one chromosome is a stretch of it
            # that part of the cohort has gained or lost, and is left out
            # how many components are structure is decided where the control-region PCs' count is: at the edge of the noise
            # bulk (pcselect). A component inside the bulk is noise of this cohort's genomes, and taking it out of them
            # would make their regions look quieter than those of the next genome
            above = pcselect.mp_select(S, *Ec.shape).n_pc
            share = np.array([[float((Wt[tab.chrom[core] == c, j] ** 2).sum()) for c in np.unique(tab.chrom[core])] for j in range(min(above, len(S)))])
            lives = share.max(1) / np.maximum(share.sum(1), 1e-30) > MAX_CHROM_SHARE if above else np.zeros(0, bool)
            keep = np.flatnonzero(~lives)[:K]
            dropped = [int(j) + 1 for j in np.flatnonzero(lives) if j < (keep[-1] if len(keep) else 0)]
            T = U[:, keep] * S[keep]
            lam2 = np.maximum((S[keep] ** 2) / n, 1e-12)
            V = np.zeros((m, len(keep)))
            for j in np.flatnonzero(np.isfinite(a)):
                r = use[:, j]
                lim = CLIP * sd[j] if np.isfinite(sd[j]) and sd[j] > 0 else 0.5
                b = _lstsq_rows(np.c_[np.ones(r.sum()), T[r]], np.clip(E[r, j], -lim, lim))
                if b.any():
                    a[j] += b[0]
                    V[j] = b[1:]
                    resid[j] = _mad(E[r, j] - b[0] - T[r] @ b[1:])
                    # the spread a genome outside the cohort will show: the residuals' own, wider by what the fit took from
                    # them and by the loadings' own errors (a region learned on few genomes is known less well)
                    nj, pj = int(r.sum()), len(b)
                    sd[j] = resid[j] * np.sqrt((nj + pj) / (nj - pj))
        typical = np.nanmedian(sd[A]) if np.isfinite(sd[A]).any() else np.nan
        bad = A & np.isfinite(sd) & (sd > MAX_REGION_SD * typical)
        sd_use = np.where(bad, np.nan, sd)
        # a second X: what it reads of the first, per region, and how it follows the modes
        u, x2 = np.full(m, np.nan), np.zeros(V.shape[1])
        if X.any() and x_from == "one X" and two_x.sum() >= 10:
            r1 = Z[two_x][:, X] - a[X][None, :] - (T[two_x] @ V[X].T if V.shape[1] else 0.0)
            ok2 = ~mask[two_x][:, X]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                u[X] = np.exp(np.nanmedian(np.where(ok2, r1, np.nan), axis=0)) - 1.0
            if V.shape[1]:
                wx = 1.0 / np.where(np.isfinite(sd[X]), sd[X], np.inf) ** 2
                R = np.nansum(wx[None, :] * (np.exp(r1) - 1.0) / u[X][None, :], axis=1) / np.nansum(np.where(np.isfinite(r1), wx[None, :], 0.0), axis=1)
                good = np.isfinite(R) & (R > 0.5) & (R < 1.5)
                if good.sum() > V.shape[1] + 10:
                    x2 = _lstsq_rows(np.c_[np.ones(good.sum()), T[two_x][good]], np.log(R[good]))[1:]
        elif X.any():
            u[X] = 1.0
        shared, n_places = _shared(tab, a, resid, V, T, E, use & A[None, :], np.isfinite(sd_use))
        model = Model(list(names), a, sd_use, V, lam2, u, x2, {}, dict(arms), n, dict(info, x_scale=x_from, components=int(V.shape[1]),
                                                                                    regions_not_used=int(bad.sum()), components_above_noise=int(above),
                                                                                    components_on_one_chromosome=dropped, places_measured=n_places),
                      local=shared, gc=dict(gc or {}))
        if rnd == rounds:
            break
        # read every genome with what is known so far: what is off its whole number stays out of the next fit
        index = np.arange(m)
        new = np.zeros((n, m), bool)
        for i in range(n):
            rd = read(model, names, np.where(np.isfinite(Y[i]), Y[i], np.nan), tau=tau, tab=tab, index=np.where(np.isfinite(Y[i]), index, -1))
            if rd is not None and rd.off is not None:
                new[i] = rd.off
        mask = new
        say(f"[karyotype] round {rnd + 1}: {int(mask.any(1).sum())} of {n} genomes hold a chromosome or a stretch off its whole number; left out of the fit there")
    return model


def _shared(tab: Table, a: np.ndarray, resid: np.ndarray, V: np.ndarray, T: np.ndarray, E: np.ndarray, use: np.ndarray, ok: np.ndarray,
            min_places: int = 20, min_genomes: int = 30) -> tuple[float, int]:
    """What the pieces of one place share beyond the modes: the SD (log) of a level common to them in one genome. For every
    place of two or more pieces, the spread over the cohort of the pieces' weighted mean is set against what the pieces'
    own spreads predict for independent pieces; the excess, per place, is the shared part's variance, and the median over
    the places is kept (zero with fewer than `min_places` places to measure it on). `E`: genomes x regions, less the
    efficiencies before the modes were fitted; `use`: the entries that count; `ok`: the regions in use."""
    est = []
    for c in tab.chromosomes():
        idx = np.flatnonzero(ok & (tab.chrom == c) & (tab.start >= 0) & np.isfinite(resid) & (resid > 0))
        if len(idx) < 2:
            continue
        idx = idx[np.argsort(tab.mid[idx], kind="stable")]
        first, _, _, nreg = loci(tab.start[idx], tab.end[idx], np.zeros(len(idx)), resid[idx])
        for f, k in zip(first, nreg):
            if k < 2:
                continue
            js = idx[f:f + k]
            rows = use[:, js].all(1)
            if rows.sum() < min_genomes:
                continue
            R = E[rows][:, js] - (T[rows] @ V[js].T if V.shape[1] else 0.0)
            R = R - np.median(R, axis=0)
            w = 1.0 / resid[js] ** 2
            mean = (R * w).sum(1) / w.sum()
            share = 1.0 - float(np.sum(w * w)) / float(w.sum()) ** 2
            est.append((float(_mad(mean)) ** 2 - 1.0 / float(w.sum())) / share)
    if len(est) < min_places:
        return 0.0, len(est)
    return float(min(np.sqrt(max(np.median(est), 0.0)), MAX_LOCAL)), len(est)


def _spread_fit(off: np.ndarray, se: np.ndarray, bins: int = 12, per_bin: int = 150) -> tuple[float, float]:
    """(factor, floor): the cohort's spread of a level as factor^2 * se^2 + floor^2, fitted to the robust spread of the
    levels in bins of their error. The factor is at least one and the floor at least zero."""
    o = np.argsort(se)
    off, se = off[o], se[o]
    nb = int(min(bins, len(off) // per_bin))
    if nb < 2:
        return float(max(_mad(off / se), 1.0)), 0.0
    cut = np.linspace(0, len(off), nb + 1).astype(int)
    e = np.array([np.median(se[i:j]) ** 2 for i, j in zip(cut[:-1], cut[1:])])
    v = np.array([_mad(off[i:j]) ** 2 for i, j in zip(cut[:-1], cut[1:])])
    good = v > 0
    e, v = e[good], v[good]
    if len(v) < 2:
        return float(max(_mad(off / se), 1.0)), 0.0
    # relative errors: (m2 * e + s2) / v = 1 in every bin
    m2, s2 = np.linalg.lstsq(np.c_[e / v, 1.0 / v], np.ones(len(v)), rcond=None)[0]
    if s2 < 0:
        m2, s2 = float(np.sum(e / v) / np.sum((e / v) ** 2)), 0.0
    if m2 < 1:
        m2, s2 = 1.0, float(max(np.sum((1.0 - e / v) / v) / np.sum(1.0 / v ** 2), 0.0))
    return float(np.sqrt(m2)), float(np.sqrt(s2))


def _floor_for(off: np.ndarray, se: np.ndarray, factor: float) -> float:
    """The floor that brings the robust SD of off / hypot(factor * se, floor) to one (zero where it is no more than one without)."""
    spread = lambda fl: float(_mad(off / np.hypot(factor * se, fl)))
    if spread(0.0) <= 1.0:
        return 0.0
    lo, hi = 0.0, float(10 * _mad(off) + 1e-6)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if spread(mid) > 1.0 else (lo, mid)
    return 0.5 * (lo + hi)


def calibrate(model: Model, readings: list[Reading | None]) -> tuple[dict[str, float], dict[str, float], dict]:
    """How far a level's error is from its regions' noise, measured on readings made before any such widening (`Chrom.se0`).

    Over the autosomes of all the genomes read, the robust spread of the levels' distances from their whole numbers is
    fitted as factor^2 * se^2 + floor^2: chromosomes differ several-fold in how many regions they hold, so the two parts
    separate. A chromosome whose z is then still wider than one is widened by its own factor. The sex chromosomes are taken
    by their whole number (`chrX:1`, `chrX:2`, `chrY:1`), with the autosomes' factor and a floor of their own: a second X
    is read on another scale, and how much of the first it reads differs a little from one cell line to the next.
    Returns (phi, floor, a summary): `Model.phi`, `Model.floor` (in copies at the reference level: two on an autosome,
    one on a sex chromosome)."""
    pairs: dict[str, list] = {}
    for r in readings:
        if r is None:
            continue
        for c, ch in r.chromosomes.items():
            if ch.se0 is None or ch.off is None or ch.status == "not read" or not np.isfinite(ch.se0) or ch.se0 <= 0:
                continue
            pairs.setdefault(c if c in AUTOSOMES else f"{c}:{ch.whole}", []).append((ch.off, ch.se0))
    auto = [np.array(v) for k, v in pairs.items() if k in AUTOSOMES]
    every = np.concatenate(auto) if auto else np.zeros((0, 2))
    m, s = _spread_fit(every[:, 0], every[:, 1]) if len(every) >= MIN_CALIBRATE else (1.0, 0.0)
    phi, floor = {}, {}
    for key in sorted(pairs, key=_natural):
        v = np.array(pairs[key])
        if len(v) < MIN_CALIBRATE:
            continue
        if key in AUTOSOMES:
            f = float(max(_mad(v[:, 0] / np.hypot(m * v[:, 1], s)), 1.0))
            phi[key], floor[key] = m * f, s * f
        elif key in ("chrX:1", "chrX:2", "chrY:1"):
            phi[key], floor[key] = m, _floor_for(v[:, 0], v[:, 1], m) / int(key[-1])
    return phi, floor, dict(factor=round(m, 3), floor=round(s, 5), levels=int(len(every)))


def assemble(vectors: list) -> tuple[list[str], np.ndarray]:
    """Per-genome (names, values) or None -> the names in first-seen order and a genomes x regions matrix, NaN where a
    genome lacks a region (a lighter fetch, an estimate written before a region existed)."""
    names: list[str] = []
    col: dict[str, int] = {}
    cache: dict[int, np.ndarray] = {}
    idx = []
    for v in vectors:
        if v is None:
            idx.append(None)
            continue
        key = hash(tuple(v[0]))
        if key not in cache:
            for n in v[0]:
                if n not in col:
                    col[n] = len(names)
                    names.append(n)
            cache[key] = np.array([col[n] for n in v[0]])
        idx.append(cache[key])
    Y = np.full((len(vectors), len(names)), np.nan, np.float32)
    for i, (v, ix) in enumerate(zip(vectors, idx)):
        if v is not None:
            Y[i, ix] = v[1]
    return names, Y


def cohort(vectors: list, arms: dict[str, int], model: Model | None = None, rules: dict | None = None, fit_own: bool | None = None,
           log=None, gc: dict[str, float] | None = None) -> tuple[list, Model | None, dict]:
    """Read a cohort's genomes. `vectors`: per genome, `gather`'s result. A cohort of at least `min_fit` genomes learns
    its own model (unless `fit_own` is False); a smaller one is read against `model` (a saved one, or the bundle's), and
    is not read at all without one. Returns (a Reading or None per genome, the model used, a summary)."""
    say = log or (lambda *a, **k: None)
    rules = rules or {}
    tau, level_z, min_off = float(rules.get("tau", TAU)), float(rules.get("level_z", LEVEL_Z)), float(rules.get("min_off", MIN_OFF))
    have = [v for v in vectors if v is not None]
    info = dict(n=len(vectors), n_with_regions=len(have))
    if not have:
        return [None] * len(vectors), None, dict(info, model="none", why="no estimate holds single-copy regions")
    names, Y = assemble(vectors)
    n_fit = int(rules.get("min_fit", MIN_FIT))
    own = fit_own if fit_own is not None else model is None or len(have) >= n_fit
    if own and len(have) < n_fit:
        if model is None:
            say(f"[karyotype] {len(have)} genomes, fewer than the {n_fit} a model of the regions needs, and no saved model: chromosomes are not read "
                "(give one with --karyotype-model)")
            return [None] * len(vectors), None, dict(info, model="none", why=f"fewer than {n_fit} genomes and no saved model")
        own = False
    if own:
        present = np.isfinite(Y).any(1)
        model = fit(names, Y[present].astype(float), arms, k_max=int(rules.get("k_max", K_MAX)), tau=tau, log=say, gc=gc)
        info.update(model="this cohort's", **model.info)
    else:
        known = set(model.names)
        unknown = sum(1 for n in names if n not in known)
        if unknown == len(names):
            say(f"[karyotype] {len(have)} genomes, fewer than the {n_fit} a model of the regions needs, and the saved model knows none of their {len(names)} regions "
                "(counts made with another bundle's regions): chromosomes are not read")
            return [None] * len(vectors), None, dict(info, model="none", why="the saved model knows none of the regions")
        info.update(model="saved", model_n=model.n, regions_unknown_to_model=unknown)
        if unknown:
            say(f"[karyotype] {unknown} of the cohort's {len(names)} regions are not in the saved model and are left out")
        if fit_own is None:
            model = repin_second_x(model, names, Y, int(rules.get("min_repin", MIN_REPIN)), say)
    tab = table(model.names, model.arms)
    pos = {n: i for i, n in enumerate(names)}
    col = np.array([pos.get(n, -1) for n in model.names])

    def read_all():
        out = []
        for i in range(len(vectors)):
            if vectors[i] is None:
                out.append(None)
                continue
            y = Y[i].astype(float)
            index = np.where((col >= 0) & np.isfinite(y[np.maximum(col, 0)]), col, -1)
            out.append(read(model, names, y, tau=tau, level_z=level_z, min_off=min_off, tab=tab, index=index))
        return out
    readings = read_all()
    if own:
        model.phi, model.floor, cal = calibrate(model, readings)
        model.info["error"] = cal
        info["error"] = cal
        say(f"[karyotype] a level's error, measured on the cohort's {cal['levels']:,} autosomal levels: {cal['factor']:.2f} times its regions' noise, "
            f"with a floor of {cal['floor']:.4f} copies; what the pieces of one place share: SD {model.local:.4f} (log)")
        readings = read_all()                                   # the statuses rest on z, and z on the calibration
    rs = [r for r in readings if r is not None]
    st = {k: sum(1 for r in rs if r.status == k) for k in ("settled", "fractional", "uncertain", "not read")}
    info.update(read=len(rs), components=model.k, regions=len(model.names), **{f"n_{k.replace(' ', '_')}": v for k, v in st.items()},
                with_event=sum(1 for r in rs if r.events))
    say(f"[karyotype] {len(rs)} genomes read against {info['model']} model ({len(model.names)} regions, {model.k} components): "
        f"{st['settled']} settled, {st['fractional']} with a chromosome or a stretch between whole numbers, {st['uncertain']} uncertain; "
        f"{sum(1 for r in rs if any(e.whole for e in r.events))} hold a whole chromosome or stretch more or less than two in every cell")
    return readings, model, info


def repin_second_x(model: Model, names: list[str], Y: np.ndarray, min_genomes: int = MIN_REPIN, log=None) -> Model:
    """A saved model read by a cohort too small to learn its own: what a second X reads of the first is a property of the
    DNA (0.95 in lymphoblastoid lines, whose inactive X replicates late; other tissues and libraries differ by a few
    percent), and one number can be measured on few genomes. With at least `min_genomes` two-X genomes, the model's
    factors are scaled to the cohort's median (a copy of the model is returned); with fewer the model stands, and the log
    says on whose scale a second X is read."""
    say = log or (lambda *a, **k: None)
    tab = table(model.names, model.arms)
    X = (tab.kind == "X") & np.isfinite(model.a) & np.isfinite(model.sd) & np.isfinite(model.u)
    A = (tab.kind == "A") & np.isfinite(model.a)
    pos = {n: i for i, n in enumerate(names)}
    col = np.array([pos.get(n, -1) for n in model.names])
    if not X.any() or (col[X] < 0).all():
        return model
    E = np.where(col[None, :] >= 0, np.asarray(Y, float)[:, np.maximum(col, 0)], np.nan) - model.a[None, :]
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lev = np.nanmedian(E[:, A], axis=1)
        r1 = E[:, X] - lev[:, None] + np.log(2.0)                               # one X at 0
        reading = np.exp(np.nanmedian(r1, axis=1))
    u_model = float(np.nanmedian(model.u[X]))
    two = np.isfinite(reading) & (reading > TWO[0] * (1 + u_model) / 2) & (reading < TWO[1] * (1 + u_model) / 2)
    if two.sum() < min_genomes:
        if np.isfinite(reading).any() and (reading > 1.5).any():
            say(f"[karyotype] a second X is read on the saved model's scale ({u_model:.3f} of the first, as measured in the model's cohort): in DNA of another "
                f"tissue or library it can read a few percent off, and a cohort of at least {min_genomes} two-X genomes measures it itself")
        return model
    u_here = float(np.median(reading[two])) - 1.0
    out = replace(model, u=np.where(np.isfinite(model.u), model.u * u_here / u_model, np.nan), phi=dict(model.phi), floor=dict(model.floor),
                  info=dict(model.info, second_x=dict(model=round(u_model, 4), cohort=round(u_here, 4), n=int(two.sum()))))
    say(f"[karyotype] a second X reads {u_here:.3f} of the first in the {int(two.sum())} two-X genomes here ({u_model:.3f} in the saved model's cohort): the saved factors are scaled to it")
    return out


def columns(rd: Reading | None, chroms: list[str]) -> dict:
    """A genome's row: the karyotype as written, its status, the sex chromosomes as letters, every chromosome's copies
    and z (the main level's distance from its whole number, in SDs), and the events with their levels."""
    if rd is None:
        return {"karyotype": None, "karyotype.status": "not read"}
    out = {"karyotype": rd.karyotype() or None, "karyotype.status": rd.status, "sex_chromosomes": rd.complement or None,
           "karyotype.events": ";".join(f"{e.label()}|{e.copies:.3f}|z{e.z:+.1f}" for e in rd.events) or "none",
           "karyotype.noise": round(rd.noise, 3), "karyotype.regions": rd.regions,
           "karyotype.gc_tilt": None if rd.tilt is None else round(rd.gc_tilt, 4),
           "karyotype.gc_tilt_z": None if rd.gc_tilt_z is None else round(rd.gc_tilt_z, 2)}
    if rd.note:
        out["karyotype.note"] = rd.note
    for c in chroms:
        ch = rd.chromosomes.get(c)
        ok = ch is not None and np.isfinite(ch.copies)
        out[f"{c}.copies"] = round(ch.copies, 3) if ok else None
        out[f"{c}.z"] = round(ch.z, 2) if ok and ch.z is not None else None
    return out


def long_table(samples: list[str], readings: list) -> list[dict]:
    """One row per chromosome of every genome, and one per stretch where a chromosome holds more than one level."""
    rows = []
    for s, rd in zip(samples, readings):
        if rd is None:
            continue
        ev = {(e.chrom, e.start, e.end): e for e in rd.events}
        for c in sorted(rd.chromosomes, key=_natural):
            ch = rd.chromosomes[c]
            if not np.isfinite(ch.copies):
                rows.append(dict(sample=s, chrom=c, kind="chromosome", regions=ch.regions, status=ch.status))
                continue
            rows.append(dict(sample=s, chrom=c, kind="chromosome", span="whole", start=ch.segments[0].start, end=ch.segments[-1].end, regions=ch.regions,
                             loci=sum(g.loci for g in ch.segments), copies=round(ch.copies, 3), se=round(ch.se, 4), whole=ch.whole, off=round(ch.off, 3),
                             z=round(ch.z, 2), status=ch.status, mean=round(ch.mean, 3), stretches=len(ch.segments),
                             p=round(ch.arms["p"][0], 3) if "p" in ch.arms else None, q=round(ch.arms["q"][0], 3) if "q" in ch.arms else None))
            if len(ch.segments) > 1:
                for g in ch.segments:
                    e = ev.get((g.chrom, g.start, g.end))
                    rows.append(dict(sample=s, chrom=c, kind="stretch", span=g.span, start=g.start, end=g.end, regions=g.regions, loci=g.loci,
                                     copies=round(g.copies, 3), se=round(g.se, 4), z=None if e is None else round(e.z, 2),
                                     status="" if e is None else "every cell" if e.whole else "part of the cells"))
    return rows


LONG_COLUMNS = ("sample", "chrom", "kind", "span", "start", "end", "regions", "loci", "copies", "se", "whole", "off", "z", "status", "mean", "stretches", "p", "q")
