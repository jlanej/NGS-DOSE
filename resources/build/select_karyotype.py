#!/usr/bin/env python3
"""Select the karyotype windows: single-copy sequence spread along every chromosome arm, for reading chromosomes in copies.

The control regions were chosen to make a denominator and a GC curve: rich in rare GC strata, wherever on the
autosomes those are. Read as chromosomes they are uneven - 89 regions on chr1, 5 on chr19, 2 on chr22, none on 18p -
and each costs a fetch more than it counts: a region of 12.5 kb lies in 1.3 CRAM slices of about 40 kb, so a quarter
of the bases fetched for it are counted. The karyotype windows are chosen for the other purpose:

  1. the same clean sequence as the controls (the complement of the NGS-PCA exclusion set: 10x SV blacklist, GEM
     mappability < 1, DGV variants, segmental duplications), in runs of at least --min-piece bp, less the regions
     the bundle already has; a run longer than --max-piece is cut into equal pieces, so that a region stays a unit
     that one copy-number variant can move and the chain along the chromosome can outvote;
  2. pieces are grouped into windows of at most --window bp, and a window is worth its clean bases per CRAM slice
     it would be fetched in ((span + slice) / slice slices): a window that fills its slices with countable bases
     costs a third of what scattered regions of the same total length do;
  3. every chromosome is given --per-chromosome clean bases (the precision of a chromosome is set by the bases
     counted on it, so every chromosome is read as precisely as every other), split between its arms, an arm that
     lacks the sequence making it up on the other; within an arm the windows are spread: the arm is cut into as
     many equal stretches as windows are wanted, the best window of each stretch is taken, and what is still
     missing is added from the best windows left;
  4. chrX outside the pseudoautosomal regions and the X-transposed region, chrY in X-degenerate sequence, as for
     the controls' known-truth regions.

  5. what sequence looks like and how it reads are two things: the selection is made --over times as large as
     wanted, counted in real genomes, and the pieces that read badly there (--drop: a BED of them, kept among the
     bundle's build inputs with the reason for each) are left out; each arm is then trimmed back to what was wanted,
     its windows with the fewest clean bases per slice going first.

Windows are numbered along each arm; every fourth is of tier 1, the rest of tier 2: tier 1 alone is a lighter set
that is still spread along every arm (`make_control_sets.py` writes the sets). On chrX and chrY every window is of
tier 1: the sex chromosomes are what a screen is most often for, and their structures (an isochromosome of the long
arm, an isodicentric Y) show only along the chromosome.

Output: a BED of pieces (name column `test:karyotype`, `test:karyotype.chrX`, `test:karyotype.chrY`) and a table of
the pieces with their window, arm and tier (--table).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ngsdose.io import FastaIndex  # noqa: E402
from select_controls import L_MAX, X_EXCLUDE_GRCH38, Y_EXCLUDE_GRCH38, complement_runs, read_exclusions  # noqa: E402

LABEL = {"chrX": "test:karyotype.chrX", "chrY": "test:karyotype.chrY"}
# GRCh38 chrY beyond the controls' exclusions: the heterochromatic long arm (the exclusion set leaves stray runs in it)
Y_HETEROCHROMATIN_GRCH38 = [(26_600_000, 57_227_415)]


def subtract(runs, taken):
    """Runs less the intervals in `taken` (both sorted lists of (start, end))."""
    out, j = [], 0
    taken = sorted(taken)
    for s, e in runs:
        while j < len(taken) and taken[j][1] <= s:
            j += 1
        pos, k = s, j
        while k < len(taken) and taken[k][0] < e:
            if taken[k][0] > pos:
                out.append((pos, taken[k][0]))
            pos = max(pos, taken[k][1])
            k += 1
        if pos < e:
            out.append((pos, e))
    return out


def pieces(runs, min_piece, max_piece):
    out = []
    for s, e in runs:
        if e - s < min_piece:
            continue
        n = max(1, -(-(e - s) // max_piece))
        cut = [s + (e - s) * i // n for i in range(n + 1)]
        out += [(a, b) for a, b in zip(cut[:-1], cut[1:])]
    return out


def windows(pcs, window, slice_bp):
    """Candidate windows: from every piece, the pieces that follow it within `window` bp. (score, first piece, last piece)"""
    out = []
    for i in range(len(pcs)):
        j = i
        while j + 1 < len(pcs) and pcs[j + 1][1] - pcs[i][0] <= window:
            j += 1
        clean = sum(b - a for a, b in pcs[i:j + 1])
        span = pcs[j][1] - pcs[i][0]
        out.append((clean / ((span + slice_bp) / slice_bp), i, j, clean))
    return out


def choose(pcs, target, window, slice_bp):
    """Windows of an arm, spread along it, holding about `target` clean bases; returns (first piece, last piece) per window.

    The arm is cut into n equal stretches and the best window that starts in each is taken, for the smallest n that
    reaches the target; an arm too sparse for that takes the best windows left as well, until the target or the end."""
    if len(pcs) < 3:
        return []
    cand = windows(pcs, window, slice_bp)
    lo, hi = pcs[0][0], pcs[-1][1]

    def spread(n):
        edges = np.linspace(lo, hi, n + 1)
        used, chosen, total = np.zeros(len(pcs), bool), [], 0
        for a, b in zip(edges[:-1], edges[1:]):
            for c in sorted((c for c in cand if a <= pcs[c[1]][0] < b), reverse=True):
                if not used[c[1]:c[2] + 1].any():
                    used[c[1]:c[2] + 1] = True
                    chosen.append((c[1], c[2]))
                    total += c[3]
                    break
        return chosen, used, total
    chosen, used, total = spread(2)
    for n in range(3, max(4, len(pcs))):
        if total >= target:
            break
        got = spread(n)
        if got[2] > total:
            chosen, used, total = got
    for c in sorted(cand, reverse=True):                          # an arm too sparse: the best of what is left
        if total >= target:
            break
        if not used[c[1]:c[2] + 1].any():
            used[c[1]:c[2] + 1] = True
            chosen.append((c[1], c[2]))
            total += c[3]
    return sorted(chosen)


def arm_targets(avail: dict, total: float) -> dict:
    """Clean bases wanted on each arm: half to each; an arm that lacks the sequence is made up on the other."""
    t = {k: min(avail[k], total / 2) for k in "pq"}
    for k, o in (("p", "q"), ("q", "p")):
        t[k] = min(avail[k], t[k] + max(0.0, total / 2 - t[o]))
    return t


def clean(window) -> int:
    return sum(e - s for s, e in window)


def density(window, slice_bp: int) -> float:
    """Clean bases per CRAM slice the window is fetched in."""
    return clean(window) / ((window[-1][1] - window[0][0] + slice_bp) / slice_bp)


def thinnest(windows, slice_bp: int):
    return min(windows, key=lambda w: density(w, slice_bp))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--exclude", required=True, help="NGS-PCA exclusion BED(.gz)")
    ap.add_argument("--controls", required=True, help="the bundle's regions so far (BED): the windows keep clear of them")
    ap.add_argument("--arms", required=True, help="JSON: per chromosome, where the short arm ends")
    ap.add_argument("--out", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--per-chromosome", type=int, default=900_000)
    ap.add_argument("--chrY", type=int, default=600_000, dest="per_y")
    ap.add_argument("--min-piece", type=int, default=5000)
    ap.add_argument("--max-piece", type=int, default=20000)
    ap.add_argument("--window", type=int, default=100_000)
    ap.add_argument("--slice", type=int, default=40_000, help="what a CRAM slice spans at 30-40x (10,000 reads)")
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--over", type=float, default=1.2, help="select this many times the clean bases wanted before --drop and the trim (1: no trim)")
    ap.add_argument("--drop", help="BED of pieces to leave out: the ones that read badly in real genomes")
    a = ap.parse_args()

    fa = FastaIndex(a.reference)
    arms = json.load(open(a.arms))
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    excl = read_exclusions(a.exclude, chroms)
    have = {c: [] for c in chroms}
    for line in open(a.controls):
        p = line.rstrip("\n").split("\t")
        if p and p[0] in have:
            have[p[0]].append((int(p[1]), int(p[2])))
    drop = set()
    if a.drop:
        for line in open(a.drop):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and not line.startswith("#"):
                drop.add((p[0], int(p[1]), int(p[2])))
    rows, summary = [], []
    for c in chroms:
        extra = {"chrX": X_EXCLUDE_GRCH38, "chrY": Y_EXCLUDE_GRCH38 + Y_HETEROCHROMATIN_GRCH38}.get(c, [])
        runs = subtract(complement_runs(excl[c], fa.length(c), extra), have[c])
        pcs = []
        for s, e in pieces(runs, a.min_piece, a.max_piece):
            if s < a.flank + L_MAX or e + a.flank + L_MAX > fa.length(c):
                continue
            if "N" in fa.fetch(c, s - L_MAX, e + L_MAX).upper():
                continue
            pcs.append((s, e))
        b = arms[c]
        by_arm = {"p": [x for x in pcs if (x[0] + x[1]) / 2 < b], "q": [x for x in pcs if (x[0] + x[1]) / 2 >= b]}
        avail = {k: sum(e - s for s, e in v) for k, v in by_arm.items()}
        total = a.per_y if c == "chrY" else a.per_chromosome
        want, over = arm_targets(avail, total), arm_targets(avail, total * a.over)
        for k in "pq":
            # more than wanted, then what read badly is left out, then the arm is trimmed back to what was wanted
            wins = [[x for x in by_arm[k][i:j + 1] if (c, x[0], x[1]) not in drop] for i, j in choose(by_arm[k], over[k], a.window, a.slice)]
            wins = [w for w in wins if w]
            while len(wins) > 2 and sum(map(clean, wins)) - clean(thinnest(wins, a.slice)) >= want[k]:
                wins.remove(thinnest(wins, a.slice))
            got = 0
            for w, pcs_w in enumerate(wins):
                tier = 1 if c in ("chrX", "chrY") or w % 4 == 1 or len(wins) < 2 else 2
                for s, e in pcs_w:
                    rows.append((c, s, e, LABEL.get(c, "test:karyotype"), f"{c[3:]}{k}.{w + 1}", k, tier))
                    got += e - s
            summary.append((c, k, len(by_arm[k]), avail[k], len(wins), got))
    rows.sort(key=lambda r: (chroms.index(r[0]), r[1]))
    with open(a.out, "w") as fh:
        for c, s, e, name, *_ in rows:
            fh.write(f"{c}\t{s}\t{e}\t{name}\n")
    with open(a.table, "w") as fh:
        fh.write("chrom\tstart\tend\twindow\tarm\ttier\n")
        for c, s, e, _, w, k, tier in rows:
            fh.write(f"{c}\t{s}\t{e}\t{w}\t{k}\t{tier}\n")
    print("arm    pieces  clean kb available   windows  clean kb chosen", file=sys.stderr)
    for c, k, n, av, nw, got in summary:
        if n:
            print(f"{c[3:]:>2}{k}   {n:6d}  {av / 1e3:10.0f}   {nw:7d}  {got / 1e3:8.0f}", file=sys.stderr)
    tot = sum(r[2] - r[1] for r in rows)
    print(f"karyotype windows: {len(rows)} pieces in {len({(r[0], r[4]) for r in rows})} windows, {tot:,} bp; "
          f"tier 1 {sum(r[2] - r[1] for r in rows if r[6] == 1):,} bp", file=sys.stderr)


if __name__ == "__main__":
    main()
