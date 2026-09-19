#!/usr/bin/env python3
"""Select single-copy control regions and known-copy-number test regions for one reference build.

Controls are what every NGS-DOSE estimate is a ratio against, and what the fragment-GC curve is
fitted on, so they must be (i) diploid and single-copy in essentially everyone, (ii) cleanly
mappable, (iii) spread over the autosomes so that an aneuploid chromosome or a replication-timing
domain cannot move the denominator, and (iv) rich in the rare GC strata - the curve has to be
known at the GC of the classes being measured (45S rDNA: 55-80%), and only ~0.05% of single-copy
sequence is that GC-rich at fragment scale.

Candidate sequence is the complement of the NGS-PCA exclusion set (10x SV blacklist, GEM
mappability < 1, DGV variants, segmental duplications). Selection is deterministic:

  1. candidate runs >= --min-run bp, long runs trimmed to their central --max-run bp;
  2. greedy fill, rarest GC stratum first, until every 1% fragment-GC bin holds --per-bin
     positions or the --n-gc budget of regions is spent;
  3. an even genomic spread of further runs up to --n-control regions;
  4. from what is left: --n-test-auto held-out autosomal runs (truth: 2 copies) and
     --n-test-x chrX runs outside the PARs and the X-transposed region (truth: 1 or 2 copies);
  5. --n-test-y chrY runs in X-degenerate sequence (truth: 1 copy in a male, none in a female);
  6. one stretch each of the mitochondrial genome and of the EBV decoy: not truths but dosages
     (copies per cell), the two standard covariates of the state of a cell line.

Output BED name column: `control`, `test:auto`, `test:chrX`, `test:chrY`, `dosage:chrM`,
`dosage:chrEBV`. Only `control` regions enter the GC curve and the denominator.
"""
import argparse
import gzip
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from ngsdose.io import FastaIndex  # noqa: E402

# GRCh38: PAR1, X-transposed region, PAR2
X_EXCLUDE_GRCH38 = [(0, 2_781_479), (88_400_000, 93_500_000), (155_701_382, 156_040_895)]
# GRCh38 chrY: PAR1 and the X-transposed region; the pericentromere; PAR2. The candidate runs were
# then screened in a female (NA12878, 30x): 54 of 56 hold < 0.4x of reads (0.05x overall); the
# two that do not - a pericentromeric run at 10.93 Mb (122x) and one at 16.53 Mb (9x) - are out.
Y_EXCLUDE_GRCH38 = [(0, 6_700_000), (10_000_000, 11_800_000), (16_528_000, 16_541_000), (56_887_902, 57_227_415)]
# Whole-contig dosages. chrM 1,000-2,500 (12S and 16S rRNA, the most conserved stretch): outside
# the D-loop and its 7S DNA, outside the 4,977-bp common deletion (8.5-13.4 kb), outside the
# stretch copied into the chr1:629-634 kb NUMT (3.9-9.7 kb), which takes mitochondrial reads -
# and short of position 3,107, where the reference carries an N (the rCRS placeholder): a window
# over an N has no GC, and a region with such windows cannot be modelled like the others.
# chrEBV 102-122 kb: unique sequence away from the W repeats (12-35 kb), the terminal repeats and
# the 139-151 kb segment that the B95-8 strain (the virus LCLs are made with) has lost. A few kb
# is plenty: at hundreds of copies per cell these are the deepest regions in the file.
DOSAGE_GRCH38 = [("chrM", 1_000, 2_500, "dosage:chrM"), ("chrEBV", 102_000, 122_000, "dosage:chrEBV")]
L_MAX = 600   # the longest fragment-GC window the engine tabulates


def read_exclusions(path, chroms):
    iv = {c: [] for c in chroms}
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            p = line.split("\t")
            if p[0] in iv:
                iv[p[0]].append((int(p[1]), int(p[2])))
    for c in iv:
        iv[c].sort()
    return iv


def complement_runs(excl, length, extra=()):
    merged = []
    for s, e in sorted(list(excl) + list(extra)):
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    runs, pos = [], 0
    for s, e in merged:
        if s > pos:
            runs.append((pos, s))
        pos = max(pos, e)
    if pos < length:
        runs.append((pos, length))
    return runs


def gc_hist(seq, L):
    a = np.frombuffer(seq.encode(), dtype=np.uint8)
    gc = np.r_[0, np.cumsum((a == 71) | (a == 67))]
    n = np.r_[0, np.cumsum(~((a == 65) | (a == 67) | (a == 71) | (a == 84)))]
    m = len(a) - L
    if m <= 0:
        return np.zeros(101, np.int64), True
    w, wn = gc[L:L + m] - gc[:m], n[L:L + m] - n[:m]
    return np.bincount((w[wn == 0] * 100 + L // 2) // L, minlength=101), bool(wn.any())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--exclude", required=True, help="NGS-PCA exclusion BED(.gz)")
    ap.add_argument("--out", required=True)
    ap.add_argument("-L", type=int, default=450)
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--min-run", type=int, default=8000)
    ap.add_argument("--max-run", type=int, default=20000)
    ap.add_argument("--per-bin", type=int, default=60000)
    ap.add_argument("--n-gc", type=int, default=500)
    ap.add_argument("--n-control", type=int, default=800)
    ap.add_argument("--n-test-auto", type=int, default=80)
    ap.add_argument("--n-test-x", type=int, default=60)
    ap.add_argument("--n-test-y", type=int, default=40)
    a = ap.parse_args()

    fa = FastaIndex(a.reference)
    autos = [f"chr{i}" for i in range(1, 23)]
    excl = read_exclusions(a.exclude, autos + ["chrX", "chrY"])
    cand, hist = [], []
    for c in autos + ["chrX", "chrY"]:
        extra = {"chrX": X_EXCLUDE_GRCH38, "chrY": Y_EXCLUDE_GRCH38}.get(c, ())
        for s, e in complement_runs(excl[c], fa.length(c), extra):
            if e - s < a.min_run:
                continue
            if e - s > a.max_run:
                mid = (s + e) // 2
                s, e = mid - a.max_run // 2, mid + a.max_run // 2
            if s < a.flank + 1 or e + a.flank + 1 > fa.length(c):
                continue
            h, has_n = gc_hist(fa.fetch(c, s, e + a.L), a.L)
            if has_n:
                continue
            cand.append((c, s, e))
            hist.append(h)
    H = np.array(hist)
    is_x = np.array([c == "chrX" for c, _, _ in cand])
    is_y = np.array([c == "chrY" for c, _, _ in cand])
    auto_idx = np.where(~is_x & ~is_y)[0]
    print(f"candidate runs: {len(auto_idx)} autosomal, {is_x.sum()} chrX, {is_y.sum()} chrY", file=sys.stderr)

    chosen = np.zeros(len(cand), bool)
    have = np.zeros(101, np.int64)
    tot = H[auto_idx].sum(0)
    for g in np.argsort(tot, kind="stable"):
        if tot[g] == 0:
            continue
        while have[g] < a.per_bin and chosen.sum() < a.n_gc:
            pool = auto_idx[~chosen[auto_idx] & (H[auto_idx, g] > 0)]
            if len(pool) == 0:
                break
            j = pool[np.argmax(H[pool, g])]
            chosen[j] = True
            have += H[j]
    rest = auto_idx[~chosen[auto_idx]]
    need = max(0, a.n_control - int(chosen.sum()))
    chosen[rest[np.linspace(0, len(rest) - 1, need).astype(int)]] = True
    role = np.array(["" for _ in cand], dtype=object)
    role[chosen] = "control"
    rest = auto_idx[~chosen[auto_idx]]
    role[rest[np.linspace(0, len(rest) - 1, a.n_test_auto).astype(int)]] = "test:auto"
    xs = np.where(is_x)[0]
    role[xs[np.linspace(0, len(xs) - 1, min(a.n_test_x, len(xs))).astype(int)]] = "test:chrX"
    ys = np.where(is_y)[0]
    role[ys[np.linspace(0, len(ys) - 1, min(a.n_test_y, len(ys))).astype(int)]] = "test:chrY"

    with open(a.out, "w") as fh:
        for (c, s, e), r in zip(cand, role):
            if r:
                fh.write(f"{c}\t{s}\t{e}\t{r}\n")
        for c, s, e, r in DOSAGE_GRCH38:
            if c in fa.index and e + a.flank <= fa.length(c):
                # every window of every position, on both strands, must be free of N
                if "N" in fa.fetch(c, s - L_MAX + 1, e + L_MAX - 1).upper():
                    sys.exit(f"{c}:{s}-{e}: an N within {L_MAX} bp - some fragment-GC windows would be undefined")
                fh.write(f"{c}\t{s}\t{e}\t{r}\n")
    sel = H[role == "control"].sum(0)
    print(f"controls: {(role == 'control').sum()} regions, {sum(e - s for (c, s, e), r in zip(cand, role) if r == 'control'):,} bp; "
          f"test:auto {(role == 'test:auto').sum()}, test:chrX {(role == 'test:chrX').sum()}, test:chrY {(role == 'test:chrY').sum()}", file=sys.stderr)
    print("fragment-GC stratum: control position-strands (forward)", file=sys.stderr)
    for g in range(15, 90, 5):
        print(f"  {g:3d}-{g + 4:3d}%  {sel[g:g + 5].sum():>12,}", file=sys.stderr)


if __name__ == "__main__":
    main()
