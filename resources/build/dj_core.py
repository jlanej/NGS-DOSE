#!/usr/bin/env python3
"""Core k-mers of the distal junction: present exactly once in each of the five CHM13 acrocentric
distal junctions and nowhere else in CHM13.

Inputs are two `ngs-dose panel --report` tables for the chr21 DJ used as the unit:
  report_a: background CHM13 with only the chr21 DJ masked  -> occurrences elsewhere, other DJs included
  report_b: background CHM13 with all five DJs masked       -> occurrences outside any DJ
A unit k-mer is core if it occurs once in the unit, 4 times in (a) and 0 times in (b).
Output: BED of k-mer start positions (unit coordinates), merged into intervals.
"""
import gzip
import sys


def load(path):
    d = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            _cls, pos, n, _shared, bg = line.split("\t")
            d[int(pos)] = (int(n), int(bg))
    return d


def main():
    report_a, report_b, name, out = sys.argv[1:5]
    a, b = load(report_a), load(report_b)
    core = sorted(p for p, (n, bga) in a.items() if n == 1 and bga == 4 and b[p][1] == 0)
    runs, s = [], None
    for p in core:
        if s is None:
            s = e = p
        elif p == e + 1:
            e = p
        else:
            runs.append((s, e + 1))
            s = e = p
    runs.append((s, e + 1))
    with open(out, "w") as fh:
        for s0, e0 in runs:
            fh.write(f"{name}\t{s0}\t{e0}\n")
    print(f"{len(a)} distinct unit k-mers, {len(core)} core, {len(runs)} intervals -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
