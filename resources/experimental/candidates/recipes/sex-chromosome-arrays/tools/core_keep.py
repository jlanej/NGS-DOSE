#!/usr/bin/env python3
"""Keep BED for a positional class whose copies differ (the DJ recipe, resources/build/dj_core.py):
unit positions whose 31-mer occurs exactly once in the unit and exactly once in EVERY copy, so that
each copy contributes the same k-mers and the class count is proportional to the number of copies.

    core_keep.py UNIT.fa COPIES.fa [--k 31] > keep.bed
"""
import argparse
import collections
import sys

COMP = str.maketrans("ACGT", "TGCA")


def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name is not None:
            out[name].append(line.strip().upper())
    return {n: "".join(s) for n, s in out.items()}


def canonical(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("unit")
    ap.add_argument("copies")
    ap.add_argument("--k", type=int, default=31)
    a = ap.parse_args()
    k = a.k
    (uname, unit), = read_fasta(a.unit).items()
    copies = read_fasta(a.copies)
    counts = [collections.Counter(canonical(s[i:i + k]) for i in range(len(s) - k + 1)) for s in copies.values()]
    ucount = collections.Counter(canonical(unit[i:i + k]) for i in range(len(unit) - k + 1))
    keep = [i for i in range(len(unit) - k + 1)
            if "N" not in unit[i:i + k]
            and ucount[canonical(unit[i:i + k])] == 1
            and all(c.get(canonical(unit[i:i + k]), 0) == 1 for c in counts)]
    # merge consecutive positions into intervals
    iv = []
    for p in keep:
        if iv and iv[-1][1] == p:
            iv[-1][1] = p + 1
        else:
            iv.append([p, p + 1])
    for s, e in iv:
        print(f"{uname}\t{s}\t{e}")
    print(f"[core_keep] {len(keep)} of {len(unit) - k + 1} unit positions occur once in each of {len(copies)} copies", file=sys.stderr)


if __name__ == "__main__":
    main()
