#!/usr/bin/env python3
"""How much of a compositional class can its panel see? For each class of a panel, the share of
read-length windows drawn from the class's own source arrays that carry at least `--min-hits` of
the class's k-mers - the engine's rule for assigning a read. A class at 99% is measured; a class
at 45% is a relative measure at best, and its mass is under-read by about that factor.

    panel_recall.py PANEL.tsv.gz CLASS=arrays.fa [CLASS=arrays.fa ...]

(`build_satellite_panel.sh` leaves the arrays of every class in $WORK/satellite_panel/<class>.fa.)
"""
import argparse
import gzip
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from ngsdose.io import read_fasta  # noqa: E402

COMP = str.maketrans("ACGT", "TGCA")


def canonical(s: str) -> str:
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("panel")
    ap.add_argument("arrays", nargs="+", metavar="CLASS=FASTA")
    ap.add_argument("--read-length", type=int, default=150)
    ap.add_argument("--min-hits", type=int, default=4)
    ap.add_argument("-n", type=int, default=3000, help="windows sampled per class")
    a = ap.parse_args()
    names, kmers, k = {}, {}, 31
    with gzip.open(a.panel, "rt") as fh:
        for line in fh:
            if line.startswith("##k="):
                k = int(line[4:])
            elif line.startswith("##class"):
                f = dict(x.split("=", 1) for x in line.rstrip("\n").split("\t")[1:])
                names[f["id"]] = f["name"]
            elif not line.startswith("#"):
                p = line.split("\t")
                kmers[p[0]] = names[p[1]]
    rng = random.Random(1)
    print("class\tMb\tpanel_kmers\tread_recall\tmedian_hits_per_read")
    for spec in a.arrays:
        cls, fasta = spec.split("=", 1)
        seqs = [s.upper() for s in read_fasta(fasta).values() if len(s) >= a.read_length]
        weights = [len(s) for s in seqs]
        hits = []
        for _ in range(a.n):
            s = rng.choices(seqs, weights=weights)[0]
            i = rng.randrange(len(s) - a.read_length + 1)
            w = s[i:i + a.read_length]
            hits.append(sum(1 for j in range(a.read_length - k + 1) if kmers.get(canonical(w[j:j + k])) == cls))
        hits.sort()
        print(f"{cls}\t{sum(weights) / 1e6:.2f}\t{sum(1 for v in kmers.values() if v == cls)}\t{sum(h >= a.min_hits for h in hits) / len(hits):.3f}\t{hits[len(hits) // 2]}")


if __name__ == "__main__":
    main()
