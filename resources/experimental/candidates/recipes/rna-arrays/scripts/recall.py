#!/usr/bin/env python3
"""Read recall of a positional class, per source copy: 150-bp reads tiled every 10 bp over each
record of a copies FASTA; a read counts if it carries >= 4 of the class's own panel k-mers (the
engine judges a positional class on its own hits only, src/count.rs classify). Also reports the
k-mers per read and, for any other class of the panel, the share of reads that would be assigned
to it (cross-talk).

    recall.py PANEL.tsv.gz CLASS=COPIES.fa [CLASS=COPIES.fa ...] [--step 10] [--read-length 150]
"""
import argparse
import gzip
import statistics

COMP = str.maketrans("ACGT", "TGCA")


def canonical(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name is not None:
            out[name].append(line.strip().upper())
    return {n: "".join(s) for n, s in out.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("panel")
    ap.add_argument("copies", nargs="+", metavar="CLASS=FASTA")
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--read-length", type=int, default=150)
    ap.add_argument("--min-hits", type=int, default=4)
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
    L = a.read_length
    print("class\tcopy\tbp\treads\trecall\tmedian_kmers_per_read\tassigned_to_other_classes")
    for spec in a.copies:
        cls, fa = spec.split("=", 1)
        for name, s in read_fasta(fa).items():
            n = ok = 0
            hits_all, other = [], 0
            for i in range(0, len(s) - L + 1, a.step):
                r = s[i:i + L]
                cnt = {}
                for j in range(L - k + 1):
                    c = kmers.get(canonical(r[j:j + k]))
                    if c:
                        cnt[c] = cnt.get(c, 0) + 1
                h = cnt.get(cls, 0)
                n += 1
                ok += h >= a.min_hits
                hits_all.append(h)
                other += any(v >= a.min_hits for c, v in cnt.items() if c != cls)
            print(f"{cls}\t{name}\t{len(s)}\t{n}\t{ok / max(n, 1):.3f}\t{statistics.median(hits_all) if hits_all else 0}\t{other}")


if __name__ == "__main__":
    main()
