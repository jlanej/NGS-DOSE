#!/usr/bin/env python3
"""Read recall of a panel on given source copies: 150-bp reads tiled every 10 bp over each
source sequence, and the share carrying at least --min-hits k-mers of the named class (the
engine's rule for assigning a read).

    recall.py PANEL.tsv.gz CLASS=source.fa [CLASS=source.fa ...]

A source file's own records are reported one by one, so a class built from one strain can be
scored against the others.
"""
import argparse
import gzip

COMP = str.maketrans("ACGT", "TGCA")


def canon(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def read_fasta(path):
    name, seq, out = None, [], []
    for line in open(path):
        if line.startswith(">"):
            if name:
                out.append((name, "".join(seq).upper()))
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        out.append((name, "".join(seq).upper()))
    return out


ap = argparse.ArgumentParser()
ap.add_argument("panel")
ap.add_argument("sources", nargs="+", metavar="CLASS=FASTA")
ap.add_argument("--read-length", type=int, default=150)
ap.add_argument("--step", type=int, default=10)
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
        elif not line.startswith("#") and line.strip():
            p = line.split("\t")
            kmers[p[0]] = names[p[1]]

print("class\tsource\tbp\treads\trecall\tmedian_hits")
for spec in a.sources:
    cls, fasta = spec.split("=", 1)
    for name, s in read_fasta(fasta):
        if len(s) < a.read_length:
            continue
        hits = []
        for i in range(0, len(s) - a.read_length + 1, a.step):
            w = s[i:i + a.read_length]
            hits.append(sum(1 for j in range(a.read_length - k + 1) if kmers.get(canon(w[j:j + k])) == cls))
        hits.sort()
        rec = sum(h >= a.min_hits for h in hits) / len(hits)
        print(f"{cls}\t{name}\t{len(s)}\t{len(hits)}\t{rec:.3f}\t{hits[len(hits) // 2]}")
