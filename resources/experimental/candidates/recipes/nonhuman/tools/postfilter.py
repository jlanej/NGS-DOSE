#!/usr/bin/env python3
"""Drop from a freshly built panel every k-mer that a shipped panel already claims, and every
low-complexity k-mer.

The engine drops a k-mer held by two loaded panels from both classes, so a candidate sharing
k-mers with the shipped panels would silently change the counts of the cohort that is running.
Removing them here means loading this panel cannot alter a shipped class.

Low complexity is dropped because the unmapped bin is where two-colour poly-G reads end up: a
k-mer of at most two distinct bases, or one holding a homopolymer of 18 bp or more, matches
those reads and not the organism.

    postfilter.py IN.tsv.gz OUT.tsv.gz REPORT.tsv SHIPPED.tsv.gz [SHIPPED.tsv.gz ...]
"""
import gzip
import sys
from collections import Counter

inp, outp, repp, *shipped = sys.argv[1:]

claimed = set()
for path in shipped:
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                claimed.add(line.split("\t", 1)[0])


def low_complexity(s: str) -> bool:
    if len(set(s)) <= 2:
        return True
    run = best = 1
    for a, b in zip(s, s[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best >= 18


names, header, rows = {}, [], []
drop_shared, drop_lc = Counter(), Counter()
with gzip.open(inp, "rt") as fh:
    for line in fh:
        if line.startswith("##class"):
            f = dict(x.split("=", 1) for x in line.rstrip("\n").split("\t")[1:])
            names[f["id"]] = f["name"]
            header.append(line)
        elif line.startswith("#"):
            header.append(line)
        elif line.strip():
            kmer, cid = line.split("\t")[0], line.split("\t")[1]
            cls = names[cid]
            if kmer in claimed:
                drop_shared[cls] += 1
            elif low_complexity(kmer):
                drop_lc[cls] += 1
            else:
                rows.append(line)

kept = Counter(names[r.split("\t")[1]] for r in rows)
with gzip.open(outp, "wt") as out:
    for line in header:
        if line.startswith("##class"):
            f = dict(x.split("=", 1) for x in line.rstrip("\n").split("\t")[1:])
            f["kmers_kept"] = str(kept[f["name"]])
            line = "##class\t" + "\t".join(f"{a}={b}" for a, b in f.items()) + "\n"
        out.write(line)
    out.writelines(rows)

with open(repp, "w") as fh:
    fh.write("#class\tkmers_before\tshared_with_shipped\tlow_complexity\tkmers_kept\n")
    for cid in sorted(names, key=int):
        c = names[cid]
        fh.write(f"{c}\t{kept[c] + drop_shared[c] + drop_lc[c]}\t{drop_shared[c]}\t{drop_lc[c]}\t{kept[c]}\n")
print(open(repp).read(), end="")
