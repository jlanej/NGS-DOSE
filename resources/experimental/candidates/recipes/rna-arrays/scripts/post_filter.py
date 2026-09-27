#!/usr/bin/env python3
"""Make a candidate panel safe to load beside the shipped panels, and specific per class.

1. Remove every k-mer that any shipped panel carries. The engine drops a k-mer claimed by two
   loaded panels from BOTH, so one shared k-mer would change a shipped class's counts.
2. The group is built with one mask per assembly (the union of every class's own loci), so a
   class's k-mer could occur at ANOTHER class's masked locus without counting as background.
   Remove, for each class, the k-mers that occur in any other class's own-locus sequence.

    post_filter.py IN.panel.tsv.gz OUT.panel.tsv.gz --shipped P.tsv.gz [P ...] \
        --loci CLASS=FASTA [CLASS=FASTA ...] --summary removed.tsv
"""
import argparse
import gzip
import sys

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


def panel_kmers(path):
    ks = set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                ks.add(line.split("\t", 1)[0])
    return ks


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--shipped", nargs="+", required=True)
    ap.add_argument("--loci", nargs="+", default=[], metavar="CLASS=FASTA")
    ap.add_argument("--summary", required=True)
    a = ap.parse_args()
    header, classes, rows = [], [], []
    with gzip.open(a.inp, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##k="):
                k = int(line[4:])
                header.append(line)
            elif line.startswith("##class\t"):
                f = dict(x.split("=", 1) for x in line.split("\t")[1:])
                classes.append(f)
            elif line.startswith("#"):
                header.append(line)
            elif line:
                rows.append(line.split("\t"))
    names = [c["name"] for c in classes]
    shipped = set()
    for p in a.shipped:
        shipped |= panel_kmers(p)
    # k-mers of each class's own loci
    loci = {}
    for spec in a.loci:
        cls, fa = spec.split("=", 1)
        s = loci.setdefault(cls, set())
        for seq in read_fasta(fa).values():
            for i in range(len(seq) - k + 1):
                w = seq[i:i + k]
                if "N" not in w:
                    s.add(canonical(w))
    removed_shipped = [0] * len(classes)
    removed_other = [0] * len(classes)
    kept_rows = []
    for r in rows:
        c = int(r[1])
        if r[0] in shipped:
            removed_shipped[c] += 1
            continue
        if any(r[0] in ks for cls, ks in loci.items() if cls != names[c]):
            removed_other[c] += 1
            continue
        kept_rows.append(r)
    kept = [0] * len(classes)
    for r in kept_rows:
        kept[int(r[1])] += 1
    with gzip.open(a.out, "wt", compresslevel=6) as w:
        for line in header:
            if line.startswith("#kmer"):
                continue
            w.write(line + "\n")
        w.write(f"##filtered=post_filter.py: removed k-mers of the shipped panels ({','.join(p.split('/')[-1] for p in a.shipped)}) and k-mers found at another class's own loci\n")
        for i, c in enumerate(classes):
            c = dict(c)
            c["kmers_kept"] = str(kept[i])
            w.write("##class\t" + "\t".join(f"{key}={c[key]}" for key in ("id", "name", "kind", "length", "circular", "source", "kmers_input", "kmers_kept")) + "\n")
        w.write("#kmer\tclass\tpos\tstrand\n")
        for r in kept_rows:
            w.write("\t".join(r) + "\n")
    with open(a.summary, "w") as s:
        s.write("class\tkept_by_builder\tremoved_shared_with_shipped\tremoved_at_other_class_loci\tfinal\n")
        for i, n in enumerate(names):
            nb = sum(1 for r in rows if int(r[1]) == i)
            s.write(f"{n}\t{nb}\t{removed_shipped[i]}\t{removed_other[i]}\t{kept[i]}\n")
    print(open(a.summary).read(), file=sys.stderr, end="")


if __name__ == "__main__":
    main()
