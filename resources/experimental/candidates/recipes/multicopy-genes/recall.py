#!/usr/bin/env python3
"""Read recall of a positional panel, per source copy: 150-bp reads tiled every 10 bp over each
copy (both strands are the same test, k-mers are canonical), and the share that carries at least
--min-hits k-mers of the class - the engine's rule for assigning a read. A copy is given as
CLASS=LABEL=FASTA (one or more records, each a copy); the output is one line per record.

    recall.py PANEL.tsv.gz CLASS=LABEL=copies.fa [...] [--step 10] [--read-length 150]

Also used for leakage: CLASS=LABEL=paralog.fa gives the share of a paralog's reads that the
class would claim (a false-positive rate).
"""
import argparse
import gzip
import sys

COMP = str.maketrans("ACGTN", "TGCAN")


def canonical(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def read_fasta(path):
    recs, name, buf = [], None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf).upper()))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf).upper()))
    return recs


def load_panel(path):
    names, kmers, k = {}, {}, 31
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##k="):
                k = int(line[4:])
            elif line.startswith("##class"):
                f = dict(x.split("=", 1) for x in line.rstrip("\n").split("\t")[1:])
                names[f["id"]] = f["name"]
            elif not line.startswith("#"):
                p = line.split("\t")
                kmers[p[0]] = names[p[1]]
    return k, kmers


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("panel")
    ap.add_argument("sources", nargs="+", metavar="CLASS=LABEL=FASTA")
    ap.add_argument("--read-length", type=int, default=150)
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--min-hits", type=int, default=4)
    a = ap.parse_args()
    k, kmers = load_panel(a.panel)
    R = a.read_length
    print("class\tlabel\tcopy\tbp\treads\trecall\tmedian_hits\tother_class_reads")
    for spec in a.sources:
        cls, label, fasta = spec.split("=", 2)
        for name, seq in read_fasta(fasta):
            if len(seq) < R:
                continue
            # per-position class of each k-mer (None if not in the panel)
            kc = [kmers.get(canonical(seq[j:j + k])) if "N" not in seq[j:j + k] else None for j in range(len(seq) - k + 1)]
            hits, other = [], 0
            for i in range(0, len(seq) - R + 1, a.step):
                win = kc[i:i + R - k + 1]
                h = sum(1 for c in win if c == cls)
                hits.append(h)
                counts = {}
                for c in win:
                    if c is not None and c != cls:
                        counts[c] = counts.get(c, 0) + 1
                if any(v >= a.min_hits for v in counts.values()):
                    other += 1
            hits.sort()
            n = len(hits)
            rec = sum(h >= a.min_hits for h in hits) / n
            print(f"{cls}\t{label}\t{name}\t{len(seq)}\t{n}\t{rec:.3f}\t{hits[n // 2]}\t{other}")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
