#!/usr/bin/env python3
"""One-mismatch background filter, in two steps around the engine's panel builder.

  write:  neighbors.py write PANEL.tsv.gz OUTDIR
          every 1-substitution neighbour of every panel k-mer, packed into positional "classes"
          (one FASTA record each, neighbours separated by N) plus a manifest, so that
          `ngs-dose panel --report` counts each neighbour's occurrences in the backgrounds
          (outside the masks). Neighbours that are themselves panel k-mers are left out.
  filter: neighbors.py filter PANEL.tsv.gz OUTDIR REPORT.tsv.gz OUT.panel.tsv.gz SUMMARY.tsv
          drop every panel k-mer with a neighbour seen in a background.

A read from a non-class locus that differs from a class k-mer by one base (a paralog, or the
sample's own allele at a near-identical copy) would otherwise be counted: in the NA12878 scan,
0.1-0.4% of each class's reads came from such loci (e.g. an LTR12C at chr6:148.70 Mb for SNORD115).
"""
import gzip
import sys
from pathlib import Path

COMP = str.maketrans("ACGT", "TGCA")
SPAN = 32                      # 31-mer + one N
PER_CLASS = (1 << 23) // SPAN - 1


def canonical(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def load(path):
    head, classes, rows = [], [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##class\t"):
                classes.append(dict(x.split("=", 1) for x in line.split("\t")[1:]))
            elif line.startswith("#"):
                head.append(line)
            elif line:
                rows.append(line.split("\t"))
    return head, classes, rows


def neighbours(k):
    for i, b in enumerate(k):
        for x in "ACGT":
            if x != b:
                yield k[:i] + x + k[i + 1:]


def write(panel, outdir):
    _, _, rows = load(panel)
    own = {r[0] for r in rows}
    seen, order = set(), []
    for r in rows:
        for n in neighbours(r[0]):
            c = canonical(n)
            if c not in own and c not in seen:
                seen.add(c)
                order.append(c)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "manifest.tsv", "w") as man, open(out / "index.tsv", "w") as idx:
        for ci in range(0, len(order), PER_CLASS):
            name = f"nb{ci // PER_CLASS}"
            chunk = order[ci:ci + PER_CLASS]
            with open(out / f"{name}.fa", "w") as fa:
                fa.write(f">{name}\n" + "N".join(chunk) + "N\n")
            man.write(f"{name}\tpositional\t{name}.fa\t0\n")
            idx.write(f"{name}\t{ci}\t{len(chunk)}\n")
        with open(out / "neighbours.txt", "w") as fh:
            fh.write("\n".join(order) + "\n")
    print(f"[neighbors] {len(rows)} panel k-mers, {len(order)} distinct neighbours in {-(-len(order) // PER_CLASS)} classes", file=sys.stderr)


def filt(panel, outdir, report, outpanel, summary):
    head, classes, rows = load(panel)
    out = Path(outdir)
    order = open(out / "neighbours.txt").read().split()
    start = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(out / "index.tsv")}
    hit = set()
    with gzip.open(report, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c, pos, _cnt, _sh, bg = line.rstrip("\n").split("\t")
            if int(bg) > 0:
                p = int(pos)
                if p % SPAN == 0:   # k-mers spanning an N do not exist; pos is the neighbour's start
                    hit.add(order[start[c] + p // SPAN])
    # a neighbour shared by several neighbour records is reported once, under its first record; the
    # set `hit` holds the neighbour itself, so every panel k-mer that generates it is dropped below
    drop = {r[0] for r in rows if any(canonical(n) in hit for n in neighbours(r[0]))}
    n_before = [0] * len(classes)
    n_drop = [0] * len(classes)
    kept = []
    for r in rows:
        ci = int(r[1])
        n_before[ci] += 1
        if r[0] in drop:
            n_drop[ci] += 1
        else:
            kept.append(r)
    with gzip.open(outpanel, "wt", compresslevel=6) as w:
        for line in head:
            if not line.startswith("#kmer"):
                w.write(line + "\n")
        w.write("##filtered=neighbors.py: removed k-mers with a one-substitution neighbour anywhere in the backgrounds outside the class loci\n")
        for i, c in enumerate(classes):
            c = dict(c)
            c["kmers_kept"] = str(n_before[i] - n_drop[i])
            w.write("##class\t" + "\t".join(f"{key}={c[key]}" for key in ("id", "name", "kind", "length", "circular", "source", "kmers_input", "kmers_kept")) + "\n")
        w.write("#kmer\tclass\tpos\tstrand\n")
        for r in kept:
            w.write("\t".join(r) + "\n")
    with open(summary, "w") as s:
        s.write("class\tbefore\tremoved_one_mismatch_background\tfinal\n")
        for i, c in enumerate(classes):
            s.write(f"{c['name']}\t{n_before[i]}\t{n_drop[i]}\t{n_before[i] - n_drop[i]}\n")
    print(open(summary).read(), file=sys.stderr, end="")


if __name__ == "__main__":
    if sys.argv[1] == "write":
        write(*sys.argv[2:4])
    else:
        filt(*sys.argv[2:7])
