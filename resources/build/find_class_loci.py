#!/usr/bin/env python3
"""Loci of a genome that are homologous to a class unit, from a minimap2 PAF of unit chunks.

    chunks:  python find_class_loci.py chunk UNIT.fa 2000 1000 > chunks.fa
    map:     minimap2 -c -x asm20 -N 50 -p 0.3 --secondary=yes GENOME.fa chunks.fa > chunks.paf
    loci:    python find_class_loci.py loci chunks.paf CLASS --min-len 500 --min-id 0.90 --pad 2000 > loci.bed

The loci are exempted from the background when the panel is built (they ARE the class), and are
where an aligner can be expected to place the class's reads.
"""
import argparse
import sys


def chunk(args):
    seq = "".join(line.strip() for line in open(args.fasta) if not line.startswith(">"))
    for i in range(0, len(seq) - args.size + 1, args.step):
        print(f">c_{i}\n{seq[i:i + args.size]}")


def loci(args):
    hits = []
    for line in open(args.paf):
        p = line.split("\t")
        alen, ident = int(p[3]) - int(p[2]), int(p[9]) / int(p[10])
        if alen >= args.min_len and ident >= args.min_id:
            hits.append((p[5], max(0, int(p[7]) - args.pad), min(int(p[6]), int(p[8]) + args.pad)))
    hits.sort()
    merged = []
    for c, s, e in hits:
        if merged and merged[-1][0] == c and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([c, s, e])
    for c, s, e in merged:
        print(f"{c}\t{s}\t{e}\t{args.name}")
    print(f"{len(merged)} loci, {sum(e - s for _, s, e in merged):,} bp", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("chunk")
    c.add_argument("fasta"); c.add_argument("size", type=int); c.add_argument("step", type=int)
    c.set_defaults(fn=chunk)
    l = sub.add_parser("loci")
    l.add_argument("paf"); l.add_argument("name")
    l.add_argument("--min-len", type=int, default=500); l.add_argument("--min-id", type=float, default=0.90)
    l.add_argument("--pad", type=int, default=2000)
    l.set_defaults(fn=loci)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
