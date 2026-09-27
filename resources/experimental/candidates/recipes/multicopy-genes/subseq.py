#!/usr/bin/env python3
"""Cut copies by where a sub-unit aligns in them: `inside` writes the aligned stretches (e.g. the
HERV-K(C4) insertion of each long C4 copy), `outside` the rest of each copy, as separate records
(e.g. the C4 gene without its insertion). Stretches shorter than --min are left out.

    subseq.py inside|outside SUBUNIT.fa COPIES.fa [--min 1000] > out.fa
"""
import subprocess
import sys


def read_fasta(path):
    recs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name:
                recs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name:
        recs[name] = "".join(buf)
    return recs


def main():
    mode, sub, copies = sys.argv[1:4]
    mn = int(sys.argv[sys.argv.index("--min") + 1]) if "--min" in sys.argv else 1000
    recs = read_fasta(copies)
    paf = subprocess.run(["minimap2", "-x", "asm20", "-c", "-P", sub, copies], capture_output=True, text=True, check=True).stdout
    hits = {}
    for line in paf.splitlines():
        p = line.split("\t")
        if int(p[10]) >= mn:
            hits.setdefault(p[0], []).append((int(p[2]), int(p[3])))
    for name, seq in recs.items():
        iv = sorted(hits.get(name, []))
        if mode == "inside":
            parts = iv
        else:
            parts, pos = [], 0
            for s, e in iv:
                parts.append((pos, s))
                pos = max(pos, e)
            parts.append((pos, len(seq)))
        for s, e in parts:
            if e - s >= 150:
                print(f">{name}:{s}-{e}\n{seq[s:e]}")


if __name__ == "__main__":
    main()
