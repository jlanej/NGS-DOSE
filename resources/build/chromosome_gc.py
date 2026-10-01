#!/usr/bin/env python3
"""The GC fraction of every primary chromosome of a reference (A, C, G and T only; N left out): bundle.json's
`karyotype.gc`. In some libraries the GC-rich chromosomes read low or high together, in proportion to this
(ngsdose/karyotype.py), and each chromosome is set against what the others say.

usage: chromosome_gc.py REFERENCE.fa  > gc.json        (needs samtools and the reference's .fai)
"""
import json
import subprocess
import sys


def main(ref: str) -> dict:
    out = {}
    for c in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]:
        p = subprocess.Popen(["samtools", "faidx", ref, c], stdout=subprocess.PIPE)
        gc = at = 0
        for n, line in enumerate(p.stdout):
            if n == 0:
                continue
            u = line.upper()
            gc += u.count(b"G") + u.count(b"C")
            at += u.count(b"A") + u.count(b"T")
        if p.wait() != 0 or gc + at == 0:
            raise SystemExit(f"samtools faidx {ref} {c} failed")
        out[c] = round(gc / (gc + at), 4)
    return out


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    json.dump(main(sys.argv[1]), sys.stdout, indent=2)
    print()
