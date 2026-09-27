#!/usr/bin/env python3
"""Positions of a unit where a paralogous copy differs from it: align the copy to the unit with
minimap2 (-x asm5 -c --cs) and list, in unit coordinates (0-based), every mismatched base, every
unit base the copy lacks, and the two unit bases around a point where the copy has extra bases.
A 31-mer of the unit that covers none of these positions is also in the copy (where the
alignment runs); one that covers any of them is not.

    diffsites.py UNIT.fa COPY.fa [COPY2.fa ...] > sites.txt
"""
import re
import subprocess
import sys


def sites(unit, copy):
    out = subprocess.run(["minimap2", "-x", "asm5", "-c", "--cs", unit, copy], capture_output=True, text=True, check=True).stdout
    s = set()
    for line in out.splitlines():
        p = line.split("\t")
        cs = next(x[5:] for x in p[12:] if x.startswith("cs:Z:"))
        t = int(p[7])
        tend = int(p[8])
        for op in re.findall(r"(:\d+|\*[a-z]{2}|\+[a-z]+|-[a-z]+)", cs):
            if op[0] == ":":
                t += int(op[1:])
            elif op[0] == "*":
                s.add(t)
                t += 1
            elif op[0] == "-":
                s.update(range(t, t + len(op) - 1))
                t += len(op) - 1
            else:  # insertion in the copy between unit bases t-1 and t
                s.update((t - 1, t))
        assert t == tend, (t, tend)
    return s


def main():
    unit, copies = sys.argv[1], sys.argv[2:]
    allsites = set()
    for c in copies:
        allsites |= sites(unit, c)
    for x in sorted(allsites):
        print(x)


if __name__ == "__main__":
    main()
