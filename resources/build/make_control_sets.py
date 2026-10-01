#!/usr/bin/env python3
"""Write the bundle's control sets: which of its single-copy regions a fetch reads.

The bundle's controls.bed holds every region: the 800 controls (the denominator and the GC curve), the known-truth
regions (held-out autosomal, chrX, chrY), the dosage regions (chrM, chrEBV) and the karyotype windows. A whole-file
scan counts them all at no cost. A fetch pays for each in bytes, so it reads a set:

  all        controls.bed itself
  base       the controls, the truths and the dosage regions: what the bundle held before the karyotype windows,
             and what the fetch menu's controls row names (a fetch that does not say otherwise reads as it always did)
  lite200    200 of the controls, the truths and the dosage regions (controls.lite200.bed, kept as it is)
  karyotype  lite200's regions and every karyotype window: every chromosome read at full precision
  screen     lite200's controls, the dosage regions and the windows of tier 1: the cheapest set that reads every
             chromosome and both sex chromosomes

Each set is a BED in the bundle's format, its regions in controls.bed's order. Only all and lite200 ship a FASTA;
`ngsdose fetchplan --controls NAME` cuts the others' from controls.fa.gz when it writes a plan.

usage: make_control_sets.py BUNDLE_DIR KARYOTYPE_TABLE   (the table select_karyotype.py writes: chrom start end window arm tier)
"""
import sys
from pathlib import Path


def rows(path):
    return [tuple(line.rstrip("\n").split("\t")) for line in open(path) if line.strip() and not line.startswith("#")]


def main():
    bundle, table = Path(sys.argv[1]), sys.argv[2]
    every = rows(bundle / "controls.bed")
    lite = {r[:3] for r in rows(bundle / "controls.lite200.bed")}
    tier = {tuple(r[:3]): int(r[5]) for r in rows(table) if r[0] != "chrom"}
    window = lambda r: r[3].startswith("test:karyotype")
    missing = [r for r in every if window(r) and r[:3] not in tier]
    if missing or not all(k in {r[:3] for r in every} for k in tier):
        sys.exit(f"{table} and {bundle / 'controls.bed'} do not hold the same karyotype windows")
    sets = {
        "base": [r for r in every if not window(r)],
        "karyotype": [r for r in every if r[:3] in lite or window(r)],
        "screen": [r for r in every if (r[:3] in lite and (r[3] == "control" or r[3].startswith("dosage:"))) or (window(r) and tier[r[:3]] == 1)],
    }
    for name, rs in sets.items():
        with open(bundle / f"controls.{name}.bed", "w") as fh:
            fh.write("".join("\t".join(r) + "\n" for r in rs))
        bp = sum(int(r[2]) - int(r[1]) for r in rs)
        print(f"controls.{name}.bed: {len(rs)} regions, {bp:,} bp ({sum(1 for r in rs if r[3] == 'control')} controls, {sum(1 for r in rs if window(r))} window pieces)", file=sys.stderr)


if __name__ == "__main__":
    main()
