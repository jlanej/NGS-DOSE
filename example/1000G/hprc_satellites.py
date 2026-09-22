#!/usr/bin/env python3
"""Satellite array mass: NGS-DOSE estimates from short reads against HPRC release-2 assemblies.

The HPRC CenSat annotation of a diploid assembly gives, per haplotype, the span of every
satellite array. Summed over both haplotypes that is the sample's array mass per class - an
assembly-based truth for the compositional classes of the experimental satellite panel, for the
~200 HPRC samples that are also in the 1000 Genomes 30x cohort.

    hprc_satellites.py --estimates estimates.tsv --censat DIR --out hprc_satellites.tsv

Caveats worth keeping next to the numbers: an assembly is not a truth for an array it failed to
span - such arrays are annotated together with their gap ("GAP,HSat2"), are tallied separately,
and a sample that has any in a class is left out of that class's comparison; and the panel's
k-mers come from CHM13, so a class whose sequence differs between people, or too few of whose
reads carry the four k-mers a read needs, is under-recovered in proportion (the recall of each
class on CHM13 itself is in resources/experimental/README.md).
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

# the parser and the comparison live in the package; this script is their command line
from ngsdose.hprc import CLASS_OF, CLASSES, MAX_GAPPED, assembly_mass, compare, labels_of  # noqa: E402,F401


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--estimates", required=True, help="table from `ngsdose estimate` or `ngsdose cohort` on scan-mode counts")
    ap.add_argument("--censat", required=True, help="directory of <sample>_<hap>_hprc_r2_v1*.cenSat.bed files")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    with open(a.estimates) as fh:
        rows, stats = compare(list(csv.DictReader(fh, delimiter="\t")), a.censat)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "cls", "assembly_Mb", "assembly_gapped_Mb", "ngsdose_Mb"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n = len({r["sample"] for r in rows})
    print(f"{n} samples with both haplotype annotations and satellite estimates -> {a.out}")
    print(f"An assembly is a truth only for the arrays it spans: a sample with more than {100 * MAX_GAPPED:.0f}% of a class's annotated sequence in "
          "gap-containing arrays is left out of that class's comparison (column 'gapped' counts them).")
    print(f"{'class':8s} {'n':>4s} {'gapped':>7s} {'assembly Mb (median, range)':>30s} {'est/assembly (median)':>22s} {'SD of log ratio':>16s} {'Pearson r':>10s} {'Spearman':>9s}")
    for cls in CLASSES:
        st = stats[cls]
        if st["n"] < 2:
            print(f"{cls:8s} {st['n']:4d} {st['n_gapped']:7d}   too few samples without gaps")
            continue
        corr = f"{st['pearson']:10.3f} {st['spearman']:9.3f}" if "pearson" in st else f"{'':>10s} {'':>9s}"
        print(f"{cls:8s} {st['n']:4d} {st['n_gapped']:7d} {st['assembly_median']:12.1f} ({st['assembly_min']:.1f}-{st['assembly_max']:.1f}){'':6s} "
              f"{st['ratio_median']:22.2f} {st['sd_log']:16.3f} {corr}")


if __name__ == "__main__":
    main()
