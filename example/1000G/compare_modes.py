#!/usr/bin/env python3
"""What does the fast mode cost? Per-sample estimates from a targeted fetch against those from
the whole-file scan of the same CRAM.

    compare_modes.py --scan single_sample.scan.tsv --fetch single_sample.fetch.tsv --out modes.tsv

Both tables come from `ngsdose estimate`. The known-truth and dosage columns are made from the
same reads in both modes and must be identical; the classes differ by whatever the sinks miss.
"""
import argparse
import csv

import numpy as np

COLUMNS = ("rDNA45S.cn_single", "rDNA45S.cn_all", "rDNA5S.cn_single", "DJ.cn_single", "truth.auto", "truth.chrX", "chrM.copies", "chrEBV.copies")


def table(path):
    with open(path) as fh:
        return {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan", required=True)
    ap.add_argument("--fetch", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    scan, fetch = table(a.scan), table(a.fetch)
    both = [s for s in scan if s in fetch]
    rows = []
    for s in both:
        row = dict(sample=s)
        for c in COLUMNS:
            try:
                x, y = float(scan[s][c]), float(fetch[s][c])
                row[c] = round(y / x, 6) if x > 0 else float("nan")
            except (KeyError, ValueError):
                row[c] = float("nan")
        rows.append(row)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", *COLUMNS], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(both)} samples counted in both modes; fetch / scan per sample -> {a.out}")
    print(f"{'column':22s} {'n':>5s} {'median':>9s} {'min':>9s} {'max':>9s} {'SD of log ratio':>16s}   lowest")
    for c in COLUMNS:
        v = np.array([r[c] for r in rows], float)
        ok = np.isfinite(v)
        if not ok.any():
            continue
        worst = sorted(((r[c], r["sample"]) for r in rows if np.isfinite(r[c])))[:3]
        sd = float(np.std(np.log(v[ok]), ddof=1)) if ok.sum() > 1 else 0.0
        print(f"{c:22s} {int(ok.sum()):5d} {np.median(v[ok]):9.5f} {v[ok].min():9.5f} {v[ok].max():9.5f} {sd:16.6f}   "
              + ", ".join(f"{s} {x:.4f}" for x, s in worst))


if __name__ == "__main__":
    main()
