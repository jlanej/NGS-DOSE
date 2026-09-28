#!/usr/bin/env python3
"""Check the estimate of the test fixture (tests/data/NA12878.subsample.bam, counted in fetch mode
with the GRCh38 bundle) against the copy numbers it is known to give, as tests/test_real_data.py
does. The container smoke tests run it on what the image counted and estimated.

    ngsdose estimate counts.json.gz -t - | python3 .github/scripts/smoke_check.py [--engine 0.1.0+abc1234]
"""
import argparse
import csv
import math
import sys

# column: (low, high); the limits of tests/test_real_data.py
LIMITS = {
    "truth.auto": (1.92, 2.08),                    # held-out autosomal sequence: 2 copies
    "truth.chrX": (1.8, 2.05),                     # NA12878 is female
    "rDNA45S.cn_single": (508 * 0.94, 508 * 1.06),  # 508 from the full 37x data
    "rDNA5S.cn_single": (150, 280),
    "DJ.cn_single": (9.2, 10.8),                   # one copy per acrocentric short arm
}

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--engine", help="the engine column the counts must carry (version+build)")
a = ap.parse_args()
rows = list(csv.DictReader(sys.stdin, delimiter="\t"))
if len(rows) != 1 or rows[0].get("sample") != "NA12878":
    sys.exit(f"expected one row for NA12878, got {len(rows)}")
row, bad = rows[0], []
for col, (lo, hi) in LIMITS.items():
    try:
        v = float(row.get(col, "nan"))
    except ValueError:
        v = math.nan
    ok = lo <= v <= hi
    print(f"{col:20s} {row.get(col, 'missing'):>10s}  [{lo:g}, {hi:g}]  {'ok' if ok else 'FAIL'}")
    if not ok:
        bad.append(col)
for col in (f"{c}.status" for c in ("rDNA45S", "rDNA5S", "DJ") if f"{c}.status" in row):
    if row[col] != "ok":
        print(f"{col:20s} {row[col]}  FAIL")
        bad.append(col)
if a.engine and row.get("engine") != a.engine:
    print(f"engine {row.get('engine')}, expected {a.engine}  FAIL")
    bad.append("engine")
if bad:
    sys.exit(f"the fixture's estimate is off in {', '.join(bad)}")
