#!/usr/bin/env python3
"""The bundle's model of its single-copy regions (karyotype.model.json.gz): what a cohort too small to learn its own is
read against, and what a single genome is read against.

A model is learned on estimates (`ngsdose estimate`) of a reference cohort: each region's efficiency, its spread, the
libraries' shared modes, what a second X reads of the first, and how far a level's error is from its regions' noise
(ngsdose/karyotype.py). Directories are taken in order, and a later one's estimate of a sample replaces an earlier
one's: the cohort's estimates first, then those of the genomes counted with more regions (the karyotype windows).
A region that only part of the cohort holds is learned on that part; the modes are learned on the regions that at
least 95% of the genomes hold.

usage: karyotype_model.py -r resources/GRCh38 -o resources/GRCh38/karyotype.model.json.gz --source "..." ESTIMATES_DIR [ESTIMATES_DIR ...]
"""
import argparse
import datetime
import gzip
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from ngsdose import __version__, karyotype as K, resources
from ngsdose.tables import load_result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("estimates", nargs="+", help="directories of *.estimate.json.gz, in order of precedence (last wins)")
    ap.add_argument("-r", "--resources", required=True, help="the bundle directory (its bundle.json names the chromosomes' arms and GC)")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--source", required=True, help="one line on the cohort, kept in the model")
    ap.add_argument("-j", "--jobs", type=int, default=8)
    a = ap.parse_args()
    bundle = resources.Bundle(a.resources)
    kar = bundle.karyotype()
    if not kar:
        raise SystemExit(f"{a.resources}/bundle.json has no karyotype entry")
    known = {n for n, role in bundle.regions() if role != "dosage"}
    files: dict[str, Path] = {}
    per_dir = []
    for d in a.estimates:
        fs = sorted(Path(d).glob("*.estimate.json.gz"))
        if not fs:
            raise SystemExit(f"{d}: no *.estimate.json.gz")
        per_dir.append((d, len(fs)))
        for f in fs:
            files[f.name[:-len(".estimate.json.gz")]] = f

    def one(item):
        s, f = item
        v = K.gather(load_result(f), kar["control_names"])
        if v is None:
            return s, None
        keep = np.array([n in known for n in v[0]])             # a region the bundle no longer holds is not learned
        return s, ([n for n, k in zip(v[0], keep) if k], v[1][keep])
    with ThreadPoolExecutor(a.jobs) as ex:
        vec = list(ex.map(one, sorted(files.items())))
    samples = [s for s, v in vec if v is not None]
    vectors = [v for _, v in vec if v is not None]
    say = lambda m: print(m, file=sys.stderr)
    say(f"[model] {len(samples):,} genomes: " + "; ".join(f"{n:,} estimates in {d}" for d, n in per_dir))
    readings, model, info = K.cohort(vectors, kar["arms"], rules=kar["rules"], fit_own=True, log=say, gc=kar["gc"])
    names, Y = K.assemble(vectors)
    held = np.isfinite(Y).sum(0)
    by_name = dict(zip(names, held.tolist()))
    counts = np.array([by_name.get(n, 0) for n in model.names])
    layouts = sorted({int(c) for c in counts})
    model.info.update(source=a.source, built=str(datetime.date.today()), ngsdose=__version__,
                      genomes_per_region=dict(min=int(counts.min()), median=int(np.median(counts)), max=int(counts.max())),
                      regions_by_genomes={str(c): int((counts == c).sum()) for c in layouts} if len(layouts) <= 12 else None,
                      settled=info["n_settled"], fractional=info["n_fractional"])
    out = Path(a.out)
    with gzip.GzipFile(out, "wb", mtime=0) as fh:                # the same bytes for the same model
        fh.write(json.dumps(model.to_json(), separators=(",", ":")).encode())
    say(f"[model] {out}: {len(model.names):,} regions, {model.k} components, learned on {model.n:,} genomes ({out.stat().st_size / 1e6:.2f} MB); "
        f"regions held by {counts.min():,} to {counts.max():,} genomes")


if __name__ == "__main__":
    main()
