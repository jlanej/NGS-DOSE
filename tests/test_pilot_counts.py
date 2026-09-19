"""Regression tests on the committed 1000 Genomes pilot counts (example/1000G/pilot): the modelling
layer on real data, with no engine and no alignment files needed."""
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from ngsdose import cohort, estimate, io, resources

ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "example" / "1000G" / "pilot"
BUNDLE = resources.Bundle(ROOT / "resources" / "GRCh38")
MALES = {"NA12891", "NA19239", "HG00512", "HG00731"}

nygc = sorted((PILOT / "counts_nygc").glob("*.json.gz"))
reps = sorted((PILOT / "counts_replicates").glob("*.json.gz"))
pytestmark = pytest.mark.skipif(len(nygc) < 6, reason="pilot counts not present")


def run(files):
    panel, units, feats, anchors = io.load_panel(BUNDLE.panel), BUNDLE.units(), BUNDLE.features(), BUNDLE.anchors()
    tabs, out = {}, {}
    for f in files:
        c = io.load_counts(f)
        L = estimate.nearest_table(c, None)["l"]
        if L not in tabs:
            tabs[L] = estimate.control_region_tables(BUNDLE.controls, L)
        out[c["sample"]] = estimate.estimate_sample(c, panel, units, feats, region_tables=tabs[L], anchors=anchors, regions=BUNDLE.regions())
    return out


@pytest.fixture(scope="module")
def ny():
    return run(nygc)


def test_known_truth_in_every_sample(ny):
    for s, r in ny.items():
        t = r["truth_regions"]
        assert abs(t["auto"]["cn"] - 2) < 0.04, (s, t["auto"]["cn"])
        x = t["chrX"]["cn"]
        # female cell lines can lose an X in part of the culture (HG00732 reads 1.61 in both libraries)
        assert (0.97 < x < 1.03) if s in MALES else (1.5 < x < 2.03), (s, x)
        assert 9.2 < r["classes"]["DJ"]["cn"] < 10.4, (s, r["classes"]["DJ"]["cn"])
        y = t["chrY"]["cn"]
        assert (0.95 < y < 1.03) if s in MALES else (y < 0.01), (s, y)


def test_both_libraries_were_counted_with_the_shipped_bundle_and_carry_the_culture_covariates(ny):
    """The committed counts are the bundle's own (same panel, controls and sinks hashes as the
    files in resources/), every input was complete, and mitochondrial and EBV dosage - different
    cultures of a line differ severalfold in both - are present for every sample."""
    import hashlib
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    want = ([sha(BUNDLE.panel)], sha(BUNDLE.controls), sha(BUNDLE.sinks))
    for f in nygc + reps:
        c = io.load_counts(f)
        assert (c["panel_sha256"], c["controls_sha256"], c["sinks_sha256"]) == want, f.name
        assert c["eof_marker"] == "present" and c["mode"] == "fetch"
    older = run(reps)
    for res in (ny, older):
        for s, r in res.items():
            t = r["truth_regions"]
            assert t["chrM"]["role"] == t["chrEBV"]["role"] == "dosage"
            assert 300 < t["chrM"]["cn"] < 3000 and 10 < t["chrEBV"]["cn"] < 1000, (s, t["chrM"]["cn"], t["chrEBV"]["cn"])
    # chrY is the same in both libraries of a man, and absent in both libraries of a woman
    for s in older:
        a, b = ny[s]["truth_regions"]["chrY"]["cn"], older[s]["truth_regions"]["chrY"]["cn"]
        assert abs(a - b) < 0.03, (s, a, b)


def test_rdna_is_plausible_and_28s_reads_low(ny):
    for s, r in ny.items():
        k = r["classes"]["rDNA45S"]
        assert 150 < k["cn"] < 900, (s, k["cn"])
        assert 0.6 < k["features"]["28S"]["cn"] / k["features"]["18S"]["cn"] < 0.95
        assert 60 < r["classes"]["rDNA5S"]["cn"] < 400


def test_calibration_is_tight_and_children_lie_within_parental_sum(ny):
    cal = cohort.calibrate(list(ny.values()), "rDNA45S", anchors=BUNDLE.anchors().get("rDNA45S"))
    cn = dict(zip(cal.samples, np.exp(cal.c)))
    assert np.all(cal.c_se < 0.01)                      # robust relative SE of each sample's estimate
    for child, father, mother in (("NA12878", "NA12891", "NA12892"), ("NA19240", "NA19239", "NA19238"),
                                  ("HG00514", "HG00512", "HG00513"), ("HG00733", "HG00731", "HG00732")):
        if all(x in cn for x in (child, father, mother)):
            assert cn[child] < cn[father] + cn[mother]


@pytest.mark.skipif(len(reps) < 6, reason="replicate counts not present")
def test_independent_libraries_agree_on_known_truth_and_rank(ny):
    hg = run(reps)
    common = [s for s in hg if s in ny]
    a = np.array([ny[s]["classes"]["DJ"]["cn"] for s in common]); b = np.array([hg[s]["classes"]["DJ"]["cn"] for s in common])
    assert np.abs(np.log(b / a)).max() < 0.08
    x = np.array([ny[s]["classes"]["rDNA45S"]["cn"] for s in common]); y = np.array([hg[s]["classes"]["rDNA45S"]["cn"] for s in common])
    lr = np.log(y / x)
    # single-sample estimates on the bundle's consensus anchors, three library chemistries, no cohort
    # calibration: no systematic offset, and pair-level disagreement of a few percent (which includes
    # real drift between cultures of the same cell line)
    assert abs(lr.mean()) < 0.03, lr.mean()
    assert lr.std(ddof=1) < 0.06, lr.std(ddof=1)
    assert np.corrcoef(x, y)[0, 1] > 0.95
    # the estimator used in the literature, for contrast: a large library offset
    fx = np.array([ny[s]["classes"]["rDNA45S"]["features"]["18S"]["cn_flat"] for s in common])
    fy = np.array([hg[s]["classes"]["rDNA45S"]["features"]["18S"]["cn_flat"] for s in common])
    assert np.log(fy / fx).mean() < -0.2


def test_cli_roundtrip(tmp_path):
    est = tmp_path / "est"
    cmd = [sys.executable, "-m", "ngsdose"]
    subprocess.run(cmd + ["estimate", *map(str, nygc[:6]), "-r", str(BUNDLE.dir), "-o", str(est), "-t", str(tmp_path / "single.tsv")],
                   check=True, cwd=ROOT)
    subprocess.run(cmd + ["cohort", *map(str, sorted(est.glob("*.json.gz"))), "-r", str(BUNDLE.dir), "-t", str(tmp_path / "cohort.tsv"),
                          "--save-efficiencies", str(tmp_path / "eff.json")], check=True, cwd=ROOT)
    rows = list(csv.DictReader(open(tmp_path / "cohort.tsv"), delimiter="\t"))
    assert len(rows) == 6 and all(float(r["rDNA45S.cn"]) > 100 for r in rows)
    assert (tmp_path / "eff.json").stat().st_size > 1000
