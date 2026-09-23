"""`ngsdose report` on the committed pilot counts: the page and its tables exist, and the numbers in
them are the pilot's numbers."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "example" / "1000G" / "pilot"
pytestmark = pytest.mark.skipif(len(list((PILOT / "counts_nygc").glob("*.json.gz"))) < 6, reason="pilot counts not present")


@pytest.fixture(scope="module")
def report(tmp_path_factory):
    out = tmp_path_factory.mktemp("report")
    ped = out / "ped.txt"
    # the four pilot trios, with sex and population, in the 1000 Genomes pedigree format
    lines = ["FamilyID SampleID FatherID MotherID Sex Population Superpopulation"]
    for fam, c, f, m, cs, pop, sp in (("1463", "NA12878", "NA12891", "NA12892", 2, "CEU", "EUR"), ("Y117", "NA19240", "NA19239", "NA19238", 2, "YRI", "AFR"),
                                       ("SH032", "HG00514", "HG00512", "HG00513", 2, "CHS", "EAS"), ("PR05", "HG00733", "HG00731", "HG00732", 2, "PUR", "AMR")):
        lines += [f"{fam} {c} {f} {m} {cs} {pop} {sp}", f"{fam} {f} 0 0 1 {pop} {sp}", f"{fam} {m} 0 0 2 {pop} {sp}"]
    ped.write_text("\n".join(lines) + "\n")
    subprocess.run([sys.executable, "-m", "ngsdose", "report", "--fetch", str(PILOT / "counts_nygc"), "-p", str(ped), "--hall", str(PILOT / "hall2021_MOESM1.txt"), "--pilot", str(PILOT),
                    "--pcs", str(ROOT / "example/1000G/ngspca/svd.pcs.txt"), "-o", str(out), "-j", "2", "--as-of", "2026-09-22"], check=True, cwd=ROOT, capture_output=True)
    return out


def test_outputs_exist_and_are_reproducible(report):
    for f in ("index.html", "report.json", "data/cohort.tsv", "data/fetch.tsv", "data/transmission.tsv", "data/flags.tsv", "data/efficiencies.json"):
        assert (report / f).stat().st_size > 0, f
    html = (report / "index.html").read_text()
    for must in ("Every number and figure on this page is recomputed", 'id="chart-auto"', 'id="chart-trio"', "Hall, Turner", "report-data", "<table"):
        assert must in html
    visible = html.split('<script id="report-data"')[0]                # the prose and tables, not the code
    assert "NaN" not in visible and "None" not in visible.replace("None of", "") and "nan" not in visible.split("<style>")[1].split("</style>")[0]
    # the second run reuses the cached estimates and produces the same page
    first = html
    subprocess.run([sys.executable, "-m", "ngsdose", "report", "--fetch", str(PILOT / "counts_nygc"), "-p", str(report / "ped.txt"), "--hall", str(PILOT / "hall2021_MOESM1.txt"), "--pilot", str(PILOT),
                    "--pcs", str(ROOT / "example/1000G/ngspca/svd.pcs.txt"), "-o", str(report), "-j", "2", "--as-of", "2026-09-22"], check=True, cwd=ROOT, capture_output=True)
    assert (report / "index.html").read_text() == first


def test_the_numbers_are_the_pilots(report):
    d = json.loads((report / "report.json").read_text())
    assert d["meta"]["n"] == 12 and d["meta"]["primary_mode"] == "fetch" and d["meta"]["n_fetch"] == 12
    kt = d["known_truth"]
    assert abs(kt["auto"]["mean"] - 1.996) < 0.003 and abs(kt["chrX"]["M"]["mean"] - 0.995) < 0.01 and 9.6 < kt["DJ"]["mean"] < 9.8
    assert kt["chrY"]["M"]["n"] == 4 and kt["chrY"]["F"]["n"] == 8 and kt["chrY"]["F"]["max"] < 0.01 and kt["sex"]["mismatch"] == []
    assert d["trios"]["n_complete"] == 4 and {t["column"] for t in d["trios"]["table"]} >= {"rDNA45S.cn", "rDNA45S.18S.flat", "truth.auto", "chrM.copies"}
    h = d["hall"]
    assert h["n"] == 5 and h["flat"]["r"] > 0.97 and 1.05 < h["flat_ratio"] < 1.1 and 1.0 < h["dup_corrected_ratio"] < 1.05
    assert any("HG00732" == s and "chrX" in f for s, f in d["flags"])            # the culture that lost an X
    rep = d["replicates"]                                                          # the same twelve people on an older technology
    assert rep["n"] == 12 and rep["table"]["calibrated"]["icc"] > 0.95 and rep["table"]["flat"]["icc"] < 0.5 and rep["table"]["flat_centred"]["icc"] > 0.8
    assert len(d["samples"]) == 12 and all(s["sex_inferred"] in ("M", "F") for s in d["samples"])
    cols = (report / "data" / "cohort.tsv").read_text().splitlines()[0].split("\t")
    assert {"sample", "rDNA45S.cn", "truth.chrY", "chrM.copies", "flags", "sex_inferred"} <= set(cols)
