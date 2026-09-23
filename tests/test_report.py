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
    # every figure stands alone in a screenshot: numbered, captioned, with its source; the 18S ratio is named as the published method
    n_fig = html.count('<figure id="fig-')
    n_supp = html.count('<div class="title">Figure S')
    assert n_fig >= 5 and html.count("<figcaption>") == n_fig and html.count('<span class="src">') == n_fig
    assert all(f'<div class="title">Figure {i}. ' in html for i in range(1, n_fig - n_supp + 1))
    assert all(f'<div class="title">Figure S{i}. ' in html for i in range(1, n_supp + 1))
    assert "18S depth ratio (published" in html and "18S depth ratio (literature)" not in html
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


def test_the_coverage_qc_comparison(report):
    """NGS-PCA's QC table for the same files: a synthetic one built from the pilot's own numbers with a
    constant duplicate-flag offset and noise must come back at r > 0.99 with the offset recovered."""
    import numpy as np

    from ngsdose import report as R
    d = json.loads((report / "report.json").read_text())
    rows = [dict(s) for s in d["samples"]]
    rng = np.random.default_rng(3)
    qc = report / "sample_qc.tsv"
    with open(qc, "w") as fh:
        fh.write("SAMPLE_ID\tMEAN_AUTOSOMAL_COV\tX_COV_RATIO\tY_COV_RATIO\tINFERRED_SEX\tMTDNA_CN\tRELEASE_BATCH\n")
        for r in rows:
            fh.write(f'{r["sample"]}\t{0.88 * r["depth"]:.3f}\t{r["truth.chrX"] / 2:.4f}\t{0.85 * r["truth.chrY"] / 2:.4f}\t{r["sex_inferred"]}\t{0.92 * r["chrM.copies"] * np.exp(rng.normal(0, 0.02)):.2f}\t2504\n')
        fh.write("NA00000\t30\t0.5\t0.4\tM\t500\t698\n")                               # a sample that was not counted
    q = R.ngspca_qc_comparison(rows, qc)
    assert q["n"] == 12 and q["n_qc"] == 13 and q["sex_agree"] == q["sex_n"] == 12
    assert q["mtdna"]["r"] > 0.99 and abs(q["mtdna"]["ratio"]["median"] - 0.92) < 0.02 and q["chrX"]["r"] > 0.999 and abs(q["depth"]["ratio"]["median"] - 0.88) < 0.01
    assert [m["sample"] for m in q["mosaic_X"]] == ["HG00732"] and rows[0]["ngspca.batch"] == "2504"


def test_a_man_with_two_x_chromosomes_is_set_apart():
    """A man whose reads show two X chromosomes and a Y (47,XXY) is not a failure of the X model: the men's
    figures are those of men with one X, and he is listed on his own."""
    from ngsdose.report import known_truth
    row = lambda s, sex, x, y: {"sample": s, "sex_inferred": sex, "truth.auto": 2.0, "truth.chrX": x, "truth.chrY": y, "DJ.cn": 10.0}
    rows = [row(f"M{i}", "M", 0.99 + 0.005 * (i % 3), 1.0) for i in range(20)] + [row(f"F{i}", "F", 1.95, 0.0) for i in range(20)]
    rows.append(row("XXY", "M", 1.97, 0.97))
    sx = known_truth(rows)["sex"]
    assert [d["sample"] for d in sx["men_extra_x"]] == ["XXY"]
    assert sx["men_one_x"]["n"] == 20 and sx["men_one_x"]["max"] < 1.01
