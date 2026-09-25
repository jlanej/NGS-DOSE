"""The cohort layer end to end, on a cohort that can be built anywhere: sixty seeded subsamples of
the NA12878 fixture, named as twenty trios. Every "sample" is the same person, so there is no
true variance at all - which makes the expected answers known: calibrated copy numbers agree to
within sampling noise, coverage PCs explain nothing beyond chance, and no estimator is
"transmitted". This is the path a cohort pipeline takes (estimate -> cohort with
control PCs -> adjust on external and internal PCs -> trios with a paired comparison), which the
twelve-sample pilot is too small to exercise."""
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("NGSDOSE_BIN", ROOT / "target" / "release" / "ngs-dose"))
BAM = ROOT / "tests" / "data" / "NA12878.subsample.bam"
BUNDLE = ROOT / "resources" / "GRCh38"
N_TRIOS = 20

pytestmark = pytest.mark.skipif(not BIN.exists() or not shutil.which("samtools"), reason="needs the engine and samtools")


def rows_of(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def cohort(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("cohort")
    names = [f"S{i:03d}" for i in range(3 * N_TRIOS)]
    counts = []
    for i, s in enumerate(names):
        bam = tmp / f"{s}.bam"
        subprocess.run(["samtools", "view", "-b", "-s", f"{i + 1}.6", "-o", str(bam), str(BAM)], check=True)
        subprocess.run(["samtools", "index", "-c", str(bam)], check=True)
        out = tmp / f"{s}.json.gz"
        subprocess.run([str(BIN), "count", "-m", "fetch", "-i", str(bam), "-p", str(BUNDLE / "panel.k31.tsv.gz"), "-c", str(BUNDLE / "controls.fa.gz"),
                        "--sinks", str(BUNDLE / "sinks.bed"), "-s", s, "-@", "2", "-o", str(out)], check=True, capture_output=True)
        bam.unlink()
        counts.append(str(out))
    # 1000 Genomes style pedigree (child, father, mother in turn) and NGS-PCA style PCs (its sample-name suffix included)
    ped = tmp / "pedigree.txt"
    lines = ["FamilyID SampleID FatherID MotherID Sex Population Superpopulation"]
    for t in range(N_TRIOS):
        c, f, m = names[3 * t:3 * t + 3]
        lines += [f"F{t} {c} {f} {m} 1 POP SUP", f"F{t} {f} 0 0 1 POP SUP", f"F{t} {m} 0 0 2 POP SUP"]
    ped.write_text("\n".join(lines) + "\n")
    rng = np.random.default_rng(5)
    pcs = tmp / "svd.pcs.txt"
    pcs.write_text("SAMPLE\t" + "\t".join(f"PC{k + 1}" for k in range(8)) + "\n"
                   + "".join(f"{s}.by1000.\t" + "\t".join(f"{x:.5f}" for x in rng.normal(size=8)) + "\n" for s in names))
    run = lambda *a: subprocess.run([sys.executable, "-m", "ngsdose", *a], check=True, cwd=ROOT, capture_output=True, text=True)
    run("estimate", *counts, "-r", str(BUNDLE), "-o", str(tmp / "est"), "-t", str(tmp / "single.tsv"), "-j", "4")
    cohort_log = run("cohort", *map(str, sorted((tmp / "est").glob("*.json.gz"))), "-r", str(BUNDLE), "-t", str(tmp / "cohort.tsv")).stderr
    cols = ["rDNA45S.18S.flat", "rDNA45S.cn_single", "rDNA45S.cn", "rDNA5S.cn", "DJ.cn", "truth.auto", "chrEBV.copies", "chrM.copies"]
    # a number of PCs given outright (these random PCs come with no spectrum to choose from) ...
    adj = {"ngspca": run("adjust", str(tmp / "cohort.tsv"), "--pcs", str(pcs), "--n-pc", "8", "-c", *cols, "-o", str(tmp / "adj_ngspca.tsv")).stderr,
           "ctrlpc": run("adjust", str(tmp / "cohort.tsv"), "--n-pc", "5", "-c", *cols, "-o", str(tmp / "adj_ctrlpc.tsv")).stderr}
    # ... and the default: as many as stand above the Marchenko-Pastur edge of the control-region spectrum
    default = run("adjust", str(tmp / "cohort.tsv"), "-c", *cols, "-o", str(tmp / "adj_default.tsv")).stderr
    # a number larger than one of the PC sets holds (one knob for two PC sets is how a pipeline dies at its last step)
    clamped = run("adjust", str(tmp / "cohort.tsv"), "--n-pc", "46", "-c", "rDNA45S.cn", "-o", str(tmp / "adj_clamped.tsv")).stderr
    sweep = run("pcsweep", str(tmp / "cohort.tsv"), "-c", "rDNA45S.cn", "DJ.cn_single", "-p", str(ped), "--boot", "100", "-o", str(tmp / "sweep.tsv")).stderr
    trios = run("trios", str(tmp / "adj_ngspca.tsv"), "-p", str(ped), "-c", *[c + ".adj" for c in cols], "--compare-to", "rDNA45S.18S.flat.adj",
                "--json", str(tmp / "transmission.json")).stdout
    return dict(tmp=tmp, cols=cols, adjust_log=adj, trios=trios, cohort_log=cohort_log, default_log=default, sweep_log=sweep, clamped_log=clamped)


def test_the_same_person_sixty_times(cohort):
    rows = rows_of(cohort["tmp"] / "cohort.tsv")
    assert len(rows) == 3 * N_TRIOS
    for col, tol in (("rDNA45S.cn", 0.08), ("DJ.cn", 0.12), ("truth.auto", 0.08), ("chrM.copies", 0.12)):
        v = np.array([float(r[col]) for r in rows])
        assert np.all(np.abs(v / np.median(v) - 1) < tol), (col, v.min(), v.max())
    assert all(f"ctrlPC{k}" in rows[0] for k in range(1, 6)) and "ctrlPC_mp" in rows[0]
    assert {r["eof_marker"] for r in rows} == {"present"} and all(float(r["truth.chrY"]) < 0.05 for r in rows)


def test_adjustment_explains_no_more_than_chance(cohort):
    for which, log in cohort["adjust_log"].items():
        lines = [l for l in log.splitlines() if l.startswith("[adjust]") and "adjusted R2" in l]      # one line per column
        assert len(lines) == len(cohort["cols"]), which
        for l in lines:
            r2_adj = float(l.split("adjusted R2 ")[1].split("%")[0])
            assert "expected by chance" in l and r2_adj < 25, l
    adj = rows_of(cohort["tmp"] / "adj_ngspca.tsv")
    assert len(adj) == 3 * N_TRIOS and all(f"{c}.adj" in adj[0] for c in cohort["cols"])


def test_nothing_is_transmitted_when_there_is_nothing_to_transmit(cohort):
    import json
    res = json.loads((cohort["tmp"] / "transmission.json").read_text())
    table = [l.split("\t") for l in cohort["trios"].splitlines() if l and not l.startswith("#")]
    assert sum(1 for t in table if t[0].endswith(".adj")) >= 2 * len(cohort["cols"]) - 1          # reliabilities + paired differences
    per = res["columns"] if "columns" in res else res
    n_checked = 0
    for col, r in per.items():
        if not isinstance(r, dict) or "reliability_midparent" not in r:
            continue
        n_checked += 1
        assert r["n_trios"] == N_TRIOS
        lo, hi = r["reliability_midparent_ci95"]
        assert lo < hi and lo < 0.6, (col, lo, hi)             # an interval far above zero would be a bug: there is no signal
    assert n_checked == len(cohort["cols"])


def test_nothing_is_regressed_out_when_there_is_nothing_but_noise(cohort):
    """One person sixty times over has no structure in the control regions: the Marchenko-Pastur
    default finds no component above the noise bulk, the default adjustment therefore changes
    nothing, and the sweep - cross-validated error of the known truths, reliability of the
    classes, for every number of PCs - recommends none either."""
    rows = rows_of(cohort["tmp"] / "cohort.tsv")
    assert "Marchenko-Pastur edge" in cohort["cohort_log"] and {r["ctrlPC_mp"] for r in rows} <= {"0", "1"}
    adj = rows_of(cohort["tmp"] / "adj_default.tsv")
    if rows[0]["ctrlPC_mp"] == "0":
        assert "regressing out 0 PCs" in cohort["default_log"]
        assert all(abs(float(r["rDNA45S.cn.adj"]) - float(r["rDNA45S.cn"])) < 1e-3 for r in adj)
    sweep = rows_of(cohort["tmp"] / "sweep.tsv")
    assert {r["column"] for r in sweep} >= {"truth.auto", "truth.chrX", "DJ.cn", "rDNA45S.cn"}
    auto = sorted((int(r["n_pc"]), float(r["sd_log_robust"])) for r in sweep if r["column"] == "truth.auto")
    assert auto[0][0] == 0 and auto[-1][0] >= 5 and auto[-1][1] > 0.9 * auto[0][1]      # more PCs buy nothing here
    picks = {l.split("\t")[0]: int(l.split("\t")[7]) for l in cohort["sweep_log"].splitlines() if l.startswith(("truth.auto", "DJ.cn"))}
    assert picks and all(v <= 2 for v in picks.values()), picks


def test_a_number_of_pcs_beyond_what_the_table_holds_is_clamped_with_a_warning(cohort):
    n_ctrl = sum(1 for c in rows_of(cohort["tmp"] / "cohort.tsv")[0] if c.startswith("ctrlPC") and c != "ctrlPC_mp")
    assert n_ctrl < 46 and f"using {n_ctrl}" in cohort["clamped_log"] and "clamped" in cohort["clamped_log"]
    assert len(rows_of(cohort["tmp"] / "adj_clamped.tsv")) == 3 * N_TRIOS
