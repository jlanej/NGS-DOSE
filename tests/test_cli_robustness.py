"""The command line on the inputs a large cohort will hand it sooner or later: a truncated counts file among
good ones, a fetch that lacked sinks, a sample counted twice, a misspelt column, a sample without PCs, a
sparse column, a pedigree without populations, a table whose ids match nothing. Each is either handled
and said, or refused with a message that names it - never a traceback, never silence."""
import csv
import gzip
import json
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


def ngsdose(*args, cwd=None):
    env = {**os.environ, "PYTHONPATH": str(ROOT) + os.pathsep + os.environ.get("PYTHONPATH", "")}
    return subprocess.run([sys.executable, "-m", "ngsdose", *map(str, args)], capture_output=True, text=True, cwd=cwd, env=env)


def rows_of(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows):
    cols = list(rows[0])
    path.write_text("\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows))


# ---------------------------------------------------------------- estimate and cohort (need the engine)

@pytest.fixture(scope="module")
def estimated(tmp_path_factory):
    if not BIN.exists():
        pytest.skip("needs the engine (cargo build --release, or NGSDOSE_BIN)")
    d = tmp_path_factory.mktemp("est")
    good = d / "good.json.gz"
    subprocess.run([str(BIN), "count", "-i", str(BAM), "-p", str(BUNDLE / "panel.k31.tsv.gz"), "-c", str(BUNDLE / "controls.fa.gz"), "-m", "fetch",
                    "--sinks", str(BUNDLE / "sinks.bed"), "-@", "2", "-o", str(good)], check=True, capture_output=True)
    shutil.copy(good, d / "again.json.gz")
    (d / "trunc.json.gz").write_bytes(good.read_bytes()[:3000])
    c = json.load(gzip.open(good, "rt"))
    json.dump({**c, "sample": "NOSINK", "sinks_missing_classes": ["DJ"], "eof_marker": "unchecked"}, gzip.open(d / "nosink.json.gz", "wt"))
    fewer = [dict(x, panel_kmers=x["panel_kmers"] - 4000) if x["name"] == "rDNA45S" else x for x in c["classes"]]
    json.dump({**c, "sample": "FEWER", "classes": fewer}, gzip.open(d / "fewer.json.gz", "wt"))
    files = [d / f"{n}.json.gz" for n in ("good", "trunc", "again", "nosink", "fewer")]
    r = ngsdose("estimate", *files, "-r", BUNDLE, "-o", d / "est", "-t", d / "t.tsv", "-j", "2")
    return d, r


def test_one_bad_file_costs_only_its_own_row(estimated):
    d, r = estimated
    assert r.returncode == 1 and "Traceback" not in r.stderr
    assert f"FAILED {d / 'trunc.json.gz'}: unreadable counts file" in r.stderr and "1 of 5 counts files failed" in r.stderr
    rows = rows_of(d / "t.tsv")
    assert [x["sample"] for x in rows] == ["NA12878", "NA12878", "NOSINK", "FEWER"]
    # every estimate is named after its input, so two files of one sample do not overwrite each other
    assert sorted(p.name for p in (d / "est").iterdir()) == [f"{n}.estimate.json.gz" for n in ("again", "fewer", "good", "nosink")]
    assert "sample NA12878 is in 2 inputs" in r.stderr


def test_what_was_not_counted_in_full_is_na_and_said(estimated):
    d, r = estimated
    row = {x["sample"]: x for x in rows_of(d / "t.tsv")}
    assert row["NA12878"]["DJ.status"] == "ok" and float(row["NA12878"]["DJ.cn_single"]) > 5
    assert row["NOSINK"]["DJ.status"] == "no_sinks_in_fetch" and row["NOSINK"]["DJ.cn_single"] == "NA"
    assert "NOSINK: DJ is not estimated: no sink intervals in the fetch" in r.stderr
    assert "NOSINK: the input's end-of-file marker could not be checked" in r.stderr
    assert row["FEWER"]["rDNA45S.status"] == "panel_mismatch" and row["FEWER"]["rDNA45S.cn_single"] == "NA"
    assert "rDNA45S was not counted with the same k-mers in every file" in r.stderr
    assert "not estimated (NA in the table): DJ in 1 sample(s): no_sinks_in_fetch" in r.stderr


def test_estimates_that_would_share_a_name_are_refused(estimated, tmp_path):
    d, _ = estimated
    for sub in ("scan", "fetch"):
        (tmp_path / sub).mkdir()
        shutil.copy(d / "good.json.gz", tmp_path / sub / "HG00096.json.gz")
    r = ngsdose("estimate", tmp_path / "scan" / "HG00096.json.gz", tmp_path / "fetch" / "HG00096.json.gz", "-o", tmp_path / "est")
    assert r.returncode == 1 and "HG00096.estimate.json.gz" in r.stderr and "Traceback" not in r.stderr
    assert not (tmp_path / "est").exists()


def test_cohort_refuses_a_sample_twice_with_a_message(estimated):
    d, _ = estimated
    r = ngsdose("cohort", d / "est" / "good.estimate.json.gz", d / "est" / "again.estimate.json.gz", "-r", BUNDLE)
    assert r.returncode == 1 and "NA12878: given twice" in r.stderr and "Traceback" not in r.stderr
    (d / "bad.estimate.json.gz").write_bytes(b"\x1f\x8b\x08\x00")
    r = ngsdose("cohort", d / "est" / "good.estimate.json.gz", d / "bad.estimate.json.gz", "-r", BUNDLE)
    assert r.returncode == 1 and "bad.estimate.json.gz: unreadable estimate file" in r.stderr and "Traceback" not in r.stderr


# ---------------------------------------------------------------- adjust and pcsweep

@pytest.fixture()
def table(tmp_path):
    rng = np.random.default_rng(3)
    rows = []
    for i in range(40):
        pc = rng.normal(size=4)
        rows.append(dict(sample=f"S{i:02d}", **{"X.cn": round(float(np.exp(4 + 0.1 * pc[0] + rng.normal(0, 0.05))), 3)},
                         **{"sparse.cn": round(float(np.exp(3 + rng.normal(0, 0.1))), 3) if i < 12 else "NA"},
                         **{f"ctrlPC{k + 1}": round(float(pc[k]), 4) for k in range(4)}, ctrlPC_mp=2))
    rows[5].update(ctrlPC3="NA", ctrlPC4="NA")               # beyond the PCs used: still adjusted
    rows[7].update(ctrlPC1="NA", ctrlPC2="NA", ctrlPC3="NA", ctrlPC4="NA")   # no control residuals: NA, but kept
    path = tmp_path / "cohort.tsv"
    write_tsv(path, rows)
    return path


def test_adjust_keeps_every_row_and_says_what_it_could_not_adjust(table, tmp_path):
    out = tmp_path / "adj.tsv"
    r = ngsdose("adjust", table, "--n-pc", "2", "-c", "X.cn", "sparse.cn", "-o", out)
    assert r.returncode == 0, r.stderr
    rows = {x["sample"]: x for x in rows_of(out)}
    assert len(rows) == 40 and rows["S07"]["X.cn.adj"] == "NA" and rows["S07"]["X.cn"] != "NA"
    assert rows["S05"]["X.cn.adj"] != "NA"                     # only the two PCs used must be there
    assert "1 of 40 samples have no finite ctrlPC1..ctrlPC2" in r.stderr and "S07" in r.stderr
    assert all(x["sparse.cn.adj"] == "NA" for x in rows.values())
    assert "sparse.cn: 11 usable values (positive) among the 39 samples with PCs; 2 PCs need at least 12" in r.stderr   # S07 has none
    sd = lambda c: np.std(np.log([float(x[c]) for x in rows.values() if x[c] != "NA"]))
    assert sd("X.cn.adj") < sd("X.cn")                          # the PC effect is regressed out


def test_adjust_refuses_a_missing_column_and_a_repeated_sample(table, tmp_path):
    r = ngsdose("adjust", table, "--n-pc", "2", "-c", "X.cnn", "-o", tmp_path / "a.tsv")
    assert r.returncode == 1 and "no column X.cnn" in r.stderr and "Traceback" not in r.stderr and not (tmp_path / "a.tsv").exists()
    lines = table.read_text().splitlines()
    dup = tmp_path / "dup.tsv"
    dup.write_text("\n".join(lines + [lines[1]]) + "\n")
    r = ngsdose("adjust", dup, "--n-pc", "2", "-c", "X.cn")
    assert r.returncode == 1 and "S00 given more than once" in r.stderr and "Traceback" not in r.stderr


def test_adjust_on_external_pcs_names_the_samples_it_lacks(table, tmp_path):
    rng = np.random.default_rng(1)
    pcs = tmp_path / "svd.pcs.txt"
    pcs.write_text("SAMPLE\tPC1\tPC2\n" + "".join(f"S{i:02d}.by1000.\t{rng.normal():.4f}\t{rng.normal():.4f}\n" for i in range(5, 40)))
    out = tmp_path / "adj.tsv"
    r = ngsdose("adjust", table, "--pcs", pcs, "--n-pc", "2", "-c", "X.cn", "-o", out)
    assert r.returncode == 0, r.stderr
    rows = rows_of(out)
    assert len(rows) == 40 and [x["sample"] for x in rows if x["X.cn.adj"] == "NA"] == ["S00", "S01", "S02", "S03", "S04"]
    assert "5 of 40 samples have no PCs in" in r.stderr and "S00, S01, S02, S03, S04" in r.stderr
    # ids that match nothing: refused, with the hint that they do not match
    bad = tmp_path / "bad" / "svd.pcs.txt"
    bad.parent.mkdir()
    bad.write_text("SAMPLE\tPC1\tPC2\n" + "".join(f"NA{i}\t0.1\t0.2\n" for i in range(40)))
    for cmd in ("adjust", "pcsweep"):
        r = ngsdose(cmd, table, "--pcs", bad, "--n-pc", "2", "-c", "X.cn")
        assert r.returncode == 1 and "only 0 samples have coverage PCs" in r.stderr and "Most ids do not match" in r.stderr, r.stderr
        assert "Traceback" not in r.stderr


def test_pcsweep_refuses_a_bad_truth_spec(table):
    for spec in ("X.cn", "X.cn=ten", "X.cn=0", "=3", "nosuch=2"):
        r = ngsdose("pcsweep", table, "--n-pc", "2", "-c", "X.cn", "--truth", spec, "--boot", "0")
        assert r.returncode == 1 and "Traceback" not in r.stderr, (spec, r.stderr)
        assert ("no column nosuch" if spec.startswith("nosuch") else "expected COLUMN=VALUE") in r.stderr


# ---------------------------------------------------------------- trios

@pytest.fixture()
def family(tmp_path):
    rng = np.random.default_rng(2)
    rows, ped, pop = [], ["child father mother"], []
    for t in range(24):
        c, f, m = f"C{t}", f"F{t}", f"M{t}"
        g = "A" if t % 2 else "B"
        vf, vm = rng.normal(100, 10) + (20 if g == "A" else 0), rng.normal(100, 10) + (20 if g == "A" else 0)
        vc = (vf + vm) / 2 + rng.normal(0, 3)
        rows += [dict(sample=s, var=round(v, 3), empty="NA", few=round(v, 3) if t < 2 else "NA") for s, v in ((c, vc), (f, vf), (m, vm))]
        ped.append(f"{c} {f} {m}")
        pop += [f"{s}\t{g}" for s in (c, f, m)]
    write_tsv(tmp_path / "t.tsv", rows)
    (tmp_path / "ped.txt").write_text("\n".join(ped) + "\n")
    (tmp_path / "pop.txt").write_text("\n".join(pop) + "\n")
    return tmp_path


def test_trios_centres_within_a_population_file(family):
    d = family
    r = ngsdose("trios", d / "t.tsv", "-p", d / "ped.txt", "-c", "var", "--perm", "0", "--json", d / "a.json")
    assert r.returncode == 0 and r.stderr.count("[trios] WARNING: values not centred within population") == 1, r.stderr
    assert json.loads((d / "a.json").read_text())["var"]["population_centred"] is False
    r = ngsdose("trios", d / "t.tsv", "-p", d / "ped.txt", "--population", d / "pop.txt", "-c", "var", "--perm", "0", "--json", d / "b.json")
    assert r.returncode == 0 and "not centred" not in r.stderr, r.stderr
    b = json.loads((d / "b.json").read_text())["var"]
    assert b["population_centred"] is True and b["n_populations"] == 2
    r = ngsdose("trios", d / "t.tsv", "-p", d / "ped.txt", "--no-population-centring", "-c", "var", "--perm", "0")
    assert r.returncode == 0 and "WARNING" not in r.stderr


def test_trios_fails_when_nothing_could_be_analysed(family):
    d = family
    r = ngsdose("trios", d / "t.tsv", "-p", d / "ped.txt", "-c", "var", "nosuch")
    assert r.returncode == 1 and "no column nosuch" in r.stderr and "Traceback" not in r.stderr
    r = ngsdose("trios", d / "t.tsv", "-p", d / "ped.txt", "--no-population-centring", "-c", "empty", "few", "--json", d / "e.json")
    assert r.returncode == 1 and "none of the 2 columns could be analysed" in r.stderr
    assert "# need at least 3 complete trios" in r.stdout and "# no numeric values" in r.stdout
    assert set(json.loads((d / "e.json").read_text())) == {"empty", "few"}
    r = ngsdose("trios", d / "t.tsv", "-p", d / "ped.txt", "--no-population-centring", "-c", "var", "few")
    assert r.returncode == 0                                    # one column is enough
    (d / "noped.txt").write_text("child father mother\nC0 0 0\n")
    r = ngsdose("trios", d / "t.tsv", "-p", d / "noped.txt", "-c", "var")
    assert r.returncode == 1 and "holds no complete trio" in r.stderr and "Traceback" not in r.stderr
