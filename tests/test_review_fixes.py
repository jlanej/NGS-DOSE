"""Fixes from the review of the fetch-menu change (2026-09-28): a foreign warning no longer hangs the trio
analysis; a column a table lacks gets an NA row, not an aborted run; a --fetch-sinks path that is not
there is refused; a bare copy of the bundle's directory says what it cannot check instead of measuring a
subset-only fetch as a family; a fetch from another pipeline than the sinks' is said so; sub-option
verdicts read an unnamed sink row as every class's and never rest on no placement at all; the fetch menu
is not taken from the source checkout beside a configured bundle; estimate files are read by content."""
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from ngsdose import contract, tables
from ngsdose.estimate import SubOption

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("NGSDOSE_BIN", ROOT / "target" / "release" / "ngs-dose"))
BUNDLE = ROOT / "resources" / "GRCh38"
FIXTURE = ROOT / "tests" / "data" / "NA12878.subsample.bam"


def ngsdose(*args, env=None):
    return subprocess.run([sys.executable, "-m", "ngsdose", *map(str, args)], capture_output=True, text=True,
                          env={**os.environ, **(env or {})}, cwd=ROOT)


def test_a_foreign_warning_inside_the_trio_analysis_is_shown_and_the_command_returns():
    """The command's own warnings are printed once; any other warning raised inside (numpy's) is shown as
    usual. It used to be re-shown inside the recording context, which recorded it again into the list being
    walked, without end - so this runs in a subprocess with a deadline."""
    code = ("import warnings\nfrom ngsdose import cli, trios\n"
            "with cli._pedigree_warnings('trios'):\n"
            "    warnings.warn('ped', trios.PedigreeWarning); warnings.warn('ped', trios.PedigreeWarning)\n"
            "    warnings.warn('numpy-like', RuntimeWarning)\n"
            "print('returned')\n")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=60, cwd=ROOT)
    assert r.returncode == 0 and "returned" in r.stdout, r.stderr
    assert r.stderr.count("[trios] WARNING: ped") == 1 and "RuntimeWarning: numpy-like" in r.stderr


@pytest.fixture
def cohort(tmp_path):
    rng = np.random.default_rng(7)
    fams = [(f"C{i}", f"F{i}", f"M{i}") for i in range(12)]
    rows = []
    for c, f, m in fams:
        a, b = rng.normal(100, 10, 2)
        rows += [(c, (a + b) / 2 + rng.normal(0, 3)), (f, a), (m, b)]
    (tmp_path / "t.tsv").write_text("sample\tvar\n" + "".join(f"{s}\t{v:.3f}\n" for s, v in rows))
    (tmp_path / "ped.txt").write_text("child father mother\n" + "".join(f"{c} {f} {m}\n" for c, f, m in fams))
    return tmp_path


def test_a_column_the_table_lacks_gets_an_NA_row_and_the_others_their_values(cohort):
    d = cohort
    r = ngsdose("trios", d / "t.tsv", "-p", d / "ped.txt", "-c", "var", "nosuch", "--perm", "0", "--json", d / "o.json")
    assert r.returncode == 0 and "Traceback" not in r.stderr, r.stderr
    assert "[trios] WARNING: no column nosuch" in r.stderr
    lines = r.stdout.splitlines()
    assert lines[1].startswith("var\t12\t") and lines[2].startswith("nosuch\t0\tNA") and "# no column nosuch" in lines[2]
    o = json.loads((d / "o.json").read_text())
    assert o["var"]["n_trios"] == 12 and "no column nosuch" in o["nosuch"]["error"]
    r = ngsdose("trios", d / "t.tsv", "-p", d / "ped.txt", "-c", "nosuch", "--perm", "0")
    assert r.returncode == 1 and "none of the columns" in r.stderr and "Traceback" not in r.stderr
    r = ngsdose("trios", d / "t.tsv", "-p", d / "ped.txt", "-c", "var", "--compare-to", "nosuch", "--perm", "0")
    assert r.returncode == 1 and "--compare-to: no column nosuch" in r.stderr


def test_a_fetch_sinks_path_that_is_not_there_is_refused(tmp_path):
    r = ngsdose("estimate", tmp_path / "x.json.gz", "--fetch-sinks", tmp_path / "plan.sinks.bd", "-t", "-")
    assert r.returncode == 1 and "--fetch-sinks" in r.stderr and "no such file" in r.stderr and "Traceback" not in r.stderr


def _fetch_counts(sinks_sha, placements, contigs=(("chrX", 156_040_895),)):
    return dict(format="ngs-dose-counts/1", sample="s", mode="fetch", sinks="x.bed", sinks_sha256=sinks_sha,
                placement_bin=1000, placement_bin_compositional=10000,
                classes=[dict(name="aSatHOR", kind="compositional", reads=100)],
                contigs=[dict(name=c, len=n) for c, n in contigs], placements=placements)


def test_sub_option_verdicts_read_unnamed_rows_as_every_class_s_and_need_a_placed_bin():
    subs = {"DXZ1": SubOption("DXZ1", "aSatHOR", [("chrX", 58_000_000, 58_010_000)], "")}
    named_only = [("chrX", 58_000_000, 58_010_000, "aSatHOR")]
    with_unnamed = named_only + [("chr1", 1_000_000, 1_010_000, "")]
    c = _fetch_counts("h", [])
    # the BED is known: only the sub-option's interval names the class -> subset only; an unnamed row outside
    # it serves every class (the engine reads it so), so the family was read beyond its sub-option
    assert set(contract.subset_only(c, subs, {"h": named_only})) == {"aSatHOR"}
    assert contract.subset_only(c, subs, {"h": with_unnamed}) == {}
    # the BED is not known: the placements decide, and no placed bin decides nothing
    assert contract.subset_only(c, subs, {}) == {}
    assert contract.subset_only(_fetch_counts("h", [dict(**{"class": "aSatHOR"}, contig="*", start=0, reads=5)]), subs, {}) == {}
    inside = [dict(**{"class": "aSatHOR"}, contig="chrX", start=58_000_000, reads=50)]
    assert set(contract.subset_only(_fetch_counts("h", inside), subs, {})) == {"aSatHOR"}
    outside = inside + [dict(**{"class": "aSatHOR"}, contig="chrX", start=60_000_000, reads=50)]
    assert contract.subset_only(_fetch_counts("h", outside), subs, {}) == {}


def test_a_fetch_from_another_pipeline_than_the_sinks_is_said_so():
    learned = dict(aligner="bwa", aligner_version="0.7.15-r1140", sq_sha256="abc", sq_n=3)
    rec = lambda pg, sq="abc", n=3: dict(format="ngs-dose-counts/1", mode="fetch", pipeline=dict(pg=pg, sq_sha256=sq, sq_n=n))
    assert contract.pipeline_check(rec([dict(id="bwa", pn="bwa", vn="0.7.15")]), learned) == "ok"
    assert contract.pipeline_check(rec([dict(id="bwa.1", pn="bwa-mem2")]), learned) == "ok"          # the same aligner's successor
    why = contract.pipeline_check(rec([dict(id="DRAGEN", pn="DRAGEN", vn="4.2.7")]), learned)
    assert "dragen" in why and "not bwa" in why
    why = contract.pipeline_check(rec([dict(pn="bwa")], sq="other"), learned)
    assert "@SQ" in why
    # nothing to compare: no record, a scan, a header without @PG lines whose @SQ set is a subset (the test fixture)
    assert contract.pipeline_check(dict(mode="fetch"), learned) is None
    assert contract.pipeline_check(dict(mode="scan", pipeline=dict(pg=[dict(pn="DRAGEN")], sq_sha256="x", sq_n=3)), learned) is None
    assert contract.pipeline_check(rec([], sq="other", n=144), learned) is None
    assert contract.pipeline_check(rec([dict(pn="bwa")]), None) is None


def test_estimate_files_are_read_by_content_not_by_name(tmp_path):
    d = {"sample": "s", "a": 1}
    with gzip.open(tmp_path / "x.estimate.json", "wt") as fh:                    # gzip content, plain name
        json.dump(d, fh)
    (tmp_path / "y.estimate.json.gz").write_text(json.dumps(d))                  # plain content, gzip name
    assert tables.load_result(tmp_path / "x.estimate.json") == d
    assert tables.load_result(tmp_path / "y.estimate.json.gz") == d


def test_the_fetch_menu_is_not_taken_from_the_checkout_beside_a_configured_bundle(tmp_path, monkeypatch):
    from ngsdose import fetchplan
    monkeypatch.setenv("NGSDOSE_RESOURCES", str(tmp_path / "GRCh38"))
    with pytest.raises(ValueError, match="no fetch menu beside the resource bundle"):
        fetchplan.default_menu()
    r = ngsdose("fetchplan", "--list", env={"NGSDOSE_RESOURCES": str(tmp_path / "GRCh38")})
    assert r.returncode == 1 and "no fetch menu beside the resource bundle" in r.stderr and "Traceback" not in r.stderr
    (tmp_path / "fetch_menu.tsv").write_text("# menu\n")
    assert fetchplan.default_menu() == tmp_path / "fetch_menu.tsv"
    monkeypatch.delenv("NGSDOSE_RESOURCES")
    assert fetchplan.default_menu() == ROOT / "resources" / "fetch_menu.tsv"


@pytest.mark.skipif(not BIN.exists() or not FIXTURE.exists(), reason="needs the engine and the fixture")
def test_a_bare_copy_of_the_bundle_says_what_it_cannot_check_instead_of_measuring(tmp_path):
    """A fetch through a sinks BED that is not the bundle's own, estimated with a bare copy of the bundle's
    directory (no resources/experimental beside it): the sub-option definitions are unknown there, so the
    fetch's satellite families are unverified, with a warning saying why - not measured as if read whole.
    With the repository's layout the same fetch is judged by the definitions."""
    bed = tmp_path / "other.bed"
    bed.write_text("# a copy of the bundle's sinks with another hash\n" + (BUNDLE / "sinks.bed").read_text())
    tel = ROOT / "resources" / "experimental" / "telomere.k31.panel.tsv.gz"
    out = tmp_path / "f.json.gz"
    subprocess.run([str(BIN), "count", "-m", "fetch", "-i", str(FIXTURE), "-p", str(BUNDLE / "panel.k31.tsv.gz"), "-p", str(tel),
                    "-c", str(BUNDLE / "controls.fa.gz"), "--sinks", str(bed), "-@", "2", "-o", str(out)], check=True, capture_output=True)
    assert hashlib.sha256(bed.read_bytes()).hexdigest() != hashlib.sha256((BUNDLE / "sinks.bed").read_bytes()).hexdigest()
    bare = tmp_path / "elsewhere" / "GRCh38"
    shutil.copytree(BUNDLE, bare)
    status = lambda r: dict(zip(r.stdout.splitlines()[0].split("\t"), r.stdout.splitlines()[1].split("\t")))
    r = ngsdose("estimate", out, "--fetch-sinks", bed, "-r", bare, "-t", "-")
    assert r.returncode == 0, r.stderr
    assert "no experimental resources beside the bundle" in r.stderr and str(bare.parent / "experimental") in r.stderr
    s = status(r)
    assert s["TEL.status"] == "unverified" and s["TEL.mass_Mb"] == "NA" and s["rDNA45S.status"] == "ok"
    r = ngsdose("estimate", out, "--fetch-sinks", bed, "-r", BUNDLE, "-t", "-")
    assert r.returncode == 0 and "no experimental resources" not in r.stderr, r.stderr
    assert status(r)["TEL.status"] == "ok"
    # named outright, a directory that is not there is reported by name
    r = ngsdose("estimate", out, "-r", bare, "-t", "-", env={"NGSDOSE_EXPERIMENTAL": str(tmp_path / "nowhere")})
    assert r.returncode == 0 and str(tmp_path / "nowhere") in r.stderr
    assert status(r)["TEL.status"] == "unverified"
    # the bundle's own sinks need no definitions: a fetch through them is measured as before
    out2 = tmp_path / "g.json.gz"
    subprocess.run([str(BIN), "count", "-m", "fetch", "-i", str(FIXTURE), "-p", str(BUNDLE / "panel.k31.tsv.gz"), "-p", str(tel),
                    "-c", str(BUNDLE / "controls.fa.gz"), "--sinks", str(BUNDLE / "sinks.bed"), "-@", "2", "-o", str(out2)], check=True, capture_output=True)
    r = ngsdose("estimate", out2, "-r", bare, "-t", "-")
    assert r.returncode == 0 and status(r)["TEL.status"] == "ok"


def test_control_pcs_are_computed_on_one_set_of_regions_in_one_order():
    """Two estimates with the same number of control regions but another set or order of them (another
    controls file, another bundle revision) must not be pooled into one PCA: the vectors are matched by
    position. The hash `estimate` records of the regions' names in order tells them apart."""
    import copy

    from test_cohort_robustness import _results, _table
    thirty = _results(30, seed=7, classes=("DJ",))
    for r in thirty:
        r["control_qc"]["region_order_sha256"] = "one-order"
    mixed = copy.deepcopy(thirty)
    for i in (3, 11, 20):
        mixed[i]["control_qc"]["region_log_ratio"] = mixed[i]["control_qc"]["region_log_ratio"][::-1]
        mixed[i]["control_qc"]["region_order_sha256"] = "another-order"
    rows, _, info, log = _table(mixed, n_control_pcs=5)
    alone, _, info_alone, _ = _table([r for i, r in enumerate(thirty) if i not in (3, 11, 20)], n_control_pcs=5)
    assert info["excluded"] == ["s03", "s11", "s20"] and info["mp"] == info_alone["mp"]
    assert all("ctrlPC1" not in rows[s] for s in ("s03", "s11", "s20"))
    for s, r in alone.items():
        assert [rows[s].get(f"ctrlPC{k}") for k in range(1, 6)] == [r[f"ctrlPC{k}"] for k in range(1, 6)]
    assert "not the same regions in the same order" in log and "s03" in log and "800 control regions" in log
    # estimates written before the hash existed are grouped by their count, as before
    old = copy.deepcopy(thirty)
    for r in old:
        del r["control_qc"]["region_order_sha256"]
    rows_old, _, info_old, _ = _table(old, n_control_pcs=5)
    assert "excluded" not in info_old and all("ctrlPC1" in r for r in rows_old.values())
