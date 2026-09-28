"""Two things the estimator accepts beyond counts made with the bundle's own resources: counts made
with a named subset of the bundle's control regions (a lighter fetch), and positional classes of
experimental panels whose unit lies outside the bundle. Everything else is refused or skipped as
before, and counts made with the bundle's resources give the same results as without these options.

Most tests need no engine and no data: the counts are simulated (test_estimate_robustness._sim),
with per-region GC tables so that counts made with a subset of the controls can be written as the
engine would. The last runs the engine on the fixture and `ngsdose estimate` on what it wrote."""
import copy
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from ngsdose import contract, estimate, resources
from ngsdose.io import Panel
from test_estimate_robustness import L, _sim

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("NGSDOSE_BIN", ROOT / "target" / "release" / "ngs-dose"))
BAM = ROOT / "tests" / "data" / "NA12878.subsample.bam"
BUNDLE = ROOT / "resources" / "GRCh38"


@pytest.fixture(scope="module")
def simc():
    return _sim(seed=3)


def _subset_counts(counts, names, tables, keep):
    """The counts an engine run with only the control regions `keep` would write: the same reads,
    with the GC tables summed over those regions. The simulation keeps per-region N, not per-region
    O, so O is rebuilt here from each region's rate x its N (exact for the simulated rates)."""
    rows = [i for i, (n, role) in enumerate(names) if role == "control" and n in keep]
    obs = {r["name"]: r["obs"] for r in counts["regions"]}
    tabs = []
    for t in counts["gc_tables"]:
        N_all, O_all = np.array(t["n"], float), np.array(t["o"], float)
        rho = np.where(N_all > 0, O_all / np.maximum(N_all, 1), 0.0)
        N = tables[rows].sum(0)
        O = sum(obs[names[i][0]] * tables[i] * rho / max((tables[i] * rho).sum(), 1e-9) for i in rows)
        tabs.append(dict(l=t["l"], n=N.tolist(), o=np.round(O).tolist()))
    regs = [r for r in counts["regions"] if r.get("role", "control") != "control" or r["name"] in keep]
    return {**counts, "regions": regs, "gc_tables": tabs}


def _half(names):
    ctrl = [n for n, role in names if role == "control"]
    return frozenset(ctrl[::2])


def test_counts_made_with_a_named_subset_of_the_controls_are_estimated(simc):
    counts, names, tables, panel, units = simc
    keep = _half(names)
    lite = _subset_counts(counts, names, tables, keep)
    run = lambda c, **kw: estimate.estimate_sample(c, panel, units, region_tables=tables, regions=names, L=L, **kw)
    with pytest.raises(ValueError, match="different controls file.*no named subset"):
        run(lite)
    with pytest.raises(ValueError, match="different controls file"):
        run(lite, control_subsets={"other": keep - {next(iter(keep))}})
    full = run(counts)
    r = run(lite, control_subsets={"lite": keep, "other": frozenset(list(keep)[:3])})
    assert r["controls_subset"] == "lite" and r["controls_used"] == len(keep) == r["control_qc"]["n_regions"]
    assert len(r["control_qc"]["region_log_ratio"]) == len(keep)
    assert abs(r["classes"]["unit"]["cn"] / full["classes"]["unit"]["cn"] - 1) < 0.03
    assert abs(r["truth_regions"]["auto"]["cn"] - 2) < 0.1 and r["classes"]["sat"]["status"] == "ok"


def test_counts_made_with_all_the_controls_are_unchanged_by_the_subsets(simc):
    counts, names, tables, panel, units = simc
    run = lambda **kw: estimate.estimate_sample(counts, panel, units, region_tables=tables, regions=names, L=L, **kw)
    a, b = run(), run(control_subsets={"lite": _half(names)})
    assert a["controls_used"] == sum(role == "control" for _, role in names) and "controls_subset" not in a
    assert json.dumps(a, default=str, sort_keys=True) == json.dumps(b, default=str, sort_keys=True)


def test_a_control_outside_the_bundle_is_refused_whatever_subsets_are_named(simc):
    counts, names, tables, panel, units = simc
    bad = copy.deepcopy(counts)
    bad["regions"][0]["name"] = "chr1:1-2"
    everything = frozenset(r["name"] for r in bad["regions"] if r["role"] == "control")
    with pytest.raises(ValueError, match="control regions not shared with the bundle"):
        estimate.estimate_sample(bad, panel, units, region_tables=tables, regions=names, L=L, control_subsets={"x": everything})


def _bundle_copy(tmp_path) -> Path:
    src = ROOT / "resources" / "GRCh38"
    d = tmp_path / "GRCh38"
    d.mkdir()
    for f in ("bundle.json", "controls.fa.gz", "controls.bed"):
        (d / f).symlink_to(src / f)
    return d


def test_the_bundle_names_its_control_subsets(tmp_path):
    d = _bundle_copy(tmp_path)
    B = resources.Bundle(d)
    assert B.control_subsets() == {}
    bed = [line.split("\t") for line in (d / "controls.bed").read_text().splitlines()]
    ctrl = [p for p in bed if p[3] == "control"][:3]
    test = next(p for p in bed if p[3].startswith("test:"))
    (d / "controls.lite3.bed").write_text("".join("\t".join(p) + "\n" for p in ctrl + [test]))
    got = resources.Bundle(d).control_subsets()
    assert got == {"lite3": frozenset(f"{p[0]}:{p[1]}-{p[2]}" for p in ctrl)}
    (d / "controls.bad.bed").write_text("chr1\t1\t2\tcontrol\n")
    with pytest.raises(ValueError, match="controls.bad.bed: a control subset must name control regions of the bundle"):
        resources.Bundle(d).control_subsets()
    (d / "controls.bad.bed").unlink()
    meta = json.loads((d / "bundle.json").read_text())
    (d / "bundle.json").unlink()
    (d / "bundle.json").write_text(json.dumps({**meta, "control_subsets": {"gone": "controls.gone.bed"}}))
    with pytest.raises(FileNotFoundError, match="control subset gone"):
        resources.Bundle(d).control_subsets()


def _experimental(tmp_path, counts, pc, name="NEWPOS", seq=None, record_hash=True, k=31):
    """An experimental panel beside a units directory, as resources/experimental/candidates lays them out."""
    cand = tmp_path / "candidates"
    (cand / "units").mkdir(parents=True, exist_ok=True)
    lines = ["##ngs-dose-panel v1", f"##k={k}",
             f"##class\tid=0\tname={name}\tkind=positional\tlength={pc.length}\tcircular={int(pc.circular)}",
             "#kmer\tclass\tpos\tstrand"] + [f"{'A' * k}\t0\t{p}\t+" for p in pc.kmer_pos]
    panel_path = cand / f"{name.lower()}.k31.panel.tsv.gz"
    panel_path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))
    (cand / "units" / f"{name}.fa").write_text(f">{name}\n{seq}\n")
    c = copy.deepcopy(counts)
    c["classes"][0]["name"] = name
    c["panel_sha256"] = ["0" * 64] + ([hashlib.sha256(panel_path.read_bytes()).hexdigest()] if record_hash else [])
    return c, resources.ExperimentalUnits([cand / "units"])


def test_an_experimental_positional_class_is_estimated_from_its_unit(simc, tmp_path):
    counts, names, tables, panel, units = simc
    pc = panel.classes["unit"]
    c, ex = _experimental(tmp_path, counts, pc, seq=units["unit"])
    run = lambda x, **kw: estimate.estimate_sample(x, Panel(31, {}), {}, region_tables=tables, regions=names, L=L, **kw)
    before = run(c)                                                   # today's behaviour without experimental units
    assert before["skipped_classes"][0]["reason"] == "not in the bundle panel" and "NEWPOS" not in before["classes"]
    r = run(c, experimental=ex)
    v = r["classes"]["NEWPOS"]
    assert r["skipped_classes"] == [] and v["status"] == "experimental" and v["estimator"] == estimate.EXPERIMENTAL_ESTIMATOR
    assert v["cn_basis"] == "all" and v["cn"] == v["cn_all"] and abs(v["cn"] / 100 - 1) < 0.05
    assert np.isnan(v["cn_anchor"]) and v["n_anchor"] is None          # no anchor windows, so no anchor figure
    assert 0 < v["cn_se_rel"] < 0.05 and v["unit_source"].endswith("units/NEWPOS.fa")
    assert len(v["windows"]) == len(estimate.estimate_sample(counts, panel, units, region_tables=tables, regions=names,
                                                             L=L)["classes"]["unit"]["windows"])
    # the same class through the bundle: the same windows and all-window estimate, and its own headline rule
    ref = estimate.estimate_sample(counts, panel, units, region_tables=tables, regions=names, L=L)["classes"]["unit"]
    assert ref["cn_all"] == v["cn_all"] and ref["status"] == "ok"


def test_an_experimental_unit_that_cannot_be_used_says_why(simc, tmp_path):
    counts, names, tables, panel, units = simc
    pc = panel.classes["unit"]
    run = lambda x, ex: estimate.estimate_sample(x, Panel(31, {}), {}, region_tables=tables, regions=names, L=L, experimental=ex)
    c, ex = _experimental(tmp_path / "a", counts, pc, seq=units["unit"], record_hash=False)
    r = run(c, ex)
    assert "NEWPOS" not in r["classes"] and "is one the counts were made with" in r["skipped_classes"][0]["reason"]
    assert r["skipped_classes"][0]["reason"].startswith("not in the bundle panel; ")
    c, ex = _experimental(tmp_path / "b", counts, pc, seq=units["unit"][:-10])
    assert "has 2990 bp, the panel says 3000" in run(c, ex)["skipped_classes"][0]["reason"]
    c, ex = _experimental(tmp_path / "c", counts, pc, seq=units["unit"])
    c["classes"][0]["panel_kmers"] = len(pc.kmer_pos) - 1
    v = run(c, ex)["classes"]["NEWPOS"]
    assert v["status"] == "panel_mismatch" and "the experimental panel" in v["reason"] and np.isnan(v["cn"])
    c, ex = _experimental(tmp_path / "d", counts, pc, seq=units["unit"], k=25)
    assert "has k=25" in run(c, ex)["skipped_classes"][0]["reason"]
    # a unit of another name: the class stays skipped, and the reason says where a unit was looked for
    c, ex = _experimental(tmp_path / "e", counts, pc, seq=units["unit"])
    c["classes"][0]["name"] = "OTHER"
    assert run(c, ex)["skipped_classes"][0]["reason"] == f"not in the bundle panel; no experimental unit OTHER.fa in {ex.dirs[0]}"
    # paths are given relative to the install root when the units lie under it
    rooted = resources.ExperimentalUnits(ex.dirs, root=tmp_path / "e")
    assert rooted.missing("OTHER") == "no experimental unit OTHER.fa in candidates/units"
    c["classes"][0]["name"] = "NEWPOS"
    v = run(c, rooted)["classes"]["NEWPOS"]
    assert v["unit_source"] == "candidates/units/NEWPOS.fa"
    assert v["unit_sha256"] == hashlib.sha256(units["unit"].upper().encode()).hexdigest()


def test_experimental_units_are_looked_for_beside_the_bundle_and_in_the_environment(tmp_path, monkeypatch):
    monkeypatch.delenv(resources.EXTRA_UNITS_ENV, raising=False)
    assert resources.experimental_unit_dirs(ROOT / "resources" / "GRCh38") == [
        ROOT / "resources" / "experimental" / "candidates" / "units"]
    monkeypatch.setenv(resources.EXTRA_UNITS_ENV, f"{tmp_path / 'a'}{os.pathsep}{tmp_path / 'b'}")
    dirs = resources.Bundle(ROOT / "resources" / "GRCh38").experimental().dirs
    assert dirs[1:] == [tmp_path / "a", tmp_path / "b"]
    (tmp_path / "b").mkdir()
    (tmp_path / "b" / "X.fa.gz").write_bytes(gzip.compress(b">X\nACGT\n"))
    assert resources.ExperimentalUnits(dirs).unit_file("X") == tmp_path / "b" / "X.fa.gz"


def test_the_all_window_headline_is_a_choice_of_the_caller(simc):
    counts, names, tables, panel, units = simc
    curve = estimate.gcmodel.fit_gc_curve(counts["gc_tables"][0]["n"], counts["gc_tables"][0]["o"], L)
    cls = counts["classes"][0]
    a = estimate.estimate_positional(cls, panel.classes["unit"], units["unit"], 31, curve, 150)
    b = estimate.estimate_positional(cls, panel.classes["unit"], units["unit"], 31, curve, 150, headline="all")
    assert a["cn_basis"] == "anchor" and b["cn_basis"] == "all" and b["cn"] == a["cn_all"] == b["cn_all"]
    with pytest.raises(ValueError, match="headline"):
        estimate.estimate_positional(cls, panel.classes["unit"], units["unit"], 31, curve, 150, headline="best")


@pytest.mark.skipif(not BIN.exists(), reason="needs the engine (cargo build --release, or NGSDOSE_BIN)")
def test_the_command_line_estimates_a_lite_fetch_and_experimental_classes(tmp_path):
    """`ngsdose estimate` as a pipeline runs it: a fetch made with the shipped lighter controls file,
    a scan with a candidate panel, a file whose controls are a subset nobody named, and a positional
    class without a unit anywhere."""
    cand = ROOT / "resources" / "experimental" / "candidates" / "macrosatellites.k31.panel.tsv.gz"
    if not (BUNDLE / "controls.lite200.fa.gz").exists() or not cand.exists():
        pytest.skip("needs controls.lite200.fa.gz and the macrosatellites candidate panel")
    count = lambda *a: subprocess.run([str(BIN), "count", "-i", str(BAM), "-p", str(BUNDLE / "panel.k31.tsv.gz"), "-@", "2", *map(str, a)],
                                      check=True, capture_output=True)
    count("-c", BUNDLE / "controls.lite200.fa.gz", "-m", "fetch", "--sinks", BUNDLE / "sinks.bed", "-o", tmp_path / "lite.json.gz")
    count("-p", cand, "-c", BUNDLE / "controls.fa.gz", "-m", "scan", "-o", tmp_path / "cand.json.gz")
    c = json.load(gzip.open(tmp_path / "lite.json.gz", "rt"))
    first = next(i for i, r in enumerate(c["regions"]) if r.get("role", "control") == "control")
    json.dump({**c, "sample": "DROPPED", "regions": c["regions"][:first] + c["regions"][first + 1:]},
              gzip.open(tmp_path / "dropped.json.gz", "wt"))
    c = json.load(gzip.open(tmp_path / "cand.json.gz", "rt"))
    json.dump({**c, "sample": "CAND"}, gzip.open(tmp_path / "cand.json.gz", "wt"))       # the fixture's sample is NA12878 in both
    json.dump({**c, "sample": "NOUNIT", "classes": [dict(x, name="NOUNIT") if x["name"] == "DXZ4" else x for x in c["classes"]]},
              gzip.open(tmp_path / "nounit.json.gz", "wt"))
    # the contract check alone: the named subset is not fatal, the unnamed one is
    B = resources.Bundle(BUNDLE)
    assert not [m for lv, m in contract.issues(json.load(gzip.open(tmp_path / "lite.json.gz", "rt")), B) if lv == "fatal"]
    assert any("no named subset" in m for lv, m in contract.issues(json.load(gzip.open(tmp_path / "dropped.json.gz", "rt")), B)
               if lv == "fatal")
    env = {**os.environ, "PYTHONPATH": str(ROOT) + os.pathsep + os.environ.get("PYTHONPATH", "")}
    env.pop(resources.EXTRA_UNITS_ENV, None)
    files = [tmp_path / f"{n}.json.gz" for n in ("lite", "cand", "dropped", "nounit")]
    r = subprocess.run([sys.executable, "-m", "ngsdose", "estimate", *map(str, files), "-r", str(BUNDLE), "-t", str(tmp_path / "t.tsv"),
                        "-o", str(tmp_path / "est")], capture_output=True, text=True, env=env)
    assert r.returncode != 0 and "1 of 4 counts files failed" in r.stderr, r.stderr
    assert "FAILED" in r.stderr and "DROPPED: made with a different controls file (a subset of them that no named subset" in r.stderr
    assert "different sets of control regions (all: 2 file(s), lite200: 1 file(s))" in r.stderr
    assert "positional class DXZ4 is not estimated" not in r.stderr and "not estimated (NA in the table): DXZ4" not in r.stderr
    assert ("NOUNIT: positional class NOUNIT is not estimated: the bundle has no panel entry for it; "
            "no experimental unit NOUNIT.fa in resources/experimental/candidates/units") in r.stderr
    with open(tmp_path / "t.tsv") as fh:
        rows = {row["sample"]: row for row in csv.DictReader(fh, delimiter="\t")}
    lite, scan, nounit = rows["NA12878"], rows["CAND"], rows["NOUNIT"]
    assert lite["controls_used"] == "200" and lite["controls_subset"] == "lite200" and lite["mode"] == "fetch"
    assert scan["controls_used"] == "800" and scan["controls_subset"] == "NA"
    assert scan["DXZ4.status"] == "experimental" and scan["DXZ4.cn_single"] != "NA" and "DXZ4.cn_se_rel" in scan
    assert scan["DXZ4.cn_anchor"] == "NA" and scan["DXZ4.n_anchor"] == "NA" and scan["rDNA45S.cn_anchor"] != "NA"
    assert scan["rDNA45S.status"] == "ok" and "rDNA45S.cn_se_rel" not in scan          # bundle classes: no such column
    assert nounit["NOUNIT.status"] == "skipped: not in the bundle panel; no experimental unit NOUNIT.fa in resources/experimental/candidates/units"
    v = json.load(gzip.open(tmp_path / "est" / "cand.estimate.json.gz", "rt"))["classes"]["DXZ4"]
    assert v["unit_source"] == "resources/experimental/candidates/units/DXZ4.fa" and len(v["unit_sha256"]) == 64
