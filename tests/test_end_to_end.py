"""Engine + estimator on a simulated genome with known truth (needs the built binary and samtools)."""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from ngsdose import estimate, io

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("NGSDOSE_BIN", ROOT / "target" / "release" / "ngs-dose"))

pytestmark = pytest.mark.skipif(not BIN.exists() or shutil.which("samtools") is None,
                                reason="needs target/release/ngs-dose (cargo build --release) and samtools")


@pytest.fixture(scope="module")
def sim(tmp_path_factory):
    from simulate_bam import simulate
    d = tmp_path_factory.mktemp("sim")
    s = simulate(d)
    run = lambda *a: subprocess.run([str(BIN), *map(str, a)], check=True, cwd=d, capture_output=True, text=True)
    run("panel", "-m", "classes.tsv", "-k", "31", "-b", "background.fa", "-o", "panel.tsv.gz")
    run("controls", "-b", "controls.bed", "-T", "ref.fa", "--flank", "700", "-o", "controls.fa.gz")
    common = ["-i", "sim.bam", "-p", "panel.tsv.gz", "-c", "controls.fa.gz", "--l-grid", "100,200,300,400"]
    run("count", *common, "-m", "scan", "-@", "2", "-o", "scan.json")
    run("count", *common, "-m", "fetch", "--sinks", "sinks.bed", "-@", "2", "-o", "fetch.json.gz")
    s["scan"], s["fetch"] = io.load_counts(d / "scan.json"), io.load_counts(d / "fetch.json.gz")
    return s


def test_scan_and_fetch_agree(sim):
    a, b = sim["scan"], sim["fetch"]
    assert a["sample"] == b["sample"] == "simulated"
    assert a["gc_tables"] == b["gc_tables"]
    assert a["classes"][0]["fwd"] == b["classes"][0]["fwd"] and a["classes"][0]["rev"] == b["classes"][0]["rev"]
    assert a["ctrl_reads"] == b["ctrl_reads"] > 50_000


def test_every_class_read_is_recovered_duplicates_included(sim):
    c = sim["scan"]["classes"][0]
    assert c["dup_flagged"] > 0, "flagged duplicates must be counted, not dropped"
    assert c["reads"] >= 0.995 * sim["n_class_reads"]
    assert c["reads"] <= sim["n_class_reads"]


def test_fragment_ends_are_placed_on_the_right_strand_and_bin(sim):
    """BAM stores reverse-strand alignments reverse-complemented; the engine must undo that, or
    every reverse read is booked as a forward read one read-length away from its true 5' end."""
    c = sim["scan"]["classes"][0]
    fwd, rev = np.array(c["fwd"]), np.array(c["rev"])
    assert abs(fwd.sum() - rev.sum()) < 0.02 * fwd.sum()
    # sequencing errors remove a handful of reads and shift none: bins agree to within those losses
    assert np.abs(fwd - sim["truth_fwd"]).sum() <= 0.005 * fwd.sum()
    assert np.abs(rev - sim["truth_rev"]).sum() <= 0.005 * rev.sum()


def test_copy_number_recovered_under_gc_bias(sim):
    panel = io.load_panel(sim["dir"] / "panel.tsv.gz")
    r = estimate.estimate_sample(sim["fetch"], panel, {"unit": sim["unit"]}, window=300, min_kmers=20)
    k = r["classes"]["unit"]
    assert r["gc_L"] == 300 and r["read_length"] == 100
    assert abs(k["cn_all"] / sim["copies"] - 1) < 0.04, k["cn_all"]
    assert abs(k["cn_anchor"] / sim["copies"] - 1) < 0.05, k["cn_anchor"]
    # every window of the unit is present at the same copy number; the GC model must make them agree
    cns = np.array([w["cn"] for w in k["windows"] if w["cn"] is not None])
    assert np.abs(cns / sim["copies"] - 1).max() < 0.15, cns
    # an estimator blind to GC is measurably worse on this unit
    reads = sim["fetch"]["classes"][0]["reads"]
    naive = 2 * reads / (2 * len(sim["unit"]) * r["ctrl_rate"])
    assert abs(naive / sim["copies"] - 1) > abs(k["cn_all"] / sim["copies"] - 1)


def test_engine_position_tables_match_python(sim):
    tab = next(t for t in sim["fetch"]["gc_tables"] if t["l"] == 300)
    rt = estimate.control_region_tables(sim["dir"] / "controls.fa.gz", 300)
    assert (rt.sum(0) == np.array(tab["n"])).all()
    assert rt.shape[0] == len(sim["fetch"]["regions"])


def test_compositional_class_mass(sim):
    """A satellite-like family without a unit: diploid mass from reads binned by their own GC."""
    panel = io.load_panel(sim["dir"] / "panel.tsv.gz")
    assert panel.classes["sat"].kind == "compositional"
    r = estimate.estimate_sample(sim["fetch"], panel, {"unit": sim["unit"]}, window=300)
    m = r["classes"]["sat"]
    assert m["gc_supported_fraction"] > 0.99
    assert abs(m["mass_bp"] / sim["sat_mass_bp"] - 1) < 0.08, (m["mass_bp"], sim["sat_mass_bp"])


def test_plan_lists_the_intervals_a_fetch_reads(sim):
    """`ngs-dose plan` writes the fetch plan as BED: the control regions padded by 600 bp and the sinks,
    merged; a file cut out with samtools along it and scanned gives the fetch's counts."""
    d = sim["dir"]
    run = lambda *a: subprocess.run([str(BIN), *map(str, a)], check=True, cwd=d, capture_output=True, text=True)
    run("plan", "-c", "controls.fa.gz", "--sinks", "sinks.bed", "-i", "sim.bam", "-o", "plan.bed")
    got = [(p[0], int(p[1]), int(p[2])) for p in (l.split("\t") for l in (d / "plan.bed").read_text().splitlines())]
    want = {}
    for line in (d / "controls.bed").read_text().splitlines():
        c, s, e = line.split("\t")[:3]
        want.setdefault(c, []).append((max(0, int(s) - 600), int(e) + 600))
    for line in (d / "sinks.bed").read_text().splitlines():
        c, s, e = line.split("\t")[:3]
        want.setdefault(c, []).append((int(s), int(e)))
    expected = []
    for c in sorted(want):
        merged = []
        for s, e in sorted(want[c]):
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        expected += [(c, s, e) for s, e in merged]
    assert got == expected
    assert run("plan", "-c", "controls.fa.gz", "--sinks", "sinks.bed").stdout.splitlines() == (d / "plan.bed").read_text().splitlines()
    subprocess.run(["samtools", "view", "-b", "-M", "-L", "plan.bed", "-o", "cut.bam", "sim.bam"], check=True, cwd=d)
    subprocess.run(["samtools", "index", "-c", "cut.bam"], check=True, cwd=d)
    # the cut, counted in fetch mode with the same bundle, sinks and padding, is the whole file's fetch
    run("count", "-i", "cut.bam", "-p", "panel.tsv.gz", "-c", "controls.fa.gz", "--l-grid", "100,200,300,400", "-m", "fetch", "--sinks", "sinks.bed", "-@", "2", "-o", "cut.json.gz")
    strip = lambda c: {k: v for k, v in c.items() if k not in ("input", "elapsed_sec")}
    assert strip(io.load_counts(d / "cut.json.gz")) == strip(sim["fetch"])
    # counted in scan mode it would pass for a whole-file scan; the sink learner and evaluator see that it is not
    run("count", "-i", "cut.bam", "-p", "panel.tsv.gz", "-c", "controls.fa.gz", "--l-grid", "100,200,300,400", "-m", "scan", "-@", "2", "-o", "cut_scan.json")
    # (the simulated genome is little but controls, so its whole-file scan looks cut too; the cohort's real scans put
    # 0.33% of their primary reads in the controls, this cut puts most of them there)
    from ngsdose import sinks as S
    cut_scan = io.load_counts(d / "cut_scan.json")
    assert S.looks_cut(cut_scan)
    with pytest.raises(ValueError, match="cut along a fetch plan"):
        S.learn([cut_scan])
    assert S.learn([cut_scan], allow_cut=True)[0]
    r = subprocess.run([sys.executable, "-m", "ngsdose", "sinks", str(d / "cut_scan.json")], capture_output=True, text=True)
    assert r.returncode != 0 and "cut along a fetch plan" in r.stderr and "Traceback" not in r.stderr


def test_fetch_refuses_a_loaded_class_without_sinks(sim):
    """A panel class the sinks BED has no interval for would be counted only where its reads fall inside
    other intervals: the engine refuses, or records the gap when told to go on."""
    d = sim["dir"]
    (d / "sinks_nosat.bed").write_text("".join(l + "\n" for l in (d / "sinks.bed").read_text().splitlines() if not l.endswith("\tsat")))
    common = [str(BIN), "count", "-i", "sim.bam", "-p", "panel.tsv.gz", "-c", "controls.fa.gz", "--l-grid", "100,200,300,400", "-m", "fetch", "--sinks", "sinks_nosat.bed", "-@", "2"]
    r = subprocess.run([*common, "-o", "nosat.json"], cwd=d, capture_output=True, text=True)
    assert r.returncode != 0 and "sat" in r.stderr and "--allow-missing-sinks" in r.stderr and not (d / "nosat.json").exists()
    r = subprocess.run([*common, "--allow-missing-sinks", "-o", "nosat.json"], cwd=d, capture_output=True, text=True, check=True)
    assert "WARNING" in r.stderr
    c = io.load_counts(d / "nosat.json")
    assert c["sinks_missing_classes"] == ["sat"]
    sat = {x["name"]: x["reads"] for x in c["classes"]}["sat"]
    assert sat < 0.05 * {x["name"]: x["reads"] for x in sim["scan"]["classes"]}["sat"]          # the undercount the flag admits to
    assert "sinks_missing_classes" not in sim["fetch"]                                            # absent when nothing is missing
