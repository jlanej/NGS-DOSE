"""Engine + estimator on a simulated genome with known truth (needs the built binary and samtools)."""
import json
import os
import shutil
import subprocess
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
