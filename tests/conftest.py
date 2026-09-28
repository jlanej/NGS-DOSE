import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

BIN = Path(os.environ.get("NGSDOSE_BIN", Path(__file__).resolve().parents[1] / "target" / "release" / "ngs-dose"))


@pytest.fixture(scope="session")
def sim(tmp_path_factory):
    """The simulated genome of simulate_bam, with a panel, controls and its scan and fetch counts
    (shared by the end-to-end and the engine input tests)."""
    if not BIN.exists() or shutil.which("samtools") is None:
        pytest.skip("needs target/release/ngs-dose (cargo build --release) and samtools")
    from ngsdose import io
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
