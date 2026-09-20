"""Real 1000 Genomes reads through the engine and the estimator.

tests/data/NA12878.subsample.bam is a 2% template subsample of NA12878 (NYGC 30x, GRCh38)
restricted to the bundle's control regions and class sinks. The assertions are the claims the
method makes about itself: known-copy-number sequence comes out at its known copy number."""
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

from ngsdose import estimate, io, resources

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("NGSDOSE_BIN", ROOT / "target" / "release" / "ngs-dose"))
BAM = ROOT / "tests" / "data" / "NA12878.subsample.bam"
BUNDLE = resources.Bundle(ROOT / "resources" / "GRCh38")

pytestmark = pytest.mark.skipif(not BIN.exists(), reason="needs target/release/ngs-dose (cargo build --release)")


def count(tmp, mode):
    out = tmp / f"{mode}.json.gz"
    cmd = [str(BIN), "count", "-i", str(BAM), "-p", str(BUNDLE.panel), "-c", str(BUNDLE.controls), "-m", mode, "-@", "2", "-o", str(out)]
    if mode == "fetch":
        cmd += ["--sinks", str(BUNDLE.sinks)]
    subprocess.run(cmd, check=True, capture_output=True)
    return io.load_counts(out)


@pytest.fixture(scope="module")
def counts(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("real")
    return {m: count(tmp, m) for m in ("scan", "fetch")}


@pytest.fixture(scope="module")
def result(counts):
    c = counts["fetch"]
    L = estimate.nearest_table(c, None)["l"]
    return estimate.estimate_sample(c, io.load_panel(BUNDLE.panel), BUNDLE.units(), BUNDLE.features(),
                                    region_tables=estimate.control_region_tables(BUNDLE.controls, L))


def test_fetch_retrieves_what_scan_finds(counts):
    s, f = counts["scan"], counts["fetch"]
    assert s["sample"] == f["sample"] == "NA12878"
    assert s["gc_tables"] == f["gc_tables"]
    for a, b in zip(s["classes"], f["classes"]):
        assert a["fwd"] == b["fwd"] and a["rev"] == b["rev"], a["name"]


def test_library_properties(counts, result):
    assert counts["fetch"]["read_length_mode"] == 150
    assert 400 < result["insert_median"] < 470
    assert 0.08 < result["ctrl_dup_frac"] < 0.13
    # the duplicate flag is set far less often inside the collapsed rDNA than in single-copy DNA,
    # which is why flagged reads are counted rather than dropped (docs/DESIGN.md, finding 1)
    assert result["classes"]["rDNA45S"]["dup_flag_frac"] < 0.7 * result["ctrl_dup_frac"]


def test_known_copy_number_sequence(result):
    t = result["truth_regions"]
    assert abs(t["auto"]["cn"] - 2.0) < 0.08            # held-out autosomal sequence: 2 copies
    assert 1.8 < t["chrX"]["cn"] < 2.05                 # NA12878 is female; late-replicating Xi reads slightly low
    dj = result["classes"]["DJ"]
    assert abs(dj["cn"] - 10.0) < 0.8, dj["cn"]         # distal junction: one copy per acrocentric short arm


def test_rdna(result):
    k = result["classes"]["rDNA45S"]
    assert k["cn_basis"] == "anchor"
    assert abs(k["cn"] / 508.0 - 1) < 0.06, k["cn"]     # 508 from the full 37x data
    f = k["features"]
    # the GC-rich 28S reads low against 18S in every published 1000 Genomes analysis (~0.77)
    assert 0.65 < f["28S"]["cn"] / f["18S"]["cn"] < 0.9
    five = result["classes"]["rDNA5S"]
    assert five["cn_basis"] == "all" and 150 < five["cn"] < 280


def test_counts_do_not_depend_on_threads_or_on_co_loaded_panels(tmp_path):
    """The same reads give the same numbers with 1 or 4 threads, and a class's counts are
    unchanged when an unrelated panel (here the experimental satellite panel) is loaded as well."""
    sat = ROOT / "resources" / "experimental" / "satellites.CHM13v2.k31.panel.tsv.gz"
    base = [str(BIN), "count", "-i", str(BAM), "-c", str(BUNDLE.controls), "--sinks", str(BUNDLE.sinks), "-m", "fetch"]
    runs = {"t1": ["-p", str(BUNDLE.panel), "-@", "1"], "t4": ["-p", str(BUNDLE.panel), "-@", "4"],
            "sat": ["-p", str(BUNDLE.panel), "-p", str(sat), "-@", "4"]}
    out = {}
    for name, extra in runs.items():
        subprocess.run(base + extra + ["-o", str(tmp_path / f"{name}.json")], check=True, capture_output=True)
        out[name] = io.load_counts(tmp_path / f"{name}.json")
    strip = lambda d: {k: v for k, v in d.items() if k != "elapsed_sec"}
    assert strip(out["t1"]) == strip(out["t4"])
    with_sat = {c["name"]: c for c in out["sat"]["classes"]}
    for c in out["t4"]["classes"]:
        assert (c["fwd"], c["rev"], c["reads"]) == (with_sat[c["name"]]["fwd"], with_sat[c["name"]]["rev"], with_sat[c["name"]]["reads"])
    assert len(with_sat) == 13 and len(out["sat"]["panel_sha256"]) == 2
    # and a whole-file scan with every panel loaded - classes, placements with their census, contig tallies -
    # is the same file whatever the number of threads
    exp = ROOT / "resources" / "experimental"
    scans = []
    for threads in ("1", "4"):
        path = tmp_path / f"scan{threads}.json"
        subprocess.run([str(BIN), "count", "-i", str(BAM), "-c", str(BUNDLE.controls), "-m", "scan", "-p", str(BUNDLE.panel), "-p", str(sat),
                        "-p", str(exp / "telomere.k31.panel.tsv.gz"), "-@", threads, "-o", str(path)], check=True, capture_output=True)
        scans.append(strip(io.load_counts(path)))
    assert scans[0] == scans[1] and len(scans[0]["classes"]) == 14


def test_counts_record_what_they_were_made_with(counts):
    import hashlib
    c = counts["fetch"]
    assert c["panel_sha256"] == [hashlib.sha256(BUNDLE.panel.read_bytes()).hexdigest()]
    assert c["controls_sha256"] == hashlib.sha256(BUNDLE.controls.read_bytes()).hexdigest()
    assert c["sinks_sha256"] == hashlib.sha256(BUNDLE.sinks.read_bytes()).hexdigest()
    assert counts["scan"]["sinks_sha256"] is None
    assert c["engine_version"] and c["engine_build"]          # the commit the engine was built from ("dev" outside CI's image build)


def test_sex_chromosome_truth_and_episome_dosage(result):
    """chrY is a known truth (NA12878 is female: none), chrM and chrEBV are dosages: copies per
    cell of the mitochondrial genome and of the EBV episome this LCL was made with."""
    t = result["truth_regions"]
    assert t["chrY"]["role"] == "test" and t["chrY"]["cn"] < 0.02
    assert t["chrM"]["role"] == t["chrEBV"]["role"] == "dosage"
    assert abs(t["chrM"]["cn"] / 850.0 - 1) < 0.05, t["chrM"]["cn"]      # 850 and 71.4 from the full 37x data
    assert abs(t["chrEBV"]["cn"] / 71.4 - 1) < 0.08, t["chrEBV"]["cn"]


def test_contigs_and_end_of_file_marker_are_recorded(counts):
    s, f = counts["scan"], counts["fetch"]
    assert s["eof_marker"] == f["eof_marker"] == "present"
    # scan: every mapped primary read is on exactly one contig of the table
    assert sum(c["reads"] for c in s["contigs"]) == s["primary"] - s["unmapped"]
    # fetch sees only what it asked for, so it reports lengths alone
    assert all("reads" not in c for c in f["contigs"])
    lengths = {c["name"]: c["len"] for c in f["contigs"]}
    assert {"chr1", "chr21", "chrY", "chrM", "chrEBV"} <= set(lengths)
    assert all(lengths[n] == ln for n, ln in BUNDLE.contig_lengths().items() if n in lengths)


def test_a_file_aligned_to_another_build_is_refused(counts):
    import copy
    c = copy.deepcopy(counts["fetch"])
    next(x for x in c["contigs"] if x["name"] == "chr1")["len"] = 249250621           # hg19's chr1
    with pytest.raises(ValueError, match="different reference build"):
        estimate.estimate_sample(c, io.load_panel(BUNDLE.panel), BUNDLE.units(), contig_lengths=BUNDLE.contig_lengths())


def test_truncated_input_is_refused(tmp_path):
    """A BAM/CRAM cut short still decodes. What is lost is whatever sorted last - in a GRCh38 file
    chr21, chr22 and the unplaced contigs, which is where the rDNA is - so the estimate would be
    wrong and look fine. The engine checks the end-of-file marker before it reads anything."""
    cut = tmp_path / "cut.bam"
    cut.write_bytes(BAM.read_bytes()[:-28])                                            # the BGZF EOF block is 28 bytes
    out = tmp_path / "cut.json.gz"
    cmd = [str(BIN), "count", "-i", str(cut), "-p", str(BUNDLE.panel), "-c", str(BUNDLE.controls), "-m", "scan", "-o", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode != 0 and "truncated" in r.stderr
    assert not out.exists() and not list(tmp_path.glob("*.partial"))
    subprocess.run(cmd + ["--allow-truncated"], check=True, capture_output=True)
    assert io.load_counts(out)["eof_marker"] == "absent"


@pytest.mark.skipif(not __import__("shutil").which("samtools"), reason="needs samtools")
def test_a_reference_without_a_contig(tmp_path):
    """A reference without the EBV decoy is still GRCh38: regions that are not controls may sit on
    contigs a file lacks, and come out as missing rather than as zero. Controls may not."""
    def without(contig, name):
        sam = subprocess.run(["samtools", "view", "-h", str(BAM)], check=True, capture_output=True, text=True).stdout
        keep = "".join(l + "\n" for l in sam.splitlines() if f"\t{contig}\t" not in l + "\t" and f"SN:{contig}\t" not in l)
        bam = tmp_path / name
        subprocess.run(["samtools", "view", "-b", "-o", str(bam), "-"], input=keep, check=True, text=True, capture_output=True)
        return bam
    out = tmp_path / "noebv.json.gz"
    base = [str(BIN), "count", "-p", str(BUNDLE.panel), "-c", str(BUNDLE.controls), "-m", "scan", "-o", str(out), "-i"]
    subprocess.run(base + [str(without("chrEBV", "noebv.bam"))], check=True, capture_output=True)
    c = io.load_counts(out)
    assert [r["name"] for r in c["regions"] if r.get("absent")] == ["chrEBV:102000-122000"]
    L = estimate.nearest_table(c, None)["l"]
    res = estimate.estimate_sample(c, io.load_panel(BUNDLE.panel), BUNDLE.units(), BUNDLE.features(),
                                   region_tables=estimate.control_region_tables(BUNDLE.controls, L))
    assert "chrEBV" not in res["truth_regions"] and res["truth_regions"]["chrM"]["cn"] > 100
    r = subprocess.run(base + [str(without("chr7", "nochr7.bam"))], capture_output=True, text=True)
    assert r.returncode != 0 and "wrong reference build" in r.stderr


def test_classes_of_a_co_loaded_panel_reach_the_estimate(tmp_path):
    """A compositional class needs nothing from the bundle, so a scan made with extra panels (the
    satellite families, the telomeric repeat) is estimated in full by the ordinary `ngsdose estimate`."""
    exp = ROOT / "resources" / "experimental"
    out = tmp_path / "merged.json.gz"
    subprocess.run([str(BIN), "count", "-i", str(BAM), "-c", str(BUNDLE.controls), "-m", "scan", "-p", str(BUNDLE.panel),
                    "-p", str(exp / "satellites.CHM13v2.k31.panel.tsv.gz"), "-p", str(exp / "telomere.k31.panel.tsv.gz"),
                    "-o", str(out)], check=True, capture_output=True)
    c = io.load_counts(out)
    res = estimate.estimate_sample(c, io.load_panel(BUNDLE.panel), BUNDLE.units())
    satellites = {"HSat1A", "HSat1B", "HSat2", "HSat3", "bSat", "aSatHOR", "ACRO", "SST1", "CER", "SATR"}
    assert {"rDNA45S", "rDNA5S", "DJ", "TEL"} | satellites == set(res["classes"]) and len(c["panel_sha256"]) == 3
    assert all(res["classes"][x]["kind"] == "compositional" and res["classes"][x]["mass_Mb"] >= 0 for x in satellites | {"TEL"})
    # the telomere panel is the six canonical 31-mers of (TTAGGG)n; every class carries the histogram of the
    # share of each read's k-mers that hit, which is what a TelSeq-style threshold is applied to afterwards
    tel = next(x for x in c["classes"] if x["name"] == "TEL")
    assert tel["panel_kmers"] == 6 and len(tel["hit_frac"]) == 11 and sum(tel["hit_frac"]) == tel["reads"]


def test_a_stalled_input_ends_with_status_75(tmp_path):
    """A dead HTTPS connection does not fail, it waits. A FIFO that nobody writes to behaves the
    same way (open blocks for ever): the watchdog must end the process with EX_TEMPFAIL so that
    the caller can retry, and leave no output behind."""
    fifo = tmp_path / "stalled.bam"
    os.mkfifo(fifo)
    out = tmp_path / "stalled.json.gz"
    r = subprocess.run([str(BIN), "count", "-i", str(fifo), "-p", str(BUNDLE.panel), "-c", str(BUNDLE.controls), "-m", "scan",
                        "--stall-timeout", "2", "-o", str(out)], capture_output=True, text=True, timeout=60)
    assert r.returncode == 75 and "nothing read for 2 s" in r.stderr
    assert not list(tmp_path.glob("stalled.json*"))


def test_counts_made_before_a_truth_set_existed_are_still_usable(counts):
    """Regions are matched to the bundle by name. Counts that lack some non-control regions (made
    with an earlier bundle, before chrY or the dosage regions were added) are estimated, and the
    missing sets are simply not reported; an unknown region or a different control set is refused."""
    import copy
    L = estimate.nearest_table(counts["fetch"], None)["l"]
    tables, regions = BUNDLE.region_tables(L), BUNDLE.regions()
    args = (io.load_panel(BUNDLE.panel), BUNDLE.units(), BUNDLE.features())
    full = estimate.estimate_sample(counts["fetch"], *args, region_tables=tables, regions=regions)
    older = copy.deepcopy(counts["fetch"])
    older["regions"] = [r for r in older["regions"] if r.get("label") not in ("chrY", "chrM", "chrEBV")]
    res = estimate.estimate_sample(older, *args, region_tables=tables, regions=regions)
    assert set(res["truth_regions"]) == {"auto", "chrX"}
    assert res["classes"]["rDNA45S"]["cn"] == full["classes"]["rDNA45S"]["cn"]
    assert res["truth_regions"]["auto"]["cn"] == full["truth_regions"]["auto"]["cn"]
    with pytest.raises(ValueError, match="different controls file"):          # row-for-row matching cannot take it
        estimate.estimate_sample(older, *args, region_tables=tables)
    fewer_controls = copy.deepcopy(counts["fetch"])
    del fewer_controls["regions"][0]
    with pytest.raises(ValueError, match="different controls file"):
        estimate.estimate_sample(fewer_controls, *args, region_tables=tables, regions=regions)
    renamed = copy.deepcopy(counts["fetch"])
    renamed["regions"][-1]["name"] = "chrEBV:1-2"
    with pytest.raises(ValueError, match="different controls file"):
        estimate.estimate_sample(renamed, *args, region_tables=tables, regions=regions)


def test_placements_carry_a_census_of_every_read_in_their_bin(counts):
    """Where class reads were aligned is recorded on the 1-kb grid of `mosdepth --by 1000`, with
    every mapped primary read of the bin beside them (and how many are flagged duplicate): what a
    depth tool sees there, and how much of it is the class. Inside the sinks, a targeted fetch
    must see exactly what the whole-file scan saw."""
    s, f = counts["scan"], counts["fetch"]
    assert s["placement_bin"] == f["placement_bin"] == 1000 and s["placement_bin_compositional"] == 10000
    placed = [p for p in s["placements"] if p["contig"] != "*"]
    assert placed and all(p["start"] % 1000 == 0 and 0 <= p["dup"] <= p["all"] for p in placed)
    assert all("all" not in p for p in s["placements"] if p["contig"] == "*")
    # class reads are among the bin's reads (bar the odd unmapped mate parked at the position)
    assert sum(p["reads"] for p in placed) <= sum(p["all"] for p in placed) and all(p["reads"] <= p["all"] + 3 for p in placed)
    # rDNA bins are nearly pure, which is what makes a depth proxy conceivable at all
    big = [p for p in placed if p["class"] == "rDNA45S" and p["reads"] >= 100]
    assert len(big) >= 20 and np.median([p["reads"] / p["all"] for p in big]) > 0.9
    # scan and fetch agree bin for bin wherever a bin lies wholly inside a sink interval
    sinks = [l.split("\t") for l in open(BUNDLE.sinks)]
    inside = lambda p: any(c == p["contig"] and int(a) + 1000 <= p["start"] and p["start"] + 2000 <= int(b) for c, a, b, _ in sinks)
    key = lambda p: (p["class"], p["contig"], p["start"])
    fetch = {key(p): p for p in f["placements"]}
    n = 0
    for p in placed:
        if inside(p):
            n += 1
            assert fetch[key(p)] == p, (p, fetch.get(key(p)))
    assert n > 100


def test_sinks_are_decided_at_ten_kb_and_trimmed_on_a_finer_grid(counts):
    """A finer placement grid must make sinks tighter, not more numerous: three stray reads in
    one kilobase do not become a sink because the grid got finer."""
    import copy
    from ngsdose import sinks
    fine = counts["scan"]
    coarse = copy.deepcopy(fine)
    coarse["placement_bin"] = 10000
    merged = {}
    for p in coarse["placements"]:
        if p["contig"] != "*" and p["class"] in ("rDNA45S", "rDNA5S", "DJ"):
            k = (p["class"], p["contig"], p["start"] // 10000 * 10000)
            merged[k] = merged.get(k, 0) + p["reads"]
    coarse["placements"] = [dict(**{"class": c}, contig=g, start=st, reads=n) for (c, g, st), n in merged.items()]
    rows_f, cap_f = sinks.learn([fine])
    rows_c, cap_c = sinks.learn([coarse])
    span = lambda rows, cls: sum(e - s0 for _, s0, e, c in rows if c == cls)
    for cls in ("rDNA45S", "DJ"):
        assert span(rows_f, cls) <= span(rows_c, cls)
        # the decision is the same one, so every fine interval lies inside a coarse interval
        for g, a, b, c in rows_f:
            if c == cls:
                assert any(g == g2 and a2 <= a and b <= b2 for g2, a2, b2, c2 in rows_c if c2 == cls), (g, a, b)
    # where reads are dense (45S, even in a 2% subsample) trimming loses nothing of note
    assert cap_f["rDNA45S"]["NA12878"] > cap_c["rDNA45S"]["NA12878"] - 0.002
    # pooling scans made on different grids works
    rows_both, _ = sinks.learn([fine, coarse])
    assert span(rows_both, "rDNA45S") >= span(rows_f, "rDNA45S")
