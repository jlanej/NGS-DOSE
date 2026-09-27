"""Sub-options: named subsets of a class's learned sink intervals (resources/experimental/subsets: DXZ1 and
DYZ3 of aSatHOR, DYZ1 of HSat3, DYZ2 of HSat1B). `ngsdose estimate` reports each one's reads and mass from
the class's placements inside its intervals, in scans and in fetches that read them; a fetch that read
only such subsets of a class leaves the class itself unmeasured. The counts are simulated
(test_estimate_robustness._sim) with placements added; the shipped files are checked against the
satellite sinks and the menu."""
import copy
import math
from pathlib import Path

import pytest

from ngsdose import contract, estimate, fetchplan, sinks
from ngsdose.tables import summary_row
from test_estimate_robustness import L, _sim

ROOT = Path(__file__).resolve().parents[1]
SUBSETS = ROOT / "resources" / "experimental" / "subsets"

# sat: 5,000 reads; 1,200 of them unmapped. The placement grid of compositional classes is 10 kb.
PLACED = [("cA", 100_000, 1000), ("cA", 110_000, 500), ("cA", 500_000, 2000), ("cB", 50_000, 300), ("*", 0, 1200)]
OPTS = {"X": estimate.SubOption("X", "sat", [("cA", 100_000, 120_000)]),      # holds 1,500
        "Y": estimate.SubOption("Y", "sat", [("cB", 49_000, 60_000)]),        # the last bin of cB, clipped at 55,000: 300
        "Z": estimate.SubOption("Z", "sat", [("cA", 105_000, 115_000)])}      # no bin wholly inside: 0


@pytest.fixture(scope="module")
def simc():
    return _sim(seed=3)


def counts_of(simc, mode, placed=PLACED, sha=None):
    c = copy.deepcopy(simc[0])
    c.update(mode=mode, placement_bin=1000, placement_bin_compositional=10_000, sinks_sha256=sha,
             contigs=[dict(name="cA", len=1_000_000), dict(name="cB", len=55_000)],
             placements=[{"class": "sat", "contig": ct, "start": s, "reads": n, "all": n, "dup": 0} for ct, s, n in placed])
    assert sum(n for _, _, n in placed) == next(x["reads"] for x in c["classes"] if x["name"] == "sat")
    return c


def run(simc, counts, **kw):
    _, names, tables, panel, units = simc
    return estimate.estimate_sample(counts, panel, units, region_tables=tables, regions=names, L=L, **kw)


def test_a_scan_reports_the_class_reads_inside_each_option_and_their_mass(simc):
    r = run(simc, counts_of(simc, "scan"), sub_options=OPTS)
    per_read = r["classes"]["sat"]["mass_bp"] / 5000
    got = r["sub_options"]
    assert {k: v["reads"] for k, v in got.items()} == {"X": 1500, "Y": 300, "Z": 0}
    assert all(v["status"] == "ok" and v["parent"] == "sat" for v in got.values())
    assert got["X"]["mass_bp"] == pytest.approx(1500 * per_read) and got["X"]["mass_Mb"] == pytest.approx(1500 * per_read / 1e6)
    assert got["Z"]["mass_bp"] == 0.0 and r["classes"]["sat"]["status"] == "ok"
    row = summary_row(r)
    assert row["X.reads"] == 1500 and row["X.status"] == "ok" and row["X.mass_Mb"] == pytest.approx(round(1500 * per_read / 1e6, 4))
    # without sub-options nothing is added: the estimate is what it was
    plain = run(simc, counts_of(simc, "scan"))
    assert "sub_options" not in plain and not any(k.startswith(("X.", "Y.", "Z.")) for k in summary_row(plain))
    assert plain["classes"] == r["classes"]


def test_a_fetch_is_measured_only_where_its_sinks_bed_is_known_and_holds_the_option(simc):
    whole = [("cA", 90_000, 130_000, "sat"), ("cA", 490_000, 520_000, "sat"), ("cB", 40_000, 55_000, "sat")]
    known = {"h1": whole, "h2": whole[:2], "h3": [(c, s, e, "") for c, s, e, _ in whole]}     # h3: rows without a class serve every class
    ok = run(simc, counts_of(simc, "fetch", sha="h1"), sub_options=OPTS, known_sinks=known)["sub_options"]
    assert {k: (v["status"], v["reads"]) for k, v in ok.items()} == {"X": ("ok", 1500), "Y": ("ok", 300), "Z": ("ok", 0)}
    lack = run(simc, counts_of(simc, "fetch", sha="h2"), sub_options=OPTS, known_sinks=known)["sub_options"]
    assert lack["Y"]["status"] == "not_fetched" and math.isnan(lack["Y"]["mass_Mb"]) and lack["Y"]["reads"] is None
    assert "1 of its 1 intervals" in lack["Y"]["reason"] and lack["X"]["status"] == "ok"
    assert run(simc, counts_of(simc, "fetch", sha="h3"), sub_options=OPTS, known_sinks=known)["sub_options"]["Y"]["status"] == "ok"
    unk = run(simc, counts_of(simc, "fetch", sha="other"), sub_options=OPTS, known_sinks=known)["sub_options"]
    assert unk["X"]["status"] == "unverified" and unk["X"]["reads"] == 1500 and math.isnan(unk["X"]["mass_Mb"])
    assert "--fetch-sinks" in unk["X"]["reason"]


def test_a_class_with_options_in_a_fetch_of_an_unknown_bed_is_not_reported(simc):
    """A fetch whose sinks BED is not known may have read the class whole or only at its options (a plan of
    xy_arrays): its placements do not tell (reads land in other fetched intervals too), so the class gets no
    value, status unverified, and the options their reads without a mass. Known, the BED decides; a scan
    is whole; a class without options is not affected."""
    whole = [("cA", 90_000, 130_000, "sat"), ("cA", 490_000, 520_000, "sat"), ("cB", 40_000, 55_000, "sat")]
    r = run(simc, counts_of(simc, "fetch", sha="other"), sub_options=OPTS, known_sinks={"h1": whole})
    sat = r["classes"]["sat"]
    assert sat["status"] == "unverified" and math.isnan(sat["mass_Mb"]) and "--fetch-sinks" in sat["reason"] and "X, Y, Z" in sat["reason"]
    assert {v["status"] for v in r["sub_options"].values()} == {"unverified"} and r["sub_options"]["X"]["reads"] == 1500
    assert summary_row(r)["sat.status"] == "unverified"
    assert run(simc, counts_of(simc, "fetch", sha="h1"), sub_options=OPTS, known_sinks={"h1": whole})["classes"]["sat"]["status"] == "ok"
    assert run(simc, counts_of(simc, "scan"), sub_options=OPTS)["classes"]["sat"]["status"] == "ok"
    others = {"Q": estimate.SubOption("Q", "not_here", [("cA", 0, 10)])}
    assert run(simc, counts_of(simc, "fetch", sha="other"), sub_options=others)["classes"]["sat"]["status"] == "ok"
    assert contract.unverified_parents(counts_of(simc, "fetch", sha="other"), OPTS, {}, {"sat": "subset_only"}) == {}


def test_a_fetch_of_only_the_options_leaves_the_class_unmeasured(simc):
    only = [("cA", 100_000, 3000), ("cA", 110_000, 500), ("cB", 50_000, 300), ("*", 0, 1200)]
    two = {k: OPTS[k] for k in ("X", "Y")}
    bed = {"h": [("cA", 100_000, 120_000, "sat"), ("cB", 49_000, 60_000, "sat")]}
    for sha, known in (("h", bed), ("unknown", {})):          # said by the BED when it is known, by the placements when not
        c = counts_of(simc, "fetch", only, sha)
        assert set(contract.subset_only(c, two, known)) == {"sat"}
        r = run(simc, c, sub_options=two, known_sinks=known)
        assert r["classes"]["sat"]["status"] == "subset_only" and math.isnan(r["classes"]["sat"]["mass_Mb"])
        assert r["sub_options"]["X"]["reads"] == 3500
        assert r["sub_options"]["X"]["status"] == ("ok" if known else "unverified")
    # the option's mass is its reads at the class's mass per read (the class record's GC mix, here the same as the scan's)
    x = run(simc, counts_of(simc, "fetch", only, "h"), sub_options=two, known_sinks=bed)
    per_read = run(simc, counts_of(simc, "scan"))["classes"]["sat"]["mass_bp"] / 5000
    assert x["sub_options"]["X"]["mass_bp"] == pytest.approx(3500 * per_read)
    # a fetch with placements elsewhere read the class's sinks; a scan is never subset-only
    assert contract.subset_only(counts_of(simc, "fetch"), two) == {}
    assert contract.subset_only(counts_of(simc, "scan", only), two) == {}
    assert contract.subset_only(counts_of(simc, "fetch", only, "h"), two, {"h": bed["h"] + [("cA", 500_000, 510_000, "sat")]}) == {}


def test_options_of_classes_the_file_cannot_measure_say_why(simc):
    opts = {"W": estimate.SubOption("W", "other", [("cA", 0, 10)]), "V": estimate.SubOption("V", "unit", [("cA", 0, 10)]),
            "sat": estimate.SubOption("sat", "unit", [("cA", 0, 10)])}
    got = run(simc, counts_of(simc, "scan"), sub_options=opts)["sub_options"]
    assert got["W"]["status"] == "not_counted" and got["V"]["status"] == "not_compositional" and got["sat"]["status"] == "name_clash"
    assert all(math.isnan(v["mass_Mb"]) for v in got.values())
    c = counts_of(simc, "scan")
    del c["placement_bin"]
    assert run(simc, c, sub_options=OPTS)["sub_options"]["X"]["status"] == "no_placements"


def test_sub_option_files_are_read_and_checked(tmp_path):
    assert estimate.load_sub_options(tmp_path / "none") == {}
    (tmp_path / "A1.bed").write_text("# a comment\nc1\t10\t20\tsat\nc1\t10\t20\tsat\nc2\t0\t5\tsat\n")
    got = estimate.load_sub_options(tmp_path)
    assert list(got) == ["A1"] and got["A1"].parent == "sat" and got["A1"].intervals == [("c1", 10, 20), ("c2", 0, 5)]
    (tmp_path / "B1.bed").write_text("c1\t10\t20\tsat\nc1\t30\t40\tother\n")
    with pytest.raises(ValueError, match="names one class"):
        estimate.load_sub_options(tmp_path)
    (tmp_path / "B1.bed").write_text("c1\t10\t20\tB1\n")
    with pytest.raises(ValueError, match="named after its own class"):
        estimate.load_sub_options(tmp_path)
    (tmp_path / "B1.bed").unlink()
    assert list(estimate.sinks_by_hash([tmp_path / "A1.bed", tmp_path / "missing.bed"]).values()) == [sinks.read_bed(tmp_path / "A1.bed")]


def test_the_shipped_sub_options():
    """The four files are subsets of the satellite sinks the menu names for their class, each class is a
    compositional class of the satellite panel, and the menu offers each as a subset of that class."""
    got = estimate.load_sub_options(estimate.sub_option_dir(ROOT / "resources" / "GRCh38"))
    assert {k: (v.parent, len(v.intervals)) for k, v in got.items()} == {
        "DXZ1": ("aSatHOR", 4), "DYZ3": ("aSatHOR", 4), "DYZ1": ("HSat3", 1), "DYZ2": ("HSat1B", 78)}
    m = fetchplan.read_menu()
    for name, opt in got.items():
        row = m.options[name]
        assert row.kind == "subset" and row.parent == opt.parent and m.resolve(row.sinks).resolve() == Path(opt.path).resolve()
        parent = m.options[opt.parent]
        learned = {r[:3] for r in sinks.read_bed(m.resolve(parent.sinks)) if r[3] == opt.parent}
        assert set(opt.intervals) <= learned, name
        assert row.panel == parent.panel and parent.kind == "compositional"
        assert opt.parent in fetchplan.panel_classes(m.resolve(row.panel))
    assert all(iv[0] == "chrX" and 58_000_000 <= iv[1] < iv[2] <= 63_000_000 for iv in got["DXZ1"].intervals)
    assert all(iv[0] == "chrY" and 10_000_000 <= iv[1] < iv[2] <= 10_700_000 for iv in got["DYZ3"].intervals)
    assert got["DYZ1"].intervals == [("chrY", 56_669_000, 56_781_000)]
    assert {iv[0] for iv in got["DYZ2"].intervals} >= {"chrX", "chrY", "chr1"}          # Y-derived HSat1B lands mostly off chrY


# ---------------------------------------------------------------- the command line, on the fixture (needs the engine)

BIN = ROOT / "target" / "release" / "ngs-dose"
BAM = ROOT / "tests" / "data" / "NA12878.subsample.bam"
BUNDLE = ROOT / "resources" / "GRCh38"


def cli(*args):
    import os
    import subprocess
    import sys
    env = {**os.environ, "PYTHONPATH": str(ROOT) + os.pathsep + os.environ.get("PYTHONPATH", "")}
    return subprocess.run([sys.executable, "-m", "ngsdose", *map(str, args)], capture_output=True, text=True, env=env)


def test_estimate_of_plan_fetches_on_the_command_line(tmp_path):
    """An ordinary fetch (the bundle's panel and sinks) never counted the options' classes: the table says
    not_counted, the log does not list it. A fetch of the core_tel and xy_arrays plan read aSatHOR, HSat3 and
    HSat1B at the arrays only: estimated without the plan's BED they are unverified (no value), with it subset_only
    and the options measured."""
    import csv
    import os
    import subprocess
    bin_ = Path(os.environ.get("NGSDOSE_BIN", BIN))
    if not bin_.exists():
        pytest.skip("needs the engine (cargo build --release, or NGSDOSE_BIN)")
    plain = tmp_path / "plain.json.gz"
    subprocess.run([str(bin_), "count", "-i", str(BAM), "-p", str(BUNDLE / "panel.k31.tsv.gz"), "-c", str(BUNDLE / "controls.fa.gz"), "-m", "fetch",
                    "--sinks", str(BUNDLE / "sinks.bed"), "-@", "2", "-o", str(plain)], check=True, capture_output=True)
    r = cli("estimate", plain, "-r", BUNDLE, "-t", tmp_path / "plain.tsv")
    assert r.returncode == 0, r.stderr
    row = next(csv.DictReader(open(tmp_path / "plain.tsv"), delimiter="\t"))
    assert row["DXZ1.status"] == "not_counted" and "DXZ1" not in r.stderr and "not_counted" not in r.stderr
    p = cli("fetchplan", "--preset", "core_tel", "xy_arrays", "--controls", BUNDLE / "controls.lite200.bed", "--engine", bin_, "-o", tmp_path / "xy")
    assert p.returncode == 0, p.stderr
    flags = (tmp_path / "xy.count_flags.txt").read_text().split()
    panels = [x for q in (tmp_path / "xy.panels.txt").read_text().split() for x in ("-p", q)]
    xy = tmp_path / "xy.json.gz"
    subprocess.run([str(bin_), "count", "-i", str(BAM), *panels, "-c", (tmp_path / "xy.controls.txt").read_text().strip(), "-m", "fetch",
                    "--sinks", str(tmp_path / "xy.sinks.bed"), *flags, "-@", "2", "-o", str(xy)], check=True, capture_output=True)
    # without the BED, a class all of whose placements lie at its options is subset_only by the placements, the others unverified
    for extra, parent, sub in (((), {"unverified", "subset_only"}, "unverified"), (("--fetch-sinks", tmp_path / "xy.sinks.bed"), {"subset_only"}, "ok")):
        r = cli("estimate", xy, "-r", BUNDLE, "-t", tmp_path / "xy.tsv", *extra)
        assert r.returncode == 0, r.stderr
        row = next(csv.DictReader(open(tmp_path / "xy.tsv"), delimiter="\t"))
        got = {c: row[f"{c}.status"] for c in ("aSatHOR", "HSat3")}      # HSat1B: some of its intervals are on contigs the fixture lacks
        assert set(got.values()) <= parent and all(row[f"{c}.mass_Mb"] == "NA" for c in got), (got, extra)
        assert row["DXZ1.status"] == sub and row["DYZ1.status"] == sub
        if not extra:
            assert "unverified" in got.values() and "--fetch-sinks" in r.stderr
            assert "may have been read whole or only at the intervals of the sub-options" in r.stderr
