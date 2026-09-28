"""`ngs-dose count --classes` and `plan --classes`: counting a subset of the loaded classes. Each named
class must come out exactly as in the full load (same panels, same fetched intervals), the others must
not appear, a fetch reads and checks the sinks of the named classes only, and without the option
nothing changes. Runs the built engine on the committed fixture and on the simulated genome; skipped
without the engine (the simulated genome also needs samtools)."""
import json
import os
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("NGSDOSE_BIN", ROOT / "target" / "release" / "ngs-dose"))
BAM = ROOT / "tests" / "data" / "NA12878.subsample.bam"
B = ROOT / "resources" / "GRCh38"
E = ROOT / "resources" / "experimental"
C = E / "candidates"
# the bundle's positional classes, the telomere and satellite families (compositional, competing for
# reads), and two candidate files: positional macrosatellites and compositional coding VNTRs
PANELS = [B / "panel.k31.tsv.gz", E / "telomere.k31.panel.tsv.gz", E / "satellites.CHM13v2.k31.panel.tsv.gz",
          C / "macrosatellites.k31.panel.tsv.gz", C / "coding-vntrs.k31.panel.tsv.gz"]

pytestmark = pytest.mark.skipif(not BIN.exists(), reason="needs target/release/ngs-dose (cargo build --release)")

# what a selection may change: the classes listed and their placements, the read-level diagnostics
# (only reads that concern a named class enter them) and the record of the selection itself
SELECTION = {"classes", "placements", "classes_selected", "ambiguous_reads", "below_threshold_reads", "elapsed_sec"}
# and between two fetches with different sinks files: the file's name and hash, and the classes it lacks
SINKS_FILE = {"sinks", "sinks_sha256", "sinks_missing_classes"}


def engine(*args):
    return subprocess.run([str(BIN), *map(str, args)], capture_output=True, text=True, timeout=600)


def ok(r):
    assert r.returncode == 0, r.stderr
    return r


def panel_args(panels):
    return [a for p in panels for a in ("-p", p)]


def count(out, *extra, panels=PANELS, inp=BAM, controls=B / "controls.fa.gz"):
    ok(engine("count", "-i", inp, "-c", controls, *panel_args(panels), "-@", 2, "-o", out, *extra))
    return json.loads(Path(out).read_text())


def rows(bed, keep=None):
    """BED rows (column 4 the class), optionally only those of the classes in `keep`."""
    out = [ln for ln in Path(bed).read_text().splitlines() if ln and not ln.startswith(("#", "track"))]
    return [ln for ln in out if keep is None or ln.split("\t")[3] in keep]


def write_bed(path, lines):
    Path(path).write_text("".join(ln + "\n" for ln in lines))
    return path


def check_selection(full, sub, chosen, differ=SELECTION):
    """`sub` counted `chosen` out of the classes of `full`: the same entries and placements for them,
    nothing for the others, everything else equal."""
    order = [c["name"] for c in full["classes"] if c["name"] in chosen]
    assert sorted(order) == sorted(chosen), "every chosen class is loaded"
    assert [c["name"] for c in sub["classes"]] == order and sub["classes_selected"] == order
    assert "classes_selected" not in full
    by = lambda c: {x["name"]: x for x in c["classes"]}
    for name in order:
        assert by(sub)[name] == by(full)[name], name
    place = lambda c: [p for p in c["placements"] if p["class"] in chosen]
    assert sub["placements"] == place(full)
    assert {p["class"] for p in sub["placements"]} <= set(chosen)
    assert sub["ambiguous_reads"] <= full["ambiguous_reads"] and sub["below_threshold_reads"] <= full["below_threshold_reads"]
    assert {k: v for k, v in sub.items() if k not in differ} == {k: v for k, v in full.items() if k not in differ}


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    d = tmp_path_factory.mktemp("classes")
    sinks = write_bed(d / "sinks.all.bed", rows(B / "sinks.bed") + rows(E / "sinks.satellites.bed"))
    return {"dir": d, "sinks": sinks, "scan": count(d / "scan.full.json", "-m", "scan")}


# DJ alone; two compositional classes (they compete with every family loaded); positional and
# compositional from three files; one satellite family
SELECTIONS = [["DJ"], ["TEL", "HSat3"], ["CER", "rDNA5S", "D4Z4"], ["bSat"], ["rDNA45S", "rDNA5S", "DJ"]]


@pytest.mark.parametrize("chosen", SELECTIONS, ids=lambda c: "+".join(c))
def test_a_scan_of_named_classes_counts_them_as_the_full_load(fx, chosen):
    sub = count(fx["dir"] / f"scan.{'_'.join(chosen)}.json", "-m", "scan", "--classes", ",".join(chosen))
    check_selection(fx["scan"], sub, chosen)
    assert sum(c["reads"] for c in sub["classes"]) > 0, "the fixture holds reads of the chosen classes"


def test_the_named_order_and_repeats_do_not_matter(fx):
    sub = count(fx["dir"] / "scan.order.json", "-m", "scan", "--classes", "D4Z4, rDNA5S,CER", "--classes", "rDNA5S")
    check_selection(fx["scan"], sub, ["CER", "rDNA5S", "D4Z4"])


@pytest.mark.parametrize("unmapped", [False, True], ids=["fetch", "fetch-unmapped"])
def test_a_fetch_reads_and_counts_the_named_classes_only(fx, unmapped):
    d, chosen = fx["dir"], ["DJ", "TEL", "HSat2", "rDNA5S"]
    extra = ["-m", "fetch"] + (["--unmapped"] if unmapped else [])
    tag = "u" if unmapped else ""
    # a sinks BED that holds the chosen classes' intervals only, as `ngsdose fetchplan` writes it
    part = write_bed(d / "sinks.part.bed", rows(fx["sinks"], set(chosen)))
    # the full load cannot fetch with it unless told to go on without sinks for the other classes ...
    r = engine("count", "-i", BAM, "-c", B / "controls.fa.gz", *panel_args(PANELS), *extra, "--sinks", part, "-o", d / "x.json")
    assert r.returncode == 1 and "no interval" in r.stderr and "--classes" in r.stderr
    full = count(d / f"fetch{tag}.full.part.json", *extra, "--sinks", part, "--allow-missing-sinks")
    assert set(full["sinks_missing_classes"]) == {c["name"] for c in full["classes"]} - set(chosen)
    # ... the selection needs no flag for the classes it leaves out; it still needs one for HSat2, two of
    # whose sinks are on contigs the fixture lacks (an engine since 0.1.1 refuses a class that loses any)
    r = engine("count", "-i", BAM, "-c", B / "controls.fa.gz", *panel_args(PANELS), *extra, "--sinks", part, "--classes", ",".join(chosen), "-o", d / "x.json")
    assert r.returncode == 1 and "HSat2 lose sink intervals" in r.stderr and "no interval" not in r.stderr, r.stderr
    sub = count(d / f"fetch{tag}.sub.part.json", *extra, "--sinks", part, "--classes", ",".join(chosen), "--allow-missing-sinks")
    assert "sinks_missing_classes" not in sub
    check_selection(full, sub, chosen, SELECTION | {"sinks_missing_classes"})
    # with the whole sinks BED it reads the chosen classes' intervals only: the same counts
    whole = count(d / f"fetch{tag}.sub.all.json", *extra, "--sinks", fx["sinks"], "--classes", ",".join(chosen), "--allow-missing-sinks")
    check_selection(full, whole, chosen, SELECTION | SINKS_FILE)
    assert whole["sinks_skipped"] == full["sinks_skipped"] and set(whole["sinks_skipped"]) == {"HSat2"}
    assert sum(c["reads"] for c in sub["classes"]) > 0


def test_missing_sinks_are_checked_for_the_named_classes_only(fx):
    d = fx["dir"]
    no_tel = write_bed(d / "sinks.noTEL.bed", [ln for ln in rows(B / "sinks.bed") if ln.split("\t")[3] != "TEL"])
    panels = PANELS[:2]
    base = ["-m", "fetch", "--sinks", no_tel]
    count(d / "fetch.dj.json", *base, "--classes", "DJ,rDNA45S", panels=panels)
    r = engine("count", "-i", BAM, "-c", B / "controls.fa.gz", *panel_args(panels), *base, "--classes", "DJ,TEL", "-o", d / "x.json")
    assert r.returncode == 1
    msg = re.search(r"for the (?:selected|loaded) class\(es\) ([^:]*):", r.stderr).group(1)
    assert msg == "TEL", r.stderr
    sub = count(d / "fetch.tel.json", *base, "--classes", "DJ,TEL", "--allow-missing-sinks", panels=panels)
    assert sub["sinks_missing_classes"] == ["TEL"]


def test_a_sinks_bed_without_class_names_is_read_whole(fx):
    d = fx["dir"]
    plain = write_bed(d / "sinks.plain.bed", ["\t".join(ln.split("\t")[:3]) for ln in rows(B / "sinks.bed")])
    r = ok(engine("count", "-i", BAM, "-c", B / "controls.fa.gz", "-p", B / "panel.k31.tsv.gz", "-m", "fetch", "--sinks", plain,
                  "--classes", "DJ", "-o", d / "plain.json"))
    assert "cannot choose among its intervals" in r.stderr
    assert [c["name"] for c in json.loads((d / "plain.json").read_text())["classes"]] == ["DJ"]


@pytest.mark.parametrize("bad, needle", [("DJ,NOPE", "'NOPE', which no loaded panel defines"), ("dj", "'dj'"), (",", "names no class")])
def test_an_unknown_class_is_refused(fx, bad, needle):
    r = engine("count", "-i", BAM, "-c", B / "controls.fa.gz", *panel_args(PANELS[:2]), "-m", "scan", "--classes", bad,
               "-o", fx["dir"] / "bad.json")
    assert r.returncode == 1 and needle in r.stderr and "define rDNA45S, rDNA5S, DJ, TEL" in r.stderr, r.stderr
    assert not (fx["dir"] / "bad.json").exists()


def test_plan_keeps_the_named_classes_sinks(fx):
    d, chosen = fx["dir"], ["TEL", "HSat2"]
    part = write_bed(d / "sinks.plan.bed", rows(fx["sinks"], set(chosen)))
    plan = lambda out, sinks, *extra: ok(engine("plan", "-c", B / "controls.fa.gz", "--sinks", sinks, "-o", out, *extra))
    for extra in ([], ["-i", BAM]):
        r = plan(d / "plan.sub.bed", fx["sinks"], "--classes", ",".join(chosen), *extra)
        assert "(of TEL, HSat2)" in r.stderr
        plan(d / "plan.part.bed", part, *extra)
        plan(d / "plan.all.bed", fx["sinks"], *extra)
        sub, part_, whole = ((d / f"plan.{n}.bed").read_text() for n in ("sub", "part", "all"))
        assert sub == part_ and sub != whole and len(sub) < len(whole)
    r = engine("plan", "-c", B / "controls.fa.gz", "--sinks", fx["sinks"], "--classes", "TEL,NOPE")
    assert r.returncode == 1 and "no sink interval for NOPE" in r.stderr
    plain = write_bed(d / "plan.plain.bed", ["\t".join(ln.split("\t")[:3]) for ln in rows(B / "sinks.bed")])
    r = engine("plan", "-c", B / "controls.fa.gz", "--sinks", plain, "--classes", "TEL")
    assert r.returncode == 1 and "gives no class" in r.stderr


def test_the_simulated_genome(sim, tmp_path):
    """A positional unit and a compositional satellite with known reads: each counted alone as in the
    full load, in scan and fetch."""
    d = sim["dir"]
    common = dict(panels=[d / "panel.tsv.gz"], inp=d / "sim.bam", controls=d / "controls.fa.gz")
    grid = ["--l-grid", "100,200,300,400"]
    full = count(tmp_path / "scan.json", "-m", "scan", *grid, **common)
    assert {c["name"]: c["reads"] for c in full["classes"]}["unit"] == sim["n_class_reads"]
    for name in ("unit", "sat"):
        check_selection(full, count(tmp_path / f"scan.{name}.json", "-m", "scan", *grid, "--classes", name, **common), [name])
        part = write_bed(tmp_path / f"sinks.{name}.bed", rows(d / "sinks.bed", {name}))
        for extra in ([], ["--unmapped"]):
            fetch = ["-m", "fetch", *grid, *extra]
            ref = count(tmp_path / "f.full.json", *fetch, "--sinks", part, "--allow-missing-sinks", **common)
            sub = count(tmp_path / "f.sub.json", *fetch, "--sinks", d / "sinks.bed", "--classes", name, **common)
            check_selection(ref, sub, [name], SELECTION | SINKS_FILE)
            assert sub["classes"][0]["reads"] == {c["name"]: c["reads"] for c in full["classes"]}[name]


# keys the engine gained on the review branch, before --classes (main 7772e32 lacks them): left out of
# the comparison only when the baseline does not write them
REVIEW_KEYS = ("pipeline", "pad", "sinks_skipped")


def without_run_time(raw):
    return re.sub(rb'"elapsed_sec":[0-9.e+-]+', b"", raw)


def same_output(ours, base):
    """Byte for byte (but the run time) when the baseline writes the same keys, as a build of this
    branch before --classes does. Against an older engine (main 7772e32), the parsed JSON, keys in the
    same order, with only REVIEW_KEYS that the baseline lacks left out of ours; a stray key such as
    classes_selected still fails."""
    a, b = json.loads(ours), json.loads(base)
    a.pop("elapsed_sec", None)
    b.pop("elapsed_sec", None)
    extra = [k for k in REVIEW_KEYS if k in a and k not in b]
    if not extra:
        return without_run_time(ours) == without_run_time(base)
    for k in extra:
        del a[k]
    return list(a) == list(b) and a == b


@pytest.mark.skipif(not os.environ.get("NGSDOSE_BASELINE_BIN"),
                    reason="set NGSDOSE_BASELINE_BIN to an engine without --classes: this branch built before it, or main 7772e32")
def test_without_the_option_nothing_changes(fx, tmp_path):
    """What an engine without --classes writes (scan, fetch, fetch --unmapped): byte for byte but the run
    time against this branch built before --classes; against main 7772e32, the same JSON apart from the
    keys the review branch added (REVIEW_KEYS)."""
    base = Path(os.environ["NGSDOSE_BASELINE_BIN"])
    fetch = ["-m", "fetch", "--sinks", fx["sinks"], "--allow-missing-sinks"]    # the candidate files have no sinks yet
    for extra in (["-m", "scan"], fetch, fetch + ["--unmapped"]):
        outs = []
        for i, b in enumerate((BIN, base)):
            r = subprocess.run([str(b), "count", "-i", BAM, "-c", B / "controls.fa.gz", *panel_args(PANELS), "-@", "2",
                                "-o", tmp_path / f"{i}.json", *map(str, extra)], capture_output=True, text=True, timeout=600)
            assert r.returncode == 0, r.stderr
            outs.append((tmp_path / f"{i}.json").read_bytes())
        assert "classes_selected" not in json.loads(outs[0]), extra
        assert same_output(*outs), extra


def test_the_baseline_comparison_is_strict():
    """same_output leaves out only the review branch's keys, and only when the baseline lacks them."""
    base = {"a": 1, "classes": [1], "elapsed_sec": 1.0}

    def enc(d):
        return json.dumps(d).encode()

    assert same_output(enc({"a": 1, "classes": [1], "pad": 600, "pipeline": None, "elapsed_sec": 2.0}), enc(base))
    assert not same_output(enc({"a": 1, "classes": [1], "pad": 600, "classes_selected": ["x"]}), enc(base))
    assert not same_output(enc({"a": 1, "classes": [2], "pad": 600}), enc(base))
    assert not same_output(enc({"classes": [1], "a": 1, "pad": 600}), enc(base))      # key order
    assert not same_output(enc({"a": 1, "pad": 600, "classes": [1]}), enc(dict(base, pad=0)))
