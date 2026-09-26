"""`ngsdose sinks`: positional classes always get sinks; a compositional class only when asked for
(the telomeric repeat, which the aligner concentrates), on its own coarser grid, clipped to the contig."""
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from ngsdose import sinks


def scan(sample, tel_bins, unit_bins, length=1_000_000):
    pl = [dict(**{"class": "unit"}, contig="c1", start=s, reads=n) for s, n in unit_bins] + [dict(**{"class": "TEL"}, contig="c1", start=s, reads=n) for s, n in tel_bins]
    return dict(format="ngs-dose-counts/1", sample=sample, mode="scan", placement_bin=1000, placement_bin_compositional=10000,
                classes=[dict(name="unit", kind="positional"), dict(name="TEL", kind="compositional")],
                contigs=[dict(name="c1", len=length, reads=0)], placements=pl)


def test_a_compositional_class_gets_sinks_only_when_named(tmp_path):
    scans = [scan("a", [(990_000, 5000), (500_000, 3000), (200_000, 3)], [(100_000, 1000)]),
             scan("b", [(990_000, 4000), (500_000, 3500)], [(100_000, 900)])]
    rows, stats = sinks.learn(scans)
    assert {r[3] for r in rows} == {"unit"} and set(stats) == {"unit"}
    rows, stats = sinks.learn(scans, classes=["TEL"])
    tel = [r for r in rows if r[3] == "TEL"]
    assert [(c, s, e) for c, s, e, _ in tel] == [("c1", 499_000, 511_000), ("c1", 989_000, 1_000_000)]   # 10-kb grid, 1-kb pad, clipped to the contig
    assert stats["TEL"]["a"] > 0.999 and stats["TEL"]["b"] == 1.0                                       # the three stray reads are not a sink
    assert sinks.capture(scans[0], rows)["TEL"] == (8003, 8000)
    # the CLI
    f = tmp_path / "a.json"
    f.write_text(json.dumps(scans[0]))
    out = subprocess.run([sys.executable, "-m", "ngsdose", "sinks", str(f), "--classes", "TEL"], check=True, capture_output=True, text=True)
    assert "TEL" in out.stdout and "[sinks] TEL: 2 intervals" in out.stderr


def test_a_bin_counts_as_captured_only_when_all_of_it_is_inside():
    """[start, min(start + bin width, contig length)] must lie inside the class's (merged) intervals:
    the last bin of a contig ends at the contig's end, and a bin that only begins inside a sink is missed."""
    s = scan("a", [(990_000, 10), (500_000, 20), (300_000, 40), (0, 80)], [], length=995_000)
    bed = [("c1", 985_000, 995_000, "TEL"), ("c1", 495_000, 505_000, "TEL"), ("c1", 305_000, 310_000, "TEL"), ("c1", 310_000, 312_000, ""),
           ("c1", 0, 4_000, "TEL"), ("c1", 4_000, 10_000, "TEL")]
    # 990,000: clipped at 995,000, inside; 500,000: runs past 505,000; 300,000: starts before the sink;
    # 0: inside two touching intervals, which merge
    assert sinks.capture(s, bed)["TEL"] == (150, 90)


def test_learn_refuses_scans_of_another_reference_or_panel_and_names_it():
    a = scan("a", [], [(100_000, 1000)])
    b = scan("b", [], [(100_000, 900)], length=2_000_000)
    for pool in ([a, b], [b, a]):
        with pytest.raises(ValueError, match="different references"):
            sinks.learn(pool)
    p1, p2, p3 = dict(a, panel_sha256=["x"]), dict(scan("c", [], [(100_000, 900)]), panel_sha256=["x", "t"]), dict(a, sample="d", panel_sha256=["y"])
    assert sinks.learn([p1, p2])[0]                                 # a scan with an extra panel: the shared one is the same
    with pytest.raises(ValueError, match="different panel files"):
        sinks.learn([p1, p3])
    assert sinks.learn([p1, p3], allow_mixed_panels=True)[0]
    far = scan("e", [], [(1_000_500, 1000)])                        # a placement beyond the contig's end
    with pytest.raises(ValueError, match="beyond the contig's end"):
        sinks.learn([far])


def test_a_class_name_no_scan_has_is_an_error(tmp_path):
    scans = [scan("a", [(500_000, 3000)], [(100_000, 1000)])]
    with pytest.raises(ValueError, match="no scan has a class TELL"):
        sinks.learn(scans, classes=["TELL"])
    f = tmp_path / "a.json"
    f.write_text(json.dumps(scans[0]))
    r = subprocess.run([sys.executable, "-m", "ngsdose", "sinks", str(f), "--classes", "TELL"], capture_output=True, text=True)
    assert r.returncode != 0 and "TELL" in r.stderr and "Traceback" not in r.stderr and not r.stdout


def test_learn_reads_each_file_once_per_pass_and_gives_the_same_answer(tmp_path):
    scans = [scan("a", [(990_000, 5000), (500_000, 3000)], [(100_000, 1000)]), scan("b", [(990_000, 4000)], [(100_000, 900), (400_000, 500)])]
    paths = []
    for s in scans:
        paths.append(tmp_path / f"{s['sample']}.json.gz")
        with gzip.open(paths[-1], "wt") as fh:
            json.dump(s, fh)
    want = sinks.learn(scans, classes=["TEL"])
    assert sinks.learn(paths, classes=["TEL"]) == sinks.learn(iter(scans), classes=["TEL"]) == want


def test_evaluate_reads_bed_files_as_the_engine_does_and_refuses_a_fetch(tmp_path):
    s = scan("a", [(990_000, 5000), (500_000, 3000)], [(100_000, 1000)])
    f = tmp_path / "a.json"
    f.write_text(json.dumps(s))
    bed = tmp_path / "s.bed"
    bed.write_text("track name=sinks\n# learned\n\nc1\t99000\t101000\tunit\nc1\t989000\t1000000\n")   # 3 columns: every class
    r = subprocess.run([sys.executable, "-m", "ngsdose", "sinks", str(f), "--evaluate", str(bed)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    out = {l.split("\t")[1]: l.split("\t") for l in r.stdout.splitlines()[1:]}
    assert out["unit"][2:5] == ["1000", "1000", "1.00000"] and out["TEL"][2:5] == ["8000", "5000", "0.62500"]
    (tmp_path / "f.json").write_text(json.dumps(dict(s, mode="fetch")))
    r = subprocess.run([sys.executable, "-m", "ngsdose", "sinks", str(tmp_path / "f.json"), "--evaluate", str(bed)], capture_output=True, text=True)
    assert r.returncode != 0 and "fetch-mode counts" in r.stderr and r.stdout.strip() == "sample\tclass\tscan_reads\tcaptured\tfraction\tunmapped"
    bad = tmp_path / "bad.bed"
    bad.write_text("c1\t5\n")
    r = subprocess.run([sys.executable, "-m", "ngsdose", "sinks", str(f), "--evaluate", str(bad)], capture_output=True, text=True)
    assert r.returncode != 0 and "bad.bed:1: expected at least 3" in r.stderr and "Traceback" not in r.stderr


BUNDLE = Path(__file__).resolve().parents[1] / "resources" / "GRCh38"
BIN = Path(os.environ.get("NGSDOSE_BIN", BUNDLE.parents[1] / "target" / "release" / "ngs-dose"))


@pytest.mark.parametrize("head, ok", [("track name=s\n# c\n\n", True), ("browser position chr1:1-100\n", False),
                                      ("   \n", False), ("chr1\t 5\t10\tDJ\n", False)])
def test_read_bed_accepts_what_the_engine_accepts(tmp_path, head, ok):
    bed = tmp_path / "s.bed"
    bed.write_text(head + (BUNDLE / "sinks.bed").read_text())
    if ok:
        assert len(sinks.read_bed(bed)) == len(sinks.read_bed(BUNDLE / "sinks.bed"))
    else:
        with pytest.raises(ValueError, match="s.bed:1: "):
            sinks.read_bed(bed)
    if BIN.exists():
        r = subprocess.run([str(BIN), "plan", "-c", str(BUNDLE / "controls.fa.gz"), "--sinks", str(bed), "-o", str(tmp_path / "p.bed")],
                           capture_output=True, text=True)
        assert (r.returncode == 0) == ok, r.stderr


def test_sinks_to_stdout_leaves_stdout_open(tmp_path, capsys):
    from ngsdose.cli import main
    f = tmp_path / "a.json"
    f.write_text(json.dumps(scan("a", [], [(100_000, 1000)])))
    main(["sinks", str(f), "-o", "-"])
    assert not sys.stdout.closed
    print("after")
    assert capsys.readouterr().out == "c1\t99000\t102000\tunit\nafter\n"
