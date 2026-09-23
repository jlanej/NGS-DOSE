"""`ngsdose sinks`: positional classes always get sinks; a compositional class only when asked for
(the telomeric repeat, which the aligner concentrates), on its own coarser grid, clipped to the contig."""
import json
import subprocess
import sys

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
