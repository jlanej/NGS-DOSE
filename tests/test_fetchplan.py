"""`ngsdose fetchplan` and `ngsdose.cost`: the bytes a fetch reads, from a CRAM index, and choosing
the options of a fetch from the menu by name, preset, budget and capture target."""
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from ngsdose import contract, cost, estimate, fetchplan, sinks

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("NGSDOSE_BIN", ROOT / "target" / "release" / "ngs-dose"))

CONTIGS = [("c1", 1_000_000), ("c2", 5_000), ("c3", 5_000)]
# ref id, 1-based start, span, container offset, slice offset, slice size
CRAI = [(0, 1, 10_000, 100, 10, 1000),
        (0, 10_001, 10_000, 2000, 20, 2000),
        (1, 1, 5_000, 5000, 30, 500),                      # one multi-reference slice, listed under c2 and c3
        (2, 1, 5_000, 5000, 30, 500),
        (-1, 0, 0, 8000, 40, 300)]                         # the unmapped bin


def write_crai(path, lines):
    with gzip.open(path, "wt") as fh:
        fh.write("".join("\t".join(map(str, x)) + "\n" for x in lines))
    return path


def write_fai(path, contigs=CONTIGS):
    path.write_text("".join(f"{n}\t{ln}\t0\t60\t61\n" for n, ln in contigs))
    return path


def test_the_planner_and_the_engine_join_queries_at_the_same_gap():
    """What the plan prices is what the engine reads only while both use one rule: cost.GROUP_GAP is src/count.rs's."""
    import re
    src = (Path(__file__).resolve().parents[1] / "src" / "count.rs").read_text()
    assert int(re.search(r"pub const DEFAULT_GROUP_GAP: i64 = ([\d_]+);", src)[1].replace("_", "")) == cost.GROUP_GAP


def test_slices_are_priced_per_fetch_with_their_compression_header(tmp_path):
    """A fetch decodes every slice its interval overlaps, once, with the compression header of each container among
    them; a fetch on another contig decodes a multi-reference slice again, and the floor (bytes_once) reads it once."""
    ix = cost.CraiIndex(write_crai(tmp_path / "x.crai", CRAI), CONTIGS)
    assert ix.total == 1000 + 2000 + 500 + 300 + 10 + 20 + 30 + 40
    got = cost.component_costs(ix, {"A": [("c1", 5_000, 15_000)], "B": [("c3", 100, 200)], "C": [("c2", 0, 10), ("c9", 0, 10)],
                                    "U": [(cost.UNMAPPED, 0, 0)]})
    assert got["A"]["bytes"] == 1000 + 10 + 2000 + 20 and got["A"]["slices"] == 2       # one fetch over two slices, two containers
    assert got["B"]["bytes"] == got["C"]["bytes"] == 530 and got["C"]["absent"] == 1
    assert got["C"]["cum_bytes"] == 3030 + 530 + 530                                    # the slice B and C share: a fetch on c3 and one on c2 each decode it
    assert got["C"]["cum_bytes_once"] == 3030 + 530                                     # the floor reads it once
    assert got["U"]["bytes"] == 340 and got["union"]["bytes"] == 3030 + 530 + 530 + 340 and got["union"]["slices"] == 4
    assert got["union"]["bytes_once"] == ix.total
    # an interval ending where a slice starts does not read it; one starting at a slice's last base does
    assert ix.keys([("c1", 0, 0)])[0] == set() and ix.keys([("c1", 9_999, 10_000)])[0] == {(100, 10)}
    assert ix.keys([("c1", 10_000, 10_001)])[0] == {(2000, 20)}


def test_a_slice_under_two_fetches_is_decoded_twice_and_the_floor_once(tmp_path):
    """An engine before 0.3.0 (gap 0) merges a plan's intervals where they touch or overlap and fetches each run: two runs
    over one slice decode it (and its container's compression header) twice, one run once; the floor is the union of
    slices, each once."""
    ix = cost.CraiIndex(write_crai(tmp_path / "x.crai", CRAI), CONTIGS, gap=0)
    apart, touching, overlapping = [("c1", 100, 200), ("c1", 300, 400)], [("c1", 100, 200), ("c1", 200, 300)], [("c1", 100, 250), ("c1", 200, 300)]
    assert ix.price(apart) == 2 * 1010 and ix.price_once(apart) == 1010
    assert ix.price(touching) == ix.price(overlapping) == 1010 == ix.price_once(touching)
    assert cost.merge(apart, 0) == {"c1": ([100, 300], [200, 400])} and cost.merge(touching, 0) == {"c1": ([100], [300])}
    assert cost.merge(overlapping + [(cost.UNMAPPED, 0, 0), (cost.UNMAPPED, 0, 0)], 0) == {"c1": ([100], [300]), cost.UNMAPPED: ([0], [0])}
    # the engine since 0.3.0 reads intervals no more than 50,000 bp apart with one query: the slice they share is decoded once
    now = cost.CraiIndex(tmp_path / "x.crai", CONTIGS)
    assert now.gap == cost.GROUP_GAP == 50_000 and now.price(apart) == 1010 == now.price_once(apart)
    assert cost.merge(apart) == {"c1": ([100], [400])}
    # a gap of exactly the limit joins, one bp more does not (src/count.rs group_plan)
    assert cost.merge([("c1", 0, 10), ("c1", 50_010, 50_020), ("c1", 100_021, 100_030)]) == {"c1": ([0, 100_021], [50_020, 100_030])}
    assert cost.merge([(cost.UNMAPPED, 0, 0), (cost.UNMAPPED, 0, 0)]) == {cost.UNMAPPED: ([0], [0])}
    # what an interval adds under that rule: nothing within reach of a query whose slices it shares, the further slice when it extends one
    plan_now = cost.merge(apart)
    assert now.extra(plan_now, "c1", 600, 700) == 0 and now.extra(plan_now, "c1", 11_000, 12_000) == 2020 and now.extra(plan_now, "c2", 0, 10) == 530
    for iv in (("c1", 600, 700), ("c1", 11_000, 12_000), ("c1", 150, 160), ("c2", 0, 10), (cost.UNMAPPED, 0, 0)):
        assert now.extra(plan_now, *iv) == now.price(apart + [iv]) - now.price(apart), iv
    # over the two slices of c1: one run decodes each once; three runs decode the first twice and the second once
    assert ix.price([("c1", 100, 200), ("c1", 200, 12_000)]) == 3030
    assert ix.price([("c1", 100, 200), ("c1", 300, 400), ("c1", 11_000, 12_000)]) == 2 * 1010 + 2020
    # what an interval adds to a plan's fetches: nothing inside or touching a run without new slices, its own fetch apart, less when it bridges two runs
    plan = cost.merge(apart, 0)
    assert ix.extra(plan, "c1", 150, 160) == 0 and ix.extra(plan, "c1", 400, 500) == 0 and ix.extra(plan, "c1", 200, 300) == 1010 - 2 * 1010
    assert ix.extra(plan, "c1", 600, 700) == 1010 and ix.extra(plan, "c2", 0, 10) == 530
    assert ix.extra(plan, cost.UNMAPPED, 0, 0) == 340 and ix.extra(cost.merge([(cost.UNMAPPED, 0, 0)], 0), cost.UNMAPPED, 0, 0) == 0
    for iv in (("c1", 150, 160), ("c1", 400, 500), ("c1", 200, 300), ("c1", 600, 700), ("c2", 0, 10), (cost.UNMAPPED, 0, 0)):
        assert ix.extra(plan, *iv) == ix.price(apart + [iv]) - ix.price(apart), iv
    # the multi-reference slice is one slice within a fetch, decoded again by a second run on its contig or a run on its other contig
    assert ix.price([("c2", 0, 10), ("c2", 20, 30)]) == 2 * 530 and ix.price([("c2", 0, 10), ("c3", 0, 10)]) == 2 * 530
    assert ix.price_once([("c2", 0, 10), ("c3", 0, 10)]) == 530 and ix.price([("c2", 0, 10), ("c2", 10, 30)]) == 530


def test_a_long_slice_is_found_behind_short_ones(tmp_path):
    lines = [(0, 1, 500_000, 100, 10, 1000), (0, 1_001, 100, 2000, 10, 50), (0, 2_001, 100, 3000, 10, 50)]
    ix = cost.CraiIndex(write_crai(tmp_path / "x.crai", lines), CONTIGS)
    assert ix.keys([("c1", 300_000, 300_001)])[0] == {(100, 10)}
    assert ix.keys([("c1", 1_050, 2_050)])[0] == {(100, 10), (2000, 10), (3000, 10)}


def test_an_index_of_another_header_is_refused(tmp_path):
    p = write_crai(tmp_path / "x.crai", CRAI)
    with pytest.raises(ValueError, match="reference id 2, but the contig list has 2 contigs"):
        cost.CraiIndex(p, CONTIGS[:2])
    with pytest.raises(ValueError, match="beyond the contig's 5,000 bp"):
        cost.CraiIndex(p, [("c1", 5_000), ("c2", 5_000), ("c3", 5_000)])
    (tmp_path / "bad.crai").write_text("0\t1\t2\n")
    with pytest.raises(ValueError, match="bad.crai:1: not a .crai line"):
        cost.CraiIndex(tmp_path / "bad.crai", CONTIGS)


def test_contig_lists_come_from_a_fai_a_header_or_names(tmp_path):
    assert cost.read_contigs(write_fai(tmp_path / "r.fa.fai")) == CONTIGS
    (tmp_path / "h.sam").write_text("@HD\tVN:1.6\n" + "".join(f"@SQ\tSN:{n}\tLN:{ln}\tM5:x\n" for n, ln in CONTIGS) + "@PG\tID:bwa\n")
    assert cost.read_contigs(tmp_path / "h.sam") == CONTIGS
    (tmp_path / "n.txt").write_text("c1\nc2\n")
    assert cost.read_contigs(tmp_path / "n.txt") == [("c1", None), ("c2", None)]
    (tmp_path / "m.txt").write_text("c1\t10\nc2\n")
    with pytest.raises(ValueError, match="expected a .fai"):
        cost.read_contigs(tmp_path / "m.txt")


# ---------------------------------------------------------------- the menu and the plan

HEAD = "\t".join(fetchplan.COLUMNS) + "\tstats\n"


def panel(path, *classes):
    path.write_text("##ngs-dose-panel v1\n##k=31\n" + "".join(f"##class\tid={i}\tname={n}\tkind={k}\tlength=100\tcircular=0\n"
                                                             for i, (n, k) in enumerate(classes)) + "ACGT\t0\t0\n")


def row(name, kind, status, tier, presets, pnl, snk, stats="-"):
    return "\t".join((name, "g", kind, status, tier, presets, pnl, snk, f"measures {name}", "-", "-", stats)) + "\n"


@pytest.fixture
def menu_dir(tmp_path):
    """A menu: controls on c1; A and B (one panel) and T (its own panel, with statistics) have sinks;
    N is a candidate. The CRAM index puts A, B and T in slices of their own."""
    d = tmp_path / "res"
    d.mkdir()
    panel(d / "ab.tsv", ("A", "positional"), ("B", "positional"))
    panel(d / "t.tsv", ("T", "compositional"))
    panel(d / "n.tsv", ("N", "compositional"))
    (d / "controls.bed").write_text("c1\t50000\t51000\tcontrol\n")
    (d / "sinks.bed").write_text("c1\t100000\t101000\tA\nc2\t0\t1000\tB\nc1\t300000\t301000\tT\nc1\t400000\t401000\tT\nc3\t0\t1000\tT\n")
    t = sinks.IntervalTally([("c1", 300000, 301000, "T"), ("c1", 400000, 401000, "T"), ("c3", 0, 1000, "T")])
    for sh in ((600, 380, 20), (500, 470, 30)):         # two scans, 1000 T reads each, all inside
        pl = [dict(contig=c, start=s, reads=n, all=n, **{"class": "T"}) for (c, s), n in zip((("c1", 300000), ("c1", 400000), ("c3", 0)), sh)]
        t.add(dict(placement_bin=1000, placement_bin_compositional=1000, classes=[dict(name="T", kind="compositional")],
                   contigs=[], placements=pl))
    with open(d / "t.stats.tsv", "w") as fh:
        sinks.write_stats(fh, t.table())
    (d / "menu.tsv").write_text(
        "# a test menu\n##preset\tcore\tA and B\n" + HEAD
        + row("controls", "controls", "shipped", "A", "-", "-", "controls.bed")
        + row("unmapped", "unmapped", "shipped", "C", "-", "-", "-")
        + row("A", "positional", "shipped", "A", "core", "ab.tsv", "sinks.bed")
        + row("B", "positional", "shipped", "B", "core", "ab.tsv", "sinks.bed")
        + row("T", "compositional", "experimental: no fetch compared yet", "C", "-", "t.tsv", "sinks.bed", "t.stats.tsv")
        + row("N", "compositional", "candidate", "C", "-", "n.tsv", "-"))
    lines = [(0, 49_001, 3_000, 100, 10, 1000),          # controls
             (0, 99_001, 3_000, 2000, 10, 2000),         # A
             (1, 1, 3_000, 5000, 10, 500),               # B
             (0, 299_001, 3_000, 8000, 10, 3000),        # T
             (0, 399_001, 3_000, 12000, 10, 4000),       # T
             (2, 1, 3_000, 17000, 10, 100000),           # T, c3: a costly slice holding 2-3% of T
             (-1, 0, 0, 200000, 10, 700)]
    write_crai(d / "a.crai", lines)
    write_crai(d / "b.crai", [(r, s, sp, co, so, 3 * sz) for r, s, sp, co, so, sz in lines])
    write_fai(d / "ref.fai")
    return d


def names(plan):
    return [r.option.name for r in plan.rows]


def test_the_menu_is_read_and_checked(menu_dir, tmp_path):
    m = fetchplan.read_menu(menu_dir / "menu.tsv")
    assert list(m.options) == ["controls", "unmapped", "A", "B", "T", "N"] and m.members("core") == ["A", "B"]
    assert m.options["T"].status == "experimental" and m.options["T"].status_note == "no fetch compared yet"
    text = (menu_dir / "menu.tsv").read_text()
    for bad, msg in ((text.replace("\tshipped\tB\t", "\tsold\tB\t"), "status 'sold'"),
                     (text.replace("\tshipped\tB\t", "\tshipped\tE\t"), "tier 'E'"),
                     (text.replace("B\tg\tpositional", "A\tg\tpositional"), "a second class name 'A'"),
                     (text.replace("\tmeasures\t", "\tmeasure\t"), "lacks the column"),
                     (text.replace("##preset\tcore", "##preset\tfull"), "preset full has no member"),
                     (text.replace("controls\tg\tcontrols", "controls\tg\tunmapped"), "exactly one row of kind controls"),
                     (text.replace("T\tg\tcompositional\texperimental: no fetch compared yet\tC\t-\tt.tsv\tsinks.bed",
                                   "T\tg\tcompositional\texperimental\tC\t-\tt.tsv\t-"), "names no sinks file")):
        (tmp_path / "m.tsv").write_text(bad)
        with pytest.raises(ValueError, match=msg):
            fetchplan.read_menu(tmp_path / "m.tsv")


def test_by_name_and_preset_with_costs(menu_dir):
    m = fetchplan.read_menu(menu_dir / "menu.tsv")
    crais = [menu_dir / "a.crai", menu_dir / "b.crai"]
    p = fetchplan.make_plan(m, presets=["core"], classes=["unmapped"], crais=crais, contigs=menu_dir / "ref.fai", log=lambda _: None)
    assert names(p) == ["controls", "A", "B", "unmapped"]              # tier order; controls always
    # the median of the two indexes: slices of 1x and 3x the size, each with a 10-byte compression header
    assert [r.mb for r in p.rows] == [pytest.approx((2 * x + 10) / 1e6) for x in (1000, 2000, 500, 700)]
    assert p.rows[-1].cum_mb == pytest.approx((2 * 4200 + 40) / 1e6) and p.rows[-1].cum_pct == pytest.approx(100 * (4240 / 111270 + 12640 / 333670) / 2)
    assert all(r.cum_mb_floor == r.cum_mb for r in p.rows)             # every interval in a slice of its own: the floor is the price
    assert p.panels == ["ab.tsv"] and p.flags == ["--unmapped"] and p.bed() == [("c1", 100000, 101000, "A"), ("c2", 0, 1000, "B")]
    assert p.rows[0].intervals == [("c1", 49400, 51600)]              # padded by 600


def test_a_candidate_is_scanned_only_and_companions_need_a_flag(menu_dir):
    m = fetchplan.read_menu(menu_dir / "menu.tsv")
    said = []
    p = fetchplan.make_plan(m, classes=["A", "N"], log=said.append)
    assert names(p) == ["controls", "A"] and [o.name for o in p.scan_only] == ["N"]
    assert p.panels == ["ab.tsv"] and p.scan_panels == ["ab.tsv", "n.tsv"]
    assert p.flags == ["--allow-missing-sinks"] and any("N is a candidate" in s for s in said) and any("also define B" in s for s in said)
    with pytest.raises(ValueError, match="no class with sinks is selected"):
        fetchplan.make_plan(m, classes=["N"], log=said.append)
    with pytest.raises(ValueError, match="no option X"):
        fetchplan.make_plan(m, classes=["X"])
    with pytest.raises(ValueError, match="choose options"):
        fetchplan.make_plan(m)


def test_a_budget_adds_by_tier_and_cost(menu_dir):
    m = fetchplan.read_menu(menu_dir / "menu.tsv")
    kw = dict(crais=[menu_dir / "a.crai"], contigs=menu_dir / "ref.fai", log=lambda _: None)
    # controls 1010, A 2010 (tier A), B 510 (B), then tier C by cost: unmapped 710, T 107030
    assert names(fetchplan.make_plan(m, budget_mb=0.0036, **kw)) == ["controls", "A", "B"]
    assert names(fetchplan.make_plan(m, budget_mb=0.0043, **kw)) == ["controls", "A", "B", "unmapped"]
    assert names(fetchplan.make_plan(m, budget_mb=1, **kw)) == ["controls", "A", "B", "unmapped", "T"]
    assert names(fetchplan.make_plan(m, budget_mb=0.0036, statuses=("shipped",), **kw)) == ["controls", "A", "B"]
    assert names(fetchplan.make_plan(m, budget_mb=0.0031, **kw)) == ["controls", "A"]      # stops at B
    assert names(fetchplan.make_plan(m, budget_mb=0.0038, fill=False, **kw)) == ["controls", "A", "B"]
    assert names(fetchplan.make_plan(m, classes=["A", "T"], budget_mb=1, **kw)) == ["controls", "A", "T"]
    with pytest.raises(ValueError, match="control regions alone read"):
        fetchplan.make_plan(m, budget_mb=0.001, **kw)
    with pytest.raises(ValueError, match="needs the CRAM index"):
        fetchplan.make_plan(m, budget_mb=1)


def test_fill_goes_past_an_option_that_does_not_fit(menu_dir):
    m = fetchplan.read_menu(menu_dir / "menu.tsv")
    kw = dict(crais=[menu_dir / "a.crai"], contigs=menu_dir / "ref.fai", log=lambda _: None)
    # 1010 + 510 + 710 fit into 2.3 kB, A (2010, tier A) does not
    with pytest.raises(ValueError, match="no class with sinks"):
        fetchplan.make_plan(m, budget_mb=0.0023, **kw)                 # strict: stops at A, and a plan without a class is refused
    p = fetchplan.make_plan(m, budget_mb=0.0023, fill=True, **kw)
    assert names(p) == ["controls", "B", "unmapped"]


def test_a_capture_target_drops_the_least_yielding_intervals(menu_dir):
    m = fetchplan.read_menu(menu_dir / "menu.tsv")
    kw = dict(log=lambda _: None)
    full = fetchplan.make_plan(m, classes=["T"], **kw).rows[1]
    assert len(full.intervals) == 3 and full.expected == pytest.approx(1.0) and full.order == ""
    # all T reads are inside and every interval is pure, so the order is by share: 300k (0.55), 400k (0.425), c3 (0.025)
    p = fetchplan.make_plan(m, classes=["T"], capture=0.95, **kw).rows[1]
    assert p.intervals == [("c1", 300000, 301000), ("c1", 400000, 401000)] and p.dropped == 1 and p.order == "read"
    assert 0.95 <= p.expected < 1 and p.mb is None
    # the two intervals hold 0.98 and 0.97 of the two scans: 10th percentile 0.971, median 0.975
    p10 = fetchplan.make_plan(m, classes=["T"], capture_class={"T": 0.974}, **kw).rows[1]
    med = fetchplan.make_plan(m, classes=["T"], capture_class={"T": 0.974}, capture_stat="median", **kw).rows[1]
    assert len(p10.intervals) == 3 and len(med.intervals) == 2 and med.expected == pytest.approx(0.975)
    # with an index the costly c3 slice is no longer read, and the table says what that saved and cost
    b = fetchplan.make_plan(m, classes=["T"], capture=0.95, crais=[menu_dir / "a.crai"], contigs=menu_dir / "ref.fai", **kw)
    r = b.rows[1]
    assert r.intervals == p.intervals and r.mb == pytest.approx(7.02e-3) and r.mb_full == pytest.approx(107.03e-3)
    # kept per byte, the expected capture is a bound: the longer of all less c3's largest share (1 - 0.03) and the
    # statistics' own curve up to the first interval not kept (here both kept intervals: 0.971)
    assert r.order == "byte" and r.expected == pytest.approx(0.971)
    line = fetchplan.table_lines(b)[2].split("\t")
    got = dict(zip(fetchplan.PLAN_COLUMNS, line))
    assert got["capture_lost"] == "0.02900" and got["mb_saved"] == "0.1" and got["order"] == "per byte"
    said = []
    fetchplan.make_plan(m, classes=["A", "T"], capture=0.95, log=said.append)
    assert any("A: no statistics, so all its 1 intervals are kept" in s for s in said)
    with pytest.raises(ValueError, match="--capture-class A=0.9: no statistics"):
        fetchplan.make_plan(m, classes=["A"], capture_class={"A": 0.9})


def test_statistics_of_other_intervals_and_missing_sinks_are_refused(menu_dir, tmp_path):
    m = fetchplan.read_menu(menu_dir / "menu.tsv")
    (tmp_path / "other.bed").write_text("c1\t100000\t101000\tA\nc2\t0\t1000\tB\nc1\t300000\t301000\tT\n")
    with pytest.raises(ValueError, match=r"statistics of T are for other intervals than its sinks \(1 of 1"):
        fetchplan.make_plan(m, classes=["T"], sinks_files=[tmp_path / "other.bed"], capture=0.9)
    said = []
    p = fetchplan.make_plan(m, classes=["T"], sinks_files=[tmp_path / "other.bed"], log=said.append)   # no target: the statistics are not needed
    assert p.bed() == [("c1", 300000, 301000, "T")] and p.rows[1].expected is None and any("other intervals" in x for x in said)
    (menu_dir / "sinks.bed").unlink()
    with pytest.raises(ValueError, match="its sinks file .* does not exist"):
        fetchplan.make_plan(m, classes=["A"])
    (tmp_path / "blank.bed").write_text("c1\t0\t10\n")
    with pytest.raises(ValueError, match="no class column"):
        fetchplan.make_plan(m, classes=["A"], sinks_files=[tmp_path / "blank.bed"])
    (tmp_path / "a.bed").write_text("c1\t0\t10\tA\n")
    with pytest.raises(ValueError, match="B: no interval of the class"):
        fetchplan.make_plan(m, classes=["B"], sinks_files=[tmp_path / "a.bed"])


def run(*args, cwd=None):
    env = {**os.environ, "PYTHONPATH": str(ROOT) + os.pathsep + os.environ.get("PYTHONPATH", "")}
    return subprocess.run([sys.executable, "-m", "ngsdose", "fetchplan", *map(str, args)], capture_output=True, text=True, cwd=cwd, env=env)


def test_the_cli_writes_what_the_engine_takes(menu_dir, tmp_path):
    out = tmp_path / "p"
    r = run("--menu", menu_dir / "menu.tsv", "--classes", "B", "T", "N", "unmapped", "--capture", "0.95", "--crai", menu_dir / "a.crai",
            "--contigs", menu_dir / "ref.fai", "-o", out, "--panel-root", "/opt/res", "--engine", tmp_path / "no-engine")
    assert r.returncode == 0, r.stderr
    table = [line.split("\t") for line in r.stdout.splitlines()]
    assert table[0] == list(fetchplan.PLAN_COLUMNS) and [t[0] for t in table[1:]] == ["controls", "B", "unmapped", "T", "N"]
    assert table[-1][-1].startswith("scan only")
    bed = sinks.read_bed(f"{out}.sinks.bed")
    assert bed == [("c1", 300000, 301000, "T"), ("c1", 400000, 401000, "T"), ("c2", 0, 1000, "B")]
    assert Path(f"{out}.panels.txt").read_text() == "/opt/res/ab.tsv\n/opt/res/t.tsv\n"
    assert Path(f"{out}.scan_panels.txt").read_text() == "/opt/res/ab.tsv\n/opt/res/t.tsv\n/opt/res/n.tsv\n"
    assert Path(f"{out}.count_flags.txt").read_text() == "--unmapped\n--allow-missing-sinks\n"
    plan = Path(f"{out}.plan.tsv").read_text()
    assert "--capture 0.95" in plan and "N is a candidate" in plan and "\nT\texperimental\tC\tcompositional\t2\t1\t" in plan
    assert Path(f"{out}.controls.txt").read_text() == "/opt/res/controls.fa.gz\n"      # the FASTA of the regions costed
    r = run("--menu", menu_dir / "menu.tsv", "--list")
    assert r.returncode == 0 and "core\tA B\tA and B" in r.stdout
    r = run("--menu", menu_dir / "menu.tsv", "--classes", "A", "--capture-class", "A0.9")
    assert r.returncode != 0 and "expected CLASS=FRACTION" in r.stderr and "Traceback" not in r.stderr
    r = run("--menu", menu_dir / "menu.tsv", "--classes", "A", "--crai", menu_dir / "a.crai")
    assert r.returncode != 0 and "--crai needs --contigs" in r.stderr and "Traceback" not in r.stderr


def test_candidates_admitted_by_a_budget_are_scanned_and_the_table_says_so(menu_dir):
    m = fetchplan.read_menu(menu_dir / "menu.tsv")
    said = []
    p = fetchplan.make_plan(m, budget_mb=1, statuses=("shipped", "experimental", "candidate"), crais=[menu_dir / "a.crai"],
                            contigs=menu_dir / "ref.fai", log=said.append)
    assert "N" not in names(p) and [o.name for o in p.scan_only] == ["N"] and "n.tsv" in p.scan_panels and "n.tsv" not in p.panels
    assert fetchplan.table_lines(p)[-1].endswith("scan only: no sinks learned yet (panel in scan_panels.txt)")
    assert any("1 candidate option(s) admitted by --status" in x for x in said) and not any("are left out" in x for x in said)
    p = fetchplan.make_plan(m, budget_mb=1, crais=[menu_dir / "a.crai"], contigs=menu_dir / "ref.fai", log=lambda _: None)
    assert p.scan_only == [] and "n.tsv" not in p.scan_panels                   # default statuses: candidates left out


def test_the_plan_names_the_controls_fasta_of_the_regions_it_costed(menu_dir, tmp_path):
    m = fetchplan.read_menu(menu_dir / "menu.tsv")
    (menu_dir / "controls.fa.gz").write_bytes(b"")
    said = []
    p = fetchplan.make_plan(m, classes=["A"], log=said.append)
    assert p.controls_fasta == str(menu_dir / "controls.fa.gz") and not any("controls" in x for x in said)
    (menu_dir / "controls.lite.bed").write_text("c1\t50000\t50500\tcontrol\n")
    (menu_dir / "controls.lite.fa.gz").write_bytes(b"")                           # a set that ships its FASTA: named as it is
    p = fetchplan.make_plan(m, classes=["A"], controls=menu_dir / "controls.lite.bed", log=said.append)
    assert p.controls_fasta == str(menu_dir / "controls.lite.fa.gz") and p.rows[0].intervals == [("c1", 49400, 51100)]
    assert any("not the menu's controls.bed: the fetch must pass -c controls.lite.fa.gz" in x for x in said)
    assert not any("no controls FASTA" in x or "is cut from" in x for x in said)
    fetchplan.write(p, str(tmp_path / "q"), m)
    assert (tmp_path / "q.controls.txt").read_text() == f"{menu_dir / 'controls.lite.fa.gz'}\n"
    assert "# controls FASTA for the fetch (-c): controls.lite.fa.gz" in (tmp_path / "q.plan.tsv").read_text()
    assert (tmp_path / "q.count_flags.txt").read_text() == "--allow-missing-sinks\n"   # never a second -c: the engine refuses one
    # a set of the bundle without a FASTA of its own, taken by name: cut from the bundle's when the plan is written
    with gzip.open(menu_dir / "controls.fa.gz", "wt") as fh:
        fh.write(">c1:50000-50500 flank=2 role=control\n" + "ACGT" * 126 + "\n>c1:90000-90100 flank=2 role=control\n" + "ACGT" * 26 + "\n")
    (menu_dir / "controls.half.bed").write_text("c1\t50000\t50500\tcontrol\n")
    said.clear()
    p = fetchplan.make_plan(m, classes=["A"], controls="half", log=said.append)
    assert p.controls_bed == str(menu_dir / "controls.half.bed") and any("is cut from the bundle's controls.fa.gz" in x for x in said)
    fetchplan.write(p, str(tmp_path / "h"), m, panel_root="/opt/res")
    cut = (tmp_path / "h.controls.fa.gz").resolve()
    assert (tmp_path / "h.controls.txt").read_text() == f"{cut}\n" and "# controls FASTA for the fetch (-c): h.controls.fa.gz" in (tmp_path / "h.plan.tsv").read_text()
    with gzip.open(cut, "rt") as fh:
        assert [line[1:].split()[0] for line in fh if line.startswith(">")] == ["c1:50000-50500"]
    with pytest.raises(ValueError, match="no control set of that name"):
        fetchplan.make_plan(m, classes=["A"], controls="nosuch")
    # a BED outside the bundle, without a FASTA: still to be built from the reference
    said.clear()
    (tmp_path / "elsewhere.bed").write_text("c1\t50000\t50500\tcontrol\n")
    fetchplan.make_plan(m, classes=["A"], controls=tmp_path / "elsewhere.bed", log=said.append)
    assert any("no controls FASTA" in x and "ngs-dose controls -b" in x for x in said)
    (tmp_path / "own.bed").write_text("c1\t50000\t50500\tcontrol\n")               # outside the menu: an absolute path, even with a root
    fetchplan.write(fetchplan.make_plan(m, classes=["A"], controls=tmp_path / "own.bed", log=lambda _: None), str(tmp_path / "r"), m,
                    panel_root="/opt/res")
    assert (tmp_path / "r.controls.txt").read_text() == f"{(tmp_path / 'own.fa.gz').resolve()}\n"


def test_controls_txt_stays_under_the_menu_through_a_symlinked_bundle(menu_dir, tmp_path):
    """A bundle reached through a symlink (a site's, or a container's) resolves elsewhere: controls.txt, like panels.txt,
    names the path under the menu's directory as the menu writes it, which --panel-root maps."""
    real = tmp_path / "real"
    real.mkdir()
    for f in ("controls.bed", "ab.tsv", "sinks.bed"):
        (real / f).write_text((menu_dir / f).read_text())
    (real / "controls.fa.gz").write_bytes(b"")
    (menu_dir / "bundle").symlink_to(real)
    text = (menu_dir / "menu.tsv").read_text().replace("\t-\tcontrols.bed\t", "\t-\tbundle/controls.bed\t").replace("\tab.tsv\t", "\tbundle/ab.tsv\t")
    (menu_dir / "menu2.tsv").write_text(text)
    m = fetchplan.read_menu(menu_dir / "menu2.tsv")
    p = fetchplan.make_plan(m, classes=["A"], log=lambda _: None)
    assert p.controls_fasta == str(menu_dir / "bundle" / "controls.fa.gz")
    fetchplan.write(p, str(tmp_path / "q"), m, panel_root="/opt/res")
    assert (tmp_path / "q.panels.txt").read_text() == "/opt/res/bundle/ab.tsv\n"
    assert (tmp_path / "q.controls.txt").read_text() == "/opt/res/bundle/controls.fa.gz\n"
    fetchplan.write(p, str(tmp_path / "r"), m)
    assert (tmp_path / "r.controls.txt").read_text() == f"{menu_dir / 'bundle' / 'controls.fa.gz'}\n"      # no root: as written
    # a menu reached through a link and the controls named by their real path: under the menu by real path, so still mapped
    (real / "menu.tsv").write_text((menu_dir / "menu.tsv").read_text())
    link = tmp_path / "link"
    link.symlink_to(real)
    m = fetchplan.read_menu(link / "menu.tsv")
    q = fetchplan.make_plan(m, classes=["A"], controls=real / "controls.bed", log=lambda _: None)
    fetchplan.write(q, str(tmp_path / "s"), m, panel_root="/opt/res")
    assert (tmp_path / "s.controls.txt").read_text() == "/opt/res/controls.fa.gz\n"


def test_a_pad_other_than_the_engines_goes_to_the_count_flags(menu_dir, tmp_path):
    m = fetchplan.read_menu(menu_dir / "menu.tsv")
    p = fetchplan.make_plan(m, classes=["A"], pad=1500, log=lambda _: None)
    assert p.rows[0].intervals == [("c1", 48500, 52500)] and p.flags[0] == "--pad=1500"
    fetchplan.write(p, str(tmp_path / "q"), m)
    assert (tmp_path / "q.count_flags.txt").read_text() == "--pad=1500\n--allow-missing-sinks\n"
    assert "--pad=600" not in fetchplan.make_plan(m, classes=["A"], pad=600, log=lambda _: None).flags
    with pytest.raises(ValueError, match="below 400 bp"):
        fetchplan.make_plan(m, classes=["A"], pad=100)


def test_the_default_menu_is_found_beside_the_bundle(tmp_path, monkeypatch):
    (tmp_path / "GRCh38").mkdir()
    (tmp_path / "fetch_menu.tsv").write_text("x\n")
    monkeypatch.setenv("NGSDOSE_RESOURCES", str(tmp_path / "GRCh38"))
    assert fetchplan.default_menu() == tmp_path / "fetch_menu.tsv"
    monkeypatch.setenv("NGSDOSE_RESOURCES", str(tmp_path / "elsewhere" / "GRCh38"))
    with pytest.raises(ValueError, match="no fetch menu beside the resource bundle"):   # nothing there: not the checkout's, silently
        fetchplan.default_menu()
    monkeypatch.delenv("NGSDOSE_RESOURCES")
    assert fetchplan.default_menu() == ROOT / "resources" / "fetch_menu.tsv"         # no bundle configured: the source checkout's
    monkeypatch.setenv("NGSDOSE_RESOURCES", str(tmp_path / "elsewhere" / "GRCh38"))
    with pytest.raises(ValueError, match="pass it with --menu"):
        fetchplan.read_menu(tmp_path / "none.tsv")


def test_the_shipped_menu(tmp_path):
    """Every row of resources/fetch_menu.tsv names a panel that defines its class; every shipped class's
    sinks are in its sinks file; a plan of the core preset is what the bundle's sinks give for those classes."""
    m = fetchplan.read_menu()
    assert {"core", "core_tel", "truths", "satellites", "biobank_lite"} <= set(m.presets)
    for o in m.options.values():
        if o.panel:
            assert o.counted in fetchplan.panel_classes(m.resolve(o.panel)), o.name
        if o.status == "shipped" and o.kind in ("positional", "compositional"):
            assert any(r[3] == o.name for r in sinks.read_bed(m.resolve(o.sinks))), o.name
    core = fetchplan.make_plan(m, presets=["core_tel"], log=lambda _: None)
    want = sorted(r for r in sinks.read_bed(ROOT / "resources" / "GRCh38" / "sinks.bed") if r[3] in ("rDNA45S", "rDNA5S", "DJ", "TEL"))
    assert core.bed() == want and core.flags == [] and len(core.rows[0].intervals) == 982
    sat = m.options["HSat2"]
    if m.resolve(sat.sinks).exists():                      # staged: the statistics the menu names are for exactly these intervals
        p = fetchplan.make_plan(m, presets=["satellites"], capture=0.99, log=lambda _: None)
        got = {r.option.name: r for r in p.rows}
        assert set(got) >= {"unmapped", "HSat1A", "HSat1B", "HSat2", "HSat3", "aSatHOR", "bSat", "ACRO", "SST1", "CER", "SATR"}
        assert all(got[n].expected is not None for n in got if n not in ("controls", "unmapped"))
        assert "--unmapped" in p.flags and "--allow-missing-sinks" not in p.flags
    if BIN.exists():
        out = tmp_path / "core"
        fetchplan.write(core, str(out), m)
        r = subprocess.run([str(BIN), "plan", "-c", str(ROOT / "resources" / "GRCh38" / "controls.fa.gz"), "--sinks", f"{out}.sinks.bed",
                            "-o", str(tmp_path / "engine.bed")], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr


# ---------------------------------------------------------------- byte order, engine flags, provenance, subsets

def stat_rows(*rows):
    """Statistics rows in rank order: (start, share_median, share_max, cum_capture); p10 = median here."""
    return [dict(contig="c", start=s, end=s + 10, rank=i + 1, scans=5, share_median=m, share_p10=m, share_max=mx, cum_capture_median=cum,
                 cum_capture_p10=cum) for i, (s, m, mx, cum) in enumerate(rows)]


def test_per_byte_trimming_is_bounded_and_never_costlier_than_per_read():
    st = stat_rows((0, 0.6, 0.62, 0.6), (10, 0.3, 0.35, 0.9), (20, 0.09, 0.5, 0.99))
    nb = {("c", 0, 10): 100, ("c", 10, 20): 1000, ("c", 20, 30): 50}       # per byte: 0, 20, 10
    price = lambda ivs: sum(nb[iv] for iv in ivs)
    ivs = sorted(nb)
    # without bytes: the statistics' own curve, exact
    assert fetchplan._trim(ivs, st, 0.85, "p10") == ([("c", 0, 10), ("c", 10, 20)], 1, 0.9, 0.99, "read")
    # per byte, 0 and 20 are kept first; their bound is the larger of 0.99 - 0.35 (the largest share of 10) and 0.6 (the
    # curve up to 10, the first rank not kept): 0.64, short of 0.85, so all three are needed (1,150 bytes). Per read
    # reaches 0.85 with 1,100: that is kept, with its exact capture
    assert fetchplan._trim(ivs, st, 0.85, "p10", nb, price) == ([("c", 0, 10), ("c", 10, 20)], 1, 0.9, 0.99, "read, fewer bytes")
    kept, dropped, expected, full, order = fetchplan._trim(ivs, st, 0.62, "p10", nb, price)
    assert kept == [("c", 0, 10), ("c", 20, 30)] and order == "byte" and expected == pytest.approx(0.64) and dropped == 1
    assert fetchplan._trim(ivs, st, 0.6, "p10", nb, price)[:3] == ([("c", 0, 10)], 2, 0.6)
    assert fetchplan._trim(ivs, st, None, "p10", nb, price)[:4] == (ivs, 0, 0.99, 0.99)


def test_statistics_without_share_max_keep_the_order_per_read(menu_dir):
    old = [line for line in (menu_dir / "t.stats.tsv").read_text().splitlines()]
    cut = old[0].split("\t").index("share_max")
    (menu_dir / "t.stats.tsv").write_text("".join("\t".join(x for i, x in enumerate(line.split("\t")) if i != cut) + "\n" for line in old))
    m = fetchplan.read_menu(menu_dir / "menu.tsv")
    said = []
    p = fetchplan.make_plan(m, classes=["T"], capture=0.95, crais=[menu_dir / "a.crai"], contigs=menu_dir / "ref.fai", log=said.append)
    assert p.rows[1].order == "read" and len(p.rows[1].intervals) == 2 and any("no share_max" in x for x in said)


def test_an_interval_the_rest_of_the_plan_reads_costs_nothing():
    st = stat_rows((0, 0.6, 0.62, 0.6), (10, 0.3, 0.35, 0.9), (20, 0.09, 0.5, 0.99))
    free = frozenset({("c", 20, 30)})
    # alone, per read, 0.62 needs 0 and 10 (the curve: 0.6, then 0.9)
    assert fetchplan._trim(sorted(st_iv(st)), st, 0.62, "p10")[:3] == ([("c", 0, 10), ("c", 10, 20)], 1, 0.9)
    # with 20 free it goes first; 20 and 0 then bound 0.99 - 0.35 = 0.64, so 10 is dropped (a bound: not the statistics' own run)
    kept, dropped, expected, full, order = fetchplan._trim(sorted(st_iv(st)), st, 0.62, "p10", free=free)
    assert kept == [("c", 0, 10), ("c", 20, 30)] and dropped == 1 and expected == pytest.approx(0.64) and order == "read"
    # per byte, a free interval ranks first whatever its bytes
    nb = {("c", 0, 10): 100, ("c", 10, 20): 100, ("c", 20, 30): 10**9}
    assert fetchplan._trim(sorted(nb), st, 0.62, "p10", nb, free=free)[0] == [("c", 0, 10), ("c", 20, 30)]
    # the capture of any set: the curve for a run of the statistics' own order, else the bound
    assert fetchplan._capture_of(st, {("c", 0, 10), ("c", 10, 20)}, "p10") == pytest.approx(0.9)
    assert fetchplan._capture_of(st, {("c", 0, 10), ("c", 20, 30)}, "p10") == pytest.approx(0.64)
    assert fetchplan._capture_of(st, set(st_iv(st)), "p10") == pytest.approx(0.99)


def st_iv(st):
    return [(r["contig"], r["start"], r["end"]) for r in st]


@pytest.fixture
def shared_menu(menu_dir):
    """menu_dir plus U, a class with its own panel that has T's costly c3 interval (0.6 of its reads) and
    T's c1:400,000 interval (0.4): U needs both for any target above 0.62."""
    d = menu_dir
    panel(d / "u.tsv", ("U", "compositional"))
    ivs = [("c3", 0, 1000, "U"), ("c1", 400000, 401000, "U")]
    (d / "u.bed").write_text("".join("\t".join(map(str, x)) + "\n" for x in ivs))
    t = sinks.IntervalTally(ivs)
    for sh in ((600, 400), (620, 380)):
        pl = [dict(contig=c, start=s, reads=n, all=n, **{"class": "U"}) for (c, s, _, _), n in zip(ivs, sh)]
        t.add(dict(placement_bin=1000, placement_bin_compositional=1000, classes=[dict(name="U", kind="compositional")],
                   contigs=[], placements=pl))
    with open(d / "u.stats.tsv", "w") as fh:
        sinks.write_stats(fh, t.table())
    with open(d / "menu.tsv", "a") as fh:
        fh.write(row("U", "compositional", "experimental: no fetch compared yet", "C", "-", "u.tsv", "u.bed", "u.stats.tsv"))
    return d


def test_an_interval_another_class_keeps_is_kept_and_saves_nothing(shared_menu):
    """T alone drops its costly c3 interval at 0.95 (test_a_capture_target_drops_the_least_yielding_intervals). U keeps
    the same interval, so the fetch reads it anyway: T keeps it too, its expected capture is that of all its intervals,
    and trimming T saves the plan nothing."""
    m = fetchplan.read_menu(shared_menu / "menu.tsv")
    kw = dict(crais=[shared_menu / "a.crai", shared_menu / "b.crai"], contigs=shared_menu / "ref.fai")
    said = []
    p = fetchplan.make_plan(m, classes=["T", "U"], capture=0.95, log=said.append, **kw)
    whole = fetchplan.make_plan(m, classes=["T", "U"], log=lambda _: None, **kw)
    got = {r.option.name: r for r in p.rows}
    t, u = got["T"], got["U"]
    assert len(t.intervals) == 3 and t.dropped == 0 and t.shared == 1 and t.expected == pytest.approx(t.full)
    assert len(u.intervals) == 2 and u.dropped == 0
    assert t.mb_saved == 0 and u.mb_saved == 0 and p.rows[-1].cum_mb == pytest.approx(whole.rows[-1].cum_mb)
    assert "1 interval(s) its target alone would drop are kept: the plan's fetches read their CRAM slices anyway (for U)" in t.note
    assert any("capture targets saved 0.0 MB of the plan" in x for x in said)
    line = dict(zip(fetchplan.PLAN_COLUMNS, fetchplan.table_lines(p)[[r.option.name for r in p.rows].index("T") + 1].split("\t")))
    assert line["mb_saved"] == "0.0" and line["capture_lost"] == "0.00000" and line["dropped"] == "0"
    assert ("c3", 0, 1000, "T") in p.bed() and ("c3", 0, 1000, "U") in p.bed()
    # without indexes: the interval lies inside U's, so it is kept by coordinates alone
    q = {r.option.name: r for r in fetchplan.make_plan(m, classes=["T", "U"], capture=0.95, log=lambda _: None).rows}
    assert len(q["T"].intervals) == 3 and q["T"].shared == 1 and q["T"].mb_saved is None
    # U untrimmed (a target for T only) keeps it as well
    r = {x.option.name: x for x in fetchplan.make_plan(m, classes=["T", "U"], capture_class={"T": 0.95}, log=lambda _: None, **kw).rows}
    assert len(r["T"].intervals) == 3 and r["U"].mb_saved is None


def test_an_interval_every_class_drops_is_a_saving_of_each(shared_menu):
    """With U at 0.3 (per byte it keeps only its cheap c1 interval: 1 - 0.62, the largest c3 share, is 0.38) and T at
    0.95, both drop c3: its slice is not read, and each option's mb_saved counts it, since keeping either class whole
    would read it; the note says the options' savings need not add up."""
    m = fetchplan.read_menu(shared_menu / "menu.tsv")
    kw = dict(crais=[shared_menu / "a.crai"], contigs=shared_menu / "ref.fai")
    said = []
    p = fetchplan.make_plan(m, classes=["T"], capture=0.95, log=said.append, **kw)
    t = {r.option.name: r for r in p.rows}["T"]
    assert t.shared == 0 and t.mb_saved == pytest.approx(0.10001) and any("capture targets saved 0.1 MB of the plan" in x for x in said)
    p = fetchplan.make_plan(m, classes=["T", "U"], capture_class={"T": 0.95, "U": 0.3}, log=said.append, **kw)
    got = {r.option.name: r for r in p.rows}
    assert got["U"].intervals == [("c1", 400000, 401000)] and got["T"].intervals == [("c1", 300000, 301000), ("c1", 400000, 401000)]
    assert got["T"].mb_saved == pytest.approx(0.10001) and got["U"].mb_saved == pytest.approx(0.10001)
    assert any("need not add up to the plan's" in x for x in said)


def fake_engine(path, *options):
    """An engine whose `count --help` lists `options`, one per line as clap prints them, after a prose line
    that names --classes and --allow-missing-sinks without offering them."""
    path.write_text("#!/bin/sh\necho 'Usage: ngs-dose count [OPTIONS]'\n"
                    "echo '          a subset BED needs no --allow-missing-sinks; see --classes'\n"
                    + "".join(f"echo '  {o}'\n" for o in options))
    path.chmod(0o755)
    return path


def test_count_flags_follow_what_the_engine_takes(menu_dir, tmp_path):
    m = fetchplan.read_menu(menu_dir / "menu.tsv")
    new = fake_engine(tmp_path / "new", "-p, --panel <PANEL>", "    --classes <CLASSES>", "    --pad <PAD>", "    --allow-missing-sinks")
    mid = fake_engine(tmp_path / "mid", "    --pad <PAD>", "    --allow-missing-sinks")            # as 7772e32
    old = fake_engine(tmp_path / "old", "-p, --panel <PANEL>", "    --unmapped")                     # as fae1124
    assert fetchplan.probe_engine(new)[0] == {"--panel", "--classes", "--pad", "--allow-missing-sinks"}
    assert fetchplan.probe_engine(old)[0] == {"--panel", "--unmapped"}                               # prose is not an option
    assert fetchplan.probe_engine(tmp_path / "missing")[0] is None and fetchplan.probe_engine(None)[0] is None
    said = []
    p = fetchplan.make_plan(m, classes=["A", "T"], engine=new, log=said.append)
    assert p.flags == ["--classes=A,T"] and p.classes_flag.startswith("--classes (") and "lists --classes" in p.classes_flag
    assert any("count_flags.txt has --classes=A,T" in x and "also define B" in x for x in said)
    fetchplan.write(p, str(tmp_path / "q"), m, header="h")
    assert (tmp_path / "q.count_flags.txt").read_text() == "--classes=A,T\n"
    assert "# count flags for the classes the panels define but the plan does not select: --classes" in (tmp_path / "q.plan.tsv").read_text()
    p = fetchplan.make_plan(m, classes=["A", "T"], engine=mid, log=lambda _: None)
    assert p.flags == ["--allow-missing-sinks"] and "lists it, not --classes" in p.classes_flag
    # neither flag (fae1124): the companions would be counted and not marked; refused unless asked for, then loudly
    with pytest.raises(ValueError, match="neither --classes nor --allow-missing-sinks.*--unmarked-companions"):
        fetchplan.make_plan(m, classes=["A", "T"], engine=old, log=lambda _: None)
    said = []
    p = fetchplan.make_plan(m, classes=["A", "T"], engine=old, unmarked_companions=True, log=said.append)
    assert p.flags == [] and p.classes_flag.startswith("none (") and "--unmarked-companions" in p.classes_flag
    assert any("] WARNING: the panels loaded also define B" in x and "--fetch-sinks" in x for x in said)
    fetchplan.write(p, str(tmp_path / "u"), m, header="h")
    assert (tmp_path / "u.count_flags.txt").read_text() == ""
    assert fetchplan.make_plan(m, classes=["A", "B"], engine=old, log=lambda _: None).flags == []     # no companions: nothing to ask
    # a flag the engine does not take is refused, not left for the fetch to fail on
    with pytest.raises(ValueError, match="needs --pad=1500 .* does not take"):
        fetchplan.make_plan(m, classes=["A", "B"], pad=1500, engine=old, log=lambda _: None)
    assert fetchplan.make_plan(m, classes=["A", "B"], pad=1500, engine=new, log=lambda _: None).flags == ["--pad=1500"]
    # no engine to probe: the flag engines since 645ae55 take, and the note says fae1124 refuses it
    said = []
    p = fetchplan.make_plan(m, classes=["A"], engine=tmp_path / "missing", log=said.append)
    assert p.flags == ["--allow-missing-sinks"] and "no engine" in p.classes_flag and "fae1124 does not" in p.classes_flag
    assert any("fae1124 does not and refuses the fetch" in x for x in said)


def test_the_cli_refuses_a_plan_the_engine_cannot_run_and_records_the_override(menu_dir, tmp_path):
    old = fake_engine(tmp_path / "old", "    --unmapped")
    r = run("--menu", menu_dir / "menu.tsv", "--classes", "A", "--engine", old, "-o", tmp_path / "x")
    assert r.returncode != 0 and "neither --classes nor --allow-missing-sinks" in r.stderr and "Traceback" not in r.stderr
    assert not (tmp_path / "x.count_flags.txt").exists()
    r = run("--menu", menu_dir / "menu.tsv", "--classes", "A", "--engine", old, "--unmarked-companions", "-o", tmp_path / "x")
    assert r.returncode == 0, r.stderr
    head = (tmp_path / "x.plan.tsv").read_text()
    assert "--unmarked-companions" in head.splitlines()[0] and "count flags: none (" in head.splitlines()[0]
    assert "[fetchplan] WARNING: the panels loaded also define B" in r.stderr


def test_a_known_sinks_bed_marks_the_classes_it_does_not_give(tmp_path):
    """fae1124 counts a class a fetch's BED gives no interval without recording it; the BED, known by its
    SHA-256, tells."""
    bed = tmp_path / "plan.sinks.bed"
    bed.write_text("# h\nc1\t0\t1000\tA\n")
    known = estimate.sinks_by_hash([bed])
    sha = next(iter(known))
    counts = dict(mode="fetch", sinks="plan.sinks.bed", sinks_sha256=sha, classes=[dict(name="A"), dict(name="B")])
    assert contract.missing_from_bed(counts, known) == ["B"] and list(contract.incomplete_sinks(counts, known)) == ["B"]
    assert "fae1124 does not" in contract.incomplete_sinks(counts, known)["B"]
    assert contract.incomplete_sinks(counts) == {}                                     # without the BED nothing tells
    assert contract.missing_from_bed(dict(counts, mode="scan"), known) == []
    assert contract.missing_from_bed(dict(counts, sinks_sha256="other"), known) == []
    plain = tmp_path / "plain.bed"
    plain.write_text("c1\t0\t1000\n")                                                  # names no class: serves every class
    plain_known = estimate.sinks_by_hash([plain])
    assert contract.missing_from_bed(dict(counts, sinks_sha256=next(iter(plain_known))), plain_known) == []
    engine_said = dict(counts, sinks_missing_classes=["B"])                             # an engine that records it keeps its own words
    assert contract.incomplete_sinks(engine_said, known)["B"] == "no sink intervals in the fetch"


# The real engines: this checkout's and those in NGSDOSE_EXTRA_BINS (':'-separated), e.g. the cohort's fae1124
# and 7772e32. Each gets a plan it can run, and its fetch must say which classes it did not read in full.
ENGINES = [BIN] + [Path(x) for x in os.environ.get("NGSDOSE_EXTRA_BINS", "").split(os.pathsep) if x]
FIXTURE = ROOT / "tests" / "data" / "NA12878.subsample.bam"


@pytest.mark.parametrize("eng", ENGINES, ids=lambda e: e.name)
def test_a_plan_runs_on_the_engine_it_was_made_for(eng, tmp_path):
    if not eng.exists():
        pytest.skip(f"no engine {eng}")
    m = fetchplan.read_menu()
    kw = dict(presets=["core_tel", "xy_arrays"], controls=ROOT / "resources" / "GRCh38" / "controls.lite200.bed", engine=eng,
              log=lambda _: None)
    if not m.resolve(m.options["aSatHOR"].panel).exists():
        pytest.skip("the satellite panel is not staged")
    opts, _ = fetchplan.probe_engine(eng)
    marks = "--classes" in opts or "--allow-missing-sinks" in opts
    if not marks:                                          # fae1124
        with pytest.raises(ValueError, match="neither --classes nor --allow-missing-sinks"):
            fetchplan.make_plan(m, **kw)
        kw["unmarked_companions"] = True
    if "--pad" not in opts:
        with pytest.raises(ValueError, match="needs --pad=1500"):
            fetchplan.make_plan(m, **dict(kw, pad=1500))
    p = fetchplan.make_plan(m, **kw)
    pre = tmp_path / "xy"
    fetchplan.write(p, str(pre), m)
    flags = Path(f"{pre}.count_flags.txt").read_text().split()
    panels = [a for x in Path(f"{pre}.panels.txt").read_text().split() for a in ("-p", x)]
    out = tmp_path / "fetch.json"
    # the fixture lacks three decoy contigs that carry DYZ2 (HSat1B) sinks: an engine since 0.1.1 refuses to fetch a
    # class that loses any sink interval to them unless told to go on (the plan itself is for inputs that have them)
    lacks = ["--allow-missing-sinks"] if "--allow-missing-sinks" in opts and "--allow-missing-sinks" not in flags else []
    r = subprocess.run([str(eng), "count", "-m", "fetch", "-i", str(FIXTURE), *panels, "-c", Path(f"{pre}.controls.txt").read_text().strip(),
                        "--sinks", f"{pre}.sinks.bed", *flags, *lacks, "-o", str(out)], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr
    counts = json.loads(out.read_text())
    selected = {"rDNA45S", "rDNA5S", "DJ", "TEL", "aSatHOR", "HSat3", "HSat1B"}
    companions = set(fetchplan.panel_classes(m.resolve(m.options["aSatHOR"].panel))) - selected
    got = {c["name"] for c in counts["classes"]}
    known = estimate.sinks_by_hash([f"{pre}.sinks.bed"])
    if "--classes" in opts:
        assert flags[-1].startswith("--classes=") and got == selected
    else:
        assert got == selected | companions
        assert companions <= set(contract.incomplete_sinks(counts, known))          # the plan's BED marks them, whatever the engine said
        if marks:
            assert companions <= set(counts["sinks_missing_classes"]) <= set(contract.incomplete_sinks(counts))
        else:
            assert not companions & set(contract.incomplete_sinks(counts))            # fae1124 records nothing: only the BED tells


def test_the_provenance_header_records_every_option(menu_dir, tmp_path):
    (menu_dir / "controls.lite.bed").write_text("c1\t50000\t50500\tcontrol\n")
    base = ["--menu", menu_dir / "menu.tsv", "--preset", "core", "--classes", "T", "--crai", menu_dir / "a.crai", menu_dir / "b.crai",
            "--contigs", menu_dir / "ref.fai", "--engine", tmp_path / "none"]
    heads = []
    for extra in (["--capture", "0.95"], ["--capture", "0.95", "--capture-stat", "median"], ["--capture", "0.95", "--stats", menu_dir / "t.stats.tsv"],
                  ["--sinks", menu_dir / "sinks.bed", "--controls", menu_dir / "controls.lite.bed", "--pad", "1500"],
                  ["--budget-mb", "1", "--status", "shipped,experimental,candidate", "--fill"]):
        out = tmp_path / f"p{len(heads)}"
        r = run(*base, *extra, "-o", out)
        assert r.returncode == 0, r.stderr
        head = Path(f"{out}.sinks.bed").read_text().splitlines()[0]
        assert head == Path(f"{out}.plan.tsv").read_text().splitlines()[0]
        heads.append(head)
    assert len(set(heads)) == len(heads)
    assert "--capture-stat p10" in heads[0] and "--capture-stat median" in heads[1] and f"--stats {menu_dir / 't.stats.tsv'}" in heads[2]
    assert f"--sinks {menu_dir / 'sinks.bed'}" in heads[3] and "--controls" in heads[3] and "--pad 1500" in heads[3]
    assert "--budget-mb 1 --status shipped,experimental,candidate --fill" in heads[4]
    assert all("--crai (2 index(es): a.crai b.crai)" in h and f"--menu {menu_dir / 'menu.tsv'}" in h and "--engine" in h for h in heads)
    assert "count flags" not in heads[0]                   # A, B and T: every class of the panels loaded is selected


def test_the_expected_capture_says_whether_its_scans_were_held_out(menu_dir, tmp_path):
    m = fetchplan.read_menu(menu_dir / "menu.tsv")
    text = (menu_dir / "t.stats.tsv").read_text()
    for line, want in (("", "of unrecorded provenance for T (2 scans)"), (sinks.HELD_OUT + "yes: none\n", "held-out for T (2 scans)"),
                       (sinks.HELD_OUT + "no: in-sample\n", "in-sample for T (2 scans)")):
        (menu_dir / "t.stats.tsv").write_text(line + text)
        said = []
        p = fetchplan.make_plan(m, classes=["T"], capture=0.9, log=said.append)
        assert want in p.capture_note and any(want in x for x in said)
        fetchplan.write(p, str(tmp_path / "q"), m)
        assert want in (tmp_path / "q.plan.tsv").read_text()


def test_sinks_statistics_record_share_max_and_whether_the_scans_were_held_out(tmp_path):
    def scan(name, reads):
        return dict(format="ngs-dose-counts/1", sample=name, mode="scan", primary=1_000_000, ctrl_reads=3000, placement_bin=1000,
                    placement_bin_compositional=10_000, classes=[dict(name="TEL", kind="compositional")], contigs=[dict(name="c1", len=100_000)],
                    placements=[{"class": "TEL", "contig": "c1", "start": s, "reads": n, "all": n} for s, n in reads])
    paths = []
    for i, reads in enumerate(([(0, 900), (50_000, 100)], [(0, 700), (50_000, 300)], [(0, 950), (50_000, 50)])):
        paths.append(tmp_path / f"s{i}.json")
        paths[-1].write_text(json.dumps(scan(f"s{i}", reads)))
    (tmp_path / "s.bed").write_text("c1\t0\t10000\tTEL\nc1\t50000\t60000\tTEL\n")
    sk = lambda *a: subprocess.run([sys.executable, "-m", "ngsdose", "sinks", *map(str, a)], capture_output=True, text=True,
                                   env={**os.environ, "PYTHONPATH": str(ROOT)})
    r = sk(*paths, "--evaluate", tmp_path / "s.bed", "--held-out", "--stats", tmp_path / "h.tsv", "--stats-note", "made for a test")
    assert r.returncode == 0, r.stderr
    st = sinks.read_stats(tmp_path / "h.tsv")["TEL"]
    assert [x["share_max"] for x in st] == [pytest.approx(0.95), pytest.approx(0.3)] and sinks.stats_held_out(tmp_path / "h.tsv") is True
    assert "# made for a test\n" in (tmp_path / "h.tsv").read_text()
    assert sk(*paths, "--evaluate", tmp_path / "s.bed", "--stats", tmp_path / "u.tsv").returncode == 0
    assert sinks.stats_held_out(tmp_path / "u.tsv") is None
    assert sk(*paths, "--classes", "TEL", "-o", tmp_path / "l.bed", "--stats", tmp_path / "l.tsv").returncode == 0
    assert sinks.stats_held_out(tmp_path / "l.tsv") is False
    r = sk(*paths, "--classes", "TEL", "-o", tmp_path / "l.bed", "--held-out", "--stats", tmp_path / "l.tsv")
    assert r.returncode != 0 and "in-sample" in r.stderr


@pytest.fixture
def subset_menu(menu_dir):
    """The test menu with S, a named subset of T's intervals (c1:300,000 only)."""
    (menu_dir / "s.bed").write_text("# S: part of T\nc1\t300000\t301000\tT\n")
    with open(menu_dir / "menu.tsv", "a") as fh:
        fh.write(row("S", "subset", "experimental", "C", "-", "t.tsv", "s.bed"))
    return menu_dir


def test_a_subset_is_fetched_as_its_class_and_says_the_class_is_not_measured(subset_menu, tmp_path):
    m = fetchplan.read_menu(subset_menu / "menu.tsv")
    assert m.options["S"].parent == "T" and m.options["S"].counted == "T" and "subset of T" in "\n".join(fetchplan.list_menu(m))
    said = []
    p = fetchplan.make_plan(m, classes=["S"], crais=[subset_menu / "a.crai"], contigs=subset_menu / "ref.fai", log=said.append)
    assert p.bed() == [("c1", 300000, 301000, "T")] and p.panels == ["t.tsv"] and p.flags == []
    assert p.rows[1].mb == pytest.approx(3.01e-3) and p.rows[1].expected is None
    assert any("S is a named subset of T's intervals, fetched without the rest of T's sinks" in x and "subset_only" in x for x in said)
    both = fetchplan.make_plan(m, classes=["S", "T"], log=said.append)
    assert both.bed() == [("c1", 300000, 301000, "T"), ("c1", 400000, 401000, "T"), ("c3", 0, 1000, "T")]
    assert any("S (a named subset of T's intervals) and T are both fetched" in x for x in said)
    with pytest.raises(ValueError, match="--capture-class S: S is a named subset"):
        fetchplan.make_plan(m, classes=["S"], capture_class={"S": 0.9})
    assert fetchplan.make_plan(m, classes=["S"], capture=0.5, log=lambda _: None).rows[1].intervals == [("c1", 300000, 301000)]
    # with another pipeline's sinks the subset holds only where they hold its class
    (tmp_path / "other.bed").write_text("c1\t299000\t302000\tT\n")
    assert fetchplan.make_plan(m, classes=["S"], sinks_files=[tmp_path / "other.bed"], log=lambda _: None).bed() == [("c1", 300000, 301000, "T")]
    (tmp_path / "other.bed").write_text("c1\t400000\t401000\tT\n")
    with pytest.raises(ValueError, match="1 of its 1 intervals .* are not inside T's intervals"):
        fetchplan.make_plan(m, classes=["S"], sinks_files=[tmp_path / "other.bed"])
    # the menu's own sinks for T must hold it too (a menu whose T row names another pipeline's sinks)
    (subset_menu / "sinks.bed").write_text("c1\t100000\t101000\tA\nc2\t0\t1000\tB\nc1\t400000\t401000\tT\n")
    with pytest.raises(ValueError, match="not inside T's intervals of .*sinks.bed"):
        fetchplan.make_plan(m, classes=["S"])


def test_a_subset_row_is_checked(subset_menu, tmp_path):
    text = (subset_menu / "menu.tsv").read_text()
    for bed, change, msg in (("c1\t0\t10\tT\nc1\t20\t30\tA\n", None, "names one class"),
                             ("c1\t0\t10\tS\n", None, "a subset of itself"),
                             (None, ("\tt.tsv\ts.bed\tmeasures S\t-\t-\t-", "\tt.tsv\ts.bed\tmeasures S\t-\t-\tt.stats.tsv"), "takes no statistics"),
                             (None, ("S\tg\tsubset\texperimental", "S\tg\tsubset\tcandidate"), "cannot be a candidate"),
                             (None, ("\ts.bed\tmeasures S", "\tnone.bed\tmeasures S"), "does not exist")):
        (tmp_path / "s.bed").write_text(bed or "c1\t300000\t301000\tT\n")
        t = text.replace("\ts.bed\t", f"\t{tmp_path / 's.bed'}\t")
        if change:
            t = t.replace(change[0].replace("s.bed", str(tmp_path / "s.bed")), change[1].replace("s.bed", str(tmp_path / "s.bed")))
        (subset_menu / "bad.tsv").write_text(t)
        with pytest.raises(ValueError, match=msg):
            fetchplan.read_menu(subset_menu / "bad.tsv")


def test_a_control_set_is_taken_by_name_and_its_fasta_cut_from_the_bundles(tmp_path):
    """`--controls NAME`: one of the sets beside the bundle's controls; a set without a FASTA of its own gets one, cut
    from the bundle's controls.fa.gz record by record when the plan is written."""
    d = tmp_path / "GRCh38"
    d.mkdir()
    recs = [("c1:100-200", "role=control"), ("c1:5000-5100", "role=control"), ("c2:10-60", "role=test label=x")]
    with gzip.open(d / "controls.fa.gz", "wt") as fh:
        for n, tags in recs:
            a, b = map(int, n.split(":")[1].split("-"))
            fh.write(f">{n} flank=5 {tags}\n" + "ACGT" * ((b - a + 10) // 4) + "AC"[: (b - a + 10) % 4] + "\n")
    (d / "controls.bed").write_text("c1\t100\t200\tcontrol\nc1\t5000\t5100\tcontrol\nc2\t10\t60\ttest:x\n")
    (d / "controls.small.bed").write_text("c1\t5000\t5100\tcontrol\nc2\t10\t60\ttest:x\n")
    sets = fetchplan.control_sets(d / "controls.bed")
    assert set(sets) == {"all", "small"} and fetchplan.resolve_controls("small", d / "controls.bed") == d / "controls.small.bed"
    assert fetchplan.resolve_controls(str(d / "controls.small.bed"), d / "controls.bed") == d / "controls.small.bed"
    with pytest.raises(ValueError, match="no control set of that name"):
        fetchplan.resolve_controls("karyotype", d / "controls.bed")
    out = tmp_path / "plan.controls.fa.gz"
    assert fetchplan.subset_fasta(d / "controls.fa.gz", ["c2:10-60", "c1:5000-5100"], out) == 2
    with gzip.open(out, "rt") as fh:
        got = fh.read()
    assert [line[1:].split()[0] for line in got.splitlines() if line.startswith(">")] == ["c1:5000-5100", "c2:10-60"]      # the source's order
    first = out.read_bytes()
    fetchplan.subset_fasta(d / "controls.fa.gz", ["c1:5000-5100", "c2:10-60"], out)
    assert out.read_bytes() == first                                         # the same bytes every time
    with pytest.raises(ValueError, match="no record for 1 of the 2 regions"):
        fetchplan.subset_fasta(d / "controls.fa.gz", ["c1:5000-5100", "c9:1-2"], tmp_path / "bad.fa.gz")
    assert not (tmp_path / "bad.fa.gz").exists() and not (tmp_path / "bad.fa.gz.part").exists()
