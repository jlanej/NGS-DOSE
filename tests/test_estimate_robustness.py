"""The estimator refuses to report what the counts do not support, and says why.

No engine and no data: counts are simulated with control regions on 22 autosomes, a positional
class and a compositional one, and then altered the way a real file can differ from a complete,
matching one (a fetch without sinks, another panel, a class the bundle lacks, an aneuploid
chromosome, a retired region, a reordered controls file, a truncated file or cache)."""
import copy
import gzip
import json
import shutil
import zlib
from pathlib import Path

import numpy as np
import pytest

from ngsdose import contract, estimate, io, resources, selftest
from ngsdose.io import Panel, PanelClass
from ngsdose.tables import load_result, summary_row, write_table

ROOT = Path(__file__).resolve().parents[1]
L, R = 400, 150
REGIONS_PER_CHROM = {f"chr{i}": (2 if i == 22 else 30) for i in range(1, 23)}


def _sim(seed=0, scale=None, depth=0.1):
    """Counts, the bundle-style (name, role) list and per-region GC tables in that order."""
    rng = np.random.default_rng(seed)
    g = np.arange(101) / 100
    lam = depth * selftest._bias(g, 0.3)
    rows, regs, names = [], [], []
    for chrom, n in REGIONS_PER_CHROM.items():
        for i in range(n):
            h = np.bincount(np.clip(rng.normal(42, 5, 2 * 10000).astype(int), 15, 75), minlength=101).astype(float)
            f = (scale or {}).get(chrom, 1.0) * np.exp(rng.normal(0, 0.02))
            o = rng.poisson(h * lam * f)
            name = f"{chrom}:{i * 100000}-{i * 100000 + 10000}"
            rows.append((h, o))
            regs.append(dict(name=name, role="control", label="", len=10000, obs=int(o.sum())))
            names.append((name, "control"))
    for label, cn in (("chrM", 400.0), ("auto", 2.0)):
        h = np.bincount(np.clip(rng.normal(44, 4, 2 * 5000).astype(int), 15, 75), minlength=101).astype(float)
        o = rng.poisson(h * lam * cn / 2)
        regs.append(dict(name=f"{label}:1-5001", role="dosage" if label == "chrM" else "test", label=label, len=5000, obs=int(o.sum())))
        names.append((regs[-1]["name"], regs[-1]["role"]))
        rows.append((h, o * 0))                                            # not a control: not in the GC table
    N = sum(h for (h, _), (_, role) in zip(rows, names) if role == "control")
    O = sum(o for (_, o), (_, role) in zip(rows, names) if role == "control")
    seq = selftest._random_unit(rng, n=3000)
    pc = PanelClass("unit", "positional", len(seq), True, np.arange(len(seq)))
    unit = selftest._simulate_sample(rng, seq, pc, 31, R, L, 100.0, strength=0.3, depth_rate=depth)["classes"][0]
    unit.update(length=len(seq), panel_kmers=len(seq), circular=True)
    comp = dict(name="sat", kind="compositional", reads=5000, dup_flagged=0,
                gc_read=np.bincount(np.clip(rng.normal(40, 5, 5000).astype(int), 0, 100), minlength=101).tolist())
    counts = dict(format=io.COUNTS_FORMAT, sample="sim", mode="fetch", engine_version="sim", k=31, insert_median=L,
                  read_length_mode=R, ctrl_reads=int(O.sum()), ctrl_dup_flagged=0, ctrl_mapq0=0, regions=regs,
                  gc_tables=[dict(l=L, n=N.tolist(), o=O.tolist()), dict(l=R, n=N.tolist(), o=O.tolist())],
                  classes=[unit, comp], eof_marker="present", unmapped_fetched=False)
    return counts, names, np.array([h for h, _ in rows]), Panel(31, {"unit": pc}), {"unit": seq}


@pytest.fixture(scope="module")
def sim():
    return _sim()


def run(counts, names, tables, panel, units, **kw):
    return estimate.estimate_sample(counts, panel, units, region_tables=tables, regions=names, L=L, **kw)


def test_a_complete_matching_file_is_estimated(sim):
    r = run(*sim)
    assert r["classes"]["unit"]["status"] == r["classes"]["sat"]["status"] == "ok"
    assert abs(r["classes"]["unit"]["cn"] / 100 - 1) < 0.05 and r["classes"]["sat"]["mass_Mb"] > 0
    assert r["skipped_classes"] == [] and r["control_qc"]["flagged_chromosomes"] == []
    assert r["untestable_chromosomes"] == ["chr22"] and abs(r["truth_regions"]["auto"]["cn"] - 2) < 0.1
    assert contract.issues(sim[0]) == []
    row = summary_row(r)
    assert row["unit.status"] == "ok" and row["unit.cn_basis"] in ("anchor", "all") and row["unit.n_anchor"] >= 0
    assert row["untestable_chromosomes"] == "chr22" and row["unmapped_fetched"] is False


@pytest.mark.parametrize("field,status", [("sinks_missing_classes", "no_sinks_in_fetch"), ("sinks_skipped", "sinks_skipped")])
def test_a_class_the_fetch_did_not_read_in_full_is_not_measured(sim, field, status, tmp_path):
    counts = copy.deepcopy(sim[0])
    counts[field] = ["unit", "sat"] if field == "sinks_missing_classes" else {n: dict(intervals=2, bp=20000) for n in ("unit", "sat")}
    r = run(counts, *sim[1:])
    u, s = r["classes"]["unit"], r["classes"]["sat"]
    assert u["status"] == s["status"] == status and u["reason"] and s["reason"]
    assert np.isnan(u["cn"]) and np.isnan(u["cn_anchor"]) and np.isnan(u["cn_all"]) and np.isnan(s["mass_Mb"])
    assert u["reads"] == sim[0]["classes"][0]["reads"]                                 # the counts stay visible
    # the window layout is the one a complete sample has, with no copy number in any window
    full = run(*sim)["classes"]["unit"]["windows"]
    assert [(w["start"], w["end"]) for w in u["windows"]] == [(w["start"], w["end"]) for w in full]
    assert all(w["cn"] is None for w in u["windows"])
    assert contract.incomplete_sinks(counts).keys() == {"unit", "sat"}
    assert [lvl for lvl, _ in contract.issues(counts)] == ["warn", "warn"]
    write_table([summary_row(r)], tmp_path / "t.tsv")
    head, vals = [l.split("\t") for l in (tmp_path / "t.tsv").read_text().splitlines()]
    row = dict(zip(head, vals))
    assert row["unit.cn_single"] == row["sat.mass_Mb"] == row["unit.cn_basis"] == "NA" and row["unit.status"] == status


def test_dropped_sink_rows_without_a_class_mark_every_class(sim):
    """A sinks BED row without a class serves every class; the engine records its loss under ''."""
    counts = copy.deepcopy(sim[0])
    counts["sinks_skipped"] = {"": dict(intervals=5, bp=500000), "unit": dict(intervals=1, bp=1000)}
    r = run(counts, *sim[1:])
    assert {n: c["status"] for n, c in r["classes"].items()} == {"unit": "sinks_skipped", "sat": "sinks_skipped"}
    assert np.isnan(r["classes"]["unit"]["cn"]) and np.isnan(r["classes"]["sat"]["mass_Mb"])
    why = contract.incomplete_sinks(counts)
    assert why.keys() == {"unit", "sat"} and "without a class" in why["sat"] and "1 sink intervals" in why["unit"]
    msgs = [m for _, m in contract.issues(counts)]
    assert len(msgs) == 2 and all(m.startswith(("sim: unit ", "sim: sat ")) for m in msgs)


def test_counts_made_with_another_panel_are_refused_per_class(sim):
    for change in (dict(panel_kmers=2000), dict(length=2999), None):
        counts = copy.deepcopy(sim[0])
        if change:
            counts["classes"][0].update(change)
        else:
            counts["k"] = 25
        r = run(counts, *sim[1:])
        assert r["classes"]["unit"]["status"] == "panel_mismatch" and np.isnan(r["classes"]["unit"]["cn"]), change
        assert "different panel" in r["classes"]["unit"]["reason"]
        assert r["classes"]["sat"]["status"] == "ok"                                   # compositional: no positions used
    # files written without these fields (the simulator's, or older engines') are not second-guessed
    counts = copy.deepcopy(sim[0])
    for key in ("panel_kmers", "length"):
        del counts["classes"][0][key]
    del counts["k"]
    assert run(counts, *sim[1:])["classes"]["unit"]["status"] == "ok"


def test_a_positional_class_the_bundle_cannot_estimate_is_listed_not_dropped(sim):
    counts = copy.deepcopy(sim[0])
    extra = copy.deepcopy(counts["classes"][0])
    extra["name"] = "NEWPOS"
    counts["classes"].append(extra)
    r = run(counts, sim[1], sim[2], sim[3], sim[4])
    assert "NEWPOS" not in r["classes"]
    assert r["skipped_classes"] == [dict(name="NEWPOS", kind="positional", reads=extra["reads"], reason="not in the bundle panel")]
    assert summary_row(r)["NEWPOS.status"] == "skipped: not in the bundle panel"
    r = run(sim[0], sim[1], sim[2], sim[3], {})
    assert r["skipped_classes"][0]["reason"] == "no unit sequence in the bundle" and set(r["classes"]) == {"sat"}


@pytest.mark.parametrize("factor", [1.2, 1.5, 0.5])
def test_a_whole_chromosome_gain_or_loss_is_flagged(factor):
    """A trisomy or monosomy moves every region of the chromosome past the region trim; the
    chromosome must still be tested, on all of its regions."""
    r = run(*_sim(seed=1, scale={"chr12": factor}))
    qc = r["control_qc"]
    assert qc["flagged_chromosomes"] == ["chr12"], (factor, qc["flagged_chromosomes"])
    n = sum(REGIONS_PER_CHROM.values())
    assert qc["n_flagged_regions"] >= 30 and abs(qc["rescale"] * (1 + (factor - 1) * 30 / n) - 1) < 0.01, qc["rescale"]
    assert summary_row(r)["flagged_chromosomes"] == "chr12"


def test_a_retired_dosage_region_is_left_out_and_order_does_not_matter(sim):
    counts, names, tables, panel, units = sim
    full = run(*sim)
    # the bundle no longer has the chrM region: old counts still estimate, without it
    keep = [i for i, (n, _) in enumerate(names) if not n.startswith("chrM")]
    r = run(counts, [names[i] for i in keep], tables[keep], panel, units)
    assert "chrM" not in r["truth_regions"] and r["regions_not_in_bundle"] == ["chrM:1-5001"]
    assert r["control_qc"] == full["control_qc"] and r["classes"]["unit"]["cn"] == full["classes"]["unit"]["cn"]
    # a control the bundle lacks is another controls file
    bad = copy.deepcopy(counts)
    bad["regions"][0]["name"] = "chr1:1-2"
    with pytest.raises(ValueError, match="different controls file"):
        run(bad, names, tables, panel, units)
    # counts written in another region order give the same vector, in the bundle's order
    shuffled = copy.deepcopy(counts)
    np.random.default_rng(1).shuffle(shuffled["regions"])
    r = run(shuffled, names, tables, panel, units)
    assert r["control_qc"]["region_log_ratio"] == full["control_qc"]["region_log_ratio"]
    assert r["control_qc"]["region_order_sha256"] == full["control_qc"]["region_order_sha256"]


def test_eof_and_panel_hash_are_reported_by_the_contract(sim, tmp_path):
    counts = copy.deepcopy(sim[0])
    counts["eof_marker"] = "unchecked"
    assert any("could not be checked" in m for lvl, m in contract.issues(counts) if lvl == "warn")
    B = resources.Bundle(ROOT / "resources" / "GRCh38")
    real = io.load_panel(B.panel)
    ok = dict(format=io.COUNTS_FORMAT, sample="x", k=real.k, panel_sha256=sorted(contract.panel_hashes(B.panel))[:1],
              regions=[dict(name=n, role=role) for n, role in B.regions()], contigs=[],
              classes=[dict(name=n, kind="positional", length=pc.length, panel_kmers=len(pc.kmer_pos), reads=1)
                       for n, pc in real.classes.items() if pc.kind == "positional"])
    assert contract.issues(ok, B) == []
    other = dict(ok, panel_sha256=["0" * 64])
    assert [m for _, m in contract.issues(other, B)] == ["x: none of the panels it was counted with is the bundle's panel file (panel.k31.tsv.gz)"]
    older = copy.deepcopy(ok)
    older["classes"][0]["panel_kmers"] -= 1
    assert contract.panel_mismatch(older, B).keys() == {older["classes"][0]["name"]}


def test_unreadable_counts_files_name_their_path(tmp_path, sim):
    good = tmp_path / "good.json.gz"
    with gzip.open(good, "wt") as fh:
        json.dump(sim[0], fh)
    plain = tmp_path / "plain.json.gz"                          # not compressed, whatever the name says
    plain.write_text(json.dumps(sim[0]))
    assert contract.load(plain)["sample"] == contract.load(good)["sample"] == "sim"
    cut = tmp_path / "cut.json.gz"
    cut.write_bytes(good.read_bytes()[:len(good.read_bytes()) // 2])
    notjson = tmp_path / "list.json"
    notjson.write_text("[1, 2]")
    other = tmp_path / "other.json"
    other.write_text(json.dumps({"format": "something-else"}))
    raw = bytearray(good.read_bytes())                          # damaged inside the deflate stream: zlib.error
    raw[len(raw) // 3:len(raw) // 3 + 400] = bytes(b ^ 0x5A for b in raw[len(raw) // 3:len(raw) // 3 + 400])
    damaged = tmp_path / "damaged.json.gz"
    damaged.write_bytes(bytes(raw))
    for bad in (cut, damaged, notjson, other, tmp_path / "absent.json"):
        with pytest.raises(contract.CountsError, match=bad.name) as e:
            contract.load(bad)
        if bad is damaged:
            assert isinstance(e.value.__cause__, zlib.error)
    assert issubclass(contract.CountsError, ValueError)
    with pytest.raises(ValueError, match=damaged.name):             # an estimate file damaged the same way
        load_result(damaged)


def test_region_table_cache_survives_a_truncated_file(tmp_path, monkeypatch):
    monkeypatch.setenv("NGSDOSE_CACHE", str(tmp_path))
    B = resources.Bundle(ROOT / "resources" / "GRCh38")
    first = B.region_tables(300)
    (f,) = tmp_path.glob("region_tables.*.npy")
    assert not list(tmp_path.glob("*.tmp.npy"))
    for content in (f.read_bytes()[:300000], b"", np.zeros((3, 101)).tobytes()):
        f.write_bytes(content)
        assert np.array_equal(B.region_tables(300), first)
    np.save(f, np.zeros((5, 101)))                                   # a valid file of the wrong shape
    assert np.array_equal(B.region_tables(300), first) and np.array_equal(np.load(f), first)


def _bundle_copy(tmp_path, **meta):
    src = ROOT / "resources" / "GRCh38"
    d = tmp_path / "b"
    d.mkdir(parents=True)
    m = json.loads((src / "bundle.json").read_text())
    m.update(meta)
    (d / "bundle.json").write_text(json.dumps(m))
    for key in ("panel", "controls"):
        shutil.copy(src / m[key], d / m[key])
    return resources.Bundle(d)


def test_a_bundle_missing_a_named_file_says_so(tmp_path):
    B = _bundle_copy(tmp_path)
    with pytest.raises(FileNotFoundError, match="anchors.json named in bundle.json is missing"):
        B.anchors()
    with pytest.raises(FileNotFoundError, match="the unit of .* named in bundle.json is missing"):
        B.units()
    assert _bundle_copy(tmp_path / "x", anchors=None).anchors() == {}
    B = resources.Bundle(ROOT / "resources" / "GRCh38")
    B.check_units(io.load_panel(B.panel))
    with pytest.raises(ValueError, match="rDNA5S: no unit"):
        B.check_units(io.load_panel(B.panel), {k: v for k, v in B.units().items() if k != "rDNA5S"})


def test_a_misspelt_region_role_is_refused(tmp_path):
    B = _bundle_copy(tmp_path)
    for header, msg in ((">chr1:100-200 flank=10 role=contrl\n", "role=contrl"), (">chr1:100-200 flank=10 role=test\n", "role=test")):
        with gzip.open(B.controls, "wt") as fh:
            fh.write(">chr1:1-50 flank=10 role=control\n" + "A" * 69 + "\n" + header + "A" * 120 + "\n")
        B = resources.Bundle(B.dir)
        with pytest.raises(ValueError, match=msg):
            B.regions()
