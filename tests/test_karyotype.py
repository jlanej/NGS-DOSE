"""Chromosomes in copies (karyotype.py), on cohorts whose truth is known: regions with efficiencies of their own, two
library modes that move whole chromosomes, counting noise, and genomes that gain or lose a chromosome, an arm or a
stretch in every cell or in part of them."""
import gzip
import json
import subprocess
import sys

import numpy as np
import pytest

from ngsdose import cohort, karyotype as K, resources

BUNDLE = resources.Bundle()
ARMS = BUNDLE.karyotype()["arms"]
GC = BUNDLE.karyotype()["gc"]
LENGTHS = BUNDLE.contig_lengths()
ACRO = ("chr13", "chr14", "chr15", "chr21", "chr22")
U2 = 0.95                    # what a second X reads of the first in the simulated cohorts


def layout(per_arm=10, pieces=1):
    """Regions of 12 kb, `per_arm` places evenly along each arm of every autosome and of chrX and along chrY's long arm;
    with `pieces` > 1 every place is that many regions 3 kb apart (a window)."""
    names = []
    for c in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]:
        b, n = ARMS[c], LENGTHS[c]
        arms = [(b + 3_000_000, min(n, 27_000_000 if c == "chrY" else n) - 3_000_000)]
        if c not in ACRO and c != "chrY":
            arms.insert(0, (3_000_000, b - 3_000_000))
        for lo, hi in arms:
            for s in np.linspace(lo, hi, per_arm).astype(int):
                names += [f"{c}:{s + k * 15_000}-{s + k * 15_000 + 12_000}" for k in range(pieces)]
    return names


def cohort_of(rng, names, n=160, noise=0.03, men=None):
    """log(observed / expected) of `n` genomes: a level each, the regions' efficiencies, two modes whose loadings have a
    part shared by the chromosome (so that a mode moves whole chromosomes) and a part of the region's own, and noise.
    Every second genome is a man unless `men` says otherwise."""
    tab = K.table(names, ARMS)
    m = len(names)
    a = rng.normal(0, 0.05, m)
    by_chrom = {c: rng.normal(0, 0.012, 2) for c in sorted(set(tab.chrom))}
    V = np.array([by_chrom[c] for c in tab.chrom]) + rng.normal(0, 0.012, (m, 2))
    T = rng.normal(0, 1, (n, 2))
    level, counting = rng.normal(0, 0.1, (n, 1)), rng.normal(0, noise, (n, m))
    men = np.arange(n) % 2 == 0 if men is None else np.asarray(men, bool)
    counting[np.ix_(~men, tab.kind == "X")] /= np.sqrt(1 + U2)          # counting noise: twice the reads, less of it
    Y = level + a[None, :] + T @ V.T + counting
    for i in range(n):
        Y[i] = sex(Y[i], tab, 1 if men[i] else 2, 1 if men[i] else 0)
    return tab, Y, men


def sex(y, tab, nx, ny):
    """Give a genome `nx` X and `ny` Y chromosomes: a second X reads U2 of the first."""
    y = y.copy()
    y[tab.kind == "X"] += np.log(max(1 + (nx - 1) * U2, 1e-3) / 2) if nx else np.log(1e-3)
    y[tab.kind == "Y"] += np.log(ny / 2) if ny else np.log(1.5e-3)
    return y


def change(y, where, share, sign=+1):
    """A copy gained (or lost) in `share` of the cells over the regions `where` of an autosome."""
    y = y.copy()
    y[where] += np.log(1 + sign * share / 2)
    return y


@pytest.fixture(scope="module")
def world():
    rng = np.random.default_rng(5)
    names = layout()
    tab, Y, men = cohort_of(rng, names)
    model = K.fit(names, Y, ARMS, gc=GC)
    readings = [K.read(model, names, y) for y in Y]
    model.phi, model.floor, _ = K.calibrate(model, readings)
    return dict(rng=rng, names=names, tab=tab, Y=Y, men=men, model=model)


def test_a_cohort_of_plain_genomes_reads_plain(world):
    model, names, Y, men = world["model"], world["names"], world["Y"], world["men"]
    rd = [K.read(model, names, y) for y in Y]
    assert [r.karyotype() for r in rd] == ["46,XY" if k else "46,XX" for k in men]
    assert all(r.status == "settled" and not r.events for r in rd)
    # the second X's factor is learned, and a chromosome's copies scatter as its regions' noise says
    assert abs(np.nanmedian(model.u) - U2) < 0.01 and model.k >= 2
    c8 = np.array([r.chromosomes["chr8"].copies for r in rd])
    assert abs(np.median(c8) - 2) < 0.003 and 0.006 < c8.std() < 0.02
    z = np.array([r.chromosomes[c].z for r in rd for c in K.AUTOSOMES])
    assert 0.8 < np.std(z) < 1.25 and np.abs(z).max() < 5
    assert max(abs(r.x - (1 if k else 2)) for r, k in zip(rd, men)) < 0.08 and max(abs(r.y - (1 if k else 0)) for r, k in zip(rd, men)) < 0.05


def test_a_trisomy_and_a_monosomy_are_whole_numbers(world):
    model, names, tab, Y = world["model"], world["names"], world["tab"], world["Y"]
    r = K.read(model, names, change(Y[1], tab.chrom == "chr21", 1.0))
    assert r.karyotype() == "47,XX,+21" and r.status == "settled"
    (e,) = r.events
    assert e.whole and e.span == "whole" and abs(e.copies - 3) < 0.06 and e.z > 20
    # the other chromosomes are not moved by it, although chr21's regions took part in nothing
    assert all(abs(r.chromosomes[c].copies - 2) < 0.05 for c in K.AUTOSOMES if c != "chr21")
    r = K.read(model, names, change(Y[0], tab.chrom == "chr18", 1.0, sign=-1))
    assert r.karyotype() == "45,XY,-18" and abs(r.chromosomes["chr18"].copies - 1) < 0.04
    # the largest chromosome in every cell, and another in a tenth of them: neither biases the other
    y = change(change(Y[1], tab.chrom == "chr1", 1.0), tab.chrom == "chr9", 0.2)
    r = K.read(model, names, y)
    assert r.karyotype() == "47,XX,+1,+9[0.20]" or r.karyotype().startswith("47,XX,+1,+9[0.")
    assert abs(r.chromosomes["chr1"].copies - 3) < 0.05 and abs(r.chromosomes["chr9"].copies - 2.2) < 0.05


def test_a_change_in_part_of_the_cells_keeps_its_level(world):
    model, names, tab, Y = world["model"], world["names"], world["tab"], world["Y"]
    for share in (0.15, 0.3, 0.6):
        r = K.read(model, names, change(Y[0], tab.chrom == "chr8", share))
        ch = r.chromosomes["chr8"]
        assert r.status == "fractional" and ch.status == "fractional" and abs(ch.copies - (2 + share)) < 0.04
        assert ch.whole == (3 if share > 0.5 else 2) and abs(ch.off - (ch.copies - ch.whole)) < 1e-9
        (e,) = r.events
        assert not e.whole and e.label() == f"+8[{e.delta:.2f}]" and r.karyotype() == f"46,XY,{e.label()}"
    # too little to tell from noise: the level is still reported, the status settled
    r = K.read(model, names, change(Y[0], tab.chrom == "chr8", 0.01))
    assert r.status == "settled" and not r.events and abs(r.chromosomes["chr8"].copies - 2.01) < 0.05


def test_a_level_between_two_grid_steps_is_one_level():
    """The chain's levels lie on a grid of 0.05 copies. Where a chromosome is known to a few thousandths of a copy, a level
    between two steps fits each half of the chromosome better on its nearer step, and the chain cut it in two (in 8 of 60
    genomes here, at every share tried). The cut is judged again on the measured levels, at the chain's own price: a
    gain of the whole chromosome in part of the cells is one level, read at its share; a gain of one arm is still an arm."""
    rng = np.random.default_rng(13)
    names = layout(per_arm=60)
    tab, Y, men = cohort_of(rng, names, n=120, noise=0.015)
    model = K.fit(names, Y, ARMS, gc=GC)
    for share in (0.07, 0.12, 0.17):
        for i in range(60):
            r = K.read(model, names, change(Y[i], tab.chrom == "chr1", share))
            (e,) = [e for e in r.events if e.chrom == "chr1"]
            assert e.span == "whole" and len(r.chromosomes["chr1"].segments) == 1 and abs(e.delta - share) < 4 * e.se, (share, i, r.karyotype())
    q = (tab.chrom == "chr1") & (tab.arm == "q")
    for i in range(20):
        (e,) = [e for e in K.read(model, names, change(Y[i], q, 0.12)).events if e.chrom == "chr1"]
        assert e.span == "q" and abs(e.delta - 0.12) < 4 * e.se


def test_an_arm_and_a_stretch_are_told_from_the_chromosome(world):
    model, names, tab, Y = world["model"], world["names"], world["tab"], world["Y"]
    q = (tab.chrom == "chr5") & (tab.arm == "q")
    r = K.read(model, names, change(Y[0], q, 1.0))
    (e,) = r.events
    assert e.span == "q" and e.whole and r.karyotype() == "46,XY,+5q" and r.chromosomes["chr5"].segments[0].span == "p"
    assert abs(r.chromosomes["chr5"].segments[0].copies - 2) < 0.05 and len(r.chromosomes["chr5"].segments) == 2
    # the last six places of chr2: a stretch that reaches the end of the long arm, in 40% of the cells
    on = np.flatnonzero(tab.chrom == "chr2")[-6:]
    r = K.read(model, names, change(Y[1], on, 0.4))
    (e,) = r.events
    assert e.span == "qter" and not e.whole and e.regions == 6 and abs(e.copies - 2.4) < 0.06
    assert abs(e.start - (tab.end[on[0] - 1] + tab.start[on[0]]) // 2) < 1 and e.end == tab.end[on[-1]]
    assert r.karyotype().startswith("46,XX,+2(") and r.karyotype().endswith("Mb)[%.2f]" % e.delta) and r.chromosomes["chr2"].whole == 2
    # an isochromosome of the long arm of X in a woman: one short arm, three long arms, two chromosomes
    y = Y[1].copy()
    xp, xq = (tab.chrom == "chrX") & (tab.arm == "p"), (tab.chrom == "chrX") & (tab.arm == "q")
    y[xp] += np.log(1 / (1 + U2))
    y[xq] += np.log((1 + 2 * U2) / (1 + U2))
    r = K.read(model, names, y)
    assert r.complement == "XX" and sorted(e.label() for e in r.events) == ["+Xq", "-Xp"] and r.karyotype() == "46,XX,-Xp,+Xq"


def test_one_region_or_one_window_is_not_a_stretch():
    rng = np.random.default_rng(8)
    names = layout(per_arm=8, pieces=3)
    tab, Y, men = cohort_of(rng, names, n=120)
    model = K.fit(names, Y, ARMS, gc=GC)
    # a deletion of one copy over a whole window (three pieces), and a region with no copy left: neither is an event,
    # and the chromosome's level stays where it was
    y = Y[0].copy()
    win = np.flatnonzero(tab.chrom == "chr3")[6:9]
    y[win] += np.log(0.5)
    y[np.flatnonzero(tab.chrom == "chr6")[4]] += np.log(0.02)
    before, r = K.read(model, names, Y[0]), K.read(model, names, y)
    assert not r.events and r.karyotype() == before.karyotype() == "46,XY"
    assert abs(r.chromosomes["chr3"].copies - before.chromosomes["chr3"].copies) < 0.02
    # the places are the windows: three pieces each
    first, xl, sl, nreg = K.loci(tab.start[tab.chrom == "chr3"], tab.end[tab.chrom == "chr3"], np.ones(48), np.full(48, 0.03))
    assert len(first) == 16 and set(nreg) == {3} and abs(sl[0] - 0.03 / np.sqrt(3)) < 1e-9
    # five windows in a row are a stretch
    r = K.read(model, names, change(Y[0], np.flatnonzero(tab.chrom == "chr3")[:15], 1.0, sign=-1))
    (e,) = r.events
    assert e.span == "pter" and e.whole and e.regions == 15 and abs(e.copies - 1) < 0.05
    # three windows are too few to be a stretch and too many to pass for noise: the chromosome is not settled, and the note says why
    r = K.read(model, names, change(Y[0], np.flatnonzero(tab.chrom == "chr3")[:9], 1.0, sign=-1))
    ch = r.chromosomes["chr3"]
    assert not r.events and (ch.far, ch.places) == (3, 16) and ch.status == "uncertain" and r.status == "uncertain"
    assert "chr3: 3 of its 16 places" in r.note and K.columns(r, ["chr3"])["karyotype.note"] == r.note
    # what the pieces of a place share is not averaged away with them
    _, _, sl, _ = K.loci(tab.start[tab.chrom == "chr3"], tab.end[tab.chrom == "chr3"], np.ones(48), np.full(48, 0.03), local=0.02)
    assert abs(sl[0] - np.sqrt(0.03 ** 2 / 3 + 0.02 ** 2 * (1 - 1 / 3))) < 1e-9


def test_the_sex_chromosomes_are_counted_on_their_own_scales(world):
    model, names, tab, Y = world["model"], world["names"], world["tab"], world["Y"]
    base = Y[0] - sex(np.zeros(len(names)), tab, 1, 1)                     # a man's regions without his sex chromosomes
    want = {(1, 0): "45,X", (2, 1): "47,XXY", (1, 2): "47,XYY", (3, 0): "47,XXX", (2, 2): "48,XXYY", (2, 0): "46,XX", (1, 1): "46,XY"}
    for (nx, ny), kary in want.items():
        r = K.read(model, names, sex(base, tab, nx, ny))
        assert r.karyotype() == kary and r.status == "settled" and abs(r.x - nx) < 0.06 and abs(r.y - ny) < 0.06, (nx, ny, r.karyotype())
    # an X lost in a fifth of a woman's cells, a Y in a third of a man's: the nearest whole numbers and the shares
    y = base.copy()
    y[tab.kind == "X"] += np.log((1 + 0.8 * U2) / 2)
    y[tab.kind == "Y"] += np.log(1.5e-3)
    r = K.read(model, names, y)
    assert r.complement == "XX" and r.status == "fractional" and abs(r.x - 1.8) < 0.04 and r.karyotype() == f"46,XX,-X[{2 - r.x:.2f}]"
    y = base.copy()
    y[tab.kind == "X"] += np.log(1 / 2)
    y[tab.kind == "Y"] += np.log(0.67 / 2)
    r = K.read(model, names, y)
    assert r.complement == "XY" and abs(r.y - 0.67) < 0.03 and r.karyotype() == f"46,XY,-Y[{1 - r.y:.2f}]" and not r.note
    # DNA of a man in a woman's sample: X below two and Y above none by the same amount
    y = base.copy()
    y[tab.kind == "X"] += np.log((0.2 * 1 + 0.8 * (1 + U2)) / 2)
    y[tab.kind == "Y"] += np.log(0.2 / 2)
    r = K.read(model, names, y)
    assert r.complement == "XX" and "two" in r.note and abs(r.y - 0.2) < 0.03


def test_chromosomes_that_follow_their_gc_together_are_a_library_not_a_loss():
    """In some libraries the GC-rich chromosomes (19, 22, 17, 16) read low or high together, by the chromosome's GC and
    with no pattern inside a chromosome. Each chromosome is set against what the others of the genome say: such a
    library is read plain, its slope reported; a loss of chromosome 19 in part of its cells still stands out from the
    others' line; and without the chromosomes' GC the same library reads as a loss."""
    rng = np.random.default_rng(23)
    names = layout(pieces=4)
    tab, Y, men = cohort_of(rng, names, n=120)
    centre = np.mean([GC[c] for c in K.AUTOSOMES])
    gc = np.array([GC[c] - centre for c in tab.chrom])
    slopes = rng.normal(0, 0.07, len(Y))                       # copies per unit of GC, as libraries differ
    slopes[:3] = (-0.7, 0.6, -0.5)
    Y = Y + (slopes[:, None] / 2) * gc[None, :]
    rd, model, info = K.cohort([(names, y) for y in Y], ARMS, gc=GC)
    assert all(not r.events and r.status == "settled" for r in rd[:3]), [r.karyotype() for r in rd[:3]]
    got = np.array([r.gc_tilt / 0.1 for r in rd])
    assert np.abs(got[:3] - slopes[:3]).max() < 0.15 and np.corrcoef(got, slopes)[0, 1] > 0.8
    assert all(abs(r.chromosomes[c].copies - 2) < 0.03 for r in rd[:3] for c in ("chr19", "chr22", "chr4"))
    assert sum(1 for r in rd if r.events) <= 1                 # the cohort's ordinary libraries read plain too
    # chromosome 19 lost in 8% of the cells of the library tilted the same way: the others' line does not explain it
    r = K.read(model, names, change(Y[2], tab.chrom == "chr19", 0.08, sign=-1))
    (e,) = r.events
    assert e.chrom == "chr19" and e.span == "whole" and not e.whole and abs(e.delta + 0.08) < 0.025
    # a chromosome's error holds the fit's share: chromosome 19 lies at the end of the line
    assert rd[5].chromosomes["chr19"].se0 > 1.1 * rd[5].chromosomes["chr18"].se0
    # without the chromosomes' GC the tilted libraries read as losses and gains
    plain, _, _ = K.cohort([(names, y) for y in Y], ARMS)
    assert plain[0].tilt is None and any(e.chrom == "chr19" for e in plain[0].events) and any(e.chrom == "chr19" for e in plain[1].events)


def test_a_gain_that_many_of_the_cohort_carry_is_not_learned_away():
    """Chromosome 12 gained in part of the cells of a fifth of the cohort (as in lymphoblastoid lines): it is no component
    of the model, its regions keep their efficiencies, and every carrier is read at its share."""
    rng = np.random.default_rng(11)
    names = layout()
    tab, Y, men = cohort_of(rng, names, n=200)
    carriers = np.arange(0, 200, 5)
    shares = rng.uniform(0.15, 0.9, len(carriers))
    for i, f in zip(carriers, shares):
        Y[i] = change(Y[i], tab.chrom == "chr12", f)
    said = []
    rd, model, info = K.cohort([(names, y) for y in Y], ARMS, log=said.append, gc=GC)
    got = np.array([rd[i].chromosomes["chr12"].copies - 2 for i in carriers])
    assert np.abs(got - shares).max() < 0.07 and abs(np.mean(got - shares)) < 0.015
    assert all(any(e.chrom == "chr12" and e.span == "whole" for e in rd[i].events) for i in carriers)
    others = [i for i in range(200) if i not in set(carriers)]
    assert all(not rd[i].events for i in others) and abs(np.median([rd[i].chromosomes["chr12"].copies for i in others]) - 2) < 0.004
    assert info["model"] == "this cohort's" and info["with_event"] == len(carriers) and any("round 1" in s for s in said)


def test_a_saved_model_reads_a_genome_alone(world, tmp_path):
    model, names, tab, Y = world["model"], world["names"], world["tab"], world["Y"]
    again = K.Model.from_json(json.loads(json.dumps(model.to_json())))
    y = change(Y[3], tab.chrom == "chr13", 1.0)
    a, b = K.read(model, names, y), K.read(again, names, y)
    assert a.karyotype() == b.karyotype() == "47,XX,+13" and abs(a.chromosomes["chr13"].copies - b.chromosomes["chr13"].copies) < 0.002
    # one genome, the saved model: read; without a model: not read, and it is said why
    rd, used, info = K.cohort([(names, y)], ARMS, model=again)
    assert rd[0].karyotype() == "47,XX,+13" and info["model"] == "saved" and used is again
    said = []
    rd, used, info = K.cohort([(names, y)], ARMS, log=said.append)
    assert rd == [None] and used is None and "no saved model" in info["why"] and "not read" in said[0]
    # a genome with fewer regions (a lighter fetch): the same reading, a little less sure; regions the model does not know are left out
    keep = np.arange(len(names)) % 2 == 0
    r = K.read(again, [n for n, k in zip(names, keep) if k] + ["chr1:5-9"], np.r_[y[keep], 0.3])
    assert r.karyotype() == "47,XX,+13" and r.regions == keep.sum() and r.chromosomes["chr13"].se > a.chromosomes["chr13"].se
    assert K.read(again, names[:10], y[:10]) is None                       # too few regions to read anything
    with pytest.raises(ValueError, match="not a karyotype model"):
        K.Model.from_json({"format": "something else"})


def test_gather_takes_new_and_old_estimates():
    new = dict(single_copy=dict(names=["chr1:100-200", "chrX:5-9"], log_ratio=[0.01, None]))
    names, lr = K.gather(new)
    assert names == ["chr1:100-200", "chrX:5-9"] and lr[0] == 0.01 and np.isnan(lr[1])
    old = dict(control_qc=dict(region_log_ratio=[0.0, 0.1]), truth_regions=dict(chrX=dict(cn=1.94), chrY=dict(cn=0.003), auto=dict(cn=2.0)))
    names, lr = K.gather(old, ["chr1:100-200", "chr2:100-200"])
    assert names == ["chr1:100-200", "chr2:100-200", "chrX:set", "chrY:set"] and abs(lr[2] - np.log(0.97)) < 1e-9
    assert K.gather(old) is None and K.gather(old, ["chr1:100-200"]) is None and K.gather({}) is None
    t = K.table(names, ARMS)
    assert list(t.kind) == ["A", "A", "X", "Y"] and t.start[2] == -1 and t.arm[2] == "" and t.arm[0] == "p"
    # genomes with different regions fill one matrix
    nm, Y = K.assemble([(["a:1-2", "b:1-2"], np.array([1.0, 2.0])), None, (["b:1-2", "c:1-2"], np.array([3.0, 4.0]))])
    assert nm == ["a:1-2", "b:1-2", "c:1-2"] and np.isnan(Y[1]).all() and Y[2, 1] == 3 and np.isnan(Y[2, 0])


def test_estimates_before_regions_were_kept_singly_give_the_chromosomes_levels():
    """chrX and chrY as one pooled value each (estimates written before 0.3.0): counted, on the same scales."""
    rng = np.random.default_rng(13)
    names = [n for n in layout() if not n.startswith(("chrX", "chrY"))]
    tab, Y, men = cohort_of(rng, names, n=120)
    lev = np.median(Y, axis=1)
    X = lev + np.where(men, np.log(1 / 2), np.log((1 + U2) / 2)) + rng.normal(0, 0.004, 120)
    Yc = lev + np.where(men, np.log(1 / 2), np.log(1.5e-3)) + rng.normal(0, 0.006, 120)
    X[3] = lev[3] + np.log((1 + 0.7 * U2) / 2)                               # a woman who lost an X in 30% of her cells
    X[4], Yc[4] = lev[4] + np.log((1 + U2) / 2), lev[4] + np.log(1 / 2)     # a man with two X
    vec = [(names + [K.aggregate_name("chrX"), K.aggregate_name("chrY")], np.r_[Y[i], X[i], Yc[i]]) for i in range(120)]
    rd, model, info = K.cohort(vec, ARMS, gc=GC)
    assert rd[4].karyotype() == "47,XXY" and rd[3].karyotype() == f"46,XX,-X[{2 - rd[3].x:.2f}]" and abs(rd[3].x - 1.7) < 0.03
    assert [r.complement for r in rd[5:15]] == ["XX" if i % 2 else "XY" for i in range(5, 15)]
    assert len(rd[0].chromosomes["chrX"].segments) == 1 and rd[0].chromosomes["chrX"].regions == 1


def test_cohort_table_writes_the_karyotype_columns(world):
    names, tab, Y = world["names"], world["tab"], world["Y"]
    Y = Y.copy()
    Y[2] = change(Y[2], tab.chrom == "chr21", 1.0)
    results = [dict(sample=f"s{i}", mode="scan", engine_version="x", depth_equiv=30.0, read_length=150, insert_median=400, gc_L=400, ctrl_dup_frac=0.1,
                    gc_rel={}, gc_curve_max_se=0.01, control_qc=None, truth_regions={}, classes={},
                    single_copy=dict(names=names, log_ratio=[round(float(v), 4) for v in Y[i]])) for i in range(len(Y))]
    out: dict = {}
    rows, _, info = cohort.cohort_table(results, {}, n_control_pcs=0, karyotype=dict(arms=ARMS, gc=GC), karyotypes=out)
    assert rows[2]["karyotype"] == "47,XY,+21" and rows[2]["karyotype.status"] == "settled" and rows[2]["sex_chromosomes"] == "XY"
    label, copies, z = rows[2]["karyotype.events"].split("|")
    assert label == "+21" and abs(float(copies) - 3) < 0.06 and z.startswith("z+") and abs(rows[2]["chr21.copies"] - 3) < 0.06 and rows[2]["chr21.z"] is not None
    assert rows[1]["karyotype"] == "46,XX" and rows[1]["karyotype.events"] == "none" and abs(rows[1]["chrX.copies"] - 2) < 0.05
    assert info["karyotype"]["read"] == len(Y) and info["karyotype"]["with_event"] == 1 and out["model"].n == len(Y)
    long = K.long_table(out["samples"], out["readings"])
    assert sum(1 for r in long if r["sample"] == "s2" and r["kind"] == "chromosome") == 24 and {r["kind"] for r in long} <= {"chromosome", "stretch"}
    assert next(r for r in long if r["sample"] == "s2" and r["chrom"] == "chr21")["whole"] == 3
    # no bundle entry, or --no-karyotype: no such columns
    rows, _, info = cohort.cohort_table(results[:5], {}, n_control_pcs=0)
    assert "karyotype" not in rows[0] and "karyotype" not in info


def test_the_cohort_command_reads_chromosomes(world, tmp_path):
    names, tab, Y = world["names"], world["tab"], world["Y"]
    Y = Y.copy()
    Y[1] = change(Y[1], (tab.chrom == "chr5") & (tab.arm == "q"), 0.5)
    files = []
    for i in range(60):
        r = dict(sample=f"s{i}", mode="scan", engine_version="x", depth_equiv=30.0, read_length=150, insert_median=400, gc_L=400, ctrl_dup_frac=0.1,
                 gc_rel={}, gc_curve_max_se=0.01, control_qc=None, truth_regions={}, classes={},
                 single_copy=dict(names=names, log_ratio=[round(float(v), 4) for v in Y[i]]))
        f = tmp_path / f"s{i}.estimate.json.gz"
        with gzip.open(f, "wt") as fh:
            json.dump(r, fh)
        files.append(str(f))
    run = lambda *args: subprocess.run([sys.executable, "-m", "ngsdose", "cohort", *args], capture_output=True, text=True)
    table, model, long = tmp_path / "t.tsv", tmp_path / "m.json.gz", tmp_path / "k.tsv"
    p = run(*files, "-t", str(table), "--control-pcs", "0", "--save-karyotype-model", str(model), "--karyotype-table", str(long))
    assert p.returncode == 0 and "[karyotype] 60 genomes read against this cohort's model" in p.stderr, p.stderr
    rows = [dict(zip(table.read_text().splitlines()[0].split("\t"), line.split("\t"))) for line in table.read_text().splitlines()[1:]]
    assert rows[1]["karyotype"].startswith("46,XX,+5q[0.") and rows[1]["karyotype.status"] == "fractional" and rows[0]["karyotype"] == "46,XY"
    lines = [line.split("\t") for line in long.read_text().splitlines()]
    assert lines[0] == list(K.LONG_COLUMNS) and sum(1 for x in lines if x[0] == "s1" and x[2] == "stretch") == 2
    # the saved model reads three genomes alone, as it read them in the cohort
    p = run(*files[:3], "-t", str(tmp_path / "t3.tsv"), "--control-pcs", "0", "--karyotype-model", str(model))
    assert p.returncode == 0 and "saved model" in p.stderr, p.stderr
    three = [line.split("\t") for line in (tmp_path / "t3.tsv").read_text().splitlines()]
    col = three[0].index("karyotype")
    assert [x[col] for x in three[1:]] == [rows[i]["karyotype"] for i in range(3)]
    # three genomes and no model: said, and no karyotype columns
    p = run(*files[:3], "-t", str(tmp_path / "t0.tsv"), "--control-pcs", "0")
    assert p.returncode == 0 and "chromosomes are not read" in p.stderr and "karyotype" not in (tmp_path / "t0.tsv").read_text().splitlines()[0].split("\t")
    p = run(*files[:3], "-t", str(tmp_path / "t1.tsv"), "--control-pcs", "0", "--no-karyotype")
    assert p.returncode == 0 and "karyotype" not in p.stderr


def test_the_bundle_names_the_arms_of_every_chromosome():
    k = BUNDLE.karyotype()
    assert set(k["arms"]) == {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    assert all(0 < v < LENGTHS[c] for c, v in k["arms"].items())
    assert len(k["control_names"]) == sum(1 for _, role in BUNDLE.regions() if role == "control")
    t = K.table([n for n, role in BUNDLE.regions() if role != "dosage"], k["arms"])
    assert set(t.kind) == {"A", "X", "Y"} and set(t.arm) == {"p", "q"} and (t.arm[t.chrom == "chr21"] == "q").all()


def test_a_small_cohort_measures_its_own_second_x(world):
    """A saved model from DNA whose second X reads 0.95 of the first, and fifteen women whose second X reads the same
    as the first: against the model as it is they would all read above two; the factor is measured on them instead."""
    model, names, tab = world["model"], world["names"], world["tab"]
    rng = np.random.default_rng(31)
    base = world["Y"][1::2][:15] - sex(np.zeros(len(names)), tab, 2, 0)[None, :]        # fifteen women's regions without their sex chromosomes
    women = np.array([b + rng.normal(0, 0.002, len(names)) for b in base])
    women[:, tab.kind == "X"] += np.log(2 / 2)                                         # two X, the second read in full
    women[:, tab.kind == "Y"] += np.log(1.5e-3)
    vec = [(names, y) for y in women]
    said = []
    rd, used, info = K.cohort(vec, ARMS, model=model, log=said.append)
    assert abs(used.info["second_x"]["cohort"] - 1.0) < 0.02 and abs(used.info["second_x"]["model"] - U2) < 0.01 and any("scaled to it" in x for x in said)
    assert [r.karyotype() for r in rd] == ["46,XX"] * 15 and max(abs(r.x - 2) for r in rd) < 0.05
    # the model as it is (--karyotype-model: a genome's reading must not depend on the others): every one of them reads high
    rd, used, _ = K.cohort(vec, ARMS, model=model, fit_own=False)
    assert used is model and min(r.x for r in rd) > 2.03
    # too few to measure it: the model stands, and it is said on whose scale a second X is read
    said.clear()
    rd, used, _ = K.cohort(vec[:4], ARMS, model=model, log=said.append)
    assert used is model and any("saved model's scale" in x for x in said)
