"""The candidate panels (resources/experimental/candidates) can be loaded in any combination next to
the shipped panels without changing a shipped count, and their reads are recognised.

(a) no candidate k-mer is in a shipped panel, none is in two candidate panels, and no class name is
    used twice (the engine drops a k-mer held by two loaded panels from both, and refuses a class name
    defined twice); no k-mer is a microsatellite or telomere-repeat variant (periodic with a period of
    1-6 bp, at most 2 mismatches: such k-mers are in every genome's reads although absent from the
    reference); every class has its unit where the estimator looks for it (units/<class>.fa or .fa.gz),
    no circular unit ends by repeating its first bases, and the manifest and the fetch menu list exactly
    these classes as candidates;
(b) a scan of the fixture with the shipped panels alone and with every candidate panel added gives the
    same entry, placements, controls and regions for every shipped class;
(c) each candidate panel loads next to the shipped panels on its own, without a shared-k-mer message;
(e) reads cut from each class's unit, added to the fixture as reads without a coordinate, are counted
    as that class.
With the GRCh38 analysis set and CHM13v2.0 in NGSDOSE_REF_DIR (default work/ref), the paralog-pair
classes are also checked against the reference: no OPN1LW k-mer in any MW gene, no OPN1MW k-mer in the LW
gene, no CDY1 k-mer at the CDY2 loci and none of CDY2 at the CDY1 loci (a build with one mask for all
classes let such k-mers through), and every OPN1MW k-mer in all five MW genes; and each circular unit cut
from an assembly is one period of its array (the base after it starts the next unit).
NGSDOSE_EXTRA_BINS (paths separated by ':') runs (b), (c) and (e) with other engine builds as well, such
as the one counting a cohort."""
import gzip
import json
import os
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

from ngsdose import fetchplan, io

ROOT = Path(__file__).resolve().parents[1]
CAND = ROOT / "resources" / "experimental" / "candidates"
SHIPPED = [ROOT / "resources" / "GRCh38" / "panel.k31.tsv.gz", ROOT / "resources" / "experimental" / "satellites.CHM13v2.k31.panel.tsv.gz",
           ROOT / "resources" / "experimental" / "telomere.k31.panel.tsv.gz"]
CONTROLS = ROOT / "resources" / "GRCh38" / "controls.fa.gz"
BAM = ROOT / "tests" / "data" / "NA12878.subsample.bam"
BIN = Path(os.environ.get("NGSDOSE_BIN", ROOT / "target" / "release" / "ngs-dose"))
ENGINES = [BIN] + [Path(p) for p in os.environ.get("NGSDOSE_EXTRA_BINS", "").split(os.pathsep) if p]
PANEL_FIELDS = {"panel", "panel_sha256", "classes", "placements", "elapsed_sec", "below_threshold_reads", "ambiguous_reads",
                "pipeline", "engine_build", "engine_version"}
COMP = str.maketrans("ACGT", "TGCA")


def candidate_panels():
    return sorted(CAND.glob("*.k31.panel.tsv.gz"))


@lru_cache(maxsize=None)
def read_panel(path):
    """(class definitions in id order, k-mers as bytes, class ids)"""
    classes, kmers, cls = [], [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##class\t"):
                classes.append(dict(f.split("=", 1) for f in line.rstrip("\n").split("\t")[1:]))
            elif not line.startswith("#"):
                k, c, _ = line.split("\t", 2)
                kmers.append(k)
                cls.append(int(c))
    return classes, np.array(kmers, dtype="S31"), np.array(cls, dtype=np.int32)


def read_fasta(path):
    recs, name = {}, None
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                recs[name] = []
            else:
                recs[name].append(line.strip().upper())
    return {n: "".join(s) for n, s in recs.items()}


def unit_file(name):
    """units/<name>.fa or units/<name>.fa.gz, exactly one of them (an index such as a stray .fai is not a unit)"""
    hits = [f for f in (CAND / "units" / f"{name}.fa", CAND / "units" / f"{name}.fa.gz") if f.is_file()]
    return hits[0] if len(hits) == 1 else None


def test_the_panels_share_no_kmer_with_each_other_or_the_shipped_ones():
    panels = candidate_panels()
    assert len(panels) >= 7
    shipped = np.unique(np.concatenate([read_panel(p)[1] for p in SHIPPED]))
    names = [c["name"] for p in SHIPPED for c in read_panel(p)[0]]
    seen = []
    for p in panels:
        classes, kmers, cls = read_panel(p)
        assert len(np.unique(kmers)) == len(kmers), p.name
        assert not np.isin(kmers, shipped).any(), f"{p.name} holds k-mers of a shipped panel"
        n = np.bincount(cls, minlength=len(classes))
        for i, c in enumerate(classes):
            assert int(c["kmers_kept"]) == n[i] > 0, (p.name, c["name"])
            assert int(c["id"]) == i
        names += [c["name"] for c in classes]
        seen.append(kmers)
    allk = np.concatenate(seen)
    assert len(np.unique(allk)) == len(allk), "a k-mer is in two candidate panels: the engine would drop it when both are loaded"
    assert len(names) == len(set(names)) <= 255, "class names must be unique over every panel a run may load"


def periodic(kmers, max_period=6, max_mismatch=2):
    m = np.frombuffer(kmers.tobytes(), dtype=np.uint8).reshape(len(kmers), -1)
    return np.logical_or.reduce([(m[:, :-p] != m[:, p:]).sum(axis=1) <= max_mismatch for p in range(1, max_period + 1)])


def test_no_panel_holds_a_simple_repeat_kmer():
    assert periodic(np.array([b"CCCTAACCCTAACCCTAGCCCTAACCCTAAC", b"ATGTATATGTATATGTATATGTATATGTATA"], dtype="S31")).all()
    for p in candidate_panels():
        classes, kmers, cls = read_panel(p)
        hit = periodic(kmers)
        assert not hit.any(), f"{p.name}: periodic k-mers in {sorted({classes[c]['name'] for c in cls[hit]})}"


def test_every_class_has_its_unit_and_is_listed_as_a_candidate():
    defined = {}
    for p in candidate_panels():
        for c in read_panel(p)[0]:
            defined[c["name"]] = (p.name, c)
            f = unit_file(c["name"])
            assert f is not None, f"no unit for {c['name']} in {CAND / 'units'}"
            if c["kind"] == "positional":                  # what the estimator reads (resources.ExperimentalUnits)
                seqs = read_fasta(f)
                assert len(seqs) == 1 and len(next(iter(seqs.values()))) == int(c["length"]), c["name"]
    rows = [line.rstrip("\n").split("\t") for line in open(CAND / "candidates.tsv") if not line.startswith("#")]
    rows = [dict(zip(rows[0], r)) for r in rows[1:]]
    assert {r["class"] for r in rows} == set(defined) and len(rows) == len(defined)
    for r in rows:
        panel, c = defined[r["class"]]
        assert (r["panel"], r["kind"], int(r["kmers"])) == (panel, c["kind"], int(c["kmers_kept"])), r["class"]
        assert r["status"] == "candidate: no sinks yet" and r["tier"] in ("A", "B", "C", "D")
    menu = fetchplan.read_menu()
    for name, (panel, c) in defined.items():
        o = menu.options.get(name)
        assert o is not None, f"{name} is not in resources/fetch_menu.tsv"
        assert o.status == "candidate" and not o.sinks and o.kind == c["kind"], name
        assert menu.resolve(o.panel) == CAND / panel and name in fetchplan.panel_classes(menu.resolve(o.panel))


def test_circular_units_do_not_repeat_their_first_bases_at_the_end():
    """A circular unit is one period: if its last bases repeat its first, the circular join makes k-mers
    that occur in no genome (TSPY, OPN1, CT45 and RNU2 did, by 21, 22, 16 and 6 bp)."""
    for p in candidate_panels():
        for c in read_panel(p)[0]:
            if c["kind"] != "positional" or c["circular"] != "1":
                continue
            seq = next(iter(read_fasta(unit_file(c["name"])).values()))
            over = [n for n in range(6, len(seq) // 2) if seq[-n:] == seq[:n]]
            assert not over, f"{c['name']}: its last {over[-1]} bp repeat its first"


REF = Path(os.environ.get("NGSDOSE_REF_DIR", ROOT / "work" / "ref"))
FASTA = {"GRCh38": REF / "GRCh38_full_analysis_set_plus_decoy_hla.fa", "CHM13": REF / "chm13v2.0.fa"}
# the other gene type's copies (whole gene +-500 bp; for CDY the other family's loci with their flanks), 1-based
PARALOG_COPIES = {
    "OPN1LW": [("GRCh38", "chrX:154182096-154197361"), ("GRCh38", "chrX:154219256-154233786"), ("GRCh38", "chrX:154257082-154271590"),
               ("CHM13", "chrX:152455758-152471023"), ("CHM13", "chrX:152492894-152507424")],
    "OPN1MW": [("GRCh38", "chrX:154143743-154159532"), ("CHM13", "chrX:152417396-152433185")],
    "CDY1": [("GRCh38", "chrY:17876446-17881220"), ("GRCh38", "chrY:18024787-18029561"),
             ("CHM13", "chrY:18783004-18787778"), ("CHM13", "chrY:18931262-18936036")],
    "CDY2": [("GRCh38", "chrY:24044265-24049016"), ("GRCh38", "chrY:25621115-25625866"),
             ("CHM13", "chrY:24813544-24818295"), ("CHM13", "chrY:26433508-26438259")],
}


# circular units cut from an assembly: (assembly, contig, first base). The base after the unit must start
# the next unit of the array, i.e. the unit's first 31-mer recurs exactly one unit length on
CIRCULAR_SOURCES = {
    "DXZ4": ("CHM13", "chrX", 114191115), "CT47": ("CHM13", "chrX", 119264861), "RS447": ("GRCh38", "chr4", 9215000),
    "MSR5p": ("GRCh38", "chr5", 17518000), "FLJ40296": ("GRCh38", "chr13", 57143000), "ZAV": ("GRCh38", "chr9", 113064500),
    "REXO1L1": ("GRCh38", "chr8", 85760231), "KIV2": ("GRCh38", "chr6", 160617277), "DEFA1A3": ("GRCh38", "chr8", 6976263),
    "TSPY": ("GRCh38", "chrY", 9462691), "OPN1": ("GRCh38", "chrX", 154182596), "GAGE": ("GRCh38", "chrX", 49532177),
    "CT45": ("GRCh38", "chrX", 135829229), "TDNA1Q23": ("GRCh38", "chr1", 161447001),
}


def _faidx(asm, region):
    out = subprocess.run(["samtools", "faidx", str(FASTA[asm]), region], capture_output=True, text=True, check=True).stdout
    return "".join(out.split("\n")[1:]).upper()


@pytest.mark.skipif(shutil.which("samtools") is None or not all(f.exists() for f in FASTA.values()),
                    reason="needs samtools and the reference assemblies (NGSDOSE_REF_DIR)")
def test_circular_units_are_one_period_of_their_array():
    circular = {c["name"] for p in candidate_panels() for c in read_panel(p)[0] if c["kind"] == "positional" and c["circular"] == "1"}
    assert circular - set(CIRCULAR_SOURCES) <= {"RNU2", "D4Z4", "PHIX"}      # GenBank units
    for name, (asm, contig, start) in CIRCULAR_SOURCES.items():
        seq = next(iter(read_fasta(unit_file(name)).values()))
        n = len(seq)
        assert _faidx(asm, f"{contig}:{start}-{start + n - 1}") == seq, f"{name} is not {asm} {contig}:{start}+{n}"
        assert _faidx(asm, f"{contig}:{start + n}-{start + n + 30}") == seq[:31], f"{name}: the next unit does not start after it"


def _region_kmers(asm, region, k=31):
    out = subprocess.run(["samtools", "faidx", str(FASTA[asm]), region], capture_output=True, text=True, check=True).stdout
    s = "".join(out.split("\n")[1:]).upper()
    return {_canon(s[i:i + k]) for i in range(len(s) - k + 1)}


def _class_kmers(name):
    for p in candidate_panels():
        classes, kmers, cls = read_panel(p)
        for i, c in enumerate(classes):
            if c["name"] == name:
                return {x.decode() for x in kmers[cls == i]}
    raise KeyError(name)


@pytest.mark.skipif(shutil.which("samtools") is None or not all(f.exists() for f in FASTA.values()),
                    reason="needs samtools and the reference assemblies (NGSDOSE_REF_DIR)")
def test_paralog_pair_classes_have_no_kmer_in_the_other_copies():
    for name, regions in PARALOG_COPIES.items():
        own = _class_kmers(name)
        for asm, region in regions:
            hit = own & _region_kmers(asm, region)
            assert not hit, f"{name}: {len(hit)} k-mers also in {asm} {region}"
    mw = _class_kmers("OPN1MW")                         # the MW genes themselves: every k-mer in each of them
    for asm, region in PARALOG_COPIES["OPN1LW"]:
        assert mw <= _region_kmers(asm, region), f"OPN1MW k-mers missing from {asm} {region}"


def _scan(engine, bam, tmp, tag, extra=()):
    out = tmp / f"{engine.parent.name}.{engine.name}.{tag}.json.gz"
    cmd = [str(engine), "count", "-m", "scan", "-i", str(bam), "-c", str(CONTROLS), "-@", "2", "-o", str(out)]
    for p in list(SHIPPED) + list(extra):
        cmd += ["-p", str(p)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "more than one panel" not in r.stderr and "WARNING" not in r.stderr, r.stderr
    return io.load_counts(out)


def _same_shipped(a, b):
    shipped = {c["name"] for c in a["classes"]}
    got = {c["name"]: c for c in b["classes"]}
    for c in a["classes"]:
        assert got[c["name"]] == c, c["name"]
    key = lambda p: json.dumps(p, sort_keys=True)
    assert sorted(map(key, a["placements"])) == sorted(key(p) for p in b["placements"] if p["class"] in shipped)
    for k in set(a) | set(b):
        if k not in PANEL_FIELDS:
            assert a.get(k) == b.get(k), k


engines = pytest.mark.parametrize("engine", ENGINES, ids=lambda p: p.name)
needs_engine = pytest.mark.skipif(not all(e.exists() for e in ENGINES), reason="needs the engine (cargo build --release)")


@needs_engine
@engines
def test_loading_every_candidate_panel_changes_no_shipped_count(engine, tmp_path):
    a = _scan(engine, BAM, tmp_path, "shipped")
    b = _scan(engine, BAM, tmp_path, "all", candidate_panels())
    _same_shipped(a, b)
    assert len(b["classes"]) == len(a["classes"]) + sum(len(read_panel(p)[0]) for p in candidate_panels())


@needs_engine
@engines
def test_each_candidate_panel_loads_on_its_own(engine, tmp_path):
    a = _scan(engine, BAM, tmp_path, "shipped")
    for p in candidate_panels():
        b = _scan(engine, BAM, tmp_path, p.name.split(".")[0], [p])
        assert {c["name"] for c in b["classes"]} - {c["name"] for c in a["classes"]} == {c["name"] for c in read_panel(p)[0]}
        _same_shipped(a, b)


def _canon(s):
    r = s.translate(COMP)[::-1]
    return min(s, r)


def _reads_of(name, kind, kmers, per_class=5, rl=150, k=31):
    """Reads cut from the class's unit that carry >= 20 of its k-mers (>= 4 for classes built from
    sparse diagnostic sites), alternately reverse-complemented."""
    own = np.sort(kmers)
    need = 20 if len(own) > 3000 else 4
    out = []
    for seq in read_fasta(unit_file(name)).values():
        step = max(1, (len(seq) - rl) // 400)
        for i in range(0, max(len(seq) - rl + 1, 0), step):
            r = seq[i:i + rl]
            if "N" in r:
                continue
            q = np.array([_canon(r[j:j + k]) for j in range(rl - k + 1)], dtype="S31")
            idx = np.searchsorted(own, q).clip(0, len(own) - 1)
            if (own[idx] == q).sum() >= need:
                out.append(r if len(out) % 2 == 0 else r.translate(COMP)[::-1])
                if len(out) == per_class:
                    return out
    return out


@pytest.mark.skipif(shutil.which("samtools") is None, reason="needs samtools")
@needs_engine
@engines
def test_reads_from_each_unit_are_counted_as_their_class(engine, tmp_path):
    reads = {}
    for p in candidate_panels():
        classes, kmers, cls = read_panel(p)
        for i, c in enumerate(classes):
            reads[c["name"]] = _reads_of(c["name"], c["kind"], kmers[cls == i])
            assert len(reads[c["name"]]) == 5, c["name"]
    sam = tmp_path / "sim.sam"
    with open(sam, "w") as fh:
        subprocess.run(["samtools", "view", "-h", str(BAM)], stdout=fh, check=True)
        for name, rs in reads.items():
            for j, r in enumerate(rs):
                fh.write(f"sim_{name}_{j}\t4\t*\t0\t0\t*\t*\t0\t0\t{r}\t*\n")
    bam = tmp_path / "sim.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)
    base = {c["name"]: c["reads"] for c in _scan(engine, BAM, tmp_path, "fixture", candidate_panels())["classes"]}
    sim_all = _scan(engine, bam, tmp_path, "sim.all", candidate_panels())
    got = {c["name"]: c["reads"] for c in sim_all["classes"]}
    short = {n: got[n] - base[n] for n in reads if got[n] - base[n] < len(reads[n])}
    assert not short, f"classes that did not get their reads: {short}"
    _same_shipped(_scan(engine, bam, tmp_path, "sim.shipped"), sim_all)
