"""The shipped GRCh38 bundle is internally consistent. A cohort run is only as good as these files,
and nothing else would notice a truncated panel or a controls FASTA that no longer matches its BED."""
import gzip
import struct
import subprocess
from pathlib import Path

import pytest

from ngsdose import io, resources

ROOT = Path(__file__).resolve().parents[1]
B = resources.Bundle(ROOT / "resources" / "GRCh38")
# sink intervals of GRCh38-v1 that run past the end of their contig (40,553 bp in all). Clipping them
# changes sinks_sha256, which every fetch counts file records, so it waits for a sinks.bed revision
# with a sinks_history entry in bundle.json.
OVERLONG_SINKS = {
    ("chr17_GL000205v2_random", 0, 191000, "DJ"), ("chr22_KI270733v1_random", 109000, 181000, "rDNA45S"),
    ("chrUn_GL000195v1", 0, 191000, "DJ"), ("chrUn_GL000219v1", 159000, 181000, "DJ"),
    ("chrUn_GL000220v1", 99000, 171000, "rDNA45S"), ("chrUn_GL000224v1", 119000, 181000, "DJ"),
    ("chrUn_JTFH01001783v1_decoy", 0, 11000, "DJ"), ("chrUn_JTFH01001847v1_decoy", 0, 11000, "DJ"),
}


def test_every_file_named_by_the_bundle_exists():
    for key in ("panel", "controls", "controls_bed", "sinks", "features", "anchors"):
        assert (B.dir / B.meta[key]).stat().st_size > 0, key
    for cls, rel in B.meta["units"].items():
        assert (B.dir / rel).exists(), cls


def test_panel_matches_units_and_is_well_formed():
    panel, units = io.load_panel(B.panel), B.units()
    assert panel.k == B.meta["k"] == 31
    assert set(panel.classes) == set(units) == {"rDNA45S", "rDNA5S", "DJ"}
    for name, pc in panel.classes.items():
        assert pc.kind == "positional" and pc.length == len(units[name])
        assert len(pc.kmer_pos) == len(set(pc.kmer_pos.tolist())), "one k-mer per unit position"
        assert pc.kmer_pos.min() >= 0 and pc.kmer_pos.max() < pc.length
    assert panel.classes["rDNA45S"].circular and panel.classes["rDNA5S"].circular and not panel.classes["DJ"].circular
    # every panel k-mer really is the unit's sequence at the position it claims
    comp = str.maketrans("ACGT", "TGCA")
    seen = 0
    with gzip.open(B.panel, "rt") as fh:
        ids = {}
        for line in fh:
            if line.startswith("##class\t"):
                kv = dict(f.split("=", 1) for f in line.rstrip("\n").split("\t")[1:])
                ids[kv["id"]] = kv["name"]
            elif not line.startswith("#"):
                kmer, cid, pos, strand = line.rstrip("\n").split("\t")
                u = units[ids[cid]]
                ref = (u + u[:30])[int(pos):int(pos) + 31]
                assert kmer == (ref if strand == "+" else ref.translate(comp)[::-1])
                seen += 1
                if seen >= 5000:
                    break
    assert seen == 5000


def test_controls_fasta_matches_its_bed():
    bed = [l.rstrip("\n").split("\t") for l in open(B.dir / B.meta["controls_bed"])]
    names, lengths = [], []
    with gzip.open(B.controls, "rt") as fh:
        n = 0
        for line in fh:
            if line.startswith(">"):
                if names:
                    lengths.append(n)
                names.append(line[1:].split())
                n = 0
            else:
                n += len(line.strip())
        lengths.append(n)
    assert len(names) == len(bed) == 982
    roles = {}
    for (c, s, e, role), (name, *tags), ln in zip(bed, names, lengths):
        assert name == f"{c}:{s}-{e}"
        t = dict(x.split("=") for x in tags)
        assert ln == int(e) - int(s) + 2 * int(t["flank"]) and int(t["flank"]) >= 600
        assert (t["role"] + (":" + t["label"] if "label" in t else "")) == role
        roles[role] = roles.get(role, 0) + 1
    assert roles == {"control": 800, "test:auto": 80, "test:chrX": 60, "test:chrY": 40, "dosage:chrM": 1, "dosage:chrEBV": 1}
    # controls are on the autosomes only, and spread over all of them
    assert {c for c, _, _, r in bed if r == "control"} == {f"chr{i}" for i in range(1, 23)}
    # every contig a region sits on has its GRCh38 length on record: what `estimate` tells builds apart by
    assert {c for c, *_ in bed} - {"chrEBV"} <= set(B.contig_lengths())


def test_every_region_position_has_a_fragment_gc_window():
    """No N within reach of any region: at the longest window the engine tabulates, every region
    holds exactly two position-strands per base. (GRCh38's chrM has an N at 3,107; a region over
    it reads 18% high unless the estimator compensates - it does, but the bundle should not need it.)"""
    from ngsdose import estimate
    names = [l[1:].split()[0] for l in gzip.open(B.controls, "rt") if l.startswith(">")]
    lens = [int(n.rsplit(":", 1)[1].split("-")[1]) - int(n.rsplit(":", 1)[1].split("-")[0]) for n in names]
    tables = estimate.control_region_tables(B.controls, 600)
    assert [int(x) for x in tables.sum(1)] == [2 * ln for ln in lens]


def test_features_anchors_and_sinks_are_inside_their_coordinate_systems():
    units = B.units()
    for cls, feats in B.features().items():
        for name, s, e in feats:
            assert 0 <= s < e <= len(units[cls]), (cls, name)
    anchors = B.anchors()
    assert set(anchors) == {"rDNA45S"}
    span = sum(e - s for s, e in anchors["rDNA45S"])
    assert 2000 <= span <= 10000 and all(0 <= s < e <= len(units["rDNA45S"]) for s, e in anchors["rDNA45S"])
    # the anchors include the gene itself, not only spacer
    assert any(3657 <= s and e <= 5526 for s, e in anchors["rDNA45S"])
    classes = set()
    lengths = B.meta["contig_lengths"]
    for line in open(B.sinks):
        c, s, e, cls = line.rstrip("\n").split("\t")
        assert int(s) < int(e) and (c not in lengths or int(e) <= lengths[c]), line
        classes.add(cls)
    assert classes == {"rDNA45S", "rDNA5S", "DJ", "TEL"}                 # the telomeric repeat is fetchable; the satellites are not
    assert B.meta["expected_copies"] == {"DJ": 10}


def bam_contig_lengths(path):
    """The @SQ names and lengths of a BAM header (BGZF reads as gzip)."""
    with gzip.open(path, "rb") as fh:
        assert fh.read(4) == b"BAM\1"
        l_text, = struct.unpack("<i", fh.read(4))
        fh.read(l_text)
        n_ref, = struct.unpack("<i", fh.read(4))
        out = {}
        for _ in range(n_ref):
            l_name, = struct.unpack("<i", fh.read(4))
            name = fh.read(l_name)[:-1].decode()
            out[name], = struct.unpack("<i", fh.read(4))
    return out


def sink_rows():
    rows = []
    for line in open(B.sinks):
        c, s, e, cls = line.rstrip("\n").split("\t")
        rows.append((c, int(s), int(e), cls))
    return rows


def contig_lengths():
    """GRCh38_full_analysis_set_plus_decoy_hla contig lengths: the fixture's header (an NYGC CRAM's, cut
    to the contigs the bundle's regions and sinks sit on), which agrees with bundle.json's primary contigs."""
    lengths = bam_contig_lengths(ROOT / "tests" / "data" / "NA12878.subsample.bam")
    for c, n in B.meta["contig_lengths"].items():
        assert lengths.get(c, n) == n, c
    return {**lengths, **B.meta["contig_lengths"]}


def test_sinks_lie_on_contigs_of_known_length():
    lengths = contig_lengths()
    rows = sink_rows()
    assert {c for c, *_ in rows} <= set(lengths), "a sink contig the fixture's header does not name: rebuild the fixture (tests/data/make_fixture.sh)"
    for c, s, e, cls in rows:
        assert 0 <= s < e and s < lengths[c], (c, s, e, cls)
    # no interval runs past its contig's end beyond the eight known ones, and the telomeric rows (clipped when learned) never do
    over = {r for r in rows if r[2] > lengths[r[0]]}
    assert over <= OVERLONG_SINKS, over - OVERLONG_SINKS


@pytest.mark.xfail(strict=True, reason="8 non-TEL sink intervals run past their contig end (OVERLONG_SINKS); clipping them changes sinks_sha256 and waits for a sinks.bed revision")
def test_sinks_end_within_their_contigs():
    lengths = contig_lengths()
    assert [r for r in sink_rows() if r[2] > lengths[r[0]]] == []


def test_shell_scripts_parse():
    for sh in list((ROOT / "resources" / "build").glob("*.sh")) + list((ROOT / "tests").rglob("*.sh")):
        subprocess.run(["bash", "-n", str(sh)], check=True)
