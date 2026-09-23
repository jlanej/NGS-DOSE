"""The comparison of satellite estimates with HPRC assemblies (example/1000G/hprc_satellites.py):
an assembly is a truth only for the arrays it spans. The first version dropped arrays annotated
together with a gap ("GAP,HSat2") without a word, and HSat2 looked over-estimated by half."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "example" / "1000G"))
import hprc_satellites as h  # noqa: E402


def test_annotation_labels():
    cases = {"HSat3": "HSat3", "GAP,HSat3": "HSat3", "active_hor(S1C1H1L)": "aSatHOR", "hor(S3C9H3-B)": "aSatHOR",
             "cenSat(SST1,SST1v)": "SST1", "cenSat(ACRO1,COMP-subunit_ACRO_rnd-1_family-2,COMP-subunit_ACRO_rnd-1_family-3)": "ACRO",
             "bSat(COMP-subunit_LSAU-BSAT_rnd-1_family-1,LSAU)": "bSat", "GAP,cenSat(CER)": "CER", "cenSat(SATR1,SATR2)": "SATR",
             "ct": None, "mon": None, "dhor(S4C20H7d)": None, "rDNA": None, "GAP": None}
    for name, want in cases.items():
        labels = h.labels_of(name)
        assert next((h.CLASS_OF[x] for x in labels if x in h.CLASS_OF), None) == want, (name, labels)
    assert "GAP" in h.labels_of("GAP,HSat2") and "GAP" not in h.labels_of("cenSat(SST1,SST1v)")


def test_gap_containing_arrays_are_kept_apart(tmp_path):
    for hap, lines in (("hap1", ["c1\t0\t30000000\tHSat2\t0\t.", "c1\t40000000\t57000000\tGAP,HSat2\t0\t.", "c2\t0\t80000000\tHSat3\t0\t."]),
                       ("hap2", ["track name=x", "c1\t0\t14000000\tHSat2\t0\t.", "c2\t0\t7000000\tHSat3\t0\t.", "c2\t9000000\t9100000\tGAP,HSat3\t0\t."])):
        (tmp_path / f"S1_{hap}_hprc_r2_v1.cenSat.bed").write_text("\n".join(lines) + "\n")
    mass, gapped, n_files = h.assembly_mass("S1", str(tmp_path))
    assert n_files == 2
    assert mass["HSat2"] == 44_000_000 and gapped["HSat2"] == 17_000_000      # a third of the class sits in an array nobody spanned
    assert mass["HSat3"] == 87_000_000 and gapped["HSat3"] == 100_000         # immaterial: 0.1%
    assert gapped["HSat2"] > h.MAX_GAPPED * (mass["HSat2"] + gapped["HSat2"]) and gapped["HSat3"] < h.MAX_GAPPED * mass["HSat3"]


def test_rdna_the_assemblies_hold(tmp_path):
    """The rDNA an assembly holds: records less than 1 kb apart on one contig are one stretch, a gap is
    not sequence, and a stray rDNA-like fragment elsewhere still counts toward the total."""
    from ngsdose.hprc import rdna_in_assembly
    for hap, lines in (("hap1", ["c1\t20000\t420000\trDNA\t0\t.", "c1\t420500\t520000\trDNA\t0\t.", "c1\t520000\t530000\tcenSat(ACRO1,COMP)\t0\t.",
                                 "c1\t36200000\t36203000\trDNA\t0\t.", "c2\t0\t2000000\tGAP,rDNA\t0\t."]),
                       ("hap2", ["track name=x", "c3\t28000\t489000\trDNA\t0\t.", "c3\t489000\t489100\tGAP\t0\t."])):
        (tmp_path / f"S1_{hap}_hprc_r2_v1.cenSat.bed").write_text("\n".join(lines) + "\n")
    total, longest, n_files = rdna_in_assembly("S1", str(tmp_path))
    assert n_files == 2
    assert longest == 500_000                                   # 20,000-520,000: the 500-bp break is not a gap between arrays
    assert total == 500_000 + 3_000 + 461_000                   # the 2-Mb gap record is not counted
    assert rdna_in_assembly("S2", str(tmp_path)) == (0, 0, 0)
