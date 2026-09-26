"""ngsdose.hprc: what the HPRC CenSat annotations of a person hold, read once whether plain or gzipped."""


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
    assert total == 499_500 + 3_000 + 461_000                   # the rDNA records' own bp; the 2-Mb gap record is not counted
    assert rdna_in_assembly("S2", str(tmp_path)) == (0, 0, 0)


def test_gzipped_annotations_read_the_same(tmp_path):
    """The results repository keeps the annotations gzipped; a haplotype present in both forms counts once."""
    import gzip

    from ngsdose.hprc import annotation_files, assembly_mass, rdna_in_assembly
    lines = {"hap1": "c1\t0\t30000000\tHSat2\t0\t.\nc1\t20000\t420000\trDNA\t0\t.\n", "hap2": "c2\t0\t7000000\tHSat3\t0\t.\nc3\t28000\t489000\trDNA\t0\t.\n"}
    plain, packed = tmp_path / "plain", tmp_path / "packed"
    plain.mkdir(), packed.mkdir()
    for hap, text in lines.items():
        (plain / f"S1_{hap}_hprc_r2_v1.cenSat.bed").write_text(text)
        with gzip.open(packed / f"S1_{hap}_hprc_r2_v1.cenSat.bed.gz", "wt") as fh:
            fh.write(text)
    assert assembly_mass("S1", str(packed)) == assembly_mass("S1", str(plain))
    assert rdna_in_assembly("S1", str(packed)) == rdna_in_assembly("S1", str(plain)) == (400_000 + 461_000, 461_000, 2)
    (packed / "S1_hap1_hprc_r2_v1.cenSat.bed").write_text(lines["hap1"])     # both forms of one haplotype
    assert len(annotation_files("S1", str(packed))) == 2 and assembly_mass("S1", str(packed))[2] == 2


def test_a_gap_record_ends_an_rdna_stretch(tmp_path):
    """rDNA, a 100-bp placeholder GAP, rDNA: two stretches, however close; the gap is not sequence."""
    from ngsdose.hprc import rdna_in_assembly
    (tmp_path / "S1_hap1_hprc_r2_v1.cenSat.bed").write_text("c1\t0\t100000\trDNA\t0\t.\nc1\t100000\t100100\tGAP\t0\t.\n"
                                                             "c1\t100100\t200100\trDNA\t0\t.\n")
    (tmp_path / "S1_hap2_hprc_r2_v1.cenSat.bed").write_text("c1\t0\t50000\trDNA\t0\t.\n")
    assert rdna_in_assembly("S1", str(tmp_path)) == (250_000, 100_000, 2)


def test_family_labels_follow_the_panel_build():
    """As build_satellite_panel.sh: the first family inside cenSat(...), ACRO/SST1/SATR by prefix, CER exactly."""
    from ngsdose.hprc import CLASS_OF, labels_of
    cls = lambda name: next((CLASS_OF[x] for x in labels_of(name) if x in CLASS_OF), None)
    assert cls("cenSat(SATR1v)") == cls("cenSat(SATR1v,SATR2)") == cls("cenSat(SATR1,SATR2)") == "SATR"
    assert cls("cenSat(SST1v)") == "SST1" and cls("cenSat(ACRO1,COMP-subunit_ACRO_rnd-1_family-2)") == "ACRO"
    assert cls("cenSat(HSAT5v1,SST1,SST1v)") is None                # named after another family first, as the build leaves it out
    assert cls("cenSat(CERv)") is None and cls("GAP") is None and cls("active_hor(S1C1H1L)") == "aSatHOR"
    assert CLASS_OF.get("SATR1v") == "SATR" and CLASS_OF.get("ct") is None


def test_a_standalone_gap_inside_an_array_is_gapped(tmp_path):
    """An active_hor, GAP, active_hor sequence: the GAP's span is gapped aSatHOR, not closed array; a
    placeholder (100 bp) is counted, and leaves the sample out only when asked; SATR1v counts as SATR;
    labels no class takes are collected when asked."""
    import collections

    from ngsdose.hprc import assembly_mass, compare
    (tmp_path / "S1_hap1_hprc_r2_v1.cenSat.bed").write_text(
        "c1\t0\t1000000\tactive_hor(S1C7H1L)\t0\t.\nc1\t1000000\t1100000\tGAP\t0\t.\nc1\t1100000\t2000000\tactive_hor(S1C7H1L)\t0\t.\n"
        "c1\t2000000\t2000100\tGAP\t0\t.\nc1\t2000100\t2010100\tcenSat(SATR1v,SATR2)\t0\t.\nc2\t0\t5000\tcenSat(FOO1)\t0\t.\n")
    (tmp_path / "S1_hap2_hprc_r2_v1.cenSat.bed").write_text("c1\t0\t2000000\tactive_hor(S1C7H1L)\t0\t.\nc1\t2000000\t2010000\tcenSat(SATR1)\t0\t.\n")
    unknown = collections.Counter()
    mass, gapped, n = assembly_mass("S1", str(tmp_path), unknown)
    assert n == 2 and mass["aSatHOR"] == 3_900_000 and mass["SATR"] == 20_000
    assert gapped["aSatHOR"] == 100_000 + 50 and gapped["SATR"] == 50      # the placeholder touches both classes: split
    assert unknown == {"FOO1": 5000}
    est = [{"sample": "S1", "aSatHOR.mass_Mb": "3.9", "SATR.mass_Mb": "0.02"}]
    rows, st = compare(est, str(tmp_path))
    assert {r["cls"]: r["assembly_gapped_Mb"] for r in rows} == {"aSatHOR": 0.1, "SATR": 0.0}
    assert st["aSatHOR"]["n_gapped"] == 1                              # 0.1 of 4.0 Mb is more than MAX_GAPPED
    assert st["SATR"]["n_gapped"] == 0 and st["SATR"]["n_unsized_gap"] == 1
    assert compare(est, str(tmp_path), exclude_unsized=True)[1]["SATR"]["n_gapped"] == 1
