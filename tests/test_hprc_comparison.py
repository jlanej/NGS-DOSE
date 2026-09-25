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
    assert total == 500_000 + 3_000 + 461_000                   # the 2-Mb gap record is not counted
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
