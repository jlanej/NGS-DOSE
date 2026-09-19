# Experimental resources

`satellites.CHM13v2.k31.panel.tsv.gz` — compositional classes HSat1A, HSat1B, HSat2, HSat3, bSat
(β-satellite) and aSatHOR (α-satellite higher-order repeats), 1.1 M k-mers, built by
`resources/build/build_satellite_panel.sh` from the CHM13 CenSat annotation.

Use it in **scan mode only** (`ngs-dose count -m scan -p this.panel.tsv.gz -c …/controls.fa.gz`):
satellite reads are placed all over a GRCh38 alignment (3–54% on decoy and unplaced contigs,
the rest on centromere models and other primary-assembly satellite), and no sinks have been
learned. The estimate is diploid array mass, `CLASS.mass_Mb`.

Status: it runs (NA12878, 1 min 50 s) and returns masses within a factor of two of twice
CHM13's haploid content in every class — HSat1A 19.7 Mb, HSat2 49.7, HSat3 53.4, bSat 11.3,
aSatHOR 150.1 — and 1.9 Mb of HSat1B in this female sample, which is right for a family that
lives mostly on Yq.

A first comparison with assemblies of the same people (HPRC release 2, CenSat annotation of
both haplotypes summed; `example/1000G/hprc_satellites.py`), two samples:

| class | HG02258 (ACB, male): assembly Mb / NGS-DOSE Mb | HG01884 (ACB, female) |
| --- | --- | --- |
| HSat1A | 28.4 / 26.5 | 24.9 / 24.2 |
| HSat1B | 11.8 / 10.6 | 2.2 / 1.8 |
| HSat2 | 44.0 / 67.0 | 30.5 / 52.8 |
| HSat3 | 87.4 / 85.1 | 71.8 / 70.3 |
| bSat | 15.6 / 10.7 | 18.1 / 13.8 |
| aSatHOR | 135.2 / 130.8 | 148.8 / 153.2 |

HSat3, HSat1A and the α-satellite HORs are within 7% of the assembly in both samples, HSat1B
within 10% in the male and 17% in the female (who has 2 Mb of it); HSat2 reads 50–75% high and
β-satellite 25–30% low. Two samples cannot say whether the estimates
*track* the assemblies across people, which is the question that matters for association work;
200 samples of the 1000 Genomes cohort have HPRC assemblies, and
`example/1000G/04_hprc_satellites.sh` makes the comparison (each assembly also carries ~30 GAP
annotations in its satellite arrays, so the truth is itself a lower bound). Until then: HSat2
and bSat should not be used, and the other four are promising, not validated.
