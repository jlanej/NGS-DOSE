# GRCh38 resource bundle (GRCh38-v1)

Everything NGS-DOSE needs that is specific to a reference build. Contig names are chr-prefixed;
the bundle was built on the 1000 Genomes analysis set (`GRCh38_full_analysis_set_plus_decoy_hla.fa`)
and works with any GRCh38 that carries the same primary and unplaced/random contigs.

| file | what | how it was made |
| --- | --- | --- |
| `panel.k31.tsv.gz` | class-diagnostic 31-mers: `rDNA45S` (31,827), `rDNA5S` (1,945), `DJ` (169,808) | `ngs-dose panel`; k-mers absent from GRCh38 and T2T-CHM13v2.0 outside `build_inputs/*.class_loci.bed` |
| `units/` | unit sequences of the positional classes | GenBank KY962518.1 (45S), X12811.1 (5S); CHM13v2.0 chr21:2,708,299-3,108,298 (distal junction) |
| `features.tsv` | named features in unit coordinates (18S, 5.8S, 28S, spacers) | KY962518.1 feature table; 5S gene located by sequence |
| `controls.bed`, `controls.fa.gz` | 800 control regions (10.1 Mb); 80 held-out autosomal, 60 chrX and 40 chrY known-truth regions; one dosage region each on chrM and chrEBV; all with 1 kb flanks | `resources/build/select_controls.py` on the complement of NGS-PCA's exclusion set, then `ngs-dose controls` |
| `sinks.bed` | where the aligner puts each class's reads (fetch mode retrieves only these) | `ngsdose sinks` on whole-file scans of NA12878 (CEU, female) and HG02258 (ACB, male), NYGC pipeline (bwa-mem 0.7.15, ALT-aware); 80 intervals, 3.3 Mb, capturing ≥ 99.94% of 45S, 99.98% of 5S and 99.77% of DJ reads in both |
| `anchors.json` | the 45S windows that set the absolute level | windows on which two sequencing chemistries agree, from the pilot's replicate pairs (`example/1000G/pilot/evaluate_pilot.py --write-anchors`) |
| `build_inputs/` | class loci masked in the background genomes, the DJ core positions, the class manifest | see `resources/build/build_grch38_bundle.sh` |

`bundle.json` also lists the GRCh38 lengths of the primary contigs; `ngsdose estimate` compares
them with the lengths recorded in a counts file and refuses a file aligned to another build.

`resources/build/build_grch38_bundle.sh` regenerates the panel and the controls from public
inputs. Sinks and anchors need sequencing data and are produced by the commands named above.

## The distal junction class

The DJ is the ~400 kb of sequence immediately distal to the rDNA array on each acrocentric short
arm. The `DJ` class uses CHM13 chr21's copy as the coordinate system and keeps only *core*
k-mers: those that occur exactly once in each of the five CHM13 distal junctions (chr13, 14, 15,
21, 22) and nowhere else in CHM13 or, outside DJ-like loci, in GRCh38. Its expected diploid
copy number is therefore 10, which makes it a known-truth control for the k-mer path on
multi-copy acrocentric sequence. (`bundle.json` records this as `expected_copies`.)

## When to rebuild what

- **Different aligner, decoy set or ALT handling** (DRAGEN, GRCh38 without decoys, …): re-learn
  `sinks.bed` from `ngs-dose count -m scan` on a handful of samples; the panel and controls stay.
- **Different library chemistry**: re-learn window efficiencies with `ngsdose cohort`; check the
  known-truth columns; treat absolute values with caution until anchors have been confirmed for
  that chemistry (DESIGN.md §7).
- **Different reference build** (GRCh37, CHM13): new `controls.*` and `sinks.bed`; the panel is
  reference-independent as long as the build's own class loci are added to the background masks.
