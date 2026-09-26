# GRCh38 resource bundle (GRCh38-v1)

Everything NGS-DOSE needs that is specific to a reference build. Contig names are chr-prefixed;
the bundle was built on the 1000 Genomes analysis set (`GRCh38_full_analysis_set_plus_decoy_hla.fa`)
and works with any GRCh38 that carries the same primary and unplaced/random contigs. Two more
things depend on the analysis set's extra contigs. The chrEBV dosage region needs chrEBV, which the
analysis set carries but plain UCSC/GRC GRCh38 does not. Two of the 59 DJ sink intervals are on
hs38d1 decoy contigs (chrUn_JTFH01001783v1_decoy, chrUn_JTFH01001847v1_decoy). Where those contigs
are missing, `count` records the dosage region as absent, and a fetch leaves the two intervals out
with a warning and records them in the counts file (`sinks_skipped`); `ngsdose estimate` then
reports DJ as not estimated (status `sinks_skipped`, NA in the table) rather than a value that may
be low. For such a reference, re-learn the sinks (see "When to rebuild what" below).

| file | what | how it was made |
| --- | --- | --- |
| `panel.k31.tsv.gz` | class-diagnostic 31-mers: `rDNA45S` (31,827), `rDNA5S` (1,945), `DJ` (169,808) | `ngs-dose panel`; k-mers absent from GRCh38 and T2T-CHM13v2.0 outside `build_inputs/*.class_loci.bed` |
| `units/` | unit sequences of the positional classes | GenBank KY962518.1 (45S), X12811.1 (5S); CHM13v2.0 chr21:2,708,299-3,108,298 (distal junction) |
| `features.tsv` | named features in unit coordinates (18S, 5.8S, 28S, spacers) | KY962518.1 feature table; 5S gene located by sequence |
| `controls.bed`, `controls.fa.gz` | 800 control regions (10.1 Mb); 80 held-out autosomal, 60 chrX and 40 chrY known-truth regions; one dosage region each on chrM and chrEBV; all with 1 kb flanks | `resources/build/select_controls.py` on the complement of NGS-PCA's exclusion set, then `ngs-dose controls` |
| `sinks.bed` | where the aligner puts each class's reads (fetch mode retrieves only these): 80 intervals for the positional classes (`rDNA45S` 19, `rDNA5S` 2, `DJ` 59; 3.34 Mb as written, of which 40.6 kb run past the ends of eight short contigs and are clipped when fetched) and, since 2026-09-22, 63 intervals (810 kb) for the telomeric repeat (`TEL`), so a fetch with `resources/experimental/telomere.k31.panel.tsv.gz` loaded measures it. Share of each class's reads inside its sinks over the 1,748 cohort scans counted by 2026-09-25 (mean / 1st percentile / minimum): `rDNA45S` 99.95 / 99.92 / 99.89%, `rDNA5S` 99.98 / 99.96 / 99.93%, `DJ` 99.76 / 99.68 / 99.65%, `TEL` 99.86 / 99.68 / 99.39% (15 of the 1,748 below 99.67%) | positional classes: `ngsdose sinks` on whole-file scans of NA12878 (CEU, female) and HG02258 (ACB, male), NYGC pipeline (bwa-mem 0.7.15, ALT-aware), on a 10-kb placement grid (scans now record 1 kb, on which the cohort run re-learns them tighter); they held ≥ 99.94% of 45S, 99.98% of 5S and 99.77% of DJ reads in those two, and 99.96 / 99.99 / 99.75% in HG01884, which they never saw. `TEL`: `ngsdose sinks --classes TEL` on the first 372 cohort scans (`bundle.json`, `sinks_learned_from.TEL`) |
| `anchors.json` | the 45S windows that set the absolute level | windows on which the pilot's three sequencing chemistries agree (NovaSeq against HiSeq 2500 and HiSeq 2000), from its replicate pairs (NGS-DOSE-1000G: `python pilot/evaluate_pilot.py --write-anchors`) |
| `build_inputs/` | class loci masked in the background genomes, the DJ core positions, the class manifest | see `resources/build/build_grch38_bundle.sh` |

`bundle.json` also lists the GRCh38 lengths of the primary contigs; `ngsdose estimate` compares
them with the lengths recorded in a counts file and refuses a file aligned to another build.

The capture figures are the share of each class's scan reads whose 10-kb (compositional) or 1-kb
(positional) placement bin lies wholly inside the class's sink intervals, as `ngsdose sinks
<scans> --evaluate resources/GRCh38/sinks.bed` reports it per sample; a read in the unmapped bin
counts as missed, although a scan classifies it (a fetch reads the unmapped bin only with
`ngs-dose count --unmapped`). In these scans the unmapped bin holds none of the `TEL` or `rDNA5S`
reads and at most 4e-5 of the `DJ` and 1.3e-6 of the `rDNA45S` reads of any sample. Regenerate the
figures when `sinks.bed` changes.

`resources/build/build_grch38_bundle.sh` rebuilds the reference-derived files (`units/`,
`build_inputs/`, the panel and the controls) from public inputs, each checked against the sha256
the shipped bundle was built from (`VERIFY=0` accepts other copies). It writes to
`work/bundle/GRCh38` and refuses to write over this directory unless `FORCE=1`; its
`build_manifest.tsv` lists every input and output with its sha256 and says which outputs are
identical to the shipped ones (a rebuild from the recorded inputs reproduces all of them byte for
byte). The NGS-PCA exclusion BED is downloaded at a pinned commit unless `EXCLUDE` names a local
copy. Sinks and anchors need sequencing data and are produced by the commands named above;
`bundle.json`, `features.tsv` and this README are written by hand. `tests/test_doc_numbers.py`
checks the interval, region and k-mer counts quoted here against the files.

## The distal junction class

The DJ is the ~400 kb of sequence immediately distal to the rDNA array on each acrocentric short
arm. The `DJ` class uses CHM13 chr21's copy as the coordinate system and keeps only *core*
k-mers: those that occur once in chr21's distal junction and four times in total across the other
four CHM13 distal junctions (chr13, 14, 15, 22), and nowhere else in CHM13 or, outside DJ-like
loci, in GRCh38 (`resources/build/dj_core.py`). Its expected diploid copy number is therefore 10,
which makes it a known-truth control for the k-mer path on multi-copy acrocentric sequence.
(`bundle.json` records this as `expected_copies`.) The rule tests the total, not one copy in each
junction: 778 of the 169,808 panel k-mers (0.46%, clustered in about 58 1-kb stretches of the
unit) are split unevenly among the five, e.g. (0,1,1,1,2), so their dosage follows duplications
and deletions of single junctions. The 10 holds for CHM13 in aggregate.

## When to rebuild what

- **Different aligner, decoy set or ALT handling** (DRAGEN, GRCh38 without decoys, …): re-learn
  `sinks.bed` with `ngsdose sinks` from whole-file scans (`ngs-dose count -m scan`) of a handful of
  that pipeline's samples, and check its capture on other scans with `ngsdose sinks --evaluate`; the
  panel and controls stay. Sinks are specific to the aligner and the reference: under DRAGEN 4.x
  with an alt-masked reference, the one genome checked had most of its 45S and DJ reads, and 64–90%
  of its HSat1A, HSat1B, β-satellite, ACRO and TEL reads, fully unmapped, but almost none of its
  HSat2 or α-satellite HOR reads. A scan classifies unmapped reads; a fetch reads them only with
  `count --unmapped`. A fetch refuses a loaded class the sinks say nothing about, or whose intervals
  all lie on contigs the input lacks (`--allow-missing-sinks` overrides and records it in
  `sinks_missing_classes`), so a panel added later needs its sinks learned first
  (`ngsdose sinks --classes`). Intervals on absent contigs are left out and recorded as
  `sinks_skipped`, and `ngsdose estimate` leaves such a class unestimated.
- **Different library chemistry**: re-learn window efficiencies with `ngsdose cohort`; check the
  known-truth columns; treat absolute values with caution until anchors have been confirmed for
  that chemistry (DESIGN.md §7).
- **Different reference build** (GRCh37, CHM13): new `controls.*` and `sinks.bed`; the panel is
  reference-independent as long as the build's own class loci are added to the background masks.
