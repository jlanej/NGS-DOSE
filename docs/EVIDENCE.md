# The case that it works

What would convince someone who doubts that rDNA copy number can be measured from short-read
whole-genome sequencing? Not a model, and not agreement with ourselves. This page collects the
evidence that exists so far, in the order a sceptic would ask for it, with the number, what it
rules out, and where it comes from. Everything here is recomputed by scripts from committed
data: the [cohort page](https://jlanej.github.io/NGS-DOSE-1000G/) (`ngsdose report`, 735 of
3,202 genomes and 149 of 602 trios as of 2026-09-23) and the
[pilot report](../example/1000G/pilot/pilot_report.md) (`evaluate_pilot.py`, 12 genomes each
sequenced twice). The numbers below are those of that date; the pages carry the current ones.

![the evidence](evidence.png)

## 1. It reads known copy numbers correctly, in every genome

Every sample carries sequence whose copy number is not in question, measured by exactly the
code that measures the rDNA: 80 held-out autosomal regions (two copies), 60 chrX and 40 chrY
regions (one or two, one or none, by sex), and the *distal junction*, a 400-kb sequence present
once on each of the ten acrocentric short arms — multi-copy, paralogous, acrocentric, i.e.
the kind of sequence the rDNA is.

| known quantity | reads | n |
| --- | --- | --- |
| held-out autosomal sequence (2) | 1.999 ± 0.008 | 735 |
| chrX in men (1) / in women with an intact culture (2) | 0.991 ± 0.006 / 1.933 ± 0.025 | 350 / 370 |
| chrY in men with an intact Y (1) / in women (0) | 0.975 ± 0.014 / 0.007 at most | 343 / 385 |
| distal junction (10) | 9.66, robust SD 0.15 | 614 |

Sex read from the X and Y agrees with the pedigree in 735 of 735: the highest X in a man is
1.00, the lowest in a woman with an intact culture 1.85. Fifteen women and seven men read below
those lines because their cell line has lost an X or a Y in part of its cells, and an independent
coverage pipeline reads the same people the same way (finding 6). **What it rules out:** that the
fragment-GC model, the k-mer assignment or the control regions are wrong in a way that shows.
*Panel a.*

## 2. A ten-copy paralog steps in whole copies, and the steps are inherited

Relative to the cohort's level, the distal junction sits at whole numbers: 4 people at −2,
16 at −1, 614 at 0, 21 at +1, with a robust SD of 0.15 copies around each (51 people sit between
steps). A step is a structural variant of an acrocentric short arm; a two-copy loss is the
signature of a Robertsonian translocation (about one person in a thousand carries one). HG00651
and her daughter HG00652 both read −2, and both also carry 10–30% less of every satellite family
of the acrocentric short arms — ACRO, SST1, β-satellite, HSat3, CER — while the pan-centromeric
α-satellite that every chromosome carries is unchanged: the arms are missing, not just the
junction. A second family, HG01204 and his child HG01206, also reads −2 in both, but only the
child carries less of the acrocentric satellites (ACRO 0.77 and SST1 0.71 of the cohort's median,
against the father's 1.11 and 1.12): not every two-copy loss of the junction takes whole arms
with it. Where a carrier parent and a child were both counted, the step was transmitted in
7 of 22 — fewer than the half a heterozygous variant passes on, but within chance at 22 pairs
(two-sided p = 0.13) — and one child, HG01517, carries a step that neither parent has (+0.77
copies): a new structural variant, or a change in part of its cell line. **What it rules out:**
that multi-copy acrocentric sequence cannot be resolved to a single copy by this path. *Panel b.*

## 3. The measured rDNA variation is inherited

A child's dosage is the average of the parents' plus segregation; measurement error is not
inherited. In 149 trios the midparent slope for the calibrated 45S estimate is 1.05 ± 0.07,
i.e. a reliability of 1 (95% interval 0.93–1.21); the interval allows a measurement error of
at most 6.4% of a person's value. Held-out autosomal sequence, which has no true variance, reads
−0.07 (−0.32 to 0.19); the cell line's mitochondrial content, which is not in the nuclear genome,
0.06 (−0.19 to 0.29). **What it rules out:** that the between-person variation is a property of
the sample preparation or the culture rather than the genome. *Panel c.*

Three honest notes. First, the design helps: every one of the 149 children is among the 698
related genomes the 1000 Genomes 30× release added to its original 2,504, and 287 of their 298
parents are among the 2,504, so no batch is shared within a family to imitate inheritance. The
same design means a difference between the batches is a difference between the generations:
the children read 4.8% more 45S than their parents, which is generation or batch, and a batch
that reads the children on a slightly different scale moves the midparent slope with it. With
the children first put on their parents' scale, the 45S reads 1.03 (0.94–1.13); for HSat2 and
ACRO, whose children spread 1.15× and 0.86× as much as their parents, the two readings part
(1.23 against 1.07, and 0.86 against 0.99). The "Mendelian" reliability, 1.5 − Var(child −
midparent) / Var(parent), turns out to be R + 1 + ρ/2 − s² with s that ratio of spreads, so on
this design it measures the batch; the page no longer shows it. Second, the spousal correlation
after centring within population is 0.17 (−0.01 to 0.34), down from 0.36 at 42 trios: spouses
share no DNA, so it is either structure shared within couples or chance, and it barely moves a
reliability near 1. Third, because people differ in 45S copy number by 22% (CV) and any
competent estimator errs by a few percent, *every* estimator has a reliability near 1 within one
pipeline: the trios establish that the variation is real, not which estimator measures it best
(the paired difference between the calibrated estimate and the published 18S ratio is +0.03,
−0.01 to +0.07). Ranking estimators takes the next two findings. The page's section 3.6 ("Class
by class") gives every metric in these terms.

## 4. The same person, sequenced twice, years apart, on different instruments

Twelve pilot genomes have an independent older library (HiSeq 2500 2×126 or HiSeq 2000 2×100,
2012–15) beside their NovaSeq 2×150 library (2019): different chemistry, read length, insert
size, depth, aligner, and a GC response that is the reverse of NovaSeq's. Test–retest
reliability (ICC) across the two technologies:

| estimator | ICC | within-person CV | offset between technologies |
| --- | --- | --- | --- |
| calibrated, anchors chosen out of sample | **0.978** | 2.8% | +2% |
| 18S depth ratio, as the literature computes it | 0.19 | 22.6% | −27% |
| the same, after removing its offset (a batch correction) | 0.87 | 7.3% | — |

The two DNA batches come from different cultures of each cell line, so these are upper bounds
on the measurement error. **What it rules out:** that the number is a property of the library.
This is the finding that separates the calibrated estimate from a depth ratio. *Panel d.*

## 5. What the model removes is the library, not the person

Even within one chemistry, libraries differ in GC bias: across the cohort the rate at which
65%-GC fragments were sequenced, relative to each library's mean, runs from 1.18 to 1.29
(middle 80%). The 18S depth ratio divided by the calibrated estimate of the same sample follows
that bias with **r = 0.84** (0.82–0.86); the same 18S region under the fragment-GC model,
r = 0.08 (0.01–0.15). Within the cohort the effect is a few percent (SD of the log ratio 0.033)
— small next to how people differ, which is why finding 3 cannot see it — but it is systematic,
and across technologies (finding 4) it is 27%. **What it rules out:** that the GC model is
decoration. *Panel e.*

## 6. Another pipeline, the same files

Hall, Turner & Queitsch (2021) published rDNA copy number for 2,419 of these genomes, from the
same CRAMs, as the 18S depth relative to chromosome 1 with duplicate-flagged reads excluded.
On the 551 samples shared so far: r = 0.984. Their values are 1.08× ours; re-applying their
duplicate exclusion to our counts brings this to 1.03 — most of the offset is the duplicate
flag, which is set for 5.2% of rDNA reads but 8.9% of single-copy reads (the collapsed rDNA
hides duplicates from the marker), by an amount that differs between samples (0.43–0.82 of the
control rate). **What it rules out:** that the number is idiosyncratic to this code — and it
explains the one difference. *Panel f.*

A second independent pipeline covers the dosage path: NGS-PCA's per-sample QC of the same CRAMs
(mosdepth coverage, duplicate-flagged reads excluded) gives mitochondrial copies per cell that
agree with ours at r = 0.995 on 735 genomes (theirs 0.92× ours, again the duplicate flag: a
16.6-kb genome at thousands-fold depth saturates the positions a duplicate marker distinguishes,
and the ratio falls with depth), a chrX coverage ratio that agrees at r = 0.9998, sex in 735 of
735, and the same fifteen women and seven men with partial loss of an X or a Y. The page's
section 3.7 carries it (`--qc`).

## 7. Against long-read assemblies

Assemblies collapse the rDNA, so they are no truth for it; but 54 of the 735 genomes belong to
people with an HPRC release-2 assembly, whose CenSat annotation of both haplotypes gives the size
of every satellite array — the same kind of sequence, measured by the same k-mer machinery. Per
genome the two measurements agree closely for most families: robust SD of the log ratio 3.1% for
α-satellite HORs and CER, 3.5% ACRO, 5.5% β-satellite, 5.9% HSat1B, 7.4% HSat1A and 7.5% HSat3.
β-satellite and CER sit on a line below equality because the panel sees only part of them (their
k-mer recall); SST1 and SATR sit further below, because the HPRC annotation labels several times
more sequence as those families than the CHM13 annotation the panels were built from. Across people only HSat1B (r = 0.99) and ACRO (0.96) keep r ≥ 0.95,
because r also depends on how much people differ: α-satellite HOR mass differs between people by
4%, so its r is 0.71. When only six people could be compared, seven families passed that mark;
six people flattered the correlations. All 16 outlying comparisons (in HSat1A, HSat1B, HSat2,
HSat3 and ACRO) are genomes whose assembly holds less than the reads show, as where part of an
array is missing from an assembly without a marked gap. HSat2 compares poorly (r = 0.30 in the
34 assemblies that close its arrays), though it is inherited as faithfully as the other classes
(R = 1.07 with the children on their parents' scale): what it measures is heritable, but the
assemblies do not yet confirm that it is the mass of HSat2. For SST1 and SATR the two
measurements disagree per genome by more than people differ. **What it rules out:** that the
k-mer path only appears to work because nothing independent has been held against it.
*Panel g.* The page's section 3.8 carries the figures, one per family.

## 8. The one-minute fetch equals the whole-file scan

Fetch mode reads the control regions and the few intervals where the aligner puts class reads —
about 0.5 GB of a 15-GB CRAM, a minute over the network, six seconds from disk. For all 735
genomes counted both ways it returns 0.9997 of the scan's 45S estimate (range 0.9992–0.9999),
0.9999 of the 5S and 0.998 of the distal junction; the lowest 45S sink capture in any scan is
99.91%. Run separately on the fetch counts alone, the whole pipeline — calibration, cohort layer
and trio test — gives the same 45S reliability. **What it rules out:** that a biobank would need
the whole files. *Panel h.*

## And PCs?

Coverage principal components (NGS-PCA's, or the internal ones from the control regions) are
regressed out of the estimates, with the number chosen at the Marchenko–Pastur edge of the
noise bulk and checked against the known truths and the trios. On this cohort they matter
little for the rDNA, and the reason is worth stating: measurement error is at most a few
percent of a person's value and the variation between people is 22%, so the technical share of
the 45S variance is about 2% — the 16 control PCs above the edge remove 4% against 2.2% expected
by chance. Where there *is* technical variance they find it: the same 16 PCs remove 34% of the
variance of the held-out autosomal estimate (which has nothing but error to remove), 22% of the
mitochondrial and 20% of the EBV dosage. The cross-validated sweep picks 5 PCs for the autosomal
control and none for the rDNA. PCs are insurance for small effects on a noisy trait; for rDNA
copy number the GC model and the calibration do the work.

## What is not yet shown

- No orthogonal assay of rDNA copy number exists for these samples; the absolute scale rests
  on unit windows where three Illumina chemistries agree, and one ddPCR value (CHM13).
- 5S copy number: reliability 0.78 with an interval from 0.46 to 1.11 at 149 trios — undecided.
- In these trios generation and batch go together (finding 3): whether children differ from
  their parents in level or spread cannot be told from the sequencing batch.
- All 735 genomes are one chemistry and one pipeline; DRAGEN alignments and other chemistries
  are untested.
- Every sample is a lymphoblastoid cell line.

## Reproduce

```bash
ngsdose report --scan counts_scan/ --fetch counts_fetch/ -p pedigree.txt --hall hall2021.txt --pilot example/1000G/pilot --censat hprc_censat/ --qc ngspca_sample_qc.tsv -o docs/
python example/1000G/evidence_figure.py --report docs/report.json --pilot example/1000G/pilot -o evidence.png
```
