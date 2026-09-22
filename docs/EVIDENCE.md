# The case that it works

What would convince someone who doubts that rDNA copy number can be measured from short-read
whole-genome sequencing? Not a model, and not agreement with ourselves. This page collects the
evidence that exists so far, in the order a sceptic would ask for it, with the number, what it
rules out, and where it comes from. Everything here is recomputed by scripts from committed
data: the [cohort page](https://jlanej.github.io/NGS-DOSE-1000G/) (`ngsdose report`, the first
333 of 3,202 genomes as of 2026-09-22) and the [pilot report](../example/1000G/pilot/pilot_report.md)
(`evaluate_pilot.py`, 12 genomes each sequenced twice). The numbers below are those of that
date; the pages carry the current ones.

![the evidence](evidence.png)

## 1. It reads known copy numbers correctly, in every genome

Every sample carries sequence whose copy number is not in question, measured by exactly the
code that measures the rDNA: 80 held-out autosomal regions (two copies), 60 chrX and 40 chrY
regions (one or two, one or none, by sex), and the *distal junction*, a 400-kb sequence present
once on each of the ten acrocentric short arms — multi-copy, paralogous, acrocentric, i.e.
the kind of sequence the rDNA is.

| known quantity | reads | n |
| --- | --- | --- |
| held-out autosomal sequence (2) | 1.998 ± 0.008 | 333 |
| chrX in men (1) / in women (2) | 0.990 ± 0.006 / 1.932 ± 0.028 | 157 / 171 |
| chrY in men (1) / in women (0) | 0.975 ± 0.018 / 0.004 at most | 152 / 176 |
| distal junction (10) | 9.62 ± 0.15 | 287 |

Sex read from the X and Y agrees with the pedigree in 333 of 333: the highest X in a man is
1.00, the lowest in a woman with an intact culture 1.85. **What it rules out:** that the
fragment-GC model, the k-mer assignment or the control regions are wrong in a way that shows.
*Panel a.*

## 2. A ten-copy paralog steps in whole copies, and the steps are inherited

Relative to the cohort's level, the distal junction sits at whole numbers: 2 people at −2,
8 at −1, 287 at 0, 6 at +1, with a robust SD of 0.15 copies around each. A step is a
structural variant of an acrocentric short arm; a two-copy loss is the signature of a
Robertsonian translocation (about one person in a thousand carries one). HG00651 and her
daughter HG00652 both read −2, and both also carry 15–30% less of every satellite family of the
acrocentric short arms — ACRO, SST1, β-satellite, HSat3, CER — while the pan-centromeric
α-satellite that every chromosome carries is unchanged: the arms are missing, not just the
junction. Where a carrier parent and a child were both counted, the step was transmitted in
3 of 6 (half is the Mendelian expectation), and no child carries a step that neither parent
has. **What it rules out:** that multi-copy acrocentric sequence cannot be resolved to a single
copy by this path. *Panel b.*

## 3. The measured rDNA variation is inherited

A child's dosage is the average of the parents' plus segregation; measurement error is not
inherited. In 42 trios the midparent slope for the calibrated 45S estimate is 1.16 ± 0.13,
i.e. a reliability of 1 (95% interval 0.90–1.56); the interval allows a measurement error of
at most 8.6% of a person's value. Held-out autosomal sequence, which has no true variance, reads
−0.05 (−0.56 to 0.32); the cell line's mitochondrial content, which is not in the nuclear genome,
0.12 (−0.27 to 0.39). **What it rules out:** that the between-person variation is a property of
the sample preparation or the culture rather than the genome. *Panel c.*

Two honest notes. First, the spousal correlation after centring within population is 0.36 ±
0.15 — spouses share no DNA, so this is either batch structure shared within families or
chance at n = 42; it barely moves a reliability near 1, and the full cohort will settle it.
Second, because people differ in 45S copy number by 23% (CV) and any competent estimator errs
by a few percent, *every* estimator has a reliability near 1 within one pipeline: the trios
establish that the variation is real, not which estimator measures it best (the paired
difference between the calibrated estimate and the literature's 18S ratio is +0.03,
−0.03 to +0.12). Ranking estimators takes the next two findings.

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
that bias with **r = 0.88**; the same 18S region under the fragment-GC model, r = 0.18. Within
the cohort the effect is a few percent (SD of the log ratio 0.033) — small next to how people
differ, which is why finding 3 cannot see it — but it is systematic, and across technologies
(finding 4) it is 27%. **What it rules out:** that the GC model is decoration. *Panel e.*

## 6. Another pipeline, the same files

Hall, Turner & Queitsch (2021) published rDNA copy number for 2,419 of these genomes, from the
same CRAMs, as the 18S depth relative to chromosome 1 with duplicate-flagged reads excluded.
On the 268 samples shared so far: r = 0.984. Their values are 1.07× ours; re-applying their
duplicate exclusion to our counts brings this to 1.03 — most of the offset is the duplicate
flag, which is set for 5.1% of rDNA reads but 8.6% of single-copy reads (the collapsed rDNA
hides duplicates from the marker), by an amount that differs between samples (0.43–0.82 of the
control rate). **What it rules out:** that the number is idiosyncratic to this code — and it
explains the one difference. *Panel f.*

## 7. Against long-read assemblies

Assemblies collapse the rDNA, so they are no truth for it; but six of the 333 genomes have HPRC
release-2 assemblies, whose CenSat annotations give the size of every satellite array — the same
kind of sequence, measured by the same k-mer machinery. Seven of ten families track the assembly
across the six people with r ≥ 0.97 (HSat1A, HSat1B, HSat3, β-satellite, α-satellite HORs, ACRO,
CER; the last three by a constant factor, their k-mer recall). **What it rules out:** that the
k-mer path only appears to work because nothing independent has been held against it. *Panel g.*

## 8. The one-minute fetch equals the whole-file scan

Fetch mode reads the control regions and the few intervals where the aligner puts class reads —
about 0.5 GB of a 15-GB CRAM, a minute over the network, six seconds from disk. For all 333
genomes counted both ways it returns 0.9997 of the scan's 45S estimate (range 0.9993–0.9999),
0.9999 of the 5S and 0.998 of the distal junction; the lowest sink capture in any scan is
99.92%. **What it rules out:** that a biobank would need the whole files. *Panel h.*

## And PCs?

Coverage principal components (NGS-PCA's, or the internal ones from the control regions) are
regressed out of the estimates, with the number chosen at the Marchenko–Pastur edge of the
noise bulk and checked against the known truths and the trios. On this cohort they matter
little for the rDNA, and the reason is worth stating: measurement error is at most a few
percent of a person's value and the variation between people is 23%, so the technical share of
the 45S variance is about 2% — nine control PCs remove 4% against 2.7% expected by chance.
Where there *is* technical variance they find it: the same nine PCs remove 36% of the variance
of the held-out autosomal estimate (which has nothing but error to remove), 19% of the
mitochondrial and 15% of the EBV dosage. The sweep against the known truths picks 8 PCs for
the autosomal control and none for the rDNA. PCs are insurance for small effects on a noisy
trait; for rDNA copy number the GC model and the calibration do the work.

## What is not yet shown

- No orthogonal assay of rDNA copy number exists for these samples; the absolute scale rests
  on unit windows where three Illumina chemistries agree, and one ddPCR value (CHM13).
- 5S copy number: reliability 0.02 with an interval from −0.8 to 1.1 at 42 trios — undecided.
- The spousal correlation of 0.36 (finding 3).
- All 333 genomes are one chemistry and one pipeline; DRAGEN alignments and other chemistries
  are untested.
- Every sample is a lymphoblastoid cell line.

## Reproduce

```bash
ngsdose report --scan counts_scan/ --fetch counts_fetch/ -p pedigree.txt --hall hall2021.txt --censat hprc_censat/ -o docs/
python example/1000G/evidence_figure.py --report docs/report.json --pilot example/1000G/pilot -o evidence.png
```
