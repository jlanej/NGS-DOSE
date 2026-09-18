# NGS-DOSE

Sequence-class dosage from WGS: copy number and array mass for rDNA, satellite,
macrosatellite and telomere classes.

**NGS-DOSE is not a CNV caller.** It measures how much of a multi-copy sequence class an
individual carries — the fraction of the genome that variant callers and CNV callers
explicitly mask. [NGS-PCA](https://github.com/jlanej/NGS-PCA) computes its denominator
from the bins it retains; NGS-DOSE measures the bins it excludes.

## What it does

For each class, `2 * observed_bases / autosomal_depth` gives diploid array mass, and
dividing by the unit length gives diploid copy number. Whether a class gets copies or
megabases is a property of the class, not a preference: rDNA, 5S, SST1, D4Z4 and DXZ4
have a defined unit; HSat1/2/3 do not, and reporting "copies" for them would invent
precision that does not exist.

Depth extraction is delegated to [mosdepth](https://github.com/brentp/mosdepth). What
NGS-DOSE owns is the part that decides whether the numbers mean anything:

- **Denominator** — NGS-PCA's `AUTO_HQ_median`, the per-sample median over bins that
  survive the exclusion set, so it is single-copy by construction.
- **GC** — bias curve fitted from the sample's own retained bins and interpolated; the
  rDNA unit is 58.1% GC against 40.9% genome-wide, where the curve is steepest.
- **Batch** — residualised on NGS-PCA coverage components rather than batch labels. The
  PC basis excludes satellite and segmental duplications, so it is disjoint from every
  target class: it can absorb library structure but not the dosage being measured.
- **Validation** — *transmission reliability*: copy number is inherited additively with
  a midparent coefficient of exactly 1, so the midparent–offspring regression slope
  estimates the reliability of the measurement rather than its heritability. Run per
  class, it is a direct test of which classes are measurable at all.

Sampling precision at 30× is 0.03–0.38%, while r = 0.015 reaches significance at
n = 127,000. The error budget is essentially all systematic, which is why the list
above is the method and the counting is not.

## Status

Design only. `docs/DESIGN.md` is the specification: estimator, class table, the four
confounders, the recipe, and the validation plan. A working reference implementation
with a data-free self-test exists and can be brought in as a starting point.

Not yet done: the k-mer backend is specified but unimplemented, and nothing has been run
against real WGS. First result that would justify a release is the 5S negative control
on 1000 Genomes 30×.

## Provenance and credit

The method design, statistics and reference implementation were developed by **Claude
Opus 5 (Anthropic)** in a September 2026 working session with **@jlanej**, as an
offshoot of a review of the acrocentric short arms. The naming, the scoping and the
decision to build it are shared; the errors are worth attributing to the machine that
made them until a human has checked each one.

It builds directly on prior work that should be cited ahead of this repository:
NGS-PCA for the denominator and the coverage components, mosdepth for depth extraction,
the T2T-CHM13 CenSat annotation for the class definitions, and Rodríguez-Algarra,
Evans & Rakyan (*Cell Genomics* 4:100562, 2024) for the demonstration that rDNA copy
number carries real phenotypic signal.

Two honesty notes. The transmission-reliability framing is a rederivation of standard
midparent regression, and no literature search has been done to confirm it is novel in
this application. And an AI system is not an author under prevailing journal policy — if
this becomes a paper, the appropriate form is a contributions or acknowledgements
statement describing what was machine-generated, with human authors taking
responsibility for verification.
