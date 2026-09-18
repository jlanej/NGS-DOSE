# NGS-DOSE design specification

Specification for sequence-class dosage estimation from WGS. Written as the design
note for `repeatcn.py` in jlanej/acrocentric_research; that module is the reference
implementation this repository is intended to supersede. Section 6's command lines use
the old module name.

A minimal, honest framing for measuring how much of a multi-copy sequence class an
individual carries, from short-read WGS. Reference implementation: `analysis/repeatcn.py` in jlanej/acrocentric_research,
validated by `python repeatcn.py selftest`, which requires no data.

This is the one measurement short reads support well inside the acrocentric short arms
and other satellite compartments: only 4.3% of short-arm sequence is short-read
accessible for genotyping (report section 4), but *dosage* of a class does not require
any base to be uniquely placed. It is also a measurement with demonstrated phenotypic
associations - rDNA copy number associates with blood cell subtype counts and renal
function in UK Biobank (Rodriguez-Algarra et al., *Cell Genomics* 4:100562, 2024).

## 1. What it measures, and what it does not

Measured: total dosage of a class in the individual, as diploid copies (classes with a
defined unit) or diploid megabases (classes without one).

Not measured: which chromosome the copies sit on, array structure or order, unit
sequence variation, or activity. Per-chromosome apportionment is a separate and harder
problem - see the IGS-barcode deconvolution proposal in report section 22.5.

## 2. The estimator

Let `d` be depth against the haploid reference, so a "30x genome" has `d = 30` and a
diploid single-copy locus also reads `d = 30`. For a class `c` with observed base
count `B_c`:

```
diploid array mass   M_c  = 2 * B_c / d_ctrl              bp
diploid copy number  CN_c = 2 * B_c / (d_ctrl * L_unit)   copies
```

`d_ctrl` is the single-copy autosomal depth, GC-corrected at the class's GC.

**Use NGS-PCA's `autosomal.median.txt` for it.** That file reports `AUTO_HQ_median`,
the per-sample median depth over exactly the bins the PCA retained - autosomal, less
the exclusion set (SV blacklist, 50-mer low mappability, DGV, segmental duplications) -
together with `N_BINS`. It is a single-copy denominator by construction, since the
regions that are not reliably single-copy in an individual are precisely the ones
excluded, and it is already computed as a by-product of the normalisation, with no
extra pass over the data. A plain genome-wide mean depth is the wrong quantity: it
includes duplicated, CNV-variable and unmappable bins, so it varies between individuals
for reasons that have nothing to do with the class being measured.

Three operational notes carried from the NGS-PCA documentation. Values are not floored,
so a failed or empty sample reads as zero rather than as missing - drop those
explicitly, as `load_auto_hq_median()` does. The file is only written by a run that
reads the mosdepth files, not by `-matrix` and not by a run reusing a cached
`tmp.mat.ser.gz` without `-overwrite`. And sample names are the mosdepth file name
minus `regions.bed.gz`, which leaves a trailing dot unless the run set `-sampleSuffix`.

The median and the PCs are complementary rather than redundant, and it is worth being
explicit about why: normalisation divides each sample by this same median before the
SVD, so the components are invariant to a sample's overall depth scale and describe
only the *shape* of its coverage. The median sets the scale; the PCs correct the shape.

**k-mer counting uses the same formulae.** A read of length `R` yields `R - k + 1`
k-mers, so counting in k-mer space multiplies numerator and denominator by the same
factor, which cancels - provided the control is counted the same way. This is why the
choice between a depth backend and a k-mer backend is an engineering choice, not a
statistical one.

**Copies or megabases is a property of the class, not a preference.** rDNA, 5S, D4Z4,
DXZ4 and SST1 have a defined repeating unit, so copy number is meaningful. HSat1, HSat2
and HSat3 have higher-order structure that varies between arrays and between
individuals; reporting "copies" for them invents precision that does not exist. Report
megabases. `repeatcn.py` enforces this from the class table.

## 3. Which classes are ripe

From `data/repeat_cn_targets.tsv`, built from the CHM13v2.0 CenSat v2.1 and
composite-repeat annotations. Reference mass is haploid CHM13; sampling SE is the
Poisson limit at 30x with 150 bp reads.

| Class | Report as | Unit | CHM13 haploid | Arrays / chrom | Sampling SE | Status |
| --- | --- | --- | --- | --- | --- | --- |
| alphaSat HOR | copies (monomers) | 171 bp | 70.33 Mb | 120 / 24 | 0.03% | ripe in aggregate; per-chromosome is hard |
| HSat3 | megabases | - | 69.33 Mb | 189 / 17 | 0.03% | ripe, unvalidated, highly variable |
| HSat2 | megabases | - | 28.71 Mb | 61 / 8 | 0.04% | ripe, unvalidated |
| HSat1B | megabases | - | 15.32 Mb | 64 / 8 | 0.06% | ripe, unvalidated |
| alphaSat monomeric | copies | 171 bp | 13.49 Mb | 191 / 23 | 0.06% | ripe |
| HSat1A | megabases | - | 13.39 Mb | 24 / 8 | 0.06% | ripe; 42% of 13p |
| **rDNA 45S** | copies | 44,838 bp | 9.93 Mb = **221 copies** | 5 / 5 | 0.07% | **ripe and validated** |
| beta-satellite | copies | 68 bp | 8.61 Mb | 184 / 17 | 0.08% | ripe; acrocentric-enriched |
| LSAU-BSAT | megabases | - | 7.04 Mb | 84 / 16 | 0.08% | ripe |
| SST1 | copies | ~1,400 bp | 2.57 Mb | 68 / 13 | 0.14% | ripe; NBL2, cancer hypomethylation |
| alphaSat divergent HOR | copies | 171 bp | 1.86 Mb | 103 / 22 | 0.16% | ripe |
| ACRO repeat | megabases | - | 1.50 Mb | 30 / 12 | 0.18% | ripe; acrocentric-specific |
| gamma-satellite | copies | 220 bp | 0.65 Mb | 64 / 14 | 0.28% | ripe |
| **rDNA 5S** | copies | ~2,213 bp | 0.39 Mb = 174 copies | 48 / 13 | 0.36% | ripe, and the best **negative control** |
| **telomere** | kb per genome | 6 bp (TTAGGG) | 0.35 Mb in reference | 24 / 10 | 0.38% | established (TelSeq); use the k-mer backend - reference telomeres are not sample telomeres |
| D4Z4 | copies | 3,300 bp | needs panel | 4q35, 10q26 | - | ripe and clinically important (FSHD) |
| DXZ4 | copies | 3,000 bp | needs panel | Xq23 | - | ripe |
| mtDNA | copies | 16,569 bp | needs panel | - | - | established; see section 4a |

The 221 haploid rDNA copies the table derives from the annotation match the published
CHM13 rDNA model count exactly, which is a useful arithmetic check on the pipeline.

D4Z4, DXZ4 and mtDNA are absent from the satellite tracks used here and need a panel
built from their own consensus before they can be estimated.

## 4. The four things that ruin it

### (a) A shared denominator

Both numerator and denominator are depths. Two classes normalised by the *same* control
depth inherit a correlation from the denominator's own noise: in logs the shared term
contributes variance `sigma_n^2` to the covariance, giving

```
r_induced = sigma_n^2 / sqrt( (sigma_c1^2 + sigma_n^2) * (sigma_c2^2 + sigma_n^2) )
```

With biological log-SDs of 0.25 and 0.30, a denominator noise of 0.05 manufactures
`r = 0.032` - about twice the genome-wide-significance threshold at n = 127,000. The
self-test reproduces this and shows it vanishing under a split denominator.

**Remedy:** when correlating two classes with each other, normalise them by disjoint
halves of the control windows (`estimate --denominator split`). The rDNA-against-mtDNA
comparison is the one most exposed to this, because both are conventionally computed as
ratios to the same nuclear depth.

**How much of this survives the high-quality median?** The statistical part, none of
it. The median is taken over the retained bins, so its own sampling noise is set by
`N_BINS`: at 2.4 million retained 1 kb bins it is 0.016%, and even at 15,000 bins after
aggressive `-sampleEvery` subsampling it is 0.20%. Feeding those into the formula gives
induced correlations of 5 x 10<sup>-7</sup> and 5 x 10<sup>-5</sup> - three to four
orders of magnitude below the 0.0153 that reaches significance at n = 127,000.
`median_sampling_se()` computes this, and the `estimate` command prints it.

So with a properly built denominator the shared-denominator risk is **systematic, not
statistical**: what remains is bias that moves the median and the class together, and
that is exactly what the coverage PCs in (c) are for. The split denominator is then
belt-and-braces rather than essential - worth keeping for a headline
class-against-class correlation such as rDNA against mtDNA, and not worth the
complexity elsewhere.

### (b) GC bias

The 47S rDNA unit is **58.1% GC** (measured from KY962518.1) against ~40.9% for the
nuclear genome and 44.4% for mtDNA. Satellite classes sit far from the genomic mean in
both directions, which is exactly where a library's GC-bias curve is steepest, and the
curve differs between PCR-free and PCR-based preparations.

`gc_curve()` fits the curve from the sample's own single-copy windows and interpolates
between knots. Interpolation is not cosmetic: with 2.5% GC bins and a realistic bias
curve, nearest-bin lookup leaves a 14% error at rDNA's GC, which interpolation reduces
to under 2% (self-test check 3).

### (c) Batch - estimate it from the data, do not rely on the labels

The UK Biobank analysis had to adjust explicitly for sequencing centre, alongside
assessment centre, genetic principal components, telomere length, age and sex. Recorded
labels are worth including, but they are a poor proxy for the thing that actually
varies: library prep chemistry, extraction method, instrument drift, flowcell, and the
GC-bias curve that follows from all of them. Much of that is unlabelled, continuous, and
changes within a nominal batch.

The better approach is to measure technical structure directly from the coverage data.
[NGS-PCA](https://github.com/jlanej/NGS-PCA) does this: mosdepth coverage in fixed-width
bins (1 kb is the recommended starting point), autosomal bins overlapping an exclusion
set removed (SV blacklist, 50-mer low mappability, DGV variants, segmental
duplications), log2 fold change within each sample relative to that sample's median bin
coverage, each bin then centred to a median of zero across samples, and a randomized SVD
(Halko et al. 2011) for the leading components. Because the centring is by median rather
than mean, this is a robust variant of PCA, and methods text should say "PCA on
median-centered log2 fold change of depth" rather than "PCA". It has been run at
142,000-sample scale.

Regressing a class estimate on those PCs is, in our reading, the right default for this
whole family of measurements, for a reason specific to how the basis is built: **the
excluded regions are exactly the repetitive ones, so the PCs are computed on sequence
disjoint from every class in the target table.** The components can absorb library,
batch and bias structure; they cannot absorb the dosage of the class being measured.
That is a much weaker assumption than the usual worry about over-adjustment.

The precedent matters here. This approach has already been applied to correct batch
effects in TelSeq leukocyte telomere length estimates and to improve depth-based bin
estimates by accounting for the PCs. Telomere length from WGS is structurally the same
measurement as rDNA copy number - a depth or k-mer ratio for a tandem repeat class with
no unique anchors - so the transfer to rDNA and satellite classes is a short step, not a
speculative one.

`repeatcn.py` implements `coverage_pcs()` (exact SVD, for modest cohorts and for the
self-test), `load_pcs()` for NGS-PCA's `svd.pcs.txt`, and `adjust_for_pcs()`, which
residualises log estimates on the components and reports the variance removed. Self-test
check 5 simulates a cohort where batch affects both the class estimate and the
phenotype: the true association is r = +0.09, the raw association is r = **-0.14** - the
wrong sign - and PC adjustment recovers r = +0.12.

### (d) Cell composition, for blood traits

A change in the proportion of a cell type changes the measured average over the sample.
This is well described for mtDNA copy number, where a rising neutrophil fraction lowers
the apparent copy number because neutrophils carry fewer mitochondria - the reverse of
the causal direction one might infer. Adjust for measured cell counts, and test the
composition ratios (neutrophil-to-lymphocyte and so on) as outcomes in their own right.

## 5. Precision is not the limiting factor - and that is the point

Sampling standard errors at 30x run from 0.03% to 0.36% across the whole table. At
n = 127,000 a correlation of r = 0.015 - 0.02% of variance - reaches genome-wide
significance. So essentially the entire error budget is systematic, and the study is
powered to detect systematic error as biology. Depth, callers and estimator details
barely matter; the denominator, GC and batch decisions are the experiment.

## 6. Recipe

Two mosdepth passes over each CRAM: one in fixed bins, which yields the denominator,
the technical covariates and the GC curve at once, and one over the target intervals,
which yields the numerators.

```bash
# 1. targets, once per reference
python analysis/repeatcn.py targets --annot data --out data
#    -> data/repeat_cn_targets.bed, data/repeat_cn_targets.tsv

# 2. per sample: fixed 1 kb bins, and the target intervals
mosdepth -n -t 1 --by 1000 --fasta chm13v2.0.fa sample.by1000 sample.cram
mosdepth --by data/repeat_cn_targets.bed --no-per-base --fast-mode sample sample.cram

# 3. cohort: NGS-PCA over the binned coverage (https://github.com/jlanej/NGS-PCA).
#    Produces svd.pcs.txt AND autosomal.median.txt in one run.
apptainer run ngs-pca.sif -input mosdepth_dir/ -outputDir ngsPCA/ \
    -numPC 100 -threads 24 -bedExclude ngs_pca_exclude...bed.gz -iters 10

# 4. per sample: copy number, denominator taken from the high-quality median
python analysis/repeatcn.py estimate --regions sample.regions.bed.gz \
    --targets data/repeat_cn_targets.tsv \
    --auto-hq-median ngsPCA/autosomal.median.txt --sample sample.cram.by1000. \
    --out sample.cn.tsv

# 5. cohort: residualise each class on the coverage PCs
python analysis/repeatcn.py adjust --estimates rDNA45S.cohort.tsv \
    --pcs ngsPCA/svd.pcs.txt --n-pc 20 --out rDNA45S.adjusted.tsv
```

For the GC correction, pass `--control` as well: a two-column table of GC and depth per
window. The natural source is the bins NGS-PCA retained (`svd.bins.txt`), whose GC is a
static per-reference annotation computed once with `bedtools nuc`, and whose depths are
already in the `by1000` output. Without it the estimator still runs, but reports
uncorrected dosage, which is only safe for classes near the genomic mean GC.

The PCs, the median and the GC curve are all reusable across every class in the table
and across any other depth-derived phenotype - telomere length included - so the
marginal cost of adding a class is one interval set, not another pass over the data.

## 7. Before believing a number

1. **Sub-region concordance.** Estimate rDNA separately over 18S, 5.8S, 28S and the
   IGS. They are the same molecule at the same copy number, so disagreement measures
   bias. `concordance()` returns the CV; report it rather than averaging it away.
2. **The 5S negative control.** 5S rDNA is tandem, multi-copy and on a different
   chromosome with different GC. An association that tracks 5S as well as 45S is
   technical.
3. **Covariate regression.** Residualise on coverage PCs (section 4c), and additionally
   on depth, duplicate rate, insert size and GC-curve summary statistics. Anything that
   survives is a candidate; anything that does not was never biology. Report how much
   estimate variance the PCs removed - if it is large, say so, because it bounds how
   much technical signal the raw estimate carried.
   Two diagnostics are worth running alongside: whether an association strengthens or
   weakens as PCs are added (a real one should plateau, a technical one should decay),
   and whether the leading PCs themselves associate with the phenotype.
4. **Replicate concordance.** Technical duplicates where available. In a pedigree,
   parent-offspring pairs are better still: acrocentric arms are transmitted intact
   (report section 23), so a parent and child share array lengths exactly, and any
   discordance beyond the de novo rate is measurement error.
5. **Denominator sanity.** Check `N_BINS` in `autosomal.median.txt` is what you
   expect and constant across samples, and that no retained sample has
   `AUTO_HQ_median` of zero. A sample whose bin count differs has had different
   sequence excluded, and its denominator is not comparable.
6. **Calibrate against CHM13**, whose true values for every class in the table are
   known by construction.

## 8. Validating on trios

The 602 complete trios in the 1000 Genomes 30x cohort are not just a sanity check. They
measure the one thing that is otherwise unmeasurable without an orthogonal assay: **how
much of your estimate is real.**

### The slope is the reliability

Array copy number is a physical DNA quantity, inherited additively with no dominance and
no environmental component: a child's diploid CN is one paternal haplotype plus one
maternal haplotype. So

```
E[T_child | parents] = (T_father + T_mother) / 2
```

holds *exactly*, with a coefficient of exactly 1. This is a stronger situation than
ordinary quantitative genetics, where the midparent slope estimates heritability and a
slope below 1 is biology. Here a slope below 1 can only be measurement error, and the
regression of the observed child value on the observed midparent value has expectation

```
slope = var_true / (var_true + var_err) = reliability
```

With the observed parental variance `V`, that single regression gives the whole
decomposition: `var_true = R * V`, `var_err = (1 - R) * V`. `trio_reliability()` returns
these as true and error CVs.

### What 602 trios can resolve

The standard error of the slope runs from 0.041 to 0.057 across plausible reliabilities,
so the smallest departure from a perfect assay detectable at two standard errors is
**R = 0.918** - about an 8% error CV for a class whose biological CV is 25%, like rDNA.
That is comfortably better than the precision any of these estimates needs.

| true R | implied error CV | SE(slope) | departure from 1 |
| --- | --- | --- | --- |
| 0.98 | 3.6% | 0.042 | 0.5 SE |
| 0.95 | 5.7% | 0.043 | 1.2 SE |
| 0.90 | 8.3% | 0.045 | 2.2 SE |
| 0.75 | 14.4% | 0.049 | 5.1 SE |
| 0.50 | 25.0% | 0.054 | 9.3 SE |

### Run it per class, and it becomes the ripeness test

This is the most useful thing the trios do. Run the regression for every class in the
target table: a class whose slope is near 1 is being measured well, and a class whose
slope is near 0 is not being measured at all - its variance is noise. The "ripe" column
in section 3 is currently an argument from array size and sampling precision; the trios
turn it into a measurement. Run 5S alongside as the negative control, and expect a
*high* slope there too, since 5S copy number is as heritable as any other array; 5S is
a control for association artefacts, not for assay noise.

### The caveat that inverts the answer

Reliability is inflated by measurement error **shared within a family**. If a trio was
extracted, libraried and sequenced together, part of the error is common to all three,
it enters the covariance, and the slope rises. Working this through, with a fraction
*p* of the error variance shared:

```
slope = (var_true + 2*p*var_err) / (var_true + (1+p)*var_err)
```

At p = 0 this is the reliability; at p = 1 it is **exactly 1, however bad the assay
is**. For a genuine 10% error variance, 50% sharing reports R = 0.957 instead of 0.909,
and full sharing reports 1.000.

1000 Genomes trios were, as a rule, processed together. So the raw number is an upper
bound, and the informative comparison is before against after `adjust_for_pcs`: if
reliability falls once coverage PCs are removed, the raw figure was inflated by shared
batch, and the adjusted one is the honest estimate. `repeatcn.py trios --pcs` prints
both and flags the drop. A permuted-family null is printed alongside, and should sit at
zero; if it does not, the estimate carries structure shared across unrelated samples.

```bash
python analysis/repeatcn.py trios --estimates rDNA45S.cohort.tsv \
    --pedigree 1kGP.3202_samples.pedigree --pcs ngsPCA/svd.pcs.txt --n-pc 20
```

Mutation, incidentally, does not bias the slope - it is mean-zero and independent of
the parents, so it adds variance to the child only and appears as residual in excess of
the Mendelian prediction. That excess is reported, and for rDNA it is a quantity worth
looking at in its own right.

### What trios cannot give you

Absolute calibration. The trios validate precision and internal consistency; they say
nothing about whether your 221-copy estimate is really 221 copies. For scale you need an
orthogonal measurement: CHM13 itself, where the true value of every class in the table
is known by construction, and the overlap between this cohort and assembled samples -
many HPRC genomes are 1000 Genomes samples, so their assemblies give per-class masses
directly. Note that rDNA specifically collapses in most assemblies, so assembly-derived
rDNA counts are a weak truth; for rDNA, CHM13 and ddPCR are the calibrators.

Two further things this cohort offers: 26 populations, so ancestry differences in array
size are measurable but also need genetic PCs included separately from coverage PCs; and
NGS-PCA already ships a worked 1000G high-coverage example, so the PCs and
`autosomal.median.txt` for this exact cohort may not need recomputing.

## 9. When to use a k-mer backend instead

Use depth over reference intervals when the class has trustworthy annotation in the
reference you aligned to. Use k-mer counting when it does not - a class absent from the
reference, a class whose arrays are collapsed in the assembly, or when you want to avoid
mapping bias entirely. Build a class-diagnostic canonical k-mer set (present in the
target consensus, rare in the rest of the genome), count it in the reads with
`meryl`, `KMC` or `jellyfish`, and apply the section 2 formulae to the counts with a
control k-mer set counted the same way. The k-mer route is the honest option for
HSat2/HSat3, where reference array lengths are one person's and mapping is unreliable.
