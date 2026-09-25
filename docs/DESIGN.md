# NGS-DOSE design

How much of a multi-copy sequence class does a person carry — 45S rDNA units, 5S units, the
acrocentric distal junction, satellite arrays — measured from the short-read WGS that biobanks
already hold. This document is the specification of the method as implemented, and the record of
the measurements on real data that decided each part of it. Where a statement was measured here
it says on what; where it is still a plan it says so.

Contents: [1 problem](#1-the-problem) · [2 what real data showed](#2-what-real-data-showed) ·
[3 estimator](#3-the-estimator) · [4 counting engine](#4-the-counting-engine) ·
[5 panels](#5-k-mer-panels) · [6 controls and the GC model](#6-controls-and-the-fragment-gc-model) ·
[7 window calibration](#7-window-efficiencies-and-the-anchor) · [8 known-truth controls](#8-known-truth-controls) ·
[9 cohort layer](#9-the-cohort-layer) · [10 transmission](#10-transmission-reliability) ·
[11 performance](#11-performance) · [12 limits](#12-what-it-does-not-do-and-known-limits) ·
[13 next](#13-what-comes-next) · [14 prior work](#14-relation-to-prior-work) · [15 audit](#15-assumption-audit-before-the-cohort-run)

## 1. The problem

A multi-copy class has no address in a linear reference. GRCh38 holds about eight copies of the
45S transcribed region (two loci on chr21p, `chrUn_GL000220v1`, `chr22_KI270733v1_random`), 17
of the ~200 5S units, and fragments of everything else; an individual's 300–600 rDNA units are
aligned onto those few copies at MAPQ 0, in whatever proportions the aligner's tie-breaking
produces. Per-base genotypes there are meaningless. *Dosage* is not: every read of the class is
still in the file, and counting them against single-copy sequence gives the amount of the class
in the genome.

That measurement is wanted at biobank scale — rDNA copy number has reported associations with
blood-cell traits, renal function and body mass in UK Biobank — and it is precise: sampling
error at 30× is a fraction of a percent. The published estimators are correspondingly simple
(reads in the reference's 18S copies over autosomal reads; depth over a unit consensus over
chr1 depth). What limits them is systematic error that differs between samples: library GC
behaviour, the duplicate flag, which reference copies the reads landed on, sequencing provider.
At n = 10⁵ a systematic effect explaining 0.02% of variance is genome-wide significant, so the
error budget that matters is entirely systematic. NGS-DOSE is built around removing those terms
one at a time, and around controls whose true copy number is known so that the removal can be
checked rather than asserted.

Measured: total diploid dosage of a class. Not measured: which chromosome the copies are on,
array structure, unit sequence variants, or activity.

## 2. What real data showed

The first implementation of this project was a specification (depth over annotated intervals
via mosdepth, NGS-PCA's median as denominator, GC curve from binned depth). Before rewriting it,
the assumptions were tested on 1000 Genomes 30× CRAMs (NYGC; TruSeq PCR-free, NovaSeq 2×150,
bwa-mem 0.7.15 to the GRCh38 analysis set). Most of them did not survive, and one bug was only
found because the rewrite is tested against simulated truth. Findings 1-5 were measured on
NA12878 unless stated; finding 6 on the twelve-sample pilot
([report](https://github.com/jlanej/NGS-DOSE-1000G/blob/main/pilot/pilot_report.md), reproducible with
[`run_pilot.sh`](https://github.com/jlanej/NGS-DOSE-1000G/blob/main/pilot/run_pilot.sh); both in NGS-DOSE-1000G).

**Finding 1 — the duplicate flag is not neutral.** MarkDuplicates keys on the mapped position of
both mates. Inside the collapsed rDNA, mates of true duplicates scatter across paralogs and
escape the flag. Flagged fraction: 10.8% in single-copy controls, 3.6% in the GL000220 rDNA
copy, 5.7% over all rDNA reads. Any tool that drops flagged reads — `samtools depth`, mosdepth
and therefore NGS-PCA's `AUTO_HQ_median` — deflates the denominator more than the numerator, by
a factor (1 − d_rDNA)/(1 − d_ctrl) that is specific to the library (duplicate rates in this
cohort run from 6% to 11%). Applying that factor to our uncorrected 18S estimate reproduces the
published value for the same CRAM (Hall et al. 2021: 301; ours 276.5 × 1.057 = 292), and
accounts for five of the eight percentage points by which their values exceed ours across the
shared pilot samples ([pilot report](https://github.com/jlanej/NGS-DOSE-1000G/blob/main/pilot/pilot_report.md) §4).
*Decision: count every primary read in numerator and denominator alike, and compute the
denominator in the same pass under the same rules.* NGS-PCA's median remains a cross-check
(32.74 dup-excluded against our 37.38 dup-included for NA12878, the expected ratio).

**Finding 2 — coverage follows the GC of the fragment, not of the read.** Fragment 5′ ends were
counted in single-copy regions and their density predicted from the GC of a window of length L
downstream; L was scored by Poisson deviance explained on held-out regions:

| L (bp) | 50 | 100 | 150 | 200 | 300 | 400 | 434 (insert median) | 500 | 600 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| deviance explained | 0.067 | 0.177 | 0.208 | 0.222 | 0.235 | 0.241 | 0.242 | 0.242 | 0.240 |

This re-derives Benjamini & Speed (2012) on current PCR-free NovaSeq data: the read-length
window leaves a seventh of the explainable signal on the table, and the optimum is flat around
the insert size. *Decision: count fragment ends, model their rate against fragment-scale GC,
take L from the sample's own insert-size distribution.* Counting ends rather than depth also
makes the estimator independent of read length.

**Finding 3 — this library over-represents GC-rich fragments.** The fitted rate relative to the
genome mean is 0.96 at 35% GC, 1.17 at 60%, 1.28 at 70%, and falls again above 78%. The shape
differs between samples (rate at 65% GC: 1.19–1.35 across the twelve pilot samples). A within-region
(fixed-effects) fit gives the same curve as the pooled fit, so this is a fragment-level effect
and not replication timing leaking in through isochores. A second, read-scale GC term
conditional on fragment GC moves rates by at most ±5% inside the range single-copy DNA can
calibrate, so it was left out.

**Finding 4 — rDNA's GC-rich segments drop out far beyond anything the genome-wide curve
predicts.** Along one unit, after the fragment-GC model: 18S (56% GC) 484, moderate-GC
intergenic spacer 505, 28S (69%) 374, ITS2 (83%) 355, with 250-bp windows down to half the 18S
level. Every unit carries each segment once, so this is technical, and single-copy sequence of
the same fragment GC is *over*-represented in the same library. It is a property of these
particular sequences (extreme local GC, stable secondary structure). The 28S/18S ratio is 0.77
± 0.05 across the 2,419 samples of Hall et al. — the dropout is universal and its size varies by
sample. *Decision: no genome-wide model can reconcile the windows; estimate a per-window
efficiency from the cohort and pin the scale on windows that can be trusted (section 7).*

**Finding 5 — placement is not something to rely on, but it is learnable.** Classifying every
read of the CRAM by k-mers, wherever it was aligned, shows where class reads actually go: 99.94%
of 45S reads sit in 17 intervals (the four rDNA loci, a genuine rDNA fragment at chr21:8.99 Mb,
and a dozen GC-rich look-alike loci such as chr2:32.91 Mb that attract 28S reads), 99.98% of 5S
reads in two, 99.8% of distal-junction reads in 56. A second scan, of a male of different
ancestry (HG02258, ACB), finds the same: those intervals hold 99.94 / 99.98 / 99.74% of his
reads, and the shipped table is the union of the two (80 intervals, 3.3 Mb). *Decision: two modes. `scan` reads
everything and is placement-independent; `fetch` retrieves only the learned sinks plus the
controls, and gives the same counts up to that capture fraction in a fraction of the I/O.*

**Finding 6 — the dropout is a property of sequence × chemistry, and that decides where the
scale may be anchored.** Every pilot sample also exists as an independent, years-older library
(HGSVC, nine samples: HiSeq 2500, 2×126, inserts ~570 bp, ~78×; Illumina Platinum, the CEU
trio: HiSeq 2000, 2×101, inserts ~300 bp, ~54×). Their GC behaviour is the reverse of NYGC's —
rate at 65% GC 0.48–0.56 and 0.84–0.88 against 1.19–1.35 — so a plain depth ratio differs
between libraries of the same person by 27–35% on average. Three observations:

- *The GC model transfers where the sequence is unremarkable.* Held-out autosomal sequence reads
  2.00 in both generations, the distal junction 9.7 and 9.9, and on a subset of 45S windows the
  two libraries agree to within a few percent although their GC corrections differ twofold.
- *Elsewhere they disagree deterministically.* The per-window log ratio between the libraries,
  each under its own GC model, runs from −1.4 to +0.3 over the moderate-GC windows, yet for a
  given window it is the same in every individual (10–90% spread across twelve pairs 0.15). At
  100-bp resolution the deficits are strand-asymmetric and one fragment wide — forward ends lost
  upstream of an element, reverse ends downstream (unit 29.4–30.7 kb is the clearest) — the
  signature of fragments that fail because they *contain* something, which mean fragment GC
  does not see. Which elements fail depends on the chemistry: on the HiSeq libraries most of
  the Alu- and microsatellite-rich spacer is inside such a zone, on NovaSeq the GC-extreme
  ETS/ITS/28S segments are. The full Benjamini–Speed model λ(fragment length, fragment GC) was
  prototyped to see whether the long, broad HiSeq inserts were the cause; it changed the
  estimates by under 1% and was not adopted. Nor is the effect regional: a within-region fit
  reproduces each library's pooled curve.
- *So an a-priori rule cannot choose the anchor.* "Fragment GC 40–60%" selects mostly clean
  windows on NovaSeq and mostly affected ones on HiSeq 2500, and estimates anchored that way
  differ by 15% between the generations.

*Decision: window efficiencies are learned per library type, and the absolute level is set on
consensus windows on which different chemistries have been shown to agree (section 7).*

**Three hypotheses that data rejected** are worth keeping on record, because each would have
been a plausible design: that the control GC curve is replication timing in disguise (the
within-region fit says no, for both library generations); that reads escape to loci outside the
sinks under a different alignment pipeline (a probe of all 2,554 unplaced, random and decoy
contigs plus the unmapped bin of an HGSVC CRAM found 9 of 1.5 million 45S reads outside the
NYGC-learned sinks); and the fragment-length model above.

**A bug worth recording.** BAM stores a reverse-strand alignment reverse-complemented. The first
engine classified the stored sequence and therefore called every aligned read "forward", booking
each reverse read one read-length away from its true 5′ end. Totals were unaffected, which is
why no real-data check caught it; the simulated end-to-end test, which knows the true strand and
bin of every fragment end, failed immediately. That test is in CI (`tests/test_end_to_end.py`).

## 3. The estimator

Let λ(g) be the expected number of fragment 5′ ends per position-strand in diploid single-copy
DNA whose fragment-length window has GC content g. For a class present in C copies per diploid
genome, a unit position p on strand s receives (C/2)·λ(g_s(p)) ends, so for any window w of the
unit

```
C_w = 2 · obs_w / exp_w        exp_w = Σ over position-strands in w of λ(g_s(p))
```

with the fragment window taken downstream of a forward 5′ end, `[p, p+L)`, and upstream of a
reverse one, `[p−L+1, p]`; units in tandem arrays are treated as circular. A position-strand
enters both sums only if (i) a read starting there would hold at least 20 panel k-mers, so its
recovery is certain, and (ii) its g lies inside the supported range of the control curve.
Everything else — the 28S 3′ end that is shadowed by pseudogene copies, windows above 82% GC —
is masked from numerator and denominator alike rather than extrapolated.

Per class the single-sample report gives `cn_anchor` (ratio of sums over anchor windows, section
7; the headline when the unit has any), `cn_all`, the median window, every window, named features
(18S, 28S, …), and `*.flat` versions with no GC model so that the effect of the model is always
visible. For a class without a unit (satellite families, section 5) the same rate gives diploid
mass, `M = Σ_g T[g] / λ_R(g)`, with reads binned by their own GC and λ_R the read-scale curve.

The denominator is not a separate quantity: it is the level of λ. It is made robust by
comparing every control region's count with its expectation under the fitted curve and dropping
regions (CNV) and whole chromosomes (aneuploidy, common in cell lines) that depart by more than
30% and 4% respectively — *and* by more than 5 standard errors, because at low depth noise alone
crosses any fixed threshold and trimming noise is biased. (That second condition was added when
a 0.75× test fixture read 3% low.)

## 4. The counting engine

`ngs-dose count` (Rust, htslib) turns one BAM/CRAM into a counts file of 70 kB (fetch) to 240 kB (scan
with the experimental panels); nothing downstream
touches the alignment again, so models can be revised without re-reading a biobank.

- **Rules**, identical for controls and classes: primary records only; duplicate flag ignored;
  no MAPQ filter (class reads are MAPQ 0 by construction; 0.2% of control reads are below 20);
  the unit of counting is the fragment 5′ end with soft clips restored.
- **Controls**: ends are tallied per position-strand over the control regions; at the end of the
  pass the insert-size median is known and the engine tabulates N[g] (position-strands) and O[g]
  (ends) in 1% GC bins for a grid of window lengths (100–600 bp and the modal read length),
  plus per-region totals. Control sequences travel with the bundle (region ± 1 kb), so no
  reference FASTA is needed for the GC tables, and the Python layer recomputes N[g]
  independently — the two implementations are required to agree exactly in the tests.
- **Classes**: every canonical 31-mer of every read is tested against the panel — first against a
  one-hash bitset (~1.5% false positives, a few nanoseconds), then, if set, against the hash
  map. Classification is therefore exact; an earlier strided screen was not (section 15). A class
  is judged on its own hits (at least 4 panel k-mers), so its counts do not depend on which
  other classes are loaded: positional classes are independent of one another and of the
  families, and only compositional families compete, winner takes all, a read split more evenly
  than 5:1 being left unassigned. A read of a positional class is placed on the unit by the
  median diagonal of its hits, oriented by majority, corrected for the BAM strand convention,
  and booked in a 50-bp bin per strand. Where class reads were *aligned* is recorded too, on
  the 1-kb grid of `mosdepth --by 1000` (10 kb for the satellite families), together with a
  census of each such bin: every mapped primary read in it, and how many of those carry the
  duplicate flag. In scan mode that histogram is what sinks are learned from; the census is
  what says how much of a bin's depth is the class, which is what any depth-based shortcut
  would have to rest on (section 13). It costs nothing measurable: a scan of NA12878 takes
  104 s with it and 105 s without, and the counts file grows from 143 to 230 kB. Several panels can be
  loaded together (`-p` repeated); a k-mer claimed by two panels is dropped from both.
- **Provenance**: every counts file records the SHA-256 of the panel(s), controls and sinks it
  was made with, and `ngsdose estimate` warns when a cohort mixes them. Regions are matched to
  the bundle by name: counts whose *control* regions are not exactly the bundle's are refused
  (they are the denominator and the GC curve), while known-truth or dosage regions that a
  later bundle revision added are simply missing for counts made before it — so extending the
  truth sets never sends anyone back to the CRAMs. Output is byte-reproducible (same input,
  any thread count, macOS or Linux).
- **Input integrity**: before reading anything the engine asks htslib whether the file ends in
  its end-of-file marker, locally or over HTTPS, and refuses one that does not: a truncated
  BAM or CRAM decodes without error, and what it has lost is whatever sorted last — chr21, chr22
  and the unplaced contigs, that is, the rDNA. The lengths of the contigs it used go into the
  counts file, and `ngsdose estimate` refuses a file aligned to a build other than the
  bundle's (hg19 has the same contig names; the controls and sinks are coordinates). Output is
  written under a temporary name and renamed, so a killed job leaves nothing a resumed run
  would take for a result. In scan mode the mapped primary reads of every contig are recorded
  as well — a cheap cross-check on the region-based dosages (it is what exposed the `N` in
  chrM; section 15).
- **Modes**: `scan` streams the file (htslib decoding threads, a reader thread doing control
  bookkeeping, a worker pool classifying batches of packed sequence). `fetch` takes the merged
  union of padded control regions and class sinks, hands intervals to worker threads that each
  hold their own indexed reader, and counts a record in the interval containing its start so
  nothing is counted twice. An interval's counts are committed only once it has been read to its
  end, and a failed interval is retried with a fresh connection, so transient network errors
  neither lose nor double-count reads. CRAM decoding is restricted to the fields used (no
  qualities, names or tags). Remote `https://` inputs work through htslib; the engine points
  the bundled libcurl at the system CA store. A dead connection does not fail, it waits, and
  the thread that would retry is the one that is blocked: a watchdog therefore ends the process
  with exit status 75 when nothing has been read for `--stall-timeout` seconds (default 300),
  and the cohort script tries the sample again (three attempts, then the scheduler's `timeout`
  as a backstop).

## 5. k-mer panels

`ngs-dose panel` builds the panel from class FASTAs and background genomes.

A **positional** class has a unit consensus; a retained k-mer occurs exactly once in the unit, so
a hit is a coordinate. A **compositional** class is a family without a stable unit; a k-mer only
says "this read is class X". In both cases a k-mer is kept only if it never occurs in the
background genomes — GRCh38 (analysis set, with decoys) and T2T-CHM13v2.0 — outside the
intervals where the class legitimately lives (`resources/GRCh38/build_inputs/*.class_loci.bed`).

| class | unit | k-mers kept | notes |
| --- | --- | --- | --- |
| `rDNA45S` | KY962518.1, 44,838 bp, circular | 31,827 of 41,476 | Alu and simple-repeat spacer k-mers removed by the background filter; 99.5% of 18S and 93% of 28S read positions remain recoverable at 150 bp |
| `rDNA5S` | X12811.1, 2,231 bp, circular | 1,945 of 2,217 | 68% GC throughout, so no anchor windows: the estimate rests on the GC model |
| `DJ` | CHM13 chr21:2,708,299–3,108,298 (400 kb distal to the rDNA array) | 169,808 | restricted to the *core*: k-mers present exactly once in each of the five CHM13 acrocentric distal junctions and nowhere else |

The filter is deliberately strict. The genome holds ~30 dispersed rDNA-derived fragments
(200–700 bp, 91–99.8% identical; chr1:91.39 Mb, chrX:109.05 Mb, chr12 in CHM13, …). They are
*not* exempted, so the k-mers they share with the unit are dropped and the corresponding stretch
of the unit goes blind (1,348 of 44,838 positions at 150 bp, chiefly the 28S 3′ end). The
estimator masks those positions exactly, which is preferable to counting reads that may come
from elsewhere. One locus *is* exempted: chr21:8,986,604–8,988,749 is a 99.6%-identical 5′ETS/18S
piece of genuine rDNA that GRCh38 happens to place outside the annotated copies; without the
exemption the 18S keeps a third of its k-mers.

Dispersed sequence — the satellite families — is compositional, and it is the one kind of class
that only a whole-file scan can measure: its reads are scattered across an alignment (3–54% of a
satellite family on decoy and unplaced contigs), so there are no sinks to fetch. A positional
class can always be added later and fetched from the few places its reads land; a dispersed one
cannot, which is why a cohort that is scanned once is scanned with the **experimental** panels
under `resources/experimental/` loaded. The telomeric repeat turned out not to be dispersed: in
372 NYGC scans the aligner put 92% of its reads within 25 kb of a chromosome end and 60% into one
10-kb bin of chr5p, essentially none on unplaced or decoy contigs, and the rest at a fixed set of
interstitial loci (chr4 at 94 kb from the 4q end, chr18 at 113 kb from the 18q end, chr2:32.9 Mb,
…). Sinks learned from 40 scans captured ≥ 99.87% of it in the other 332, so the bundle carries
them (`ngsdose sinks --classes TEL`) and fetch mode measures the class when its panel is loaded.

- *Ten satellite families* from the CHM13 CenSat annotation (HSat1A, HSat1B, HSat2, HSat3,
  β-satellite, α-satellite HORs, and four smaller families of the acrocentric short arms and
  pericentromeres: ACRO1 composites, SST1, CER, SATR; 1.13 M k-mers that recur at least ten times
  in a family's CHM13 arrays and occur in no other family and nowhere outside annotated
  satellite). What a panel can see is measured by its *recall*: the share of 150-bp reads from the
  family's own CHM13 arrays that carry the four k-mers a read needs. HSat1A/1B/2/3 and the HORs:
  97–99.9%. ACRO 90%, β-satellite 69%, SST1 61%, SATR 50%, CER 42% — relative measures, under-read
  by about their recall. Left out: gamma satellite (13%), divergent HORs (10%), HSat4 (four k-mers
  survive), and monomeric α, which as a class of its own takes 17% of the HORs' k-mers with it,
  because a k-mer shared between classes is dropped from both.
- *Against assemblies of the same people* (HPRC release 2, CenSat annotation of both haplotypes
  summed; two samples so far, section 15): HSat3 0.97 and 0.98 of the assembly, HSat1A 0.93 and
  0.97, HORs 0.97 and 1.03, HSat1B 0.90 in the male and 0.83 in the female (who has 2 Mb of it:
  the family is mostly on Yq); β-satellite 0.69 and 0.76, which is its recall; ACRO 0.79 in both,
  CER 0.42 and 0.43 (its recall again): under-read, but by a stable factor. **HSat2 cannot be judged in either**: both assemblies have gaps inside
  their HSat2 arrays (17 and 9 Mb of gap-containing arrays beside 44 and 30 Mb of spanned ones),
  and the estimates — 1.52 and 1.73 of the spanned arrays, 1.09 and 1.33 with the gapped ones
  counted at their annotated size — sit where a truth that is a lower bound leaves them. An
  assembly is a truth only for the arrays it spans; NGS-DOSE-1000G's `pipeline/hprc_satellites.py` keeps the
  two apart and leaves a sample out of a class's comparison when gap-containing arrays are more
  than 2% of what is annotated. SST1 and SATR are annotated
  several times more generously in the HPRC assemblies than in CHM13, so their absolute ratios
  mean nothing. Two hundred samples of the cohort have such assemblies, and whether the estimates
  *track* them across people is the question that matters (`04_hprc_satellites.sh`).
- *The telomeric repeat* (`TEL`: the six canonical 31-mers of (TTAGGG)n, unfiltered). A relative
  measure, not a telomere length: a read is assigned with 34 bp of perfect repeat, which
  interstitial telomeric sequence also has, and an exact 31-mer is lost to one sequencing error
  where TelSeq's hexamer count is not. The per-read share of k-mers that hit is recorded for every
  class (`hit_frac`), so the threshold can be chosen, and calibrated against TelSeq, afterwards.
  Fetchable through the bundle's sinks (above); the whole-file count and the count inside the
  chromosome-end windows an NGS-TL/TelSeq-style query would retrieve agree at r = 0.9997 across
  the 372 scans.

## 6. Controls and the fragment-GC model

`resources/build/select_controls.py` chooses the controls deterministically from the complement
of NGS-PCA's exclusion set (10x SV blacklist, 100-mer mappability < 1, DGV, segmental
duplications): 800 runs of 8–20 kb, 10.1 Mb, on all autosomes. The rare GC strata are filled
first — only ~0.05% of single-copy sequence is ≥70% GC at fragment scale, and the curve has to be
known where the classes are — then an even genomic spread. 80 further autosomal runs, 60 chrX
runs (outside the PARs and the X-transposed region) and 40 chrY runs (X-degenerate sequence;
candidates screened in a female, where 54 of 56 hold under 0.4× and the two that do not are
excluded by coordinate) are held out as known-truth regions. Two more regions are dosages rather
than truths: 1.5 kb of the mitochondrial genome (12S/16S rRNA — outside the D-loop, the common
deletion and the stretch copied into the chr1 NUMT, and short of the `N` that GRCh38's chrM
carries at 3,107; section 15) and 20 kb of the EBV decoy (unique sequence
away from the W and terminal repeats and from the segment the B95-8 strain lacks). Only
`control` regions enter the GC curve, the denominator, the library statistics and the control
PCs; every other role is tallied and reported. A reference that lacks a non-control contig
(no EBV decoy, say) is accepted and the value is missing, not zero.

λ(g) is a Poisson regression of O on a natural cubic spline in g (knots every 5% GC) with offset
log N, fitted per sample. Its supported range is where the controls hold at least 2,000
position-strands per 1% GC bin — 12–82% — and outside it nothing is estimated. The range is a
property of the bundle, deliberately not of the sample: the first version defined it by how
precisely the curve was known, which made the set of usable windows shrink with depth and the
all-window estimate drift with it (section 15). How well the curve is determined is reported
per sample instead (`gc_curve_max_se`). The headline
estimate is insensitive to the size of the control set — 508.0, 508.3 and 508.2 with 800, 400
and 200 regions (513 with 100) — because moderate-GC windows sit where the curve is best
determined; GC-rich features move by 3–4%, which is one reason they are not the headline.

## 7. Window efficiencies and the anchor

After the GC model, within one unit

```
log C_iw = c_i + a_w + e_iw
```

c_i is sample i's log copy number, a_w a window efficiency shared by all samples of one library
type (findings 4 and 6), e_iw what is specific to the sample. `ngsdose cohort` fits this by
median polish, so c_i draws on every window and is robust to a run of deviant ones; its
window-based relative standard error is below 1% in every pilot sample.

The model only identifies c_i + constant. The constant is fixed by requiring median a_w = 0
over the **anchor** windows, and finding 6 says they have to be chosen by evidence rather than
by rule. The bundle ships them (`anchors.json`): the moderate-GC windows on which the NYGC
library and an older library of the same person agree (|median log ratio| < 0.06, 10–90% spread
< 0.15 across twelve pairs and three chemistries) — 13 windows, 3.25 kb, in the 18S and at 22,
27, 31 and 35–38 kb of the spacer. Two chemistries that disagree by up to fourfold elsewhere are
unlikely to share a bias where they agree. Without `anchors.json`, and for classes that have
none (5S, DJ), the rule is fragment GC 40–60%, or every usable window if there are none.

What this buys, out of sample (anchors learned with each family left out in turn): the two
library generations agree with a mean log ratio of +0.02 and a pair-to-pair SD of 0.032
(r = 0.98), against −0.31 and 0.10 (r = 0.86) for the 18S depth ratio. Without any cohort —
one sample, ratio of sums over the anchor windows — +0.006 and 0.046.

What it does not buy: a library chemistry that has not been compared with another has
unknown dropout zones. Its estimates are comparable within its batch (a_w absorbs the zones);
their absolute level is provisional. The residuals e_iw are kept: their principal components
are sample-level covariates orthogonal to c_i by construction, and a sample with a coherent run
of deviant windows is a candidate structural variant of the unit rather than a bad measurement.
Clean spacer segments differ reproducibly from one another (roughly 0.9 to 1.05 of the anchor
level, in both chemistries), which suggests that spacer segments are not all present exactly
once per unit; the gene count is what the 18S measures, and the anchors include an 18S window.

## 8. Known-truth controls

Four sequences of known copy number are measured in every sample by exactly the code paths
used for the classes, so accuracy is observed per sample rather than assumed:

- **held-out autosomal regions** (truth 2) — the alignment-position path and the GC model;
- **chrX** (truth 1 or 2) — the same, across a twofold dosage step;
- **chrY** (truth 1 or 0) — the same at the bottom of the scale, where a floor of mis-placed
  reads would show; with chrX it is also the sex check and the mosaic-loss screen;
- **the distal junction** (truth 10: one per acrocentric short arm) — the k-mer path on a
  multi-copy, paralogous, acrocentric sequence, i.e. the path rDNA takes. A Robertsonian
  translocation carrier should read 8; that is a prediction, not yet tested on a known carrier.

Pilot, twelve NYGC samples: held-out autosomal 1.996 ± 0.006; chrX 0.995 ± 0.004 in males and
1.91–1.98 in females; chrY 0.979 ± 0.004 in males and 0.002 in females; DJ 9.73 ± 0.13 (range
9.50–9.93). In the twelve older libraries: 2.007 ± 0.006, 1.002, 1.94, chrY 0.979 and 0.002,
and 9.93 ± 0.17. (chrY reads 2% low in every male and in both libraries alike, so it is a
property of those 40 regions on this reference — X-homologous reads placed on the X, or
sequence the reference Y carries and these men do not — not of the culture; what matters for
its job is that it is the same 0.98 everywhere and 0.002 where there is no Y.) Children's DJ
sits within 0.1 of the midparent in all four trios. Two things follow. The controls catch what they should: HG00732 reads 1.61 X copies in
the 2019 DNA and 1.84 in the 2015 DNA — a culture losing an X. And the NYGC libraries read both
late-replicating controls (the inactive X, the acrocentric DJ) about 3% low where the older
libraries do not, which is what a larger S-phase fraction in the source culture would do; that
is a hypothesis for the cohort, where coverage PCs and 1,600 female X values can test it (a
first look, the leading component of the control regions' residual depth against DJ and female
X in twelve samples, has the expected sign and no power).

The two dosage regions say how different the two DNA batches of a cell line are. Between the
NYGC culture and the older one, mitochondrial genomes per cell differ by a factor of 1.22 (SD of
the log ratio; 495–2,170 copies overall) and EBV episomes by 1.88 (30–440 copies) — against
1.038 for calibrated 45S. That is the scale of culture-to-culture biology that the
cross-library comparison of section 7 has inside it, which is why it is an upper bound on
measurement error and not an estimate of it. Across the twelve pairs the 45S difference does
not follow either (r = +0.34 with mitochondrial content, +0.06 with EBV load; n = 12 decides
nothing); the cohort's trios can ask the question properly (section 10).

The DJ is also measured *along* its 400 kb, and that turned out to be the sharpest check of all.
In 10-kb windows each sample is flat to 3.5% after a shared profile is removed — except for
steps of one copy in ten that are inherited: HG00512 and his daughter HG00514 both read 0.77–0.81
over the distal 20 kb (about 8 of 10 copies), NA12891 and his daughter NA12878 both read +1 copy
there and −1 near 380 kb, HG00731 and his daughter HG00733 both read −1 near 130 kb. Single-copy
changes in a ten-copy paralogous sequence, seen in parent and child independently: the k-mer
path resolves what it claims to. At cohort scale this becomes a Mendelian test on integer
states (section 13).

## 9. The cohort layer

- `ngsdose estimate` — counts → per-sample estimates and the summary table.
- `ngsdose cohort` — window calibration; writes `CLASS.cn`, its robust relative SE, the profile
  roughness and profile PCs, and the efficiency table.
- `ngsdose adjust` — residualises estimates on coverage PCs from
  [NGS-PCA](https://github.com/jlanej/NGS-PCA). The PC basis is built from autosomal bins with
  repeats, segmental duplications and low-mappability sequence excluded, so it is disjoint from
  every class measured here: it can absorb library and batch structure but cannot absorb the
  dosage itself. It *can* absorb ancestry, so an association analysis should carry genetic PCs
  separately, and the variance removed by coverage PCs should be reported, not hidden.
- `ngsdose cohort` also reports **control PCs**: principal components of the 800 control regions'
  log(observed/expected) after the GC model — a small internal version of the same idea, for
  cohorts on which NGS-PCA has not been run; `ngsdose adjust` uses them when no `--pcs` file
  is given. Like the coverage PCs they are computed on sequence disjoint from every class.
- **How many PCs.** By default, those that clear the edge of the noise bulk of the SVD they came
  from (`--n-pc mp`; any number overrides it). The textbook way to find that edge — fit the
  Marchenko–Pastur law to the spectrum, for instance by matching its median — assumes every
  entry of the matrix has the same noise variance, and coverage does not: noise falls with a
  sample's depth and varies with a bin's mappability. On simulated noise whose rows differ in SD
  by ±40% the textbook fit called 34 components where 7 were planted and 148 where 5 were, and
  on NGS-PCA's spectrum of this cohort (the top 200 singular values of a 3,200 × 142,070
  matrix) it gives 59, 66 or 80 depending on whether 100, 150 or 200 values are kept. What does
  survive unequal variances is the *shape* of the edge: the density of any such bulk vanishes
  like a square root at its top, so the j-th largest noise value lies at E − a·j^(2/3). The edge
  E is therefore fitted, with a, to the lower half of the leading singular values, iterated on
  the number of components set aside as signal; a component counts if it clears E by a margin
  of 1% (plus four residual SDs of the fit and the Tracy–Widom scale, which matters for dozens
  of samples and not for thousands). With noise that differs between rows by ±40% and between
  columns by ±30% that recovers planted components exactly and finds nothing in noise alone;
  on the 1000 Genomes spectrum it gives 38–40 components however much of the spectrum is kept
  (45–49 without the margin; the fitted edge is 1.5 times as broad as equal-variance noise
  would make it). `ngsdose cohort` makes the same choice for its control PCs, records it
  (`ctrlPC_mp`), and writes twice as many. Two limits, both measured:
  - *Bins with heavy-tailed noise variance make the count lean high.* A bin whose variance is
    several times the typical one is a component of its own — real structure in the matrix, of
    no interest. With log-normal bin SDs (σ = 0.3; 1,500 × 12,000) noise alone yields 2.8
    components per run without the margin and 1.1 with it; at σ = 0.5, 24 and 19. The margin is
    a palliative that costs nothing that was measured (components planted 3–15% above the
    largest noise value are all found at margins up to 2%), not a cure; fitting nearer the edge
    collapses when many components sit at the threshold, and a curvature term is unstable.
  - *The count is a property of the matrix; whether components that deep are usable is a
    property of the solver and the sample set.* Two NGS-PCA runs of this cohort (3,200 and
    3,202 samples) share their leading 20 PCs exactly (largest principal angle 1°) and then
    part ways: 6° at 25, 18° at 30, 40° at 46; PC 40 of one run lies 86% inside the other's
    leading subspace, PC 46 80% (NGS-PCA's own seed control shows the same). Within one run
    such components are a legitimate variance sink at n = 3,200; they are not portable between
    runs, and an analysis that has to be reproduced from a fresh SVD should not lean on them.
  So the default is principled and on the high side, and what gets reported should rest on the
  sweep.
- **Whether that number is right** is an empirical question, and the cohort carries what it
  takes to answer it: sequence of known copy number in every sample, and trios.
  `ngsdose pcsweep` regresses out 0, 1, 2, … PCs and reports, for each number, the error of the
  known truths (held-out autosomal, chrX and chrY by sex, the distal junction) and the
  transmission reliability of the classes with its paired difference from no adjustment.
  Everything is cross-validated — residual variance falls with every regressor whether it
  means anything or not: forty random PCs "explain" a third of pure noise in 120 samples, and
  out of fold make it 10% worse — and the recommendation follows the one-standard-error rule:
  the fewest PCs that do as well as the best number, to within the sampling error of "best".
  A known truth says when adjustment has stopped removing noise; reliability says when it has
  started removing signal. They need not agree with the spectrum, and where they do not, they
  win.
- Two confounders the counts cannot remove and an analysis must carry: **replication timing**
  in DNA from cycling cells — the inactive X reads 1.94 rather than 2 in female LCLs here, and
  late-replicating satellites and silent rDNA will be under-represented in the same way, by an
  amount that tracks the S-phase fraction of the culture; blood DNA is largely free of it — and
  **cell composition** for blood-derived traits.

## 10. Transmission reliability

Array dosage is inherited additively: E[T_child | parents] = (T_f + T_m)/2 exactly. With X = T
+ e, measurement error independent between individuals and ρ the observed spousal correlation,

```
midparent slope      b = R(1+ρ_T)/(1+Rρ_T)            ⇒  R = b − ρ(1 − b)
single-parent slope    = (R + ρ)/2                     ⇒  R = b_f + b_m − ρ
Mendelian variance   Var(X_c − X_mid) = σ_T²/2 + 1.5σ_e²  ⇒  R = 1.5 − D/V
```

where R = σ_T²/(σ_T² + σ_e²) is the reliability. The first version of this design used the raw
midparent slope on log values; the slope is inflated by anything that correlates spouses
(population structure in a pooled cohort), and inheritance is additive on the natural scale, not
the log scale. `ngsdose trios` therefore centres values within population, works in copies, and
reports all three estimators, the spousal correlation and a permuted-family null. All three are
inflated by error shared within a family — a trio libraried together — and the spousal
correlation after population centring is the direct test for that, since spouses share no
dosage by descent. Trios bound reliability from above; technical replicates measure it. Run per
estimator, the table answers which estimator to use: the uncorrected 18S ratio, `cn_anchor`, the
calibrated `cn`, each before and after coverage-PC adjustment (NGS-DOSE-1000G's `pipeline/02_cohort.sh`).
Confidence intervals come from resampling families, which carries the uncertainty of ρ into R;
with 602 simulated trios the 95% interval is about ±0.09. That is too wide to rank estimators
whose reliabilities differ by a few hundredths, but two estimators of the same quantity are
strongly correlated, so `ngsdose trios --compare-to` bootstraps the *difference* R_A − R_B on the
same resampled families (in simulation a gap of 0.08 is resolved at 602 trios with a 95%
interval of +0.03 to +0.14). The algebra is checked in `ngsdose selftest`: normal and skewed
(log-normal) array sizes, additive and multiplicative error, spousal correlation, and error
shared within families.

The Mendelian line also assumes the children are measured on their parents' scale. With s the
children's standard deviation over the parents', it is exactly R_mendel = R + 1 + ρ/2 − s², so it
adds to the slope only whether the children vary more or less than their parents. In the 1000
Genomes 30× release every child is among the 698 related genomes sequenced after the original
2,504, where nearly every parent is, so generation and batch go together: a batch that reads a
quantity on a scale k multiplies b, and R, by k. The cohort page therefore reports s, the
children's level against their parents', and R with the children first rescaled to their
parents' spread (b/s for b). A batch's scale moves R and leaves the rescaled R alone; new
variation arising in the children does the reverse; the two bracket the reliability. On the
first 135 trios they agree within a few hundredths for most metrics and part for HSat2 (s = 1.16,
slope above 1, which new variation cannot produce) and ACRO (s = 0.84, which it cannot either).
`trios.by_sex` asks whether the four parent–child pairings by sex differ, referring Cochran's Q
on Fisher's z to its distribution under shuffled children's sexes and swapped parental roles:
the chi-square reference is far too liberal for skewed quantities (the distal junction's
pairings: p < 0.001 by chi-square, 0.51 by the shuffles).

The cohort table carries its own controls for this analysis. `truth.auto` has no variance but
error, so its "reliability" should be nil. The EBV load and the mitochondrial content of the
culture (`chrEBV.copies`, `chrM.copies`; section 6) vary several-fold between cell lines and are
not transmitted through the nuclear genome: if they come out transmitted, members of a family
share something other than DNA — a batch, a culture protocol — and every reliability in the
table is inflated by about as much. The same two columns are the first handle on the variance
that is not transmitted *and* not measurement error: change of the array in culture. Earlier
work relates rDNA dosage to mitochondrial abundance (Gibbons et al. 2014); in LCL-derived data
that relation and a culture effect are confounded unless both are measured.

## 11. Performance

Apple M1 Max, NA12878 (15.8 GB CRAM, 758 M primary reads), bundle GRCh38-v1:

| mode | input | wall | CPU | data read |
| --- | --- | --- | --- | --- |
| scan, 10 threads | local CRAM | 1 min 40 s | 13.5 min | 15.8 GB |
| scan, 8 threads, with the experimental panels (1.1 M k-mers, 14 classes in all) | local CRAM | 1 min 40 s | 13.5 min | 15.8 GB |
| scan, the same | the same CRAM over HTTPS, home connection (~15 MB/s) | 15–18 min | 14 min | 15.8 GB, network-bound |
| fetch, 8 threads | local CRAM | 2–3 s | ~15 s | 0.48 GB (1,973 of 79,637 slices; chrM and chrEBV are 14 MB of it, chrY 5 MB) |
| fetch, 16 threads | `https://1000genomes.s3.amazonaws.com/…` from a home connection | 60–100 s | ~20 s | 0.48 GB |
| estimate | counts JSON | ~1 s | | 70–240 kB |

In fetch mode the 3,202-sample cohort is therefore about 55 hours of single-stream time and never
needs a CRAM on disk; scanned whole from staged CRAMs it is about 750 CPU-hours and 48 TB of
transfer, a few files at a time (NGS-DOSE-1000G's `pipeline/01_stage_and_dose.sh`). Peak memory is 0.6 GB
(fetch) and 1.3 GB (scan) at 8 threads.

## 12. What it does not do, and known limits

- **Absolute scale of rDNA** rests on the anchors: windows on which three Illumina chemistries
  agree after their own GC corrections. That is evidence, not proof, of being unbiased. The one
  orthogonal check so far is ddPCR on twelve lymphoblastoid lines (Potapova et al. 2025; the
  results repository's `assembly_rdna` study, `python -m report --ddpcr` in NGS-DOSE-1000G): NGS-DOSE reads about
  0.97× the assay. Relative dosage within a library type does not depend on it.
- **A new chemistry** (DNBSEQ, Element, Ultima, a future Illumina) has unknown dropout zones.
  Window efficiencies must be re-learned on it, and its absolute level is provisional until it
  has been compared with a known one on the same samples.
- **Cell-line DNA.** Late-replicating sequence is probably under-represented in DNA from cycling
  cultures (section 8), by an amount that differs between cultures; rDNA and satellites are
  mostly late-replicating. Blood DNA should be largely free of this. The cross-library
  agreement quoted here (3–5%) includes it, and real drift between cultures as well.
- **Units that are not the consensus.** Truncated or rearranged units contribute only the
  windows they contain; the window profile absorbs a population-average of that into a_w.
- **Sinks are aligner- and reference-specific.** The shipped table was learned from two NYGC
  CRAMs (bwa-mem 0.7.15, ALT-aware, GRCh38 analysis set; a CEU female and an ACB male). A probe of one bwakit-0.7.12 + postalt
  CRAM found nothing outside them on unplaced, decoy or unmapped reads, but DRAGEN, or a
  different decoy set, needs `scan` on a handful of samples and `ngsdose sinks`.
- **GRCh38 only**, chr-prefixed names. A CHM13 or GRCh37 bundle is a rebuild of the controls and
  sinks; panels are reference-independent.
- **Dispersed families are experimental** (section 5): six satellite families with a first
  comparison against assemblies, four smaller ones without, and a telomeric-repeat class that is a
  relative measure only (exact 31-mers are less tolerant of sequencing error than TelSeq's
  hexamer count; the per-read hit fractions are kept so that a threshold can be chosen later).
- **Low depth** costs precision, not accuracy: the anchor estimate is unbiased down to 0.4×
  (1.004 ± 0.007 of the full-depth value over 20 subsamples; section 15), and a 0.75× subsample
  of NA12878 returns 499 against 504 from the full data. The GC curve itself gets noisy
  (`gc_curve_max_se`), which the GC-extreme features feel first.
- **Four trios** are the only pedigree data analysed so far; section 10 is specification and
  simulation until the cohort run.

## 13. What comes next

1. All 3,202 samples of the 1000 Genomes 30× cohort (NGS-DOSE-1000G's `pipeline/`) **scanned whole**, from
   staged CRAMs, with the satellite panel loaded — and fetched as well, three seconds more on a
   local file (its `01b_dose_sample.sh`). The scan is the full-accuracy mode and the reference for
   everything cheaper: placement-independent counts, the satellite families against an assembly
   truth for the 200 samples that have HPRC assemblies (its `04_hprc_satellites.sh`), sinks
   evaluated and re-learned across 26 populations and both sexes, and fetch / scan for every
   sample (its `03_compare_modes.sh`). From the cohort: the per-estimator reliability table with
   paired comparisons and its negative-control traits (its `02_cohort.sh`), efficiencies at
   n = 3,202, adjustment by NGS-PCA's coverage PCs and by the internal control PCs, DJ and chrY
   as cohort-wide accuracy checks, comparison with Hall et al.'s table sample by sample.
   - *How far can it be simplified for a biobank?* Three levels, each judged against the scan
     in the same people, with transmission reliability (the paired bootstrap of section 10) as
     the arbiter rather than correlation with the full estimate, which shares its errors.
     (i) **Targeted fetch**, which exists: about half a gigabyte and a minute per genome, no
     file staged; the cohort says what it loses, sample by sample, and the sinks re-learned on
     the 1-kb grid say how much smaller the retrieval can be made. (ii) **A depth proxy from
     bins a cohort already has.** NGS-PCA's mosdepth run (1-kb bins, all contigs, no MAPQ
     filter, duplicate-flagged reads excluded) is kept for this cohort. In NA12878, 99.8% of
     the aligned 45S reads fall in 273 of those bins (chr21's rDNA models, KI270733, GL000220);
     the engine's census of those bins reproduces a duplicate-excluded depth computed
     independently to 3%. What a proxy gives up is known in kind — no fragment-GC model, no
     window calibration, no anchor, and a duplicate filter that removes 5.5% of rDNA reads but
     10.8% of single-copy reads in that sample, a different gap in every sample — and the scans,
     the kept mosdepth outputs and the census are what it takes to put numbers on each. Within
     one library type and one pipeline several of those terms are constants, which is the case
     for trying. (iii) **Less of everything**: fewer controls, smaller sinks, lower depth — all
     of which can be tried on the counts files alone, because they hold positions and not
     summaries. None of this needs an answer before the run; it needs the run to keep what the
     answers will be computed from, and it does.
   - *A Mendelian test on integer states.* DJ windows carry inherited ±1-copy steps (section 8).
     Across 602 trios, a step present in a child and in neither parent is either a de novo
     event or an error, and a step in a parent is transmitted half the time: that calibrates
     the k-mer path at single-copy resolution, far more sharply than any regression on totals.
   - *The S-phase hypothesis.* If late-replicating sequence is under-represented in DNA from
     cycling cultures, DJ, female X and the leading control PC should move together across the
     cohort, and adjustment should tighten DJ around 10.
   - *The running record.* NGS-DOSE-1000G's `python -m report` turns whatever counts exist into one page — what is
     measured and why, and the evidence that it works: the known truths in every sample, fetch
     against scan, the trios, the cell-line covariates, the satellites against assemblies, the
     coverage PCs and their sweep — with every number recomputed from the counts files and
     every table beside it. It is published as the run proceeds (NGS-DOSE-1000G's `pipeline/05_report.sh`;
     the results repository's GitHub Pages), partial results and flags included, so that what
     the cohort shows is on record at every stage and not only at the end.
   - *Culture or error?* Whether a child's departure from the midparent tracks the EBV load,
     the mitochondrial content or the leading control PCs of the culture it was sequenced from;
     and, with the same columns, the published relations between rDNA dosage, 5S dosage and
     mitochondrial abundance (section 14), inside one library type and against known truths.
2. Assemblies collapse rDNA, so for rDNA the external anchors remain CHM13 (ddPCR 409 ± 9) and
   HG002; HPRC assemblies are a truth for the satellite classes only (item 1).
3. Sinks for DRAGEN-aligned CRAMs (UK Biobank, All of Us).
4. An orthogonal rDNA calibration wider than the twelve ddPCR lines.

## 14. Relation to prior work

Gibbons et al. (2014, 2015) introduced rDNA dosage from WGS depth. Hall, Turner & Queitsch
(*Sci Rep* 2021) analysed the same 1000 Genomes CRAMs used here, showed how library and batch
drive such estimates, and published the per-sample table this project is checked against.
Rodriguez-Algarra, Evans & Rakyan (*Cell Genomics* 2024) measured an 18S read-count ratio in UK
Biobank and demonstrated phenotypic signal; `18S.flat` is our counterpart of their estimator.
They also report that its mean differs between the two UK Biobank sequencing centres and adjust
for centre. Raj et al. (medRxiv, 13 January 2026, doi:10.64898/2026.01.09.26343685, Calico)
estimated 45S and 5S copy number in 490,383 UK Biobank genomes and report associations with
metabolic disease, adiposity and blood traits; that much is checked against its abstract (2026-09-22).
Its method is described in a literature summary as a fragment-GC Poisson model on length-matched
control loci, which would make it the closest relative of section 6, arrived at independently.
We have not read the method in full, and its code (github.com/calico/rDNA_CN) was not public on
2026-09-22, so the method description is unverified and not cited. T2T (Nurk et al. 2022) and
CONKORD (Potapova et al.) estimate rDNA from 18S k-mers in reads against GC-matched windows.
Benjamini & Speed (*NAR* 2012) is the source of the fragment-GC model.

Long reads do not yet replace a short-read measurement. In T2T-CHM13 three of the five rDNA
arrays are model sequences, because ultra-long nanopore reads could not order their units (Nurk
et al. 2022); in 156 acrocentric short arms assembled from HiFi, ultra-long nanopore and Hi-C
data the array collapsed in every one (Lin et al., *Cell* 2026); per-array size and activity have
been measured in fifteen genomes with long reads and imaging (Potapova et al., *Cell Genomics*
2025). The HPRC release-2 assemblies of cohort members hold about a third of the rDNA their
measured copy number implies, in pieces of at most about a megabase (`python -m report --censat`,
`data/rdna_hprc.tsv`). The cohort page's section 1 sets this out with the references.

What is specific to NGS-DOSE: placement-independent read assignment with learned sinks, so that
the same counts come from a full scan or a one-minute targeted retrieval; a per-position
recoverability mask; the separation of a library GC curve from sequence-specific window
efficiencies with an explicit anchor; known-copy-number controls measured in every sample,
including a multi-copy acrocentric one; the duplicate-flag finding; transmission reliability with
its spousal-correlation correction as the criterion for choosing between estimators; and coverage
PCs from a repeat-free basis for batch. The transmission algebra is textbook quantitative
genetics applied to a trait with heritability one; no claim of novelty is made for it.

## 15. Assumption audit before the cohort run

Before committing 3,202 genomes to it, every assumption that the counts files depend on — the one
layer that cannot be revised without re-reading the CRAMs — was listed and tested, along with
the modelling assumptions that were cheap to test. What was wrong is recorded as plainly as what
held.

| assumption | test | result | consequence |
| --- | --- | --- | --- |
| The strided k-mer screen loses no class read ("lossless for min-hits ≥ stride", as the first engine's help text said) | full NA12878, stride 4 vs stride 1 | **false**: 0.04–0.06% of class reads lost (hits in short, broken runs) | replaced by an exact bitset prefilter, no slower; unit test with exactly such a read |
| A class's counts do not depend on what else is in the panel | bundle panel alone vs with the satellite panel | **false**: DJ lost 1.7% of its reads, because a read carrying β-satellite k-mers as well was called ambiguous | classes judged independently (section 4); identical counts now asserted in CI |
| Output is reproducible | same input, 1 vs 4 vs 8 threads | counts identical; placement list ordered by hash iteration on ties | total ordering; byte-identical output asserted in CI |
| Sinks learned from one CEU female transfer to other samples | whole-file scan of HG02258 (ACB, male) | hold: 99.94 / 99.98 / 99.74% of 45S / 5S / DJ reads inside them, against 99.94 / 99.98 / 99.77%; chrY takes 66 DJ reads. A third scan that the sinks never saw (HG01884, ACB, female): 99.96 / 99.99 / 99.75% | shipped sinks are the union of the first two scans (80 intervals); the cohort run re-evaluates them on 200 samples (`03_compare_modes.sh`) |
| … and to another alignment pipeline | probe of all 2,554 unplaced/decoy contigs and the unmapped bin of a bwakit 0.7.12 + postalt CRAM | hold: 9 of 1.5 M 45S reads outside | DRAGEN remains untested |
| The estimate does not depend on depth | NA12878 subsampled, 20 replicates per depth | anchor estimate unbiased: 1.009 ± 0.005 at 1.1×, 1.004 ± 0.007 at 0.37×. **All-window estimate biased**: +2.4% and +6.5%, because the GC support, defined by the curve's precision, narrowed with depth and GC-rich low-efficiency windows dropped out | support defined by control positions (section 6): 0.998 ± 0.004 and 0.996 ± 0.005; asserted in CI |
| The robust denominator is neutral when nothing is aberrant | same series | its scale factor was not exactly 1 at low depth | now the *relative* change from dropping flagged regions: exactly 1 when none is |
| Analysis choices do not move the headline | NA12878: min k-mers 5–60, window 100–1000 bp, L 300–600 bp | 45S anchor 503.6–505.2; DJ 9.75–9.87; held-out autosomal 2.001–2.002. Features move more (18S 457–489 with the k-mer threshold, 5S 202–217 with L) | none; features are diagnostics, not headlines |
| The DJ deficit (9.5–9.9 rather than 10) is not structural | DJ in 10-kb windows | uniform within a sample to 3.5%, plus inherited one-copy steps (section 8) | deficit is global, consistent with the S-phase hypothesis; steps become a cohort test |
| Trio estimators hold for skewed sizes and multiplicative error | simulation, log-normal sizes (CV 35%), 10% multiplicative error | recovered: 0.852 for a true 0.852 | in `ngsdose selftest` |
| 602 trios can rank estimators | simulation | separate CIs are ±0.09, too wide; the paired bootstrap resolves a gap of 0.08 | `ngsdose trios --compare-to` |
| A file that decodes is a whole file | fixture BAM and a 200 MB head of a real CRAM, end-of-file block removed | **false, silently**: both decode without error, and what a truncated coordinate-sorted GRCh38 file has lost is chr21, chr22 and the unplaced contigs - the rDNA | the engine checks the EOF marker before reading (local or HTTPS) and refuses; `--allow-truncated` overrides; recorded in every counts file; asserted in CI |
| The file is aligned to the bundle's build | - | nothing checked it: hg19 has the same contig names, and the controls and sinks are coordinates | contig lengths are recorded in the counts file and `estimate` refuses a mismatch; asserted in CI |
| A failed connection fails | the pilot re-count, left running overnight on a laptop that slept | **false**: two fetch processes sat for eight hours on connections that were never going to answer; per-interval retries never fire, because nothing errs | stall watchdog in the engine (exit 75 after `--stall-timeout` s without a record), retry loop in the cohort script; asserted in CI with a FIFO that nobody writes to |
| A killed job leaves nothing a resumed run would mistake for a result | - | the engine wrote its output in place | output written under a temporary name, synced, renamed |
| A class in the counts file reaches the estimate table | merged-panel scan through `ngsdose estimate` | **false**: classes absent from the bundle's own panel were skipped, so the satellite columns of a cohort scan would have been empty | compositional classes need nothing from the bundle and are always estimated; asserted in CI |
| The bottom of the scale is zero | 40 chrY regions in a female | hold: 0.003 copies in NA12878 (after excluding two candidate runs that attract X-derived reads) | chrY added as a known truth (section 8) |
| Every position of a region has a fragment-GC window | chrM region-based estimate against the scan's whole-contig read tally | **false for chrM**: GRCh38 carries an `N` at chrM:3,107; windows over it are in no GC table but reads starting there were in the observed count, and the first chrM estimate was 18% high (978 against 830) | region moved off the `N`; the bundle build refuses such regions; expected counts are scaled to all position-strands, so a custom bundle cannot repeat it; asserted in CI |
| Non-transmitted variance is measurement error | - | untestable without a handle on the state of the culture | chrM and chrEBV dosage in both modes (NA12878: 850 and 71 copies per cell), and as negative-control traits in the transmission table (section 10) |
| The cohort script works | a mock cohort: seeded subsamples of the NA12878 fixture named after real 1000 Genomes trio members (66 samples, 22 trios), through the `ngsdose` commands that the cohort scripts (now NGS-DOSE-1000G's `pipeline/`) call, with the real pedigree and NGS-PCA's real PC file | it ran - for the first time; the twelve-sample pilot is too small for control PCs, PC adjustment or the paired bootstrap. With one person sixty-six times over there is nothing to transmit, and every reliability came out at zero within its interval. It also showed `adjust` reporting that twenty PCs "explain" 30% of pure noise at n = 66 | a sixty-sample version is a CI test of the whole cohort layer (`tests/test_cohort_flow.py`); `adjust` reports the chance expectation k/(n-1) and the adjusted R² beside the raw figure |
| Counts survive a bundle revision | pilot counts made before the chrY and dosage regions existed, estimated with the new bundle | refused: regions were matched to the bundle row for row | matched by name: the control set must be identical, other regions may be missing (section 4); asserted in CI |
| The manifest is right | every URL of `00_setup.sh`'s manifest requested (HEAD) | three sample names in the 1000 Genomes sequence index end in a space, which put a space into three URLs and file names | names trimmed, file names taken from the index's path column, the manifest validated; all 6,404 CRAM and index URLs answer 200 |
| The cohort run keeps what a cheaper method would be judged on | what a depth proxy needs: the scan estimate, the cohort's kept mosdepth bins, and the composition of those bins | placements were on a 10-kb grid and held class reads only: enough to learn sinks, not to say how much of a mosdepth bin is rDNA | placements on mosdepth's 1-kb grid with a census of every read in the bin (all, duplicate-flagged); against an independent duplicate-excluded depth over chr21's rDNA bins the census is within 3%; scan and fetch agree bin for bin inside the sinks, asserted in CI |
| The Marchenko–Pastur law can be fitted to a coverage spectrum to choose the number of PCs | simulated noise with unequal variances (rows ±40%, columns ±30%) and planted components; NGS-PCA's 1000 Genomes spectrum truncated at 100, 150, 200 values | **false**: the median- or quantile-matched fit called 34 components for 7 planted and 148 for 5, and 59 / 66 / 80 on the real spectrum depending on the truncation | the edge is fitted from its universal square-root shape instead, with a 1% margin (section 9): exact on those simulations, none in that noise, 38–40 on the real spectrum at every truncation. It still leans high when bins have heavy-tailed variances (a few components in noise alone), and components beyond the leading 20–25 are not reproducible between SVD runs - so `ngsdose pcsweep` lets known truths and transmission decide; all asserted in CI |
| A known truth's error under adjustment can be read off a regression of the estimate | review of `pcsweep`: a two-valued truth (chrX: 1 or 2), PCs unrelated to it, n = 3,200 | **false**: the estimate carries the variance of the truth itself (SD of log truth 0.35), every regressor adds 0.35·√(k/n) of estimation noise out of fold, and the reported error rose 0.010 → 0.040 at 46 PCs - adjustment would always have looked harmful for chrX | the error, log(estimate / truth), is what is regressed; asserted in CI, including a PC that follows sex. (Same review: one `N_PC` knob for two PC sets of different size killed `02_cohort.sh` at its last step - now `N_PC` and `N_CTRL_PC`, and a number beyond what a table holds is clamped with a warning; the one-standard-error band used 1.25 instead of 1.65 for a MAD-based SD.) |
| The experimental satellite panel measures array mass | whole-file scans of two 1000 Genomes samples that have HPRC release-2 assemblies (HG02258, ACB male; HG01884, ACB female) | HSat3 0.97, 0.98 of the assembly; HSat1A 0.93, 0.97; α-satellite HORs 0.97, 1.03; HSat1B 0.90, 0.83; β-satellite 0.69, 0.76 - and its k-mer recall on CHM13 itself is 0.69. **My first reading of HSat2 (1.52, 1.73: "not usable") was wrong**: the comparison script dropped arrays annotated together with an assembly gap ("GAP,HSat2": 17 and 9 Mb), and an array with a gap in it is no truth | the script tallies gap-containing arrays separately and compares a class only where it has none; HSat2 is undecided until gap-free samples are compared; each class's recall is measured and documented (`resources/build/panel_recall.py`); 200 samples of the cohort have assemblies: `04_hprc_satellites.sh` |
| Every loaded class has sinks in fetch mode | fetch with a panel class the sinks BED names nowhere (a satellite or telomere panel against an older sinks file) | **false**: the class was counted only where its reads fell inside other intervals, and nothing said so | the engine refuses; `--allow-missing-sinks` continues and lists the class in `sinks_missing_classes` |
| A trio table survives a column without variance | `ngsdose trios` on a constant column (EBV in blood-derived DNA) and on one with too few complete trios | **false**: division by the parents' variance raised, and one column's error ended the run | NaN for what cannot be estimated, a note for what cannot be run, the other columns unaffected; the pedigree may also be a PLINK PED/FAM or a child/father/mother table |
| A site that cannot let the engine read CRAMs can reproduce a fetch | - | the fetch intervals (controls padded by 600 bp, the sinks, merged) existed only inside the engine | `ngs-dose plan` writes them as BED, dropping contigs the input's header lacks (samtools passes over such a region and exits 0); the cut is counted in fetch mode, which reproduces the whole-file fetch exactly, and a cut counted in scan mode is refused by `ngsdose sinks` (control ends above 5% of primary reads) |
| `ngs-dose plan` needs no input | the GRCh38 bundle without `-i` | **false**: every contig was filed under one id and the controls' overlap check faulted regions on different chromosomes against each other | an id per contig name |
| A column without a slope has no permutation p | a constant column with 10 or more trios | **false**: the null was all NaN, nothing exceeded NaN, and p came out 1/(n+1) | NaN when the slope is NaN |
| A pedigree header may start with `#` | `#kid dad mom ...`, PLINK's `#FID IID ...` | **false**: dropped as a comment, and the layout guessed from the data | a first line starting with `#` is a header when it names the columns |

Not tested, and the cohort run will not test them either: a chemistry other than Illumina's;
DRAGEN alignments; an orthogonal assay for the absolute rDNA scale beyond the twelve ddPCR lines.
