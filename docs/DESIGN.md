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
At n = 10⁵ a systematic effect explaining 0.03% of variance is genome-wide significant (at UK
Biobank's 5 × 10⁵, 0.006%), so the error budget that matters is entirely systematic. NGS-DOSE is
built around removing those terms one at a time, and around controls whose true copy number is
known so that the removal can be checked rather than asserted.

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
a factor (1 − d_rDNA)/(1 − d_ctrl) that is specific to the library (duplicate-flag rates in the
single-copy controls run from 6% to 12% across the twelve pilot samples, and from 5.4% to 23%
across the first 1,748 cohort genomes, 1st–99th percentile 5.9–17.4%). Applying that factor to
our uncorrected 18S estimate reproduces the published value for the same CRAM (Hall et al. 2021:
301; ours 276.4 × 1.057 = 292), and
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
reads, and the shipped table for the three positional classes is the union of the two (80
intervals, 3.3 Mb); since 2026-09-22 `sinks.bed` also holds the telomeric repeat's sinks (TEL, 63
intervals, 810 kb, learned from 372 cohort scans; section 5), 143 rows covering 4.09 Mb in all.
*Decision: two modes. `scan` reads everything and is placement-independent; `fetch` retrieves
only the learned sinks plus the controls, and gives the same counts up to that capture fraction
in a fraction of the I/O.*

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
sinks under another bwa-mem pipeline on the same analysis-set reference (a probe of all 2,554
unplaced, random and decoy contigs plus the unmapped bin of an HGSVC bwakit 0.7.12 + postalt CRAM
found only 9 of 1.5 million 45S reads there, and a whole-file scan of a Google bwa 0.7.17 / hs38DH
HG002 BAM found 99.93 / 99.97 / 99.74% of 45S / 5S / DJ reads inside the NYGC-learned sinks; this
says nothing about other aligners: under DRAGEN 4.x most 45S and DJ reads of the one genome
checked are unmapped, section 12); and the fragment-length model above.

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
mass, `M = Σ_g T[g] / λ_R(g) / f`, with reads binned by their own GC, λ_R the read-scale curve
(the table nearest the modal read length) and the sum over the GC bins it supports; f, the
fraction of the class's reads in those bins, is reported as `gc_supported_fraction`.

A class whose counts cannot be trusted is reported with every value NaN and a `status`, never
with a number: a fetch that lacked sinks for it (`no_sinks_in_fetch`) or left some of them out
because the file's header lacks their contigs (`sinks_skipped`), or counts made with a panel whose
k, unit length or k-mer count differs from the bundle's (`panel_mismatch`); or a satellite
family with sub-options that a fetch read only at them (`subset_only`), or through a sinks BED
`estimate` does not know, so that which of its intervals were read cannot be told (`unverified`:
aSatHOR, HSat1B and HSat3, until `--fetch-sinks` names the BED; section 4). A positional class
that the bundle has no panel entry or unit for is listed in `skipped_classes`.

The denominator is not a separate quantity: it is the level of λ. It is made robust by
comparing every control region's count with its expectation under the fitted curve and dropping
regions (CNV) and whole chromosomes (aneuploidy, common in cell lines) that depart by more than
30% and 4% respectively — *and* by more than 5 standard errors, because at low depth noise alone
crosses any fixed threshold and trimming noise is biased. (That second condition was added when
a 0.75× test fixture read 3% low.) Each chromosome is tested on its own control regions, trimmed
of outliers against its own median rather than the genome's, so that a full trisomy or monosomy,
whose every region departs from the genome, is not trimmed out of its own test; the genome-wide
level it is compared with leaves out the chromosomes already flagged. A chromosome left with
fewer than 3 control regions is reported in `untestable_chromosomes` (chr22, which has 2, in this
bundle) rather than as unflagged.

## 4. The counting engine

`ngs-dose count` (Rust, htslib) turns one BAM/CRAM into a counts file of 70 kB (fetch) to 240 kB (scan
with the experimental panels); nothing downstream
touches the alignment again, so models can be revised without re-reading a biobank.

- **Rules**, identical for controls and classes except for unmapped reads (below): primary
  records only (secondary, supplementary and QC-fail records are skipped); duplicate flag
  ignored; no MAPQ filter (class reads are MAPQ 0 by construction; 0.2% of control reads are
  below 20); the unit of counting is the fragment 5′ end with soft clips restored (a hard clip,
  whose bases the record no longer holds, is not; a soft clip behind a hard clip is restored from
  0.1.1 on - `5H10S85M` starts at pos − 10, where 0.1.0 restored only an outermost soft clip -
  which changes counts only for hard-clipped primary alignments, of which bwa `-Y` and DRAGEN's
  default output have none). Unmapped reads are k-mer classified like any
  other read and count toward class totals. One that carries its mate's position is placed at
  that position, and a fetch whose interval holds that position returns it. A fully unmapped
  read is placed under contig `*` and is read in fetch mode only with `--unmapped`. Unmapped
  reads are left out of the control, per-contig and per-bin census tallies.
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
  was made with, and `ngsdose estimate` warns when a cohort mixes them. It also records the
  pipeline that made the input, from its header (`pipeline`: each `@PG` line's id, program and
  version, and `sq_sha256`, a hash of the sorted `@SQ` name, length and M5 lines, which
  `samtools view -H` reproduces), and in fetch mode the padding used (`pad`). Sinks are specific
  to a pipeline (section 12), so this is what says which sinks a file's fetch can be trusted
  with. Regions are matched to the bundle by name: counts whose *control* regions are not
  exactly the bundle's are refused (they are the denominator and the GC curve), while
  known-truth or dosage regions that a later bundle revision added are simply missing for counts
  made before it — so extending the truth sets never sends anyone back to the CRAMs. Output is
  reproducible (same input, any thread count, macOS or Linux): identical apart from the run's
  `elapsed_sec`, input path and engine build; placements are totally ordered.
- **Input integrity**: before reading anything the engine asks htslib whether the file ends in
  its end-of-file marker, locally or over HTTPS, and refuses one that does not: a truncated
  BAM or CRAM decodes without error, and what a truncated coordinate-sorted GRCh38 analysis-set
  file loses first is the unmapped reads, then the HLA, decoy, chrEBV, alt and unplaced contigs
  (where much of the rDNA and DJ lands, e.g. `chrUn_GL000220v1`, `chr22_KI270733v1_random`), then
  chrM, chrY and chrX. A stream whose marker cannot be checked up front (a pipe, a server without
  range requests) is checked at its end in scan mode and refused then if the marker is missing;
  a CRAM stream decoded with more than one thread cannot be, and stays `unchecked`, which
  `ngsdose estimate` reports. A CRAM is refused without `-T` or a non-empty `REF_PATH`, rather than letting
  htslib fetch its reference from the EBI server, and a `-T` FASTA that does not match the CRAM
  (checked against the `@SQ` M5 of the contig that failed to decode) ends the run at once
  instead of being retried. A run in which no read starts in a control region is refused: its
  file would hold nothing to estimate from. The lengths of the contigs it used go into the
  counts file, and `ngsdose estimate` refuses a file aligned to a build other than the
  bundle's (hg19 has the same contig names; the controls and sinks are coordinates). Output is
  written under a temporary name and renamed, so a killed job leaves nothing a resumed run
  would take for a result. In scan mode the mapped primary reads of every contig are recorded
  as well — a cheap cross-check on the region-based dosages (it is what exposed the `N` in
  chrM; section 15).
- **Modes**: `scan` streams the file (htslib decoding threads, a reader thread doing control
  bookkeeping, a worker pool classifying batches of packed sequence). `fetch` takes the merged
  union of the control regions, each padded by `--pad` bp (default 600, at least 400), and the
  class sinks. With `--unmapped` it also reads the unmapped bin (htslib region `*`; median 0.15%
  of primary reads in 1,748 1000 Genomes scans, 1st–99th percentile 0.10–0.26%) as one extra
  job. It hands the intervals to worker threads that each hold their own indexed reader, and
  counts a record in the interval containing its start so nothing is counted twice. An
  interval's counts are committed only once it has been read to its end, and a failed interval
  is retried with a fresh connection, so transient network errors neither lose nor double-count
  reads. CRAM decoding is restricted to the fields used (no
  qualities, names or tags). Remote `https://` inputs work through htslib (not `s3://` or
  `gs://`: this build of htslib has no S3 or GCS plugin); the engine points the bundled libcurl
  at the system CA store. A remote input that still fails after `--retries` attempts (default 5, per
  open and per interval) ends the process with exit status 75 (`EX_TEMPFAIL`), so that the caller
  can try again later; a missing file, a refused request (404, 403) or a wrong reference exits 1 at
  once, as does a remote file with no index beside it. A lost index request (a 503 on a remote
  BAM's `.csi`, for one) is retried like any other open. The engine's own messages and the counts'
  `input` carry a URL with its query string, fragment and user information redacted; for an input
  or index URL with a query string (a signed URL) htslib's own messages, which print the URL in
  full, are turned off before the first open, and one line on stderr says so (a URL without a
  query keeps them). One transient failure still exits 1 at
  once, without a retry: a read error in the middle of a remote scan, which cannot resume. A dead connection does not fail, it waits, and the thread that would retry is
  the one that is blocked: a watchdog therefore ends the process with exit status 75 when nothing
  has been read for `--stall-timeout` seconds (default 300), and the 1000 Genomes cohort script
  (NGS-DOSE-1000G's `pipeline/01_count.sh`) tries the sample again: three attempts, each run under
  GNU `timeout` (default 5,400 s, `SAMPLE_TIMEOUT`) as a backstop where that command exists.

  A class can be fetched only if the sinks BED has intervals for it; in the shipped bundle these
  are the positional classes (rDNA45S, rDNA5S, DJ) and TEL. The satellite families have no
  shipped sinks and are measured by scan, unless sinks are learned for them from scans of the
  same pipeline (`ngsdose sinks --classes`; section 5); for NYGC bwa-mem an experimental set,
  learned from 100 cohort scans, is in `resources/experimental/sinks.satellites.bed`, and no
  fetch through it has been compared with a scan yet. Sink intervals on contigs the file's
  header lacks cannot be fetched. When the sinks BED names its classes, a fetch by an engine since
  645ae55 refuses a loaded class that keeps no interval, whether the BED lacks it or every interval
  of it is on an absent contig, and from 0.1.1 on also a class that loses any interval to such a
  contig (a lost row without a class counts against every class): `ngsdose estimate` reports such
  a class as NaN (section 3), so the fetch would be spent for nothing. With `--allow-missing-sinks`
  (also since 645ae55) it is counted: a class without any interval is recorded in
  `sinks_missing_classes`, the intervals lost are recorded per class (`sinks_skipped`: intervals
  and bp) with a warning; `--classes` (engines from 2026-09-26 on) can leave such a
  class out: it is then neither counted nor listed. fae1124, the engine of the running cohort, has
  none of the three: it does not refuse, it counts such a class only where its reads fall inside
  other intervals without saying so, so `ngsdose fetchplan` refuses to plan such a fetch for it
  (below: what a fetch costs).
  `ngs-dose plan` writes the same interval set (without the unmapped bin) as BED. Sites that must
  cut the reads out with `samtools view -M -L` then count the cut in fetch mode (section 15);
  `plan --unmapped` prints the samtools steps that add the unmapped bin to the cut.

### What a fetch costs, and how it is chosen

A fetch decodes whole CRAM slices, so its cost is the bytes of every slice whose alignment span
overlaps one of its intervals, not the bases of the intervals. `ngsdose/cost.py` computes it from
a CRAM's own index, as `ngs-dose count -m fetch` reads the file. The engine merges the plan's
intervals where they touch or overlap (`fetch_plan`) and makes one indexed fetch per run
(`read_interval`); htslib decodes a slice again for every fetch whose interval overlaps it, and
re-reads its container's compression header. `ngsdose.cost` therefore prices a plan as the sum,
over the merged runs, of the slices each overlaps, a slice under k runs counted k times, and
reports beside it the floor with every slice decoded once, which a reader that sorted the plan's
slices would reach (`cum_mb_floor`). Slices are keyed by container and slice offset, so within
one fetch a slice is counted once however many contigs list it (in HG00096, 327 slices, 190 MB,
are each listed under 2 to 329 contigs, median 6, one under 329 of which 328 are HLA contigs; the
listings are on decoy, HLA, alt and unplaced contigs); a long slice behind short ones is still
found; and each fetch adds the compression header of every container it enters (all of them
together are 109 MB, 0.7%, of HG00096's 15.74-GB file). On the 13 NYGC indexes the engine reads
2-3% more than the floor in men and 11-13% more in women for `core_tel` (546.9 against 522.7 MB,
medians), up to 27% more on the controls alone and 20-27% on `xy_arrays` in women: a woman's chrY
holds few, sparse slices, and the 40 chrY truth regions decode the same ones one by one. Decoding
each slice once (sorting the plan's slices, or one pass over a chromosome's runs) is a possible
engine improvement worth that difference; the floor says what it would save. The index numbers contigs by their place in the header, so the contig order comes
from the reference's `.fai` or `.dict`, or `samtools view -H`; a counts file's `contigs` list will
not do, as it holds only the contigs that held reads (2,168 of 3,366 for HG00096). The index must
be the CRAM's own: an index of another file of the same sample gave wrong numbers, and nothing in
the index alone shows it. `ngsdose fetchplan` uses these costs to choose what a fetch reads.

- **Where the bytes go.** On 13 NYGC bwa-mem CRAMs (median 16.5 GB), the controls with all of the
  bundle's sinks (rDNA45S, rDNA5S, DJ, TEL: the `core_tel` preset, NGS-DOSE-1000G's fetch
  configuration on an image newer than fae1124) read 546.9 MB per genome, 3.33% of the file
  (medians; 522.7 MB with every slice decoded once). What each option costs alone, and what other plans cost, is in
  [fetch_examples.md](fetch_examples.md), which `resources/build/fetch_examples.sh` writes from
  `fetchplan` on the same 13 indexes; this document does not repeat those figures. The control
  file is always read; alone, its 800 control regions cost 182.1 MB, the 180 truth regions
  48.7 MB (up to 99.1 in a woman, whose chrY slices are decoded once per region) and chrM and
  chrEBV 10.0 MB (medians over the 13, padded as fetched). The bytes are
  set by the reads the intervals share their slices with. The smallest possible cost of a class
  is its share of the primary reads times the file size; against it, rDNA45S costs 3.7 times that
  floor (a floor of 57.7 MB), HSat2 2.9, the α-satellite HORs 2.0, but TEL 116 (a floor of
  1.2 MB) and SATR 160. The ten satellite
  families are 4.79% of all reads, so no fetch of all of them can be cheap. One 12-kb pile-up
  bin, chr2:32,909,000–32,921,000 (1.3 M reads in HG00096's scan, few of any class, and inside the
  sinks of TEL, rDNA45S and eight satellite families), costs 95.4 MB alone (69.5–126.2); without
  it the lowest held-out capture (1,375 scans) would fall from 99.39% to 96.46% for TEL and from
  99.89% to 99.80% for rDNA45S. It stays in the sinks; a capture target can drop it from rDNA45S, but TEL keeps it
  at targets above about 0.979 (below).
- **The menu.** `resources/fetch_menu.tsv` lists every option: the controls, the unmapped bin, each
  class with its panel, sinks file and interval statistics, a status (`shipped`: sinks in the bundle;
  `experimental`: sinks learned from scans, not in the bundle; `candidate`: no sinks), a tier (A
  core, B shipped and measured, C experimental with a named question and truth data, D
  experimental, costly or of lower capture) and presets. A row of kind `subset` names part of a
  class's sinks as an option of its own (below). Its tiers and presets are a first assignment,
  meant to be edited. `fetchplan --classes` and `--preset` select by name; `--budget-mb` adds
  options in tier order, cheapest first within a tier, and stops at the first that does not fit,
  so a lower tier never displaces a higher one (`--fill` goes on past it); the controls count
  toward the budget. A candidate is never fetched: its panel goes to `PREFIX.scan_panels.txt`, for
  the scans that will learn its sinks. `shipped` does not mean that a fetch has been compared with
  a scan for every such class: for rDNA45S, rDNA5S and DJ the cohort's fetches return what the
  scans placed in the sinks (fetch / scan reads inside them 1.00000–1.00019 over 1,748 genomes),
  but those fetches did not load TEL, so TEL's support is, as for the satellites, the capture of
  scan placements (99.39 / 99.86% held out, min / median).
- **Capture targets.** `ngsdose sinks --evaluate BED --stats FILE` writes, per interval, its share
  of the class in each scan (median, 10th percentile and the largest in any one scan), the reads of
  any class in its placement bins, and the cumulative capture with intervals ranked by share of the
  class per read of any class; `--held-out` records that the scans were not used to learn the
  intervals, and `fetchplan` says for each class whether its expected capture is held out, in-sample
  or of unrecorded provenance. A bin is scored as in `capture` and attributed to one interval, so in
  each scan the shares add up to the capture exactly, and the full-set cumulative capture equals
  `--evaluate`'s median and 10th percentile. The bundle's sinks have such statistics from 1,375
  cohort scans that none of their intervals was learned from (`resources/GRCh38/sinks.stats.tsv`:
  the 1,748 counted by 2026-09-25, less the 372 TEL was learned from and HG02258); held-out capture,
  min / 10th percentile / median, is 99.89 / 99.93 / 99.95% for rDNA45S, 99.93 / 99.97 / 99.98% for
  rDNA5S, 99.65 / 99.72 / 99.76% for DJ and 99.39 / 99.79 / 99.86% for TEL. The satellites' come
  from the 1,648 scans held out from their sinks (`resources/experimental/sinks.satellites.stats.tsv`,
  marked `# held-out: yes`). Both files carry the largest share (`share_max`), so every class of the
  menu can be trimmed per byte. `fetchplan --capture F` keeps, per class, the
  intervals of highest yield until the capture reaches F at the 10th percentile (or the median).
  Without an index, yield is share per read and the expected capture is the curve's. With `--crai`,
  yield is share per byte of the interval's own slices, so an interval of low share on costly
  slices (a decoy slice shared with other contigs, a pile-up bin) is among the first dropped. The
  chr2:32.91 Mb pile-up bin is one of rDNA45S's first drops, but it holds 1.3% of TEL's reads in
  the median held-out scan (up to 3.4%), so TEL keeps it at targets above about 0.979 and drops it
  below that. The stats file holds quantiles of the per-read curve, not per-scan shares, so the
  per-byte expected capture is a lower bound that holds scan by scan: the larger of the full
  capture less the largest share of each dropped interval, and the per-read curve up to the first
  interval not kept. Summing median shares instead overstated HSat3 by up to 0.28. Against each
  held-out scan's own capture, at targets of 0.95–0.999, the bound was never above the truth for
  the bundle's classes, and for the satellites never by more than the file's six-digit rounding.
  Because the bound is conservative, the per-read set can reach the target with fewer bytes;
  `fetchplan` then keeps it (`order` in the plan says which). An index prices the plan as the engine
  reads it, so a plan priced with one is dearer than the floor by the slices its runs share (a few
  percent in men, 11-13% in women on `core_tel`), and an option can even make a plan cheaper by
  bridging two runs that decode the same slices. Classes share intervals and slices. An interval that another selected class keeps is read
  anyway, so every class that has it keeps it; with `--crai`, so is any interval whose slices the
  plan reads for another option, and each class is trimmed again with those slices free. A class
  can therefore keep more intervals in a larger plan: at 0.995 rDNA45S drops the pile-up bin when
  TEL is not in the plan, and keeps it, read for TEL at no extra cost, when TEL is. Per option, the plan
  reports `mb_saved`, what its trimming saved of the plan with the rest of the plan as it is
  (slices another option reads are no saving, so the options' savings need not add up to the
  plan's, which a note gives), and `capture_lost`. A sub-option is never trimmed: it is fetched
  whole. HSat1B is not trimmed at the default 10th-percentile statistic: all its sinks hold
  97.17% there, so any target above 0.9717 keeps all its intervals. With `--capture-stat median`
  its full capture is 0.990, so a target up to 0.990 does trim it.
  [fetch_examples.md](fetch_examples.md) shows trimmed plans (examples 4, 6, 11 and 13). The
  expected capture leaves out a class's reads in other options' intervals, which a fetch also
  counts.
- **Sub-options.** The sinks hold the arrays the reads are placed on, so a family's reads can be
  split by array without a new panel. `resources/experimental/subsets/` names four such parts of
  the NYGC satellite sinks, each a menu row of kind `subset` (preset `xy_arrays`): DXZ1 and DYZ3
  (the aSatHOR intervals at the chrX and chrY centromere models), DYZ1 (the Yq12 HSat3 interval)
  and DYZ2 (78 HSat1B intervals, mostly on autosomes and chrX: 66 autosomal, 7 on chrX, 3 on
  decoys, 2 on chrY; where men carry reads and women do not: 21.7% of men's HSat1B, about a
  quarter of its Y-derived part, so an index of it rather than its whole mass). A plan fetches a sub-option whole
  and writes it under its family in `PREFIX.sinks.bed`; `ngsdose estimate` reports `DXZ1.reads` and
  `DXZ1.mass_Mb` (the family's mass per read times the reads inside) for scans and fetches alike. A
  fetch of a sub-option without the rest of its family counts the family only there, so the family
  itself is marked `subset_only` and not measured. `estimate` knows a fetch's sinks BED by the
  sha256 the counts record (the bundle's, `resources/experimental/*.bed` and the sub-options' own);
  a fetchplan BED is given with `--fetch-sinks`. Without it, a fetch through a BED it does not know
  leaves the sub-options `unverified` (reads, no mass), and the families that have sub-options
  (aSatHOR, HSat1B, HSat3) `unverified` too (NaN) wherever the fetch counted them, even in a plan
  that selected no sub-option: the placements alone do not show which intervals were read. Their
  costs are in [fetch_examples.md](fetch_examples.md) (examples 7, 15 and 16).
- **Partial panels.** Fetching a subset of a panel's classes (HSat2 without the other satellites)
  loads the whole panel. `ngs-dose count --classes A,B` (the fae1124 engine predates it) counts only
  the named classes, each exactly as the full load does in scan mode (every panel is still loaded and merged,
  so k-mers shared between panels are dropped as before, and every compositional family keeps its
  k-mers when a compositional class is named, since families compete for a read), and a fetch
  then reads and checks only the named classes' sinks, so a named class's reads inside other
  classes' sinks are not read (1 of 54,767 rDNA45S reads in the fixture); the counts record
  `classes_selected`.
  When the panels a plan loads define classes it does not select, `fetchplan --engine` runs
  `ENGINE count --help` and writes what that engine takes to `PREFIX.count_flags.txt`: an engine
  with `--classes` (from 2026-09-26 on) gets `--classes=...`, and the others are neither counted nor
  listed; one with only `--allow-missing-sinks` (645ae55 up to that change, e.g. 7772e32) gets that flag,
  the counts list the others in `sinks_missing_classes`, and `ngsdose estimate` leaves them NA
  (`no_sinks_in_fetch`); fae1124, the cohort's engine, takes neither, and `fetchplan` refuses the
  plan. `--unmarked-companions` writes it anyway: the fetch then counts those classes only inside
  the plan's intervals and the counts do not say so, and only `estimate --fetch-sinks
  PREFIX.sinks.bed` marks them `no_sinks_in_fetch`. A plan of whole panel files (`core`, `core_tel`)
  has no such classes and needs no flag, so fae1124 runs it. Without an engine to probe,
  `fetchplan` writes `--allow-missing-sinks`, which fae1124 rejects as an unknown option, so that fetch fails at once rather than undercounting.
- **The controls.** They are the largest fixed cost, and a lighter set is a named subset of the
  bundle's (section 6 gives its cost and what it changes). `fetchplan --controls` costs a plan on
  it and writes the FASTA the fetch must use to `PREFIX.controls.txt`.

These costs belong to one pipeline and one file layout: another aligner puts the reads elsewhere
and another writer cuts slices differently, so a site costs its options on its own indexes
(`fetchplan --crai`, the median over several) with sinks learned from its own scans (`--sinks`).
The experimental sinks and sub-options here were learned from NYGC bwa-mem scans of the GRCh38
analysis set and are for those CRAMs only.

## 5. k-mer panels

`ngs-dose panel` builds the panel from class FASTAs and background genomes.

A **positional** class has a unit consensus; a retained k-mer occurs exactly once in the unit, so
a hit is a coordinate. A **compositional** class is a family without a stable unit; a k-mer only
says "this read is class X". For the bundle's classes (all positional), a k-mer is kept only if
it never occurs in the background genomes — GRCh38 (analysis set, with decoys) and
T2T-CHM13v2.0 — outside the intervals where the class legitimately lives
(`resources/GRCh38/build_inputs/*.class_loci.bed`). The experimental compositional panels are
built differently. A satellite k-mer must recur at least ten times in its family's CHM13 arrays,
occur in no other family, and occur nowhere in CHM13 outside CenSat-annotated satellite; CHM13 is
the only background, and GRCh38 is not screened. The TEL panel is unfiltered.

| class | unit | k-mers kept | notes |
| --- | --- | --- | --- |
| `rDNA45S` | KY962518.1, 44,838 bp, circular | 31,827 of 41,476 | Alu and simple-repeat spacer k-mers removed by the background filter; 99.5% of 18S and 93% of 28S read positions remain recoverable at 150 bp |
| `rDNA5S` | X12811.1, 2,231 bp, circular | 1,945 of 2,217 | 68% GC throughout, so no anchor windows: the estimate rests on the GC model |
| `DJ` | CHM13 chr21:2,708,299–3,108,298 (400 kb distal to the rDNA array) | 169,808 | restricted to the *core*: k-mers that occur once in the chr21 unit and five times in CHM13, every occurrence inside one of the five acrocentric distal junctions (169,277 of the 170,056 occur once in each junction); 169,808 of them survive the panel's GRCh38 filter, which excludes its 23 DJ-like loci, and the cross-class filter |

The filter is deliberately strict. The genome holds ~30 dispersed rDNA-derived fragments
(200–700 bp, 91–99.8% identical; chr1:91.39 Mb, chrX:109.05 Mb, chr12 in CHM13, …). They are
*not* exempted, so the k-mers they share with the unit are dropped and that stretch of the unit
goes blind: at 150 bp, 367 positions per strand at the 28S 3′ end (11.6–12.4 kb) yield reads
with fewer than the four panel k-mers a read needs. In total, 1,348 of 44,838 positions per
strand are blind. Most of the rest (936) lie in the spacer at 20.7–29.5 kb, where the unit
repeats itself (an internal duplication, two Alus, and CATA/CT simple repeats), so a positional
panel cannot keep those k-mers. The estimator masks more conservatively. With its default
`--min-kmers 20` it drops 2,827 position-strands per strand (18S 92.5% usable, 28S 88.3%) from
numerator and denominator alike, which is preferable to counting reads that may come from
elsewhere. One locus *is* exempted: chr21:8,986,604–8,988,749 is a 99.6%-identical 5′ETS/18S
piece of genuine rDNA that GRCh38 happens to place outside the annotated copies; without the
exemption the 18S keeps a third of its k-mers.

**A new class is added scan first.** Its panel is loaded in whole-file scans, which classify
every read wherever the aligner put it, the unmapped ones included. Its sinks are learned from
those scans (`ngsdose sinks`: a 10-kb neighbourhood becomes a sink where it holds at least 1e-5 of
the class's reads and at least 25 reads in any training scan; the sink keeps the placement bins
that hold reads, 1 kb for a positional class and 10 kb for a compositional one, padded by 1 kb).
Their capture is checked on held-out scans (`ngsdose sinks --evaluate`, with `--stats` for the
per-interval statistics that capture targets use; section 4). Only then does the class get a sinks
file in the fetch menu and become fetchable. A class's reference coordinates are never taken as
its sinks: the aligner places reads by what it can place them on, which depends on the pipeline
(alt, decoy and unplaced contigs, alt-aware or alt-masked alignment) and on the person. In one
whole-file scan of NA12878 with the candidate panels below, 58 of the 72 classes with at least 25
reads had at least 98% of them at their GRCh38 reference copies, and 14 did not: GSTT1 65%
(GRCh38 carries it only on `chr22_KI270879v1_alt`, and ALT-aware bwa moves reads off it), the 4qA
D4Z4 end 81%, CCL3L 85%, eight others 90.5–97.6%: reads placed away from the reference copies. Three
more had few or none there, which measures panel specificity, not placement: UGT2B17 (0% of 174;
deleted in NA12878), DAZ (0.9% of 109) and RBMY (0% of 47) in a female. The first build also gave
HHV7, which has no GRCh38 copy, 283 reads at chromosome ends from simple-repeat k-mers; the
periodic-k-mer filter now applied to every candidate panel (below) removed them. Among the satellites,
about 87% of male HSat1B is Y-derived (DYZ2), yet it lands on autosomes, chrX and decoys, not on
chrY; the DYZ2 sub-option (78 intervals where men carry reads and women do not) holds 21.7% of
men's HSat1B, a male-specific index rather than the whole Y-derived mass
(`resources/experimental/README.md`).
`candidates/reference_copies.GRCh38.bed` holds each candidate's reference copies on the placement
grid, so that `--evaluate` can measure on held-out scans what fetching them would lose.

To estimate a positional class, the estimator needs its unit sequence (expected counts per window)
and its panel entry (which read starts are callable). The bundle's classes have both. For another
positional class, `ngsdose estimate` takes them from experimental units
(`resources.ExperimentalUnits`: `<class>.fa` in `resources/experimental/candidates/units/` or a
directory in `NGSDOSE_EXTRA_UNITS`, one record, and the panel file beside it, which must be among
the panels whose sha256 the counts record, with the counts' k and the unit's length). Such a class
is estimated with the same windows, masks and GC model as a bundle class, but with no anchor
windows and no cohort calibration: status `experimental`, `cn` = `cn_all`, and a per-sample
`cn_se_rel` (the windows' log-SD over the square root of their number, plus the Poisson term). An
all-window estimate can sit off the anchored one by a level: in 60 cohort genomes the all-window
45S estimate is 0.935–0.978 of the anchored one (median 0.958), DJ's 0.996–1.025 (median 1.012).
Without a unit, or a custom bundle (`-r`) that carries the class, it is reported as `skipped: not in
the bundle panel`, with where a unit was looked for.

Dispersed sequence — the satellite families — is compositional. Its reads are spread over hundreds
of loci, on unplaced, unlocalized and decoy contigs as much as on the chromosomes: by family, a
median of 3–89% of its reads lie off the assembled chromosomes (medians over 200 cohort scans:
α-satellite HORs 3%, HSat2 4.5%, HSat1A 18%, HSat3 27%, SST1 39%, HSat1B 51%, β-satellite 54%, CER
and SATR 58%, ACRO 89%; rDNA45S, fetched through the shipped sinks, 45%). Those contigs can be
fetched like any other region, and in the NYGC bwa-mem CRAMs the aligner puts each family on a
small, stable set of intervals. Sinks learned with `ngsdose sinks --classes <families>` from 30
whole-file scans held ≥ 99.8% of HSat1A, HSat2, HSat3, α-satellite HORs, β-satellite, ACRO, SST1,
CER and SATR in every one of 200 held-out genomes (≥ 99.85% in two of three random draws of the 30
and the 200), in 0.2–3.5 Mb of intervals per family, 60 Mb
for the HORs. Together with HSat1B's they come to about 81–82 Mb (merged; 81.3 Mb in 1,840 intervals
on 488 contigs in the draw used here, 82.4 Mb in a second random draw of 30 and 200 scans), and
their placement bins hold about 6% of a genome's primary reads (5.0–6.8% in 30 of the held-out
genomes, half of it the HORs'), against 100% for a scan. HSat1B, mostly on Yq, is the exception:
≥ 96.6% (median 98.8%), because part of it is left fully unmapped in the 698-genome batch; it needs
`count --unmapped` and a capture calibrated on scans. These figures come from the scans' placement
histograms, with a placement bin counted as captured only when all of it lies inside a sink; no
satellite fetch has been run yet. Sinks depend on the aligner and the reference (section 12), so
they are learned where files are scanned whole: in a whole cohort, or in a biobank's scanned subset.
That is why a cohort that is scanned once is scanned with the **experimental** panels under
`resources/experimental/` loaded. No satellite sinks ship with this bundle. An experimental set
for NYGC bwa-mem sits beside the panels (`resources/experimental/sinks.satellites.bed`, with
`sinks.satellites.stats.tsv`): learned by the same rule from 100 cohort scans (25 per release
batch and inferred sex, 15 populations), 2,189 intervals, 85.1 Mb summed over the families. In the
1,648 cohort scans not used to learn it, it holds at the lowest 99.83% (SATR) to 99.99% (HSat2) of
nine families and 96.59% of HSat1B (97.92% with the unmapped bin); the medians are 99.89–100.00%,
and 99.02% for HSat1B. What fetching each family costs, alone and under a capture target, is in
[fetch_examples.md](fetch_examples.md) (examples 11 and 12). For the shipped positional classes
the fetch returns what the placements say (fetch / scan reads inside the sinks 1.00000–1.00019 in
1,748 genomes); for the satellites that has still to be shown with real fetches.

The telomeric repeat turned out not to be dispersed: in 372 NYGC scans the aligner put 92% of its
reads within 25 kb of a chromosome end and 60% into one 10-kb bin of chr5p, essentially none on
unplaced or decoy contigs, and the rest at a fixed set of interstitial loci (chr4 at 94 kb from
the 4q end, chr18 at 113 kb from the 18q end, chr2:32.9 Mb, …). Sinks learned from 30 scans (about
0.7 Mb) held ≥ 99.56% of it in each of 200 others. The set the bundle carries
(`ngsdose sinks --classes TEL`: 63 intervals, 0.81 Mb, learned from all 372) captures a median
of 99.87% and at least 99.39% (1st percentile 99.68%) in each of 1,748 cohort scans, so fetch
mode measures the class, to that capture, when its panel is loaded.

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
  summed; the 200 genomes of the 1000 Genomes cohort that have one, in the published results page,
  NGS-DOSE-1000G `python -m report --censat`): the median estimate/assembly ratio is HSat3 0.94,
  HSat1A 0.91, HSat1B 0.86, α-satellite HORs 1.02 and ACRO 0.79. β-satellite (0.74) and CER (0.43)
  read low by about their recall. Per genome the two measurements agree to a robust SD of 2–8%
  for eight families; across people, r = 1.00 for HSat1B, 0.95 for ACRO, 0.93 for β-satellite,
  0.91 for CER, 0.79 for HSat1A and 0.76 for HSat3. For the HORs r is 0.70, because people differ
  by only 5% while the two agree to 3%. HSat2 agrees poorly even in the 79 assemblies with at most
  2% of its arrays in marked gaps (median ratio 1.08, r = 0.39), so
  what the HSat2 panel measures is heritable but not yet confirmed to be HSat2 mass. An assembly
  is a truth only for the arrays it spans;
  `ngsdose/hprc.py` leaves a sample out of a class when gap-containing arrays exceed 2% of what is
  annotated. (Standalone `GAP` records next to a class's arrays count as gaps of that class.) SST1 and SATR are annotated several times more generously in the HPRC assemblies
  than in CHM13, so their absolute ratios mean nothing.
- *The telomeric repeat* (`TEL`: the six canonical 31-mers of (TTAGGG)n, unfiltered). A relative
  measure, not a telomere length: a read is assigned with 34 bp of perfect repeat, which
  interstitial telomeric sequence also has, and an exact 31-mer is lost to one sequencing error
  where TelSeq's hexamer count is not. For every class the counts keep an eleven-bin histogram
  (`hit_frac`) of the share of each assigned read's k-mers that hit the panel, in tenths, so a
  stricter threshold can be chosen, and calibrated against TelSeq, afterwards; reads below the
  scan-time minimum of 4 panel k-mers are not recorded, so it cannot be lowered. Fetchable through
  the bundle's sinks (above); the whole-file count and the count inside the
  chromosome-end windows an NGS-TL/TelSeq-style query would retrieve agree at r = 0.9997 across
  the 372 scans.

**Candidate classes** (`resources/experimental/candidates/`): 83 classes in seven panel files,
5,076,728 k-mers, built from a literature assessment of which multi-copy sequence has truth data
for 1000 Genomes samples. None has sinks; they are for the scans that will learn them. The files
are final for the cohort's remaining scans; their sha256 are listed in `candidates/README.md`.

| panel | classes | k-mers | what |
| --- | --- | --- | --- |
| macrosatellites | 12 positional | 45,362 | DXZ4, CT47, RS447, MSR5p, FLJ40296, RNU2, D4Z4, ZAV, REXO1L1, and the 4qA, 4qB and 10q D4Z4 ends |
| multicopy-genes | 27 positional | 544,098 | LPA KIV-2, C4 and its HERV-K insertion, CYP21A2, AMY1, AMY2B, SMN, SMN1, RHD, HBA, the common deletions (GSTM1, GSTT1, UGT2B17, LCE3B/C, APOBEC3B), HPR, CCL3L, DEFB, NOTCH2NL, ... |
| sex-chromosome-arrays | 12 positional, 1 compositional | 126,346 | TSPY, RBMY, DAZ, BPY2, CDY1, CDY2, DYZ19; the opsin array and its LW/MW exon 5, GAGE, CT45, SPANXB |
| rna-arrays | 5 positional | 148,026 | the 1q23 tRNA-gene array, U1 and U3 snRNA units, SNORD116, SNORD115 |
| nonhuman | 5 positional, 1 compositional | 395,211 | HHV-6A, HHV-6B, HHV-7, SMRV, phiX, EBV type 2 |
| nonhuman-myco | 1 compositional | 3,731,523 | culture Mycoplasma, five species as one class |
| coding-vntrs | 19 compositional | 86,162 | long-unit coding VNTRs (ACAN, MUC1, MUC19, FLG, NEB, ...), measured as array length |

Each group was built with the panel builder against GRCh38 (analysis set) and CHM13 at
`--max-bg 0`, with only the class's own copies masked, so a k-mer a paralog also carries is
dropped; the RNA arrays went through one more filter that also drops a k-mer one substitution
away from a background occurrence, after a first scan put 0.1–0.4% of each class's reads at single
foreign loci. The assembly step then removes from every candidate panel each k-mer that is
periodic with a period of 1–6 bp at ≤ 2 mismatches (78 k-mers, 55 of them viral), which took
the simple-repeat floor of HHV7 and HHV6A in NA12878 to 0. Of the 17 circular units, the 14 cut from
an assembly are each one exact period of their array (the first 31-mer recurs at the base after the
unit; REXO1L1 was re-cut, and TSPY, OPN1, CT45 and DEFA1A3 trimmed by 1–22 bp); of the GenBank
units, the RNU2 clone was trimmed of the 6 bp its ends repeat, D4Z4 needed no trim, and PHIX is a
whole circular genome. Two joins remain partial because of allele differences, not unit ends: 17 of
RNU2's 30 junction 31-mers are in the arrays, and 12 of DEFA1A3's in CHM13 (30 in GRCh38). The two
hs38d1 decoys that are pieces of the MUC19 and MUC6 arrays are masked as those classes' own copies.
No candidate k-mer is in a shipped panel, and none is in two candidate files, so a class has the
same k-mers whichever panels a run loads. Loading all seven changed no count of a shipped class, on
the test fixture with three engine builds (fae1124, 7772e32 and the fetch-menu build) and on a whole 30×
NA12878 CRAM, and each class counted, in all three builds, exactly the reads cut from its unit that
the classification rule predicts. With the final files, fae1124 scanned that CRAM in 119.5 s and
1.65 GB of peak memory, against 106.1 s and 1.43 GB without them (the fetch-menu build: 122.9 s and 1.55 GB
with them). fae1124 runs repeated over the review rounds on the same shared machine gave
109.6–119.5 s and 1.49–1.81 GB with them (the fetch-menu build 107.0–122.9 s, 1.48–1.71 GB). The
Mycoplasma file, 74% of the k-mers, is most of the memory and can be left out.
`candidates/candidates.tsv` gives each class's unit, recall of reads cut from its copies, tier
(A: per-sample truth for cohort members, B partial, C relative or technical, D exploratory), what
it measures, the truth available and its caveats; several are kept below 50% recall on an argument
recorded there (SMN1, whose k-mers sit only over the 32 sites where it differs from SMN2; RNU1;
the D4Z4 ends). Some classes carry a floor of off-target reads that the cohort's scans must size:
one-substitution alleles at other loci (UGT2B17, GSTM1, APOBEC3B, the chrY classes in women) and,
for EBV type 2, type-1 reads carrying a sequencing error or a minor variant at a site where the
types differ, 1.9e-5 of chrEBV records in NA12878, so type 2 is called from the EBV2/chrEBV ratio,
not from a count. From the placement
bins of the one NA12878 scan, fetching all 83 would add about 109 MB to the controls and the
bundle's sinks (38 tier-A classes: 33 MB; medians over 13 NYGC bwa-mem CRAMs, with every slice
decoded once, the floor; the engine's figure is higher, since many small bins in few slices is
where its per-run decoding costs most): planning figures,
not sinks, and from one female genome, so the chrY candidates have no bins there apart from
off-target reads (DAZ 0.27 MB) and their real sinks are not priced.

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
log N, fitted per sample on the GC table of the tabulated fragment length nearest the sample's
insert median. Its supported range is where the controls hold at least 2,000 position-strands per
1% GC bin (12–82% at L = 450, 12–83% at L = 400; 3–87% for the read-length curve used for
compositional classes), and outside it nothing is estimated. A sample whose controls cannot
support a curve at all (complete GC dropout) is refused with that reason rather than fitted. The
range depends on the bundle and L, deliberately not on depth: the first version defined it by how
precisely the curve was known, which made the set of usable windows shrink with depth and the
all-window estimate drift with it (section 15). How well the curve is determined is reported
per sample instead (`gc_curve_max_se`). The headline
estimate is insensitive to the size of the control set — 508.0, 508.3 and 508.2 with 800, 400
and 200 regions (513 with 100) — because moderate-GC windows sit where the curve is best
determined; GC-rich features move by 3–4%, which is one reason they are not the headline.

**A lighter control set.** The controls are the largest fixed cost of a fetch (section 4), so the
bundle carries a named subset, `controls.lite200.bed` / `controls.lite200.fa.gz`: 200 of the 800
control regions with all 182 truth and dosage regions, chosen from the reference alone (at least
three regions per chromosome, both of chr22's two, GC strata within each, and of 300 such draws the one whose
position-GC distribution at L = 150, 400 and 450 is closest to the 800's), its records
byte-identical to the bundle's. It reads 120.9 MB instead of 260.0 (median of 13 NYGC CRAMs; 93.1
and 231.4 with every slice decoded once). The
estimator accepts counts whose control regions are exactly the bundle's or exactly one named subset
(`controls.<name>.bed` in the bundle, or `control_subsets` in `bundle.json`), and still refuses any
other set; the result records `controls_used` and `controls_subset`. On the 1,748 cohort genomes
(the subset's tables rebuilt from each scan's per-region counts, a rebuild checked against real
engine runs on the pilot's region cuts to an SD of 0.2–0.6%), lite200 against all 800:

- 45S copy number +0.21% (median; SD 0.30%, largest 1.1%), against a cross-library SD of 3.8–3.9%
  in the pilot; 5S −0.06% (SD 0.67%); DJ +0.16% (SD 0.27%); truth.auto +0.11% (SD 0.19%); the
  satellite masses SD 0.25–0.78%, except HSat1B (+1.45%, SD 2.27%, one genome 23%).
- Transmission reliability in 385 trios: 45S 0.9993 and 0.9992, 5S 0.945 and 0.937, DJ 0.751 and
  0.752, HSat1B 0.969 and 0.985 (bootstrap intervals about ±0.1 or wider).
- The GC support narrows from 12–82% to 18–78%, so 45S's usable fraction falls from 0.852 to 0.790
  and its all-window estimate rises by 1.48% (the anchored headline does not).
- The aneuploidy test flags 20 genomes, all among the 32 the full set flags: 12 are missed.

A choice of the widest GC support instead kept a region with odd AT-rich positions and made HSat1B
unstable (SD 13%, one sample off by 2.7×). Below 200 regions the gain is small and the losses grow:
100 regions save a further 22.2 MB (median of the same 13 CRAMs, every slice decoded once: 70.9 MB against 93.1), but the chromosome test raises 14 new flags and DJ's
reliability falls to 0.731; at 50, 5S has outliers of 15% and 14 chromosomes cannot be tested.
Thinning the truth regions as well (20 per set) would save 24 MB more in HG00096 (every slice decoded
once; more as the engine reads, since a woman's 180 truth regions cost up to 99.1 MB there) and nearly doubles
truth.auto's SD (0.008 to 0.015); it is not shipped. A cohort should not mix fetches with the two
control sets: lite200 shifts the median level of each headline class by up to about 0.3% (SST1
+0.29%, 45S +0.21%), but the all-window 45S estimate and HSat1B by about 1.5%, and `ngsdose estimate` warns when its
inputs were made with different sets.

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
unlikely to share a bias where they agree. Anchor windows must also have fragment GC 40–60%.
Without `anchors.json` (`--gc-rule-anchors`; a bundle that names a missing `anchors.json` is
refused), and for classes that have none (5S, DJ), the GC rule alone picks them. If too few
remain (fewer than 5 windows present in ≥ 90% of samples in `ngsdose cohort`, or fewer than
1,000 fragment ends in the anchor windows in a single-sample estimate), the level rests on every retained window instead
(every usable window for a single sample); the single-sample table says which in
`CLASS.cn_basis`.

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

### A class with rules of its own: the distal junction

One model of every window at one copy number does not hold for the distal junction, and the
cohort run with the HPRC assemblies of 28 of its genomes showed how (NGS-DOSE-1000G,
`docs/DJ.md`). The bundle's `calibration.json` carries what follows for the class, and `ngsdose
cohort` applies it to estimates of the unit the rules were made for (`--no-class-rules`
calibrates as before; `DJ.cn_unit` is what `DJ.cn` was):

- **The level is set on the core.** Junction copies differ in parts of the unit: some lack the
  distal 22 kb, some 197–217 kb, and smaller intervals near 132, 162 and 227 kb vary too.
  `level_exclude` (0–30, 128–137, 155–170 and 190–232 kb) leaves them out of c_i and of the
  anchors; their windows keep their efficiencies and stay in the profile.
- **The scale is pinned to the cohort's mode.** The junction has ten copies in nearly everyone,
  and nine genomes whose assembly holds ten complete copies read 9.37–9.79 on the scale the
  fragment-GC model sets in the GC-rule anchors: the model reads this class 2.2% low (the
  cohort's mode is 9.78), a residual in repeat context that also shows as window efficiencies
  falling with GC (−0.34 log units per unit GC here, −0.53 for the 45S, whose shipped anchors
  correct it). With fifty genomes or more, `scale: mode` takes δ = log(mode / 10) from c_i and
  adds it to a_w, and records it in the efficiency table, so that a saved table carries the pin
  to genomes counted later. A smaller cohort is not pinned, and says so.
- **Polymorphic intervals are put on whole numbers.** Where a deletion is common the cohort's
  median genome lacks part of a copy, so the median over genomes gives the interval's windows an
  efficiency too low by that part, and every genome reads the interval too high by one factor: a
  comb of whole numbers half a copy off at 197–217 kb, nearly a whole copy at 5–15 kb.
  `cohort.comb_offset` finds the offsets that put the cohort on whole numbers and takes the
  smallest that leaves at most 15% of the genomes above their own level, since a copy can lack
  an interval but seldom holds it twice. On the cohort that is ×1.075 (5–15 kb), ×1.028
  (15–23 kb) and ×1.046 (197–217 kb), each within 0.02 of what the assemblies give; every copy
  holds 197–217 kb in 59% of genomes, one lacks it in 30%, two in 9%. A cohort too small for a
  comb takes the offsets the bundle records.
- **Whole numbers of copies are called along the unit** (`segments.py`). A genome's calibrated
  profile is read as a chain over its windows whose states are whole numbers: a window in state
  k is expected at k · f · exp(g · t), f the genome's scale, g its lean (the log change across
  the unit) and t the window's position; a deviation counts at most 3.5 SDs of the window's
  noise, a change of state costs τ = 16 log units, and the most probable chain is found for every
  f and g of a grid under Gaussian priors (SD 1.5% and 3%). Runs of fewer than 16 windows go to
  the neighbour that fits them better. One copy over about 10 kb, or two over 4 kb, can be found.
  A genome is described against ten copies where it holds ten over 40 kb or more of the core,
  otherwise against the state that holds most of it (`DJ.copies`); a gain that reaches an end of
  the unit is a partial copy, a loss that does a partial loss (`DJ.partial`), anything else a
  local event (`DJ.variants`). A genome whose level lies between two whole numbers throughout
  has two readings; `DJ.call_gap` is how many log units the second is behind, and below three
  the call is `uncertain`.
- **What the whole numbers leave is kept** (`segments.find_fractions`). Whole numbers are what a
  germ line holds. A junction lost or gained in part of the cells leaves the profile a fraction
  of a copy off them, and the chain would hide it: the nearest whole number is called and the
  difference goes into f. So the difference is measured. A genome's scale is measured without
  a prior and set against the scales of the cohort: `DJ.off` is its level less its whole
  numbers, in copies, and `DJ.off_z` the same in the cohort's robust SDs. And steps of
  fractional height are looked for in the log of each core window over its called state: a
  level per stretch is fitted by weighted least squares together with a lean and a slope on
  the windows' GC, whose priors are the cohort's own spreads of them; steps are added one at a
  time while each improves the fit by τ and is a quarter of a copy high (a pair at once, for 2τ,
  where no single step does: a stretch inside the unit); and a step is kept where its stretch
  stands four robust SDs from what the same fit finds in the cohort's other genomes. The cohort
  is the measure because its profiles wander more than counting alone allows: the contrast
  between two stretches of 50 kb has a robust SD of 0.10 to 0.20 copies over the cohort where
  counting gives 0.07 to 0.12. A call that is not uncertain is then `settled`, or `fractional`
  where the level lies three SDs from its whole number or a stretch stands off
  (`DJ.fractional`). With fewer than fifty genomes no step is looked for and the level is
  judged against the spread that a saved efficiency table carries, or the bundle's.

The lean is there because some profiles rise or fall smoothly along the unit, by 1.5% (robust
SD) and by more than 5% in one genome in eighteen, which a chain of whole numbers would break
into a step. It goes weakly with the release batch (more than 5% in 5.1% of the genomes of one
and 7.4% of the other) and with the library's GC response (|r| at most 0.08), and not with depth
or mitochondrial load; none of them accounts for it, and its cause is not known.

What this buys, on the cohort. Of 3,202 calls 3,124 are settled: 2,747 genomes hold ten copies
throughout (the polymorphic intervals aside), 106 nine, 7 eight, 23 eleven, 4 more; 207 carry a
copy that holds or lacks an end of the unit, the commonest a copy of the first 316 kb (60
genomes) and one of the first 262 kb (45). 49 calls are fractional (29 by the level, of which
23 lie below their whole number, where a normal scatter of the scales would give 9 in all; 20
by a stretch of the unit) and 29 uncertain. Against the assemblies, every one of the 7 partial
copies in the 19 genomes whose assembly resolves the junction is called, with the breakpoint
within 4 kb, and the called state equals the assembly's in 90.9% of the 5-kb blocks of the 17 of
them whose call is settled; the other two are fractional, a fraction in the reads where the
assembly, made from another culture of the line, holds a whole number. In the trios, of the 556
whose three calls are settled, 553 have no 5-kb block of the core that the parents' states do
not allow (99.6% of 33,151 blocks), and a parent explains 95% of the blocks in which a child
deviates. The common deletions, which are in the germ line beyond doubt, pass as a germline
variant should: where one parent lacks one copy of 197-217 kb and the other none, 67 of 137
children lack one.

What it does not buy. A genome counted alone has no cohort to pin its scale or to fit its
polymorphic intervals: it takes both from a saved efficiency table of its own chemistry, and
without one its level rests on the GC model, a few percent low. The shorter polymorphisms are
left out of the level but not genotyped, and in the short polymorphic intervals a deletion is
called in 8-16% of the children of two parents without it: an event of 10 kb is near the limit
of what the chain can find. A fraction is seen in a window of its own: put into genomes of the
cohort, a change in two fifths to three fifths of the cells is fractional or uncertain in half
to nine tenths of them (a fifth to a third for a stretch of 110 kb); a change in fewer cells
passes as none and one in more as a whole copy, though the level reads it at every share. What
a fraction is, the measurement cannot say: a tenth of the level's offset lies between groups of
samples handled together, and spouses' offsets go together (r = 0.20). And the trios hold
something the calls do not explain: a whole-copy step passes to 15 of 44 children (two-sided
binomial p = 0.05 against one half) and a partial copy, looked for in the child by its
breakpoint, to 29 of 70 (p = 0.19), and a father's loss of a junction to 4 of 23 where a
mother's passes to 10 of 16 (Fisher exact p = 0.007, a comparison made after the fact). The
readings leave little room (where a step was not passed on, the parent's level lies within 0.3
copies of its whole number in 28 of 29 pairs, and none of these parents is fractional), and the
events did not arise in culture at a rate the children's lines would show. The results
repository's `docs/DJ.md` has the tests.

## 8. Known-truth controls

Four sequences of known copy number are measured in every sample under the same counting rules,
fragment-GC curve and control denominator as the positional classes (the region truths by
alignment position, as whole regions without the classes' GC-support and k-mer callability
masks; DJ through the k-mer path the classes take), so accuracy is observed per sample rather than
assumed:

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
libraries do not, which is what a larger S-phase fraction in the source culture would do. The
cohort does not bear this out. In the NGS-DOSE-1000G report (all 3,202 genomes; the 1,538 women
whose line has kept both X chromosomes, chrX 1.85–2.15), DJ and female X do not move together
(r = 0.01, 95% CI −0.04 to 0.06), where the hypothesis predicts a positive r. DJ's inherited copy-number spread between people can hide a small shared effect, so this
weakens the S-phase explanation without ruling it out; the test against the leading control PC
is still to come. The assemblies say where the junction's deficit lies, whatever its cause: nine
genomes whose HPRC assembly holds ten complete junction copies read 9.37–9.79 here, so the
deficit is in the measurement's scale and not in the copies, and the cohort layer now pins the
scale of this class to its mode (section 7).

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
path resolves what it claims to. At cohort scale this became whole numbers of copies called
along the unit and a Mendelian test on them (section 7).

## 9. The cohort layer

- `ngsdose estimate` — counts → per-sample estimates and the summary table. Each counts file
  is checked first: a wrong format, a wrong build or a control set other than the bundle's makes
  that file fail (the others are still estimated, and the run exits 1); an unchecked or missing
  end-of-file marker, sinks missing or skipped for a class, or a class counted with another
  panel than the bundle's is printed as a warning, and such a class is NaN with its `status`
  (section 3).
  Each estimate is named after its input file, so two inputs with one file name are refused.
- `ngsdose cohort` — window calibration of the positional classes; writes `CLASS.cn`, its robust
  relative SE (`cn_se_rel`), the profile roughness (`profile_sd`) and profile PCs
  (`profilePC1..3`), and with `--save-efficiencies` the efficiency table, which `--efficiencies`
  applies to new samples instead of refitting (refused when its windows are not this cohort's).
  A class is calibrated on the samples that have usable windows for it; the others get NA and are
  named in a warning. A sample id given twice is refused: a cohort takes one estimate per sample.
- `ngsdose adjust` — residualises estimates on coverage PCs from
  [NGS-PCA](https://github.com/jlanej/NGS-PCA). By default the regression is on log values and
  the cohort's mean log (the geometric mean) is restored, so estimates ≤ 0 become NA; `--no-log`
  fits on the natural scale and restores the arithmetic mean. It writes the result as `COL.adj`
  and reports on stderr the share of variance the PCs remove, raw, beside the share expected by
  chance (k/(n−1)) and adjusted for the number of regressors. Every input row is kept: a sample
  without the PCs gets NA and is named on stderr, and a column with fewer than k + 10 usable
  values is left NA rather than fitted. The PC basis is built from autosomal bins with
  repeats, segmental duplications and low-mappability sequence excluded, so it is disjoint from
  every class measured here: it can absorb library and batch structure but cannot absorb the
  dosage itself. It *can* absorb ancestry, so an association analysis should carry genetic PCs
  separately, and the variance removed by coverage PCs should be reported, not hidden.
- `ngsdose cohort` also reports **control PCs**: principal components of the 800 control regions'
  log(observed/expected) after the GC model (all of them, including regions flagged by control
  QC). Before the SVD each sample's row is median-centred, then each region's column is
  median-centred across samples, and the entries are clipped to ±0.5 in natural log (a factor of
  about 1.65). The clip caps the weight of a single sample's deleted or amplified region, so such a
  region cannot dominate a component; the Marchenko–Pastur count below is taken from the singular
  values of this clipped matrix. Samples without control residuals (`--no-control-qc`) or with
  another number of regions get NA and are named. It is a small internal version of the same
  idea, for cohorts on which NGS-PCA has not been run; `ngsdose adjust` uses them when no
  `--pcs` file is given. Like the coverage PCs they are computed on sequence disjoint from every class.
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
  on the 1000 Genomes spectrum it gives 38–40 components whether 100, 150 or 200 singular values
  are kept (45–49 without the margin; 34, or 37 without the margin, if only 80 are kept); the
  fitted edge is 1.5–1.8 times as broad as equal-variance noise would make it. `ngsdose cohort`
  makes the same choice for its control PCs, records it (`ctrlPC_mp`), and writes twice as many,
  at least 20 but no more than one per five samples (none for fewer than 10 samples with control
  residuals), so that `pcsweep` can look beyond the choice. Two limits, both measured:
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
- **Whether that number is right** is an empirical question, and the cohort carries what it takes
  to answer it: sequence of known copy number in every sample, and trios. `ngsdose pcsweep`
  regresses out 0, 1, 2, … PCs and reports, for each number, the error of the known truths
  (held-out autosomal, chrX and chrY by sex, the distal junction) and the transmission reliability
  of the classes with its paired difference from no adjustment. Everything is cross-validated —
  residual variance falls with every regressor whether it means anything or not: forty random PCs
  "explain" a third of pure noise in 120 samples, and out of fold make its SD about a quarter
  worse (1.27× over 200 simulations, 10-fold) — and the recommendation follows the
  one-standard-error rule: the fewest PCs that do as well as the best number, to within the
  sampling error of "best". A known truth says when adjustment has stopped removing noise;
  reliability says when it has started removing signal. They need not agree with the spectrum, and
  where they do not, they win.
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
the log scale. `ngsdose trios` therefore centres values within population (labels from the
pedigree or from `--population FILE`; it warns when the trios' samples carry neither at least
two labelled populations of 5 or more nor one label shared by all, and records
`population_centred`), works in copies, and
reports all three estimators, the spousal correlation and a permuted-family null. All three are
inflated by error shared within a family — a trio libraried together — and the spousal
correlation after population centring is the direct test for that, since spouses share no
dosage by descent. Trios bound reliability from above; technical replicates measure it. Run per
estimator, the table answers which estimator to use: the uncorrected 18S ratio, `cn_anchor`, the
calibrated `cn`, each before and after coverage-PC adjustment (NGS-DOSE-1000G's `pipeline/02_cohort.sh`).
Confidence intervals come from resampling families (trios that share anyone, siblings or three
generations, form one family), which carries the uncertainty of ρ into R;
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
Genomes 30× release nearly every child (597 of the 602 fully sequenced trios) is among the 698 related genomes
sequenced after the original 2,504, where most parents (91%) are, so generation and batch go
together: a batch that reads a
quantity on a scale k multiplies b, and R, by k. The cohort page therefore reports s, the
children's level against their parents', and R with the children first rescaled to their
parents' spread (b/s for b). A batch's scale moves R and leaves the rescaled R alone; new
variation arising in the children does the reverse; the two bracket the reliability. On the
602 trios of the published page they agree within a few hundredths for most metrics and part for
HSat2 (s = 1.08, slope 1.03: R 1.03 from the slope, 0.96 rescaled) and ACRO (s = 0.96: 0.91 and
0.95); the 45S reads 0.95 and 0.94.
`trios.by_sex` asks whether the four parent–child pairings by sex differ, referring Cochran's Q
on Fisher's z to its distribution under shuffled children's sexes and swapped parental roles:
the chi-square reference is far too liberal for skewed quantities (on the 602 trios: the EBV
load's pairings p = 0.0005 by chi-square, 0.09 by the shuffles; the distal junction's 1 × 10⁻⁷ and
0.27; the 45S's 0.04 and 0.14). A contrast between the father-to-son and father-to-daughter
correlations counts as sex linkage on the page only when this shuffle test agrees (p ≤ 0.01): at
602 trios the 5S array's contrast reaches 3.2 standard errors on Fisher's scale with a shuffle p of
0.14, and it is on chromosome 1.

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
| fetch, 8 threads | local CRAM | 2–3 s (cohort: 5.4 s, see below) | ~15 s | ~0.51 GB (about 2,010 of 79,637 slices with the current sinks, TEL included; chrM and chrEBV about 10 MB of it, chrY about 5 MB) |
| fetch, 16 threads | `https://1000genomes.s3.amazonaws.com/…` from a home connection | 60–100 s | ~20 s | ~0.51 GB |
| estimate | counts JSON | ~1 s | | 70–340 kB |

The fetch timings were measured on one genome before the TEL sinks were added (0.47 GB, about
1,920 slices, without them); the data read is the current plan's slices in the CRAM index. In the
1000 Genomes run, fetches of staged local CRAMs with the 80-interval sinks (no TEL) took a median
5.4 s (10th–90th percentile 4.2–7.7 s; `elapsed_sec` of the 1,748 fetch files counted by
2026-09-25).

In fetch mode the 3,202-sample cohort is therefore 53–89 hours of single-stream time (60–100 s
per genome) and never needs a CRAM on disk; scanned whole it is about 750 CPU-hours and 48 TB of
transfer. The 1000 Genomes run does both on each staged CRAM, a few files at a time, and deletes
the CRAM afterwards (NGS-DOSE-1000G's `pipeline/01_stage_and_dose.sh`), so that fetch can be
checked against scan in every genome. Peak memory is 0.6 GB (fetch) and 1.3 GB (scan) at 8
threads.

## 12. What it does not do, and known limits

- **Absolute scale of rDNA** rests on the anchors: windows on which three Illumina chemistries
  agree after their own GC corrections. That is evidence, not proof, of being unbiased. The one
  orthogonal check so far is ddPCR on twelve lymphoblastoid lines (Potapova et al. 2025; the
  results repository's `assembly_rdna` study, `python -m report --ddpcr` in NGS-DOSE-1000G):
  NGS-DOSE reads about 0.96× the assay (median ratio 0.957, r = 0.95, n = 12; NGS-DOSE-1000G
  `docs/data/ddpcr.tsv`), 0.97× on the nine 1000 Genomes lines alone. Relative dosage within a
  library type does not depend on it.
- **A new chemistry** (DNBSEQ, Element, Ultima, a future Illumina) has unknown dropout zones.
  Window efficiencies must be re-learned on it, and its absolute level is provisional until it
  has been compared with a known one on the same samples.
- **Cell-line DNA.** Late-replicating sequence is probably under-represented in DNA from cycling
  cultures (section 8), by an amount that differs between cultures; rDNA and satellites are
  mostly late-replicating. Blood DNA should be largely free of this. The cross-library
  agreement quoted here (3–5%) includes it, and real drift between cultures as well.
- **Units that are not the consensus.** Truncated or rearranged units contribute only the
  windows they contain; the window profile absorbs a population-average of that into a_w.
- **Sinks are aligner- and reference-specific.** The positional classes' 80 intervals (rDNA45S,
  rDNA5S, DJ) were learned from two NYGC CRAMs (bwa-mem 0.7.15, ALT-aware, GRCh38 analysis set; a
  CEU female and an ACB male), and TEL's 63 from 372 NYGC whole-file scans (section 5). For rDNA
  and DJ, two other bwa-mem pipelines on hs38DH-type references agree: a bwakit-0.7.12 + postalt
  probe found 9 of 1.5 million 45S reads outside the sinks, on unplaced, decoy or unmapped reads,
  and a whole-file scan of a Google HG002 BAM (bwa 0.7.17) had ≥ 99.7% of each class inside them.
  TEL was not checked on another pipeline. Another aligner is another matter. The public DRAGEN
  re-analyses of the 1000 Genomes CRAMs were fetched for one genome (HG00096) with the shipped
  bundle. Under DRAGEN 3.7.6 (alt-aware graph hash, the family UK Biobank and All of Us use at
  3.7.8, which is an inference: no header of theirs has been read) the sinks held 99.96% of the
  45S and 5S reads of the NYGC scan and 98.9% of DJ's. Under DRAGEN 4.2.7 (alt-masked) they held
  9.7% of 45S and 43% of DJ: most of those reads are fully unmapped there, and the 4.2.7 unmapped
  bin also holds 64–90% of HSat1A, β-satellite, ACRO, HSat1B and TEL (but almost no HSat2 or
  α-satellite HORs). With `--unmapped`, a DRAGEN 4.4.7 fetch recovered 99.7 / 99.9 / 99.6% of 45S
  / 5S / DJ, 72 / 38 / 53% of it from the unmapped bin. A fetch of such a file is refused from 0.1.1 on
  unless `--allow-missing-sinks` (two DJ sinks are on hs38d1 decoy contigs, absent from a no-alt
  reference), which leaves those intervals out, records them in `sinks_skipped` and warns;
  `ngsdose estimate` reports those classes as NaN either way, so a pipeline fetching no-alt or
  DRAGEN inputs passes the flag only for classes it can do without, or plans with `--classes`. The counts' `pipeline` record
  (`@PG` lines, `@SQ` hash) says which pipeline a file came from. DRAGEN, or a different decoy
  set, needs its own scanned subset and `ngsdose sinks` (`--classes TEL ...` for compositional
  classes).
- **Biobank scale: which classes need a scan** (status as of 2026-09-28, the cohort run complete). A biobank scans a
  subset whole, budgeted here at 0.1–1% of CRAMs (about 500–5,000 at UK Biobank), and that
  subset is where each pipeline's sinks are learned and checked. For the rest, on an NYGC-like
  bwa-mem pipeline:
  - *always fetched*: the controls, the truth regions, chrM and chrEBV are reference regions and
    part of every fetch (identical counts in fetch and scan);
  - *fetched through the shipped sinks*: rDNA45S, rDNA5S and DJ (DJ at a steady 99.6–99.8% of its
    scan count, an offset a scanned subset can calibrate), and TEL when its panel is loaded (its
    sinks ship, but no fetch through them has yet been compared with its scan);
  - *fetchable in principle (sinks learned and checked on scans; no fetch compared yet)*
    (section 5; an experimental NYGC set learned from 100 scans is in `resources/experimental/`,
    not in the bundle): HSat1A, HSat2, HSat3, α-satellite HORs (60 Mb of intervals),
    β-satellite, ACRO, SST1, CER and SATR;
  - *fetchable with a lower capture*, on the same footing: HSat1B (all its sinks hold 97.17% at
    the 10th percentile of held-out scans), from learned sinks plus the unmapped bin (`count --unmapped`; htslib region `*`, 0.10–0.26% of primary reads in the NYGC
    CRAMs, 1st–99th percentile);
  - *scan-only*: the 83 candidate classes, until sinks are learned for them from scans that load
    their panels.

  No fetch through the satellite sinks has been compared with its scan, so the satellite (and
  TEL) figures rest on scan placements.

  What is fetched is then a matter of bytes (section 4). On the NYGC CRAMs the controls with the
  bundle's four classes read 3.33% of a file; the lite control set and capture targets (by the
  held-out statistics that ship with the sinks) read less, and the satellite families several
  times more ([fetch_examples.md](fetch_examples.md));
  `ngsdose fetchplan` chooses between them by name, preset or budget. Each biobank pipeline has to
  cost its own files: slice sizes and read placement are properties of the aligner and the writer.

  Under DRAGEN 3.7.x (UK Biobank, All of Us) the positional sinks carried over in the one genome
  tested, with DJ's offset larger (98.9%); satellites and TEL are unmeasured there. Under DRAGEN
  4.x the unmapped bin is needed for 45S, DJ, TEL and most satellite families other than HSat2
  and the HORs (and for 5S under 4.4.7, not 4.2.7), with sinks for what stays mapped. TOPMed
  and CCDG used the same functional-equivalence bwa-mem pipeline and analysis-set reference as
  NYGC, so the shipped sinks are expected to transfer to them, but a scanned subset still has to
  confirm it.
- **GRCh38 only**, chr-prefixed names. A CHM13 or GRCh37 bundle is a rebuild of the controls and
  sinks; panels are reference-independent.
- **Dispersed families are experimental** (section 5). All ten are compared with HPRC release-2
  assemblies of 200 cohort members (the NGS-DOSE-1000G page, `docs/data/satellites_hprc.tsv`; a
  class is compared only where gap-annotated arrays are at most 2% of it). The median
  estimate/assembly ratio is 1.02 for α-satellite HORs, 0.94 for HSat3 (162 samples), and 0.91 and
  0.86 for HSat1A and HSat1B. β-satellite reads at 0.74, ACRO at 0.79 and CER at 0.43. These are
  stable under-reads, set by their k-mer recall, that track the assemblies across people (r
  0.91–0.95). HSat2's median ratio is 1.08 in the 79 samples compared, but it tracks the assembly
  poorly from person to person (r = 0.39), so its mass is not yet confirmed. SST1 and SATR (0.26, 0.09) cannot be judged in absolute terms, because HPRC annotates them
  several times more generously than the CHM13 annotation the panels came from. None of the ten
  has sinks in the bundle; the experimental NYGC set (section 5) has been checked on scan
  placements only. There is also a telomeric-repeat class that is a relative
  measure only (exact 31-mers are less tolerant of sequencing error than TelSeq's hexamer count;
  for each class the counts keep an eleven-bin histogram (`hit_frac`) of the share of each
  assigned read's k-mers that hit the panel, in tenths, so that a stricter share threshold can be
  chosen later; reads below the scan-time minimum of 4 panel k-mers are not recorded, so the
  threshold cannot be lowered afterwards).
- **Low depth** costs precision, not accuracy: the anchor estimate is unbiased down to 0.4×
  (1.004 ± 0.007 of the full-depth value over 20 subsamples; section 15), and a 0.75× subsample
  of NA12878 returns 499 against 504 from the full data. The GC curve itself gets noisy
  (`gc_curve_max_se`), which the GC-extreme features feel first.
- **Trios.** Section 10 now runs on the cohort. The results page (NGS-DOSE-1000G) estimates transmission from all
  602 complete trios (the pedigree's 603rd names a parent the release never sequenced): 45S 0.95
  (0.86–1.04), 5S 1.02 (0.88–1.15), the distal junction 0.78, the satellite families 0.44–1.03, the
  negative controls −0.31 to 0.15 (insert size's −0.31 is the parents' shared batch, ρ = 0.20,
  subtracted by the spousal correction from a slope of zero). Nearly every child is in the 698-genome
  batch (597 of the 602 fully sequenced trios) and most parents (91%) in the original 2,504, so
  generation and batch largely go together. A batch that reads on a different scale moves R,
  which is why the page also reports s and R with the children rescaled to their parents'
  spread; these bracket the reliability rather than remove the batch effect.
- **Fetch costs are estimates from the index.** `ngsdose.cost` counts every slice whose span
  overlaps an interval as read whole, which is what the engine decodes; it does not model the
  bytes htslib reads ahead, retries, or the index requests. It reads CRAM indexes only (no
  `.bai`/`.csi`). A capture target's expected capture comes from the scans behind the statistics;
  kept per byte it is a lower bound (section 4). The costs and the experimental sinks are those of NYGC bwa-mem CRAMs of the GRCh38 analysis set, and
  no fetch through the satellite or TEL sinks, or through a sub-option, has yet been compared with
  a scan of the same genome.
- **Candidate classes and experimental units.** The 83 candidates have been scanned in one genome
  only (NA12878), and their sinks, capture and costs are unknown until the cohort's scans carry
  them. Estimates of positional candidates are single-sample all-window estimates without an
  anchor or calibration, and may carry a level offset of several per cent (section 5); `ngsdose
  cohort` does not calibrate them (their status is not `ok`). Loading candidate
  panels part way through a cohort changes `panel_sha256` in the counts files. Scans with and
  without the candidates can still be pooled, because one set of panels holds the other;
  `ngsdose sinks --allow-mixed-panels` is needed only when two scans loaded different candidate
  subsets, neither of which holds the other.
- **The lite control set** keeps 45S copy number within 0.3% (SD) of the full set, 5S within 0.67%
  and the satellites within 0.78%, except HSat1B (2.3%), but loses power in
  the aneuploidy test (12 of 32 flagged genomes missed) and widens HSat1B's spread (section 6). Its
  numbers come from rebuilt GC tables, not yet from fetches of whole CRAMs with the lite FASTA.

## 13. What comes next

1. All 3,202 samples of the 1000 Genomes 30× cohort (NGS-DOSE-1000G's `pipeline/`) **scanned
   whole**, from staged CRAMs, with the satellite and telomere panels loaded, and fetched as well,
   a few seconds more on a local file (its `01b_dose_sample.sh`). Done on 2026-09-28: all 3,202
   genomes scanned and fetched, 1,748 with engine build fae1124 and 1,454 with 7772e32, builds
   whose counts NGS-DOSE's tests hold byte-identical; the first 1,748 fetches hold the bundle's
   positional classes with the sinks as they were before TEL (sha256 9dd52ba1…), the other 1,454
   the telomere panel through the current sinks as well. The results page (see "The running
   record" below) covers every genome and 602 complete trios, including fetch against scan, the
   Hall et al. and ddPCR comparisons, and satellites against HPRC assemblies in 200 samples. The scan is the full-accuracy mode and the
   reference for everything cheaper: placement-independent counts, the satellite families
   against an assembly truth for the 200 samples that have HPRC assemblies (its
   `04_hprc_satellites.sh`), sinks
   evaluated and re-learned across 26 populations and both sexes, and fetch / scan for every
   sample (its `03_compare_modes.sh`). From the cohort: the per-estimator reliability table with
   paired comparisons and its negative-control traits (its `02_cohort.sh`), efficiencies at
   n = 3,202, adjustment by NGS-PCA's coverage PCs and by the internal control PCs, DJ and chrY
   as cohort-wide accuracy checks, comparison with Hall et al.'s table sample by sample.
   - *How far can it be simplified for a biobank?* Three levels, each judged against the scan
     in the same people, with transmission reliability (the paired bootstrap of section 10) as
     the arbiter rather than correlation with the full estimate, which shares its errors.
     (i) **Targeted fetch**, which exists: about half a gigabyte and a minute per genome, no file
     staged; the cohort says what it loses, sample by sample, and the sinks re-learned on the 1-kb
     grid say how much smaller the retrieval can be made. For the satellite families the scans
     already say what sinks would hold (section 5), and the sinks file learned from 100 cohort
     scans now exists, kept apart from the bundle's `sinks.bed` so that fetches without the
     satellite panel neither change their `sinks_sha256` nor read several times more
     (`resources/experimental/sinks.satellites.bed`); what remains is real fetches of about a
     hundred genomes through it, compared with their scans (`--unmapped` for HSat1B). (ii)
     **A depth proxy from bins a cohort already has.** NGS-PCA's mosdepth run (1-kb bins, all
     contigs, no MAPQ filter, duplicate-flagged reads excluded) is kept for this cohort. In
     NA12878, 99.8% of the aligned 45S reads fall in 273 of those bins (chr21's rDNA models,
     KI270733, GL000220); the engine's census of those bins reproduces a duplicate-excluded depth
     computed independently to 3%. What a proxy gives up is known in kind — no fragment-GC model,
     no window calibration, no anchor, and a duplicate filter that removes 5.5% of rDNA reads but
     10.8% of single-copy reads in that sample, a different gap in every sample — and the scans,
     the kept mosdepth outputs and the census are what it takes to put numbers on each. Within one
     library type and one pipeline several of those terms are constants, which is the case for
     trying. (iii) **Less of everything**: fewer controls, smaller sinks, lower depth, all of
     which can be tried on the counts files alone, because they hold positions and not summaries.
     The first two now have tools and first numbers (sections 4 and 6): the lite control set and
     capture targets of 0.995 take the controls and the bundle's sinks from 3.33% to 2.39% of a
     NYGC 1000 Genomes CRAM ([fetch_examples.md](fetch_examples.md), examples 2 and 13), and fetching about ten genomes with
     both control sets on the cluster would confirm the lite set on whole files. None of this
     needs an answer before the run; it needs the run to keep what the answers will be computed
     from, and it does.
   - *A Mendelian test on integer states.* DJ windows carry inherited ±1-copy steps (section 8).
     Across 602 trios, a step present in a child and in neither parent is either a de novo
     event or an error, and a step in a parent is transmitted half the time: that calibrates
     the k-mer path at single-copy resolution, far more sharply than any regression on totals.
     Done, with the calls of section 7: 553 of the 556 trios with three settled calls are
     Mendelian at every position of the core, two children carry an event neither parent has,
     and steps pass to fewer than half of the children, a father's loss least of all, which is
     now the open question.
   - *The S-phase hypothesis.* If late-replicating sequence is under-represented in DNA from
     cycling cultures, DJ, female X and the leading control PC should move together across the
     cohort, and adjustment should tighten DJ around 10. DJ and female X do not, so far
     (section 8); the control-PC test remains.
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
3. Sinks for DRAGEN-aligned CRAMs (UK Biobank, All of Us). One genome per DRAGEN version has
   been fetched so far (section 12); the public DRAGEN re-analyses of the 1000 Genomes CRAMs,
   trios included, are where DRAGEN sinks can be learned from scans and checked, with
   `--unmapped` for DRAGEN 4.x, before any biobank access.
4. An orthogonal rDNA calibration wider than the twelve ddPCR lines.
5. **New classes through the scans still to run.** About 1,450 genomes of the cohort are still to
   be scanned. Loading the final candidate panels in them (NGS-DOSE-1000G's `CANDIDATE_PANELS`;
   fae1124 loads them unchanged, with identical shipped counts) gives every candidate its scans: 30 are
   enough for `pipeline/06_learn_sinks.sh` to learn sinks from half of them and check the rest,
   and every further scan is another held-out check. A class whose capture holds gets sinks and
   statistics in its menu row (status `experimental`, from `candidate`) and can be fetched; one
   that does not stays scan-only. Several candidates have per-sample truth
   among genomes not yet counted (the Schaap 2013 HapMap macrosatellite sizes, the 12 CEU lines in
   which Kojima 2021 found SMRV, Telford 2018's HHV-6A carrier NA18999), so the rest of the run is
   also their first validation.
6. **The fetch menu's open ends.** Real fetches through the TEL and satellite sinks and the
   sub-options, compared with the scans of the same genomes (needs an image newer than fae1124);
   and the costs of a DRAGEN pipeline, from its own indexes and sinks learned
   from its own scans.

## 14. Relation to prior work

**rDNA and the distal junction from short reads.** Gibbons et al. (*Nat Commun* 2014; *PNAS*
2015) introduced rDNA dosage and mitochondrial abundance from WGS depth and reported concerted
5S–45S variation. Hall, Turner & Queitsch (*Sci Rep* 2021) estimated 18S, 28S and 5S copy number
in the high-coverage 1000 Genomes data analysed here (and in the SSC), showed how library and
centre drive such estimates, found no meaningful 5S–45S correlation, and attributed the earlier
concerted variation to mixed data quality; their per-sample table is the one this project is
checked against (section 2, finding 1). Rodriguez-Algarra, Evans & Rakyan (*Cell Genomics* 2024)
counted reads in GRCh38's 18S-analogue intervals against all reads on the numbered chromosomes,
validated this 18S Ratio against realignment (R = 0.97 in 94 genomes) and measured it in UK
Biobank; Sanger and Sanger Vanguard read higher than deCODE, and they adjust for centre.
`18S.flat` is the closest counterpart here. The same group (Rodriguez-Algarra et al., *Cell
Genomics* 2026) learned GRCh38 "rDNA-analogue regions" from where fully realigned reads of
monozygotic twins had first mapped, and from simulated reads, then extracted only those reads
from UK Biobank CRAMs to call rDNA variants. Raj et al. (medRxiv 2026) estimated 45S and 5S copy
number in 490,383 UK Biobank genomes. They located the GRCh38 loci where simulated T2T rDNA reads
map (five for 45S on chr21, chr22_KI270733v1_random and chrUn_GL000220v1; three on chr1 for 5S),
retrieved only the reads at those loci, and fitted per sample a Poisson model of fragment counts
as a function of fragment GC on 1,000 control loci matched in length to the rDNA loci, after
Benjamini & Speed; they regress out QC fields, of which sequencing provider is the strongest,
find 45S and 5S essentially uncorrelated, and report associations with metabolic disease,
adiposity and blood traits. That is the model of section 6, reached independently, and their
loci are this bundle's rDNA sinks. Rhie et al. (bioRxiv 2026, DJCounter) measure the distal
junction in 4,172 newborns and relatives, 490,416 UK Biobank genomes and all 3,202 genomes of
this cohort, with a 31-mer mode (k-mers present once in each of the five CHM13 junctions, the
rule behind this bundle's DJ core) and a GRCh38 mode that counts reads in a fixed target, joined
by a fitted linear transform; they also estimate 45S in a CHM13-based mode. They find one
junction lost (9 copies) in 2.8–3.4% of people and gained (11 or more) in 8.4–9.3%, read about 8
in karyotyped Robertsonian carrier lines, and see 1000 Genomes parents shifted slightly above
their children. CONKORD (Potapova et al., *Cell Genomics* 2025) counts 18S 31-mers in reads
against 31-mers of GC-matched windows and rescales on single-copy genes of similar GC, with ddPCR
validation. Benjamini & Speed (*NAR* 2012) is the source of the fragment-GC model.

**Long reads** do not yet replace a short-read measurement. In T2T-CHM13 three of the five rDNA
arrays are model sequences, because ultra-long nanopore reads could not order their units (Nurk
et al. 2022). Copy number and activity of individual arrays have been measured in several genomes
(Potapova et al. 2025). The HPRC release-2 assemblies of cohort members hold about a third of the
rDNA their measured copy number implies (median 36%, range 22–71% over 96 genomes), in pieces of
at most 1.7 Mb, where the average array would be about 2.2 Mb (NGS-DOSE-1000G: `python -m report
--censat`, `docs/data/rdna_hprc.tsv`).

**Assignment by sequence, and targeted retrieval.** Counting class k-mers in reads wherever the
aligner put them is established: QuicK-mer2 (Shen & Kidd 2020), run on 2,457 of the NYGC 1000
Genomes CRAMs; GeneToCN (Pajuste & Remm 2023), KILDA (Molitor et al. 2025; LPA KIV-2 in 2,459
genomes of this cohort), ctyper (Ma & Chaisson 2025) and danbing-tk (Lu & Chaisson 2021);
CONKORD and DJCounter above; and for satellites, HSat3/DYZ1 array size from subfamily 24-mers
(Altemose et al. 2014), array-specific 24-mers of DXZ1 and DYZ3 (Miga et al. 2014) and k-Seek's
simple satellites on 2,504 of the NYGC CRAMs (Said, Barbash & Clark 2024). Paralog-diagnostic
sequence goes back to singly unique nucleotides (Sudmant et al. 2010), and counting reads on one
unit rather than where the genome alignment put them to Lucotte et al. (2018). All of these read
the whole file. Targeted retrieval at biobank scale is established too: the rDNA loci of Raj et
al. and the analogue regions of Rodriguez-Algarra et al. (2026) — learned from realigned reads,
a sink in all but name — the GRCh38 target of Rhie et al., reads on the EBV decoy in 490,560 UK
Biobank and 245,394 All of Us genomes (Nyeo et al. 2026), and unmapped plus EBV-decoy reads
realigned to 31 viral genomes (Kamitaki et al. 2026); Chrisman et al. (2022) found the decoy
contig on which HHV-6 reads land. What is not in this literature is the two joined: one read-level
k-mer classifier used identically on a whole-file scan and on a retrieval of intervals learned
from scans, with each class's capture measured, so that fetch reproduces scan class by class.
For rDNA the sinks are the loci others found by simulation or BLAST; what learning them from
scans adds is the measured capture and its reach to classes nobody has fetched: in these NYGC
CRAMs, sinks learned from 30 scans hold at least 99.8% of nine satellite families in each of 200
held-out genomes (99.85% in two of three random draws; 59.6 Mb for the α-satellite HORs, 0.2–3.5
Mb for each of the others; HSat1B about 96.6%), though no satellite fetch has yet been compared
with its scan. Sinks belong to an aligner and a reference: in the one genome checked under DRAGEN
4.x with an alt-masked reference, most 45S and DJ reads, and 64–90% of the reads of HSat1A,
HSat1B, β-satellite, ACRO and the telomeric repeat, were left unmapped (almost none of HSat2 or the
α-satellite HORs), so each pipeline needs its own scanned subset (section 12).

**The rest of the estimator.**
- *GC.* The fragment-GC curve fitted per library on single-copy controls is Benjamini & Speed's
  and Raj et al.'s. QuicK-mer2, Parascopy (Prodanov & Bansal 2022), AmpliCoNE (Vegesna et al.
  2019) and Said et al. also fit a per-sample GC relation on single-copy sequence; CONKORD
  matches by GC instead, and ctyper and GeneToCN drop GC-extreme k-mers. Tessereau et al. (2014)
  noted that depth under-counts the GC-rich RNU2 unit. Not found elsewhere: window efficiencies
  shared by a library type, kept apart from the GC curve and pinned on anchor windows chosen
  where two chemistries agree (section 7).
- *Mask.* Measuring only informative positions has analogues: singly unique nucleotides, the
  unique k-mers of QuicK-mer2 and KILDA, AmpliCoNE's family-specific positions, Miga et al.'s
  per-site specificity of the X and Y arrays, and Nyeo et al.'s per-base mask of the EBV decoy,
  learned from cohort coverage and checked with synthetic reads. Here the mask is derived from
  the panel (section 3) and applied to observed and expected counts alike, per position-strand of
  a multi-copy unit.
- *Known-copy controls.* Sizing a known sequence in every sample by the class's own calculation
  is in Altemose et al. (a 49-kb control region, used to correct and to screen samples), and
  AmpliCoNE reads single-copy X-degenerate genes as a check. Here four such readouts span 0, 1, 2
  and 10 copies and stay out of the denominator. The 10 is a population mode, not a per-sample
  truth: by Rhie et al. about one person in eight carries 9 or 11 or more, so the DJ is a truth
  per sample only where its integer state is known (window profile and inheritance, section 8)
  and a cohort-level check otherwise.
- *Duplicates.* No precedent was found for finding 1 of section 2. Rhie et al. saw batch effects
  that depended on how the BAM was made when duplicate, secondary and supplementary alignments
  were counted in the GRCh38 DJ target, and removed them; their filtered modes read about 0.7
  copy above their k-mer mode, the direction finding 1 predicts, though they give no mechanism
  and filtered all three kinds at once. Here secondary and supplementary records are excluded
  too, but the duplicate flag is ignored in numerator and denominator alike, as Rodriguez-Algarra
  et al. (2024) did implicitly. QuicK-mer2 and Said et al. removed flagged reads on these same
  CRAMs, and AmpliCoNE removes duplicates; finding 1 applies to them.
- *Transmission.* Relatives have served as validation before: heritability, twins and replicates
  (Hall et al.), relative pairs within a centre (Rodriguez-Algarra et al. 2024), Mendelian
  concordance in the same 602 trios (Parascopy) and in 641 trios (ctyper), trios and twins in UK
  Biobank (Raj et al.). Choosing an adjustment by an outside criterion is also not new (Taub et
  al. 2022: Southern blot and age; Burren et al. 2024: SNP heritability). Not found: midparent
  reliability with the spousal correction and the generation-batch rescaling, compared between
  estimators by a paired bootstrap (section 10).
- *Batch.* Coverage PCs from mappable autosomal bins, with structural variants and CNVs
  excluded, as covariates of a repeat-derived measure are Taub et al.'s (TOPMed telomere length,
  about 150,000 bins of 1 kb, 200 PCs) and Burren et al.'s (UK Biobank, 300 PCs per batch of
  about 20,000 genomes); NGS-PCA's basis is of the same kind. Sequencing centre or provider as a
  covariate is in Rodriguez-Algarra et al. (2024) and Raj et al.

**Other classes and dosage regions.** TelSeq (Ding et al. 2014) counts reads carrying at least k
TTAGGG repeats against reads of 48–52% GC; `hit_frac` is kept so that TEL can be calibrated
against it. For `chrEBV.copies`, Mandage et al. (2017) give EBV genomes per cell in 1,753 1000
Genomes LCLs from a masked EBV reference, validated by qPCR (r² = 0.88), with passage explaining
18% of the variance against 72% between individuals — a per-sample comparator — and Nyeo et al.
(2026) the biobank-scale version. Read depth along chromosomes of 1000 Genomes LCLs reflects
replication timing (Koren et al. 2014), the precedent for the S-phase hypothesis of section 8.

**What is specific to NGS-DOSE**, then: one k-mer classifier on both the whole-file and the
targeted path, with sinks learned from scans and their capture measured per class and per
pipeline, so that fetch reproduces scan, including for dispersed satellite families; a
recoverability mask derived from the panel and applied per position-strand of each unit;
sequence-specific window efficiencies kept apart from the library GC curve and pinned on anchors
chosen where chemistries agree; known-copy readouts at 0, 1, 2 and (as a population mode) 10
copies through the class code paths in every sample; the duplicate-flag finding; and midparent
reliability with its spousal correction and generation rescaling as the criterion for choosing
between estimators and numbers of PCs. The fragment-GC model and coverage-PC adjustment are
borrowed. The transmission algebra is textbook quantitative genetics for an additively inherited
trait; no claim of novelty is made for it. rDNA dosage is not a trait with heritability one: a
cluster rearranges in more than 10% of meioses (Stults et al. 2008) and can change in culture,
which lowers the slope, while error shared within a family raises it (section 10).

References for this section (DOI; PMID where indexed):

- Altemose N et al. *PLoS Comput Biol* 2014;10:e1003628, doi:10.1371/journal.pcbi.1003628, PMID 24831296.
- Benjamini Y, Speed TP. *Nucleic Acids Res* 2012;40:e72, doi:10.1093/nar/gks001, PMID 22323520.
- Burren OS et al. *Nat Genet* 2024;56:1832–1840, doi:10.1038/s41588-024-01884-7, PMID 39192095.
- Chrisman BS et al. *Virol J* 2022;19:225, doi:10.1186/s12985-022-01941-9, PMID 36566197.
- Ding Z et al. *Nucleic Acids Res* 2014;42:e75, doi:10.1093/nar/gku181, PMID 24609383.
- Gibbons JG et al. *Nat Commun* 2014;5:4850, doi:10.1038/ncomms5850, PMID 25209200.
- Gibbons JG et al. *PNAS* 2015;112:2485–2490, doi:10.1073/pnas.1416878112, PMID 25583482.
- Hall AN, Turner TN, Queitsch C. *Sci Rep* 2021;11:449, doi:10.1038/s41598-020-80049-y, PMID 33432083.
- Kamitaki N et al. *Nature* 2026 (DNA virome), doi:10.1038/s41586-026-10288-y, PMID 41882355.
- Koren A et al. *Cell* 2014;159:1015–1026, doi:10.1016/j.cell.2014.10.025, PMID 25416942.
- Lu TY, Chaisson MJP. *Nat Commun* 2021;12:4250, doi:10.1038/s41467-021-24378-0, PMID 34253730.
- Lucotte EA et al. *Genetics* 2018;209:907–920, doi:10.1534/genetics.118.300826, PMID 29769284.
- Ma W, Chaisson MJP. *Nat Genet* 2025;57:2909–2919, doi:10.1038/s41588-025-02346-4, PMID 41107550.
- Mandage R et al. *PLoS One* 2017;12:e0179446, doi:10.1371/journal.pone.0179446, PMID 28654678.
- Miga KH et al. *Genome Res* 2014;24:697–707, doi:10.1101/gr.159624.113, PMID 24501022.
- Molitor C et al. *NAR Genom Bioinform* 2025;7:lqaf070, doi:10.1093/nargab/lqaf070, PMID 40453649.
- Nurk S et al. *Science* 2022;376:44–53, doi:10.1126/science.abj6987, PMID 35357919.
- Nyeo SS et al. *Nature* 2026, doi:10.1038/s41586-025-10020-2, PMID 41606327.
- Pajuste FD, Remm M. *Sci Rep* 2023;13:17765, doi:10.1038/s41598-023-44636-z, PMID 37853040.
- Potapova TA et al. *Cell Genomics* 2025;5:101031, doi:10.1016/j.xgen.2025.101031, PMID 41043432.
- Prodanov T, Bansal V. *Nat Commun* 2022;13:3221, doi:10.1038/s41467-022-30930-3, PMID 35680869.
- Raj A et al. medRxiv 2026-01-13 (preprint), doi:10.64898/2026.01.09.26343685.
- Rhie A et al. bioRxiv 2026-03-10 (preprint), doi:10.64898/2026.03.08.710242, PMID 41959047.
- Rodriguez-Algarra F, Evans DM, Rakyan VK. *Cell Genomics* 2024;4:100562, doi:10.1016/j.xgen.2024.100562, PMID 38749448.
- Rodriguez-Algarra F et al. *Cell Genomics* 2026;6:101213, doi:10.1016/j.xgen.2026.101213, PMID 41966685.
- Said I, Barbash DA, Clark AG. *Genome Biol Evol* 2024;16:evae153, doi:10.1093/gbe/evae153, PMID 39018452.
- Shen F, Kidd JM. *Genes* 2020;11:141, doi:10.3390/genes11020141, PMID 32013076.
- Stults DM et al. *Genome Res* 2008;18:13–18, doi:10.1101/gr.6858507, PMID 18025267.
- Sudmant PH et al. *Science* 2010;330:641–646, doi:10.1126/science.1197005, PMID 21030649.
- Taub MA et al. *Cell Genomics* 2022;2:100084, doi:10.1016/j.xgen.2021.100084, PMID 35530816.
- Tessereau C et al. *Nucleic Acids Res* 2014;42:9121–9130, doi:10.1093/nar/gku639, PMID 25034697.
- Vegesna R et al. *PLoS Genet* 2019;15:e1008369, doi:10.1371/journal.pgen.1008369, PMID 31525193.

## 15. Assumption audit before the cohort run

Before committing 3,202 genomes to it, every assumption that the counts files depend on — the one
layer that cannot be revised without re-reading the CRAMs — was listed and tested, along with
the modelling assumptions that were cheap to test. What was wrong is recorded as plainly as what
held. The rows from "The satellite families can only be scanned" on come from reviews made while
the run was under way: the first with 1,748 genomes counted, the second of the fetch-menu change
itself (2026-09-28). None changes a count the engine has written for the cohort (the 5′-end rule
for hard-clipped primary alignments, below, touches none of its bwa `-Y` alignments).

| assumption | test | result | consequence |
| --- | --- | --- | --- |
| The strided k-mer screen loses no class read ("lossless for min-hits ≥ stride", as the first engine's help text said) | full NA12878, stride 4 vs stride 1 | **false**: 0.04–0.06% of class reads lost (hits in short, broken runs) | replaced by an exact bitset prefilter, no slower; unit test with exactly such a read |
| A class's counts do not depend on what else is in the panel | bundle panel alone vs with the satellite panel | **false**: DJ lost 1.7% of its reads, because a read carrying β-satellite k-mers as well was called ambiguous | classes judged independently (section 4); identical counts now asserted in CI |
| Output is reproducible | same input, 1 vs 4 vs 8 threads | counts identical; placement list ordered by hash iteration on ties | total ordering; identical counts (every field but the elapsed time) for 1 and 4 threads, in scan and fetch mode, asserted in CI |
| Sinks learned from one CEU female transfer to other samples | whole-file scan of HG02258 (ACB, male) | hold: 99.94 / 99.98 / 99.74% of 45S / 5S / DJ reads inside them, against 99.94 / 99.98 / 99.77%; chrY takes 66 DJ reads. A third scan that the sinks never saw (HG01884, ACB, female): 99.96 / 99.99 / 99.75% | shipped positional sinks are the union of the first two scans (80 intervals; TEL's 63 were added on 2026-09-22 from 372 cohort scans, 143 in all); NGS-DOSE-1000G's `pipeline/03_compare_modes.sh` evaluates and re-learns them on every scanned sample: across 1,748 scans they hold ≥ 99.89 / 99.93 / 99.65% of 45S / 5S / DJ reads (the cohort's fetches so far used the 80-interval file) |
| … and to another alignment pipeline | probe of all 2,554 unplaced/decoy contigs and the unmapped bin of a bwakit 0.7.12 + postalt CRAM | hold: 9 of 1.5 M 45S reads outside | a bwa-mem pipeline only. DRAGEN, one genome per version (section 12): 3.7.6 holds for 45S and 5S and loses 1.1% of DJ; 4.x leaves most 45S and DJ reads unmapped, which `--unmapped` recovers |
| The estimate does not depend on depth | NA12878 subsampled, 20 replicates per depth | anchor estimate unbiased: 1.009 ± 0.005 at 1.1×, 1.004 ± 0.007 at 0.37×. **All-window estimate biased**: +2.4% and +6.5%, because the GC support, defined by the curve's precision, narrowed with depth and GC-rich low-efficiency windows dropped out | support defined by control positions (section 6): 0.998 ± 0.004 and 0.996 ± 0.005; asserted in CI |
| The robust denominator is neutral when nothing is aberrant | same series | its scale factor was not exactly 1 at low depth | now the *relative* change from dropping flagged regions: exactly 1 when none is |
| Analysis choices do not move the headline | NA12878: min k-mers 5–60, window 100–1000 bp, L 300–600 bp | 45S anchor 503.6–505.2; DJ 9.75–9.87; held-out autosomal 2.001–2.002. Features move more (18S 457–489 with the k-mer threshold, 5S 202–217 with L) | none; features are diagnostics, not headlines |
| The DJ deficit (9.5–9.9 rather than 10) is not structural | DJ in 10-kb windows | uniform within a sample to 3.5%, plus inherited one-copy steps (section 8) | deficit is global, consistent with the S-phase hypothesis, which the cohort has not borne out so far (section 8); steps become a cohort test |
| Trio estimators hold for skewed sizes and multiplicative error | simulation, log-normal sizes (CV 35%), 10% multiplicative error | recovered: 0.852 for a true 0.852 | in `ngsdose selftest` |
| 602 trios can rank estimators | simulation | separate CIs are ±0.09, too wide; the paired bootstrap resolves a gap of 0.08 | `ngsdose trios --compare-to` |
| A file that decodes is a whole file | fixture BAM and a 200 MB head of a real CRAM, end-of-file block removed | **false, silently**: both decode without error, and what a truncated coordinate-sorted GRCh38 analysis-set file loses first is the unmapped reads, then the HLA, decoy, chrEBV, alt and unplaced contigs (where much of the rDNA and DJ lands, e.g. chrUn_GL000220v1, chr22_KI270733v1_random), then chrM, chrY and chrX | the engine checks the EOF marker before reading (local or HTTPS) and refuses; `--allow-truncated` overrides; recorded in every counts file; asserted in CI. A stream that cannot be checked up front (a pipe, a server without range requests) is checked at its end in scan mode (not a CRAM decoded with more than one thread), and `ngsdose estimate` warns on `unchecked` |
| The file is aligned to the bundle's build | - | nothing checked it: hg19 has the same contig names, and the controls and sinks are coordinates | contig lengths are recorded in the counts file and `estimate` refuses a mismatch; asserted in CI |
| A failed connection fails | the pilot re-count, left running overnight on a laptop that slept | **false**: two fetch processes sat for eight hours on connections that were never going to answer; per-interval retries never fire, because nothing errs | stall watchdog in the engine (exit 75 after `--stall-timeout` s without a record), retry loop in the cohort script (NGS-DOSE-1000G's `pipeline/01_count.sh`); asserted in CI with a FIFO that nobody writes to. A remote input that still fails after `--retries` attempts also exits 75, and opens are retried too |
| A killed job leaves nothing a resumed run would mistake for a result | - | the engine wrote its output in place | output written under a temporary name, synced, renamed |
| A class in the counts file reaches the estimate table | merged-panel scan through `ngsdose estimate` | **false**: classes absent from the bundle's own panel were skipped, so the satellite columns of a cohort scan would have been empty | compositional classes need nothing from the bundle and are always estimated; asserted in CI |
| The bottom of the scale is zero | 40 chrY regions in a female | hold: 0.003 copies in NA12878 (after excluding two candidate runs that attract X-derived reads) | chrY added as a known truth (section 8) |
| Every position of a region has a fragment-GC window | chrM region-based estimate against the scan's whole-contig read tally | **false for chrM**: GRCh38 carries an `N` at chrM:3,107; windows over it are in no GC table but reads starting there were in the observed count, and the first chrM estimate was 18% high (978 against 830) | region moved off the `N`; the bundle build refuses such regions; expected counts are scaled to all position-strands, so a custom bundle cannot repeat it; asserted in CI |
| Non-transmitted variance is measurement error | - | untestable without a handle on the state of the culture | chrM and chrEBV dosage in both modes (NA12878: 850 and 71 copies per cell), and as negative-control traits in the transmission table (section 10) |
| The cohort script works | a mock cohort: seeded subsamples of the NA12878 fixture named after real 1000 Genomes trio members (66 samples, 22 trios), through the `ngsdose` commands that the cohort scripts (now NGS-DOSE-1000G's `pipeline/`) call, with the real pedigree and NGS-PCA's real PC file | it ran - for the first time; the twelve-sample pilot is too small for control PCs, PC adjustment or the paired bootstrap. With one person sixty-six times over there is nothing to transmit, and every reliability came out at zero within its interval. It also showed `adjust` reporting that twenty PCs "explain" 30% of pure noise at n = 66 | a sixty-sample version is a CI test of the whole cohort layer (`tests/test_cohort_flow.py`); `adjust` reports the chance expectation k/(n-1) and the adjusted R² beside the raw figure |
| Counts survive a bundle revision | pilot counts made before the chrY and dosage regions existed, estimated with the new bundle | refused: regions were matched to the bundle row for row | matched by name: the control set must be identical, other regions may be missing (section 4); asserted in CI |
| The manifest is right | every URL of NGS-DOSE-1000G's `pipeline/00_setup.sh` manifest requested (HEAD) | three sample names in the 1000 Genomes sequence index end in a space, which put a space into three URLs and file names | names trimmed, file names taken from the index's path column, the manifest validated; all 6,404 CRAM and index URLs answer 200 |
| The cohort run keeps what a cheaper method would be judged on | what a depth proxy needs: the scan estimate, the cohort's kept mosdepth bins, and the composition of those bins | placements were on a 10-kb grid and held class reads only: enough to learn sinks, not to say how much of a mosdepth bin is rDNA | placements on mosdepth's 1-kb grid with a census of every read in the bin (all, duplicate-flagged); against an independent duplicate-excluded depth over chr21's rDNA bins the census is within 3%; scan and fetch agree bin for bin inside the sinks, asserted in CI |
| The Marchenko–Pastur law can be fitted to a coverage spectrum to choose the number of PCs | simulated noise with unequal variances (rows ±40%, columns ±30%) and planted components; NGS-PCA's 1000 Genomes spectrum truncated at 100, 150, 200 values | **false**: the median- or quantile-matched fit called 34 components for 7 planted and 148 for 5, and 59 / 66 / 80 on the real spectrum depending on the truncation | the edge is fitted from its universal square-root shape instead, with a 1% margin (section 9): exact on those simulations, none in that noise, 38–40 on the real spectrum at every truncation. It still leans high when bins have heavy-tailed variances (a few components in noise alone), and components beyond the leading 20–25 are not reproducible between SVD runs - so `ngsdose pcsweep` lets known truths and transmission decide; all asserted in CI |
| A known truth's error under adjustment can be read off a regression of the estimate | review of `pcsweep`: a two-valued truth (chrX: 1 or 2), PCs unrelated to it, n = 3,200 | **false**: the estimate carries the variance of the truth itself (SD of log truth 0.35), every regressor adds 0.35·√(k/n) of estimation noise out of fold, and the reported error rose 0.010 → 0.040 at 46 PCs - adjustment would always have looked harmful for chrX | the error, log(estimate / truth), is what is regressed; asserted in CI, including a PC that follows sex. (Same review: one `N_PC` knob for two PC sets of different size killed NGS-DOSE-1000G's `pipeline/02_cohort.sh` at its last step - now `N_PC` and `N_CTRL_PC`, and a number beyond what a table holds is clamped with a warning; the one-standard-error band used 1.25 instead of 1.65 for a MAD-based SD.) |
| The experimental satellite panel measures array mass | whole-file scans of two 1000 Genomes samples that have HPRC release-2 assemblies (HG02258, ACB male; HG01884, ACB female) | HSat3 0.97, 0.98 of the assembly; HSat1A 0.93, 0.97; α-satellite HORs 0.97, 1.03; HSat1B 0.90, 0.83; β-satellite 0.69, 0.76 - and its k-mer recall on CHM13 itself is 0.69. **My first reading of HSat2 (1.52, 1.73: "not usable") was wrong**: the comparison script dropped arrays annotated together with an assembly gap ("GAP,HSat2": 17 and 9 Mb), and an array with a gap in it is no truth | the script tallies gap-containing arrays separately and compares a class only where it has none. In the cohort run so far (96 HPRC samples; NGS-DOSE-1000G `docs/data/satellites_hprc.tsv`), 43 assemblies have no gap-containing HSat2 array at all. There NGS-DOSE reads a median 1.07 of the assembly, but it follows the assembly poorly across people (r = 0.18, SD of the log ratio 0.28). HSat2 is heritable, but the assemblies do not yet confirm that it measures HSat2 mass (NGS-DOSE-1000G `docs/EVIDENCE.md`). Each class's recall is measured and documented (`resources/build/panel_recall.py`). 200 cohort members have assemblies (NGS-DOSE-1000G `pipeline/04_hprc_satellites.sh`) |
| Every loaded class has sinks in fetch mode | fetch with a panel class the sinks BED names nowhere (a satellite or telomere panel against an older sinks file) | **false**: the class was counted only where its reads fell inside other intervals, and nothing said so | the engine refuses (since 645ae55; fae1124 still undercounts silently, so `fetchplan` refuses such a plan for it); `--allow-missing-sinks` continues and lists the class in `sinks_missing_classes`, and `ngsdose estimate` reports it as NaN with status `no_sinks_in_fetch` |
| A trio table survives a column without variance | `ngsdose trios` on a constant column (EBV in blood-derived DNA) and on one with too few complete trios | **false**: division by the parents' variance raised, and one column's error ended the run | NaN for what cannot be estimated, a note for what cannot be run, the other columns unaffected; the pedigree may also be a PLINK PED/FAM, a child/father/mother table or the tab-separated 1000 Genomes `.ped` with its multi-word header (a file without a header never supplies a population; `--population` can); a column the table lacks gets an NA row with the reason, and `pcsweep` leaves it out (2026-09-28) |
| A site that cannot let the engine read CRAMs can reproduce a fetch | - | the fetch intervals (controls padded by 600 bp, the sinks, merged) existed only inside the engine | `ngs-dose plan` writes them as BED, dropping contigs the input's header lacks (samtools passes over such a region and exits 0) and reporting the sink intervals it drops per class; the cut is counted in fetch mode, which reproduces the whole-file fetch exactly, except for the unmapped bin, which `-L` never outputs (`plan --unmapped` prints the samtools steps that add it; tested equal to a whole-file `fetch --unmapped`); a cut counted in scan mode is refused by `ngsdose sinks` (control ends above 5% of primary reads) |
| `ngs-dose plan` needs no input | the GRCh38 bundle without `-i` | **false**: every contig was filed under one id and the controls' overlap check faulted regions on different chromosomes against each other | an id per contig name |
| A column without a slope has no permutation p | a constant column with 10 or more trios | **false**: the null was all NaN, nothing exceeded NaN, and p came out 1/(n+1) | NaN when the slope is NaN |
| A pedigree header may start with `#` | `#kid dad mom ...`, PLINK's `#FID IID ...` | **false**: dropped as a comment, and the layout guessed from the data | a first line starting with `#` is a header when it names the columns |
| The satellite families can only be scanned (this document said so before this revision: "no sinks to fetch") | review of the cohort's scans: sinks learned with `ngsdose sinks --classes` from 30 scans, capture in 200 others | **false** for this pipeline: ≥ 99.8% of nine families in every held-out genome (≥ 99.85% in two of three draws), HSat1B ≥ 96.6% (part of it unmapped) | section 5 rewritten; no satellite sinks ship until a real satellite fetch has been compared with its scan (section 13) |
| Sinks learned from 40 scans held ≥ 99.87% of TEL in each of the other 332 (as first recorded here and in `bundle.json`) | TEL sinks learned from 40 of the 372 learning scans (the first 40, and ten random draws), each scored on the other 332 | **false**: the lowest capture in the other 332 is 99.42–99.59% across the eleven draws (99.48% for the first 40), the median about 99.8%; 99.87% is the median in-sample capture of the 372-scan set. For context, the shipped 372-scan set over all 1,748 cohort scans: median 99.87%, minimum 99.39%, 1st percentile 99.68% | the figures quoted are the measured ones (section 5) |
| A capture score says what a fetch retrieves | review of `ngsdose sinks`: a 10-kb placement bin counted as captured when only its start fell inside a sink | **false** at the margin: a fetch reads only the part of the bin inside the sink | a bin counts only when all of it (clipped at the contig end) lies inside; TEL capture moves by at most 5.5 × 10⁻⁴ over 1,748 scans, the positional classes not at all |
| A fetch that cannot read a sink says so | sinks on contigs the file's header lacks | **false**: the intervals were skipped with a note, the class undercounted, and a class with every interval skipped was not refused | skipped intervals recorded per class (`sinks_skipped`), a class that keeps none refused unless `--allow-missing-sinks` (engines since 645ae55), and from 0.1.1 on one that loses any interval, since the estimator would not measure it and the fetch would be spent for nothing (fae1124 neither refuses nor takes the flag), and `ngsdose estimate` reports such classes as NaN (section 4) |
| Control QC flags a whole-chromosome aneuploidy | HG00096 with one chromosome's control counts scaled by 0.5–1.5 | **false** for a full trisomy (scaled 1.45 and above) or monosomy (0.5): its regions were trimmed as outliers against the genome's median, out of their own chromosome's test | regions trimmed against their own chromosome's median, and the genome level taken without flagged chromosomes (section 3); on the 1,748 counted genomes one flag is added (HG03363, a chr11 gain over 9 of its 38 regions) and chr22, with 2 control regions, is reported as untestable |
| One estimate per sample reaches the cohort | the same genome's scan and fetch estimates given to `ngsdose cohort` together | **false**: rows with one sample id were merged | a repeated sample id is refused, and `ngsdose estimate` names each estimate after its input file |
| A foreign warning inside the trio analysis is shown | `ngsdose trios` on a column whose values overflow (a numpy RuntimeWarning) | **false, and worse**: the command never returned; the warning was re-shown inside the recording context, which recorded it again into the list being walked | re-shown after the context is left; asserted in CI under a deadline |
| A bare copy of the bundle's directory is a bundle | `NGSDOSE_RESOURCES` pointing at a copy of `resources/GRCh38` alone, and a fetch of a `core_tel xy_arrays` plan | **false, silently**: the sub-option definitions and the experimental sinks are looked for beside the bundle; absent, a fetch that read aSatHOR only at DXZ1 and DYZ3 was reported as the family's mass with status ok | `ngsdose estimate` says what it cannot find (`NGSDOSE_EXPERIMENTAL` names another place) and reports a fetch through any sinks BED but the bundle's own as unverified there; asserted in CI |
| The pipeline record guards the sinks | a fetch of a CRAM aligned by another program to the same contig set (DRAGEN 4.x on alt-masked hg38) | **false**: `pipeline` was written and compared with nothing, so such a fetch through the NYGC sinks read 45S low with status ok (section 12) | the bundle records the aligner and the @SQ hash its sinks were learned under (`sinks_learned_from.pipeline`); `ngsdose estimate` warns on either difference and records `sinks_pipeline` |
| Control PCs are computed on one set of regions | estimates whose control regions differ in order but not in number | **false, silently**: samples were grouped by the length of their residual vector, and the region hash `estimate` writes was read by nothing; three reversed vectors among twelve made PC1 the "which order" axis | grouped by `region_order_sha256` (estimates without it by their count); the others are named and left NA; asserted in CI |
| A sub-option verdict rests on evidence | a fetch through a BED with an unnamed row; a family with no placed bin | **false**: `subset_only` was the one place that did not read an unnamed row as every class's, and `all([])` judged a family with no placement subset-only | unnamed rows count for the family; no placement decides nothing (unverified); asserted in CI |
| `--fetch-sinks` names a file | a mistyped path | **false, silently**: the path was passed over, the parents stayed unverified, and the warning told the user to pass `--fetch-sinks` | a path that is not a file stops the run |
| The fetch menu goes with the bundle | `NGSDOSE_RESOURCES` set with no `fetch_menu.tsv` beside it (a release install: the tarball packed no menu) | **false**: `fetchplan` silently took the source checkout's menu, whose rows name the checkout's panels and sinks, so a plan mixed two bundles; from the release alone it could not run at all | the tarball carries the menu; `fetchplan` stops rather than take another bundle's; the release smoke test plans and fetches through it |
| A column the table lacks is that column's failure | `ngsdose trios -c var nosuch` | **false** since the fetch-menu change (it had held since 2026-09-24): the run stopped with no table and no JSON | an NA row with the reason; `pcsweep` leaves the column out; `adjust` still stops, since its output would lack what was asked for |
| Signed URLs stay out of the logs | a URL with a query string; a 503 answered with retries | **false**: the engine's own messages were redacted, but htslib's (`[E::hts_open_format] Failed to open file "…?X-Amz-Signature=…"`) printed the signature on every failed open | htslib's messages are turned off for an input or index URL with a query string, said once; asserted in CI against a local server |
| A fetch that loses part of a class's sinks says so before it is spent | one sink interval on a contig the file lacks | **false**: the engine warned and fetched, `ngsdose estimate` then reported the class NaN, and no flag covered the partial case on either side | one rule (0.1.1): a class that loses any interval is refused unless `--allow-missing-sinks`, which records the loss (`sinks_skipped`) as before |
| A rule that changes counts is versioned | the 5′ end of a soft clip behind a hard clip (`5H10S85M`) | **false**: 0.1.0 restored only an outermost soft clip, the fetch-menu build restores one behind a hard clip, and both wrote `engine_version` 0.1.0 (the cohort's bwa `-Y` alignments carry no hard-clipped primaries, so its counts are unchanged) | engine and package 0.1.1; the rule is stated in section 3 and pinned by a unit test |
| The plan's price is what the fetch reads | the fixture as a CRAM under a read-logging shim; 13 NYGC indexes priced per interval | **false**: `cost.py` priced the union of slices once, while the engine makes one indexed fetch per run of intervals and htslib decodes a slice again for every run overlapping it - the 1.67-MB fixture CRAM was read as 37 MB, and `core_tel` costs 2-3% more than priced in men and 11-13% in women (chrY's few sparse slices decoded once per truth region), `xy_arrays` up to 27% | priced as the engine reads (one fetch per merged run, a slice once per run), the floor with every slice once reported beside it (`cum_mb_floor`); every documented figure regenerated (`core_tel` 546.9 MB, 3.33%) |

Not tested, and the cohort run will not test them either: a chemistry other than Illumina's;
DRAGEN alignments beyond one genome per version (section 12); an orthogonal assay for the
absolute rDNA scale beyond the twelve ddPCR lines.
