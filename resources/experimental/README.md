# Experimental resources: dispersed sequence

A positional class (rDNA, the distal junction) has a unit, its reads land in a few places, and a
targeted fetch finds them. The satellite families have no unit, and their reads are spread over
centromere models, decoys and hundreds of other contigs (about 82 Mb of intervals for the ten
families together), but within one pipeline those places are stable. In the NYGC bwa-mem alignments
to the GRCh38 analysis set, sinks learned from 30 cohort scans
(`ngsdose sinks --classes HSat1A HSat1B ...`) held at least 99.8% of `HSat1A`, `HSat2`, `HSat3`,
`aSatHOR`, `bSat`, `ACRO`, `SST1`, `CER` and `SATR` in every one of 200 other genomes (≥ 99.85% in
two of three random draws of the 30 and the 200), in 0.2–3.5 Mb of intervals per family and 59.6 Mb (the centromere
models) for `aSatHOR`. `HSat1B` is the exception: about half of its reads sit on hs38d1 decoys, up
to 3% are fully unmapped, and its 12–13 Mb of sinks held only 96.6% of it in the worst genome
(median 98.8%; 96.8% and 98.8% in a second draw). Sinks depend on the aligner and the reference
(under DRAGEN 4.x with an alt-masked reference, 64–90% of the `HSat1A`, `HSat1B`, `bSat`, `ACRO` and
telomeric reads of the one genome checked were fully unmapped, but almost none of its `HSat2` or
`aSatHOR` reads), so each pipeline has to learn them from whole-file scans of its own. The bundle
does not ship satellite sinks. For the NYGC pipeline, `sinks.satellites.bed` in this directory
holds a set learned from 100 cohort scans and checked on 1,648 others
([below](#sinkssatellitesbed--where-nygc-bwa-mem-puts-the-satellite-reads)), but no fetch through it
has yet been compared with a scan of the same genome, so these panels are still measured by a
whole-file scan. That is why a cohort that is scanned once should be scanned with them loaded
(NGS-DOSE-1000G's pipeline does; `EXTRA_PANELS` in its `config.sh`). A cohort too large to scan
would scan a subset of each pipeline with them loaded, learn the satellite sinks there, check their
capture on held-out scans (`ngsdose sinks --evaluate`), and fetch the rest with those sinks and
`-p satellites.CHM13v2.k31.panel.tsv.gz` (`ngsdose fetchplan --preset satellites`, or a subset of
the families, writes the files for it). A scan classifies fully unmapped reads; a fetch reads them
only with `ngs-dose count --unmapped`, which a class like `HSat1B` needs. The telomeric repeat
already has its sinks: the aligner concentrates its reads at the chromosome ends, the bundle's
`sinks.bed` carries the intervals, and a fetch with `-p telomere.k31.panel.tsv.gz` measures it
(`FETCH_PANELS`); the cohort's fetches so far did not load it, so no `TEL` fetch has yet been
compared with its scan.

```bash
ngs-dose count -m scan -i sample.cram -T ref.fa -c ../GRCh38/controls.fa.gz -p ../GRCh38/panel.k31.tsv.gz \
    -p satellites.CHM13v2.k31.panel.tsv.gz -p telomere.k31.panel.tsv.gz -o sample.scan.json.gz
ngsdose estimate sample.scan.json.gz          # CLASS.mass_Mb: diploid sequence mass of each family
```

Loading them changes nothing about the bundle's classes (asserted in CI), costs no measurable
time (1 min 40 s for a 30× genome either way) and adds ~130 kB to a counts file (105 → 237 kB).
`resources/build/build_satellite_panel.sh` rebuilds both panels; `resources/build/panel_recall.py`
produces the recall column below.

## `satellites.CHM13v2.k31.panel.tsv.gz` — ten families, 1.13 M k-mers

From the T2T-CHM13v2.0 CenSat annotation. A k-mer is kept if it occurs at least ten times in the
family's CHM13 arrays (arrays of at least 2 kb), in none of the other nine families, and nowhere
in CHM13 outside CenSat-annotated satellite. K-mers that a family shares with satellites that have
no class of their own are kept. These include monomeric α, divergent α HORs, γ-satellite, HSat4,
rDNA and the other CenSat families, so reads from those satellites can be counted as the panelled
family. For example, making monomeric α a class would take 17% of `aSatHOR`'s k-mers with it (see
below). *Recall* is the share of 150-bp reads drawn from the family's own CHM13 arrays that
carry the four k-mers a read needs to be assigned: what the panel can see of the genome it was
built from, and so an upper bound on what it sees of anyone else's.

| class | what | CHM13 (haploid) | k-mers | recall |
| --- | --- | --- | --- | --- |
| `HSat1A` | human satellite 1A | 13.4 Mb | 89,892 | 99.8% |
| `HSat1B` | human satellite 1B, mostly Yq | 15.3 Mb | 78,360 | 99.2% |
| `HSat2` | human satellite 2 | 28.7 Mb | 126,147 | 99.0% |
| `HSat3` | human satellite 3 | 69.3 Mb | 425,310 | 97.3% |
| `aSatHOR` | α-satellite higher-order repeats, active and inactive | 70.3 Mb | 325,529 | 99.9% |
| `bSat` | β-satellite | 8.6 Mb | 52,314 | 69% |
| `ACRO` | ACRO1 composites of the acrocentric short arms | 1.5 Mb | 13,939 | 90% |
| `SST1` | SST1 arrays (where Robertsonian translocations break) | 0.5 Mb | 6,590 | 61% |
| `CER` | centromeric repeat | 1.1 Mb | 4,062 | 42% |
| `SATR` | SATR1/2 | 0.3 Mb | 4,363 | 50% |

`HSat1A`, `HSat1B`, `HSat3` and `aSatHOR` are measured; `HSat2` has full recall but the assemblies
do not confirm it (see below); `ACRO` nearly; `bSat`, `SST1`, `CER` and `SATR` are relative
measures — comparable between people, under-read in absolute terms by about their recall
(β-satellite came out at 0.69 and 0.76 of the assembly in the comparison below, and its recall is
0.69). Left out, because they cannot be measured this way or cost more than they give: gamma
satellite (13%), divergent α HORs (10%), HSat4 (four k-mers survive); and monomeric α (45%), which
as a class of its own takes 17% of `aSatHOR`'s k-mers with it, because a k-mer shared between two
classes is dropped from both.

### Against assemblies of the same people

HPRC release 2, CenSat annotation of both haplotypes summed, for two 1000 Genomes samples scanned
whole (NGS-DOSE-1000G's `pipeline/hprc_satellites.py`). Cells are assembly Mb / NGS-DOSE Mb (ratio).

| class | HG02258 (ACB, male) | HG01884 (ACB, female) |
| --- | --- | --- |
| `HSat1A` | 28.4 / 26.5 (0.93) | 24.9 / 24.2 (0.97) |
| `HSat1B` | 11.8 / 10.6 (0.90) | 2.2 / 1.8 (0.83) |
| `HSat2` | 44.0 + 17.2 in gap-containing arrays / 67.0 | 30.5 + 9.2 in gap-containing arrays / 52.8 |
| `HSat3` | 87.4 / 85.1 (0.97) | 71.8 / 70.3 (0.98) |
| `aSatHOR` | 135.2 / 130.8 (0.97) | 148.8 / 153.2 (1.03) |
| `bSat` | 15.6 / 10.7 (0.69) | 18.1 / 13.8 (0.76) |
| `ACRO` | 3.3 / 2.6 (0.79) | 3.8 / 3.0 (0.79) |
| `CER` | 1.8 / 0.8 (0.42) | 2.0 / 0.9 (0.43) |
| `SST1` | 2.9 / 0.7 | 2.4 / 0.7 |
| `SATR` | 2.7 / 0.2 | 2.8 / 0.3 |

- HSat3, HSat1A and the HORs are within 7% of the assembly in both people; HSat1B within 10% in the
  male and 17% in the female, who has 2 Mb of it.
- The relative measures are under-read by a *stable* factor: `ACRO` 0.79 in both, `CER` 0.42 and
  0.43 (its recall), `bSat` 0.69 and 0.76 (its recall).
- **HSat2 cannot be judged here.** Both assemblies have gaps inside their HSat2 arrays, and an array
  with a gap in it is a lower bound, not a truth: the estimates are 1.52 and 1.73 of the spanned
  arrays, 1.09 and 1.33 if the gapped ones are counted at their annotated size. (An earlier version
  of the comparison dropped gap-annotated arrays silently and read this as HSat2 being
  over-estimated by 50-75%.) The script tallies such arrays separately and leaves a sample out of a
  class's comparison when they are more than 2% of what is annotated.
- `SST1` and `SATR` are annotated several times more generously in the HPRC assemblies than in the
  CHM13 annotation the panel was built from (2.4-2.9 Mb against 1.0 Mb diploid for SST1), so the
  absolute ratio means nothing; whether they track is a question for more samples.

Two samples say nothing about whether the estimates *track* the assemblies across people, which is
what association work needs. The cohort run in NGS-DOSE-1000G (`pipeline/04_hprc_satellites.sh`;
its `docs/EVIDENCE.md` section 7 and the cohort page's section 3.8) makes that comparison for the
cohort members with HPRC release-2 assemblies (200 in all; 96 counted as of 2026-09-24), leaving a
sample out of a class when gaps could hide more than 2% of the class in its assembly. Per genome,
most families agree within a robust SD of the log ratio of about 3–8%. Across people, the
correlation also depends on how much people differ. The page of 2026-09-24 shows `HSat1B`
r = 0.99, `ACRO` 0.95, β-satellite 0.94, `CER` 0.89, `HSat1A` 0.87, `HSat3` 0.79. α-satellite HOR
mass reaches only r = 0.65: its mass differs between people by about 4%. `HSat2` does not track
the assemblies with at most 2% of its arrays in marked gaps (r = 0.09 in 47; 0.18 in the 43 with
none), so its estimate is heritable but
not yet confirmed as HSat2 mass. Recomputed on the same table with the revised gap accounting,
which also counts standalone gap records next to an array, the HORs give r = 0.63 in 92 genomes,
`HSat3` 0.81 in 81 and `HSat2` 0.02 in 38; the page shows these once it is regenerated. The
numbers change as more of the cohort is counted; NGS-DOSE-1000G holds the current ones.

## `sinks.satellites.bed` — where NYGC bwa-mem puts the satellite reads

Learned from scans, not from the reference: `ngsdose sinks --classes` with the bundle's rule
(10-kb neighbourhoods holding at least 1e-5 of a family's reads and at least 25 reads in any
training scan; the placement bins with reads inside them, 10 kb for these compositional families,
merged and padded by 1 kb) over 100 whole-file scans of the 1000 Genomes cohort (engine fae1124), 25 per
release batch (2504, 698) and inferred sex, across the 15 populations scanned so far. The file's
header records the rule, the training samples and the capture below. It is for the NYGC bwa-mem
alignments to the GRCh38 analysis set only, and it is kept apart from the bundle's `sinks.bed`, so
that a fetch without the satellite panel neither reads it nor changes its `sinks_sha256`.
`sinks.satellites.stats.tsv` holds each interval's share of its family (median, 10th percentile
and the largest in any one scan, `share_max`) and the capture curve over the 1,648 cohort scans
not used for learning (`ngsdose sinks HELD_OUT --evaluate sinks.satellites.bed --held-out --stats`;
its header says `# held-out: yes`), which `ngsdose fetchplan --capture` trims by; the menu
(`../fetch_menu.tsv`) names both files. With `--crai`, `fetchplan` keeps the satellites'
intervals by share per byte of their CRAM slices, as for the bundle's classes, and its expected
capture is then a lower bound over those held-out scans.

Held-out capture over those 1,648 scans (a placement bin counts only if all of it lies inside a
sink; a fully unmapped read counts as missed):

| class | intervals | Mb | capture, min / median |
| --- | --- | --- | --- |
| `HSat1A` | 30 | 0.22 | 99.98 / 99.99% |
| `HSat1B` | 1,232 | 14.72 | 96.59 / 99.02% (97.17% at the 10th percentile) |
| `HSat2` | 53 | 1.92 | 99.99 / 100.00% |
| `HSat3` | 92 | 3.56 | 99.91 / 99.97% |
| `aSatHOR` | 248 | 59.76 | 99.85 / 99.89% |
| `bSat` | 408 | 2.41 | 99.93 / 99.96% |
| `ACRO` | 22 | 0.43 | 99.92 / 99.96% |
| `SST1` | 34 | 0.61 | 99.90 / 99.97% |
| `CER` | 49 | 1.16 | 99.95 / 99.98% |
| `SATR` | 21 | 0.35 | 99.83 / 99.93% |

`HSat1B` holds 97.92 / 99.42% with the unmapped bin (`count --unmapped`). Its sinks together never
reach a capture target of 0.99 or 0.995 at the 10th percentile (97.17%), so at the default statistic
such a target keeps all of them; with `--capture-stat median` their full capture is 0.990, so a target
up to 0.990 trims them. What fetching the families costs, alone and together, with all their intervals and with
capture targets, is in [`docs/fetch_examples.md`](../../docs/fetch_examples.md) (examples 11 and
12, on 13 NYGC bwa-mem CRAMs); all ten with the controls and the bundle's sinks come to
3,351.6 MB, 21.06% of the median CRAM, as `ngs-dose count -m fetch` reads them (3,099.5 MB with
every slice decoded once). No fetch of all ten can be cheap: they are 4.79% of all
reads, so even a CRAM in which they sat apart from every other read would need about
785 MB read for them.

- **One bin costs more than its reads.** chr2:32,909,000–32,921,000 is a pile-up that held
  1.3 M reads in HG00096's scan, few of them of any class, and it lies in the sinks of `TEL`, `rDNA45S`
  and eight satellite families. Alone it costs 95.4 MB (69.5–126.2). It is why `HSat1A`, `ACRO` and
  `SATR` fall so steeply under a capture target when nothing else in the plan reads it (`SATR`
  drops one interval of 21 and most of its bytes; example 4 of `docs/fetch_examples.md`). Taking it
  out of the bundle's sinks would cost `TEL` capture (lowest of 1,375 held-out scans: 99.39% →
  96.46%) and a little of `rDNA45S`'s (99.89% → 99.80%). It stays in. It holds 1.3% of `TEL`'s reads in the median
  held-out scan (up to 3.4%), so `TEL` keeps it at targets above about 0.979; a class whose target
  alone would drop it keeps it too when the plan reads it for another option, since its slices are
  then read anyway.
- **Per-array readouts come from the placements.** The sinks hold the arrays the reads are placed
  on, so a family's reads can be split by array without a new panel. `subsets/` names four such
  parts of these sinks, each a fetch-menu option of kind `subset` (preset `xy_arrays`; `DXZ1` and
  `DYZ3` are also in `truths`). A fetch reads a sub-option whole and counts it as its family, and
  `ngsdose estimate` reports `NAME.reads`, `NAME.mass_Mb` (the reads times the family's mass per
  read) and `NAME.status` for scans and fetches alike. A sub-option is always fetched whole: a
  capture target does not trim it. What each costs is in `docs/fetch_examples.md` (examples 7, 15 and
  16).

  | option | family, intervals | measures | truth, and what the held-out scans show |
  | --- | --- | --- | --- |
  | `DXZ1` | `aSatHOR`, 4 (4.21 Mb) at chrX's centromere model | DXZ1 array mass over all X chromosomes | HPRC r2 chrX HOR array, both haplotypes: r = 0.992 (n = 124); women read 2.02 times men |
  | `DYZ3` | `aSatHOR`, 4 (0.38 Mb) at chrY's centromere model | DYZ3 array mass | HPRC r2 chrY HOR array: r = 0.835 (60 men); 88.5% of men above the highest woman, 99.9% of women below the lowest man |
  | `DYZ1` | `HSat3`, 1 (0.11 Mb) at chrY:56.67–56.78 Mb | Yq12 HSat3 mass | sex only (the assemblies rarely close Yq12); 19.8% of men's `HSat3`, 0.04% of women's |
  | `DYZ2` | `HSat1B`, 78 (0.91 Mb): 66 on autosomes, 7 on chrX, 3 on decoys, 2 on chrY | an index of male-specific (Y-derived) `HSat1B`: the part in 78 male-only intervals, about a quarter of it | sex only; 21.7% of men's `HSat1B`, 0.3% of women's |

  DYZ2 shows why sinks are learned rather than read off the reference: about 87% of male `HSat1B` is
  Y-derived, but the 13 `HSat1B` intervals on chrY hold 0.16% of it; the reads land mostly on
  autosomes, chrX and decoys (2 of the 78 intervals are on chrY), and the 78 intervals were chosen
  on the 100 training scans as those where men carry it and women do not. They hold 21.7% of men's
  `HSat1B`, so about a quarter of its Y-derived part; the rest lands in intervals that women's reads
  also reach, and `DYZ2.mass_Mb` is an index of the male-specific mass, not the whole of it. In
  30 cohort scans (14 men, 16 women)
  `estimate` gave median masses of 3.46 Mb (men) and 6.17 Mb (women) for DXZ1, 0.79 and 0.013 Mb for
  DYZ3, 15.3 and 0.025 Mb for DYZ1, and 2.62 and 0.005 Mb for DYZ2. A fetch of a sub-option without
  the rest of its family counts the family only there: `estimate` then marks the family
  `subset_only`, which it can tell from the fetch's sinks BED (a fetchplan `PREFIX.sinks.bed` is
  passed with `--fetch-sinks`). Without it, a fetch through a BED `estimate` does not know leaves
  the sub-options `unverified`, with reads and no mass, and the three families that have
  sub-options (`aSatHOR`, `HSat1B`, `HSat3`) `unverified` and NaN wherever the fetch counted them,
  even when the plan selected no sub-option. These readouts expose a few genomes worth checking: HG01683 (inferred female, with
  DYZ3, DYZ1 and DYZ2; possibly XXY) and HG02966 (male, without Yq12 heterochromatin).
- **What is not established.** No fetch with `-p satellites.CHM13v2.k31.panel.tsv.gz` and these
  sinks, or with a sub-option, has been run and compared with the scan of the same genome; all the
  figures above are scan placements (the test fixture has no reads at these arrays, so only the
  mechanics of such a fetch have been run). For the shipped positional classes the fetch returns
  what the placements say (fetch / scan inside the sinks 1.00000–1.00019 over 1,748 genomes), which
  is the reason to expect the same here, and the reason it still has to be shown. The training set
  covers the 15 HG-prefixed populations scanned by then (the NA-prefixed samples were scanned later in the run), and
  some of the 698 batch's children in training have parents among the held-out scans.

## `telomere.k31.panel.tsv.gz` — class `TEL`, six k-mers

The six canonical 31-mers of (TTAGGG)n, unfiltered. A read is assigned with as few as four of
them — 34 bp of perfect repeat, which interstitial telomeric sequence also has — and an exact
31-mer is lost to a single sequencing error where TelSeq's hexamer count is not. So `TEL.mass_Mb`
as it stands is **a relative measure of telomeric content, not a telomere length**: NA12878 gives
46,822 reads, 27% of them with at least half their k-mers telomeric, which would be about 1 kb
per chromosome end from the strict reads and 4 kb from all of them. What makes it worth carrying
is that every class's counts come with the histogram of the share of each read's k-mers that hit
(`hit_frac`, eleven bins), so a threshold can be chosen — and calibrated against TelSeq on a few
samples — after the cohort has been scanned, not before. Like the mitochondrial and EBV dosages it
is first of all a covariate of the state of a cell line.

**Telomeric reads are not dispersed.** Across the first 372 NYGC scans of the 1000 Genomes cohort,
92% of `TEL` reads (85.8–94.9% per genome) were placed within 25 kb of a chromosome end — the 48
windows an NGS-TL/TelSeq-style targeted query retrieves — and 60% in the single bin
chr5:10,000–20,000, where bwa-mem puts pure telomeric reads; essentially none on unplaced or decoy
contigs; the remaining 8% at a fixed set of interstitial loci (chr4:190.12 Mb, chr18:80.26 Mb,
chr2:32.91 Mb, chr1:180 kb, …). Normalised by the controls, the count inside the 48 end windows
and the whole-file count agree at r = 0.9997 across genomes (SD of the log ratio 0.019), so a
targeted query measures the same thing as the scan. Sinks learned from 40 scans
(`ngsdose sinks --classes TEL`: 10-kb neighbourhoods holding ≥ 1e-5 of the class and ≥ 25 reads;
the 10-kb placement bins with reads inside them, merged and padded 1 kb)
captured a median 99.8% of the class in the other 332 (lowest 99.4–99.6% across eleven draws of
the 40); sets learned from 30 scans do no better, at 99.4–99.6% lowest in 200 others over three
draws. The bundle's `sinks.bed` now carries the set learned from all 372 (63 intervals, 810 kb),
which holds a median 99.87% (lowest 99.67%) of the class in those 372 and, across all 3,202 cohort
scans, a median of 99.86%, a 1st percentile of 99.68% and a minimum of 99.39%. A placement bin counts as captured only when all of it lies inside a sink.

Two caveats travel with the number. First, a scan classifies fully unmapped reads, but a fetch
retrieves them only with `ngs-dose count --unmapped`; an unmapped read with a mapped mate sits at
its mate's position and is fetched with its sink regardless. Across the 1,748 cohort scans not one
`TEL` read was placed in the unmapped bin (bwa-mem places every pure telomeric read somewhere, most
of them at chr5p), so nothing is lost this way here. Second, the implied length is 4.4 kb per
chromosome end counting whole reads or 1.7 kb weighting each read by its telomeric k-mer share, an
under-read that a calibration against TelSeq or NGS-TL on the same genomes will size.

## `candidates/` — new classes, for scans only

83 further classes in seven panel files (macrosatellites, multi-copy genes, sex-chromosome
arrays, RNA-gene arrays, viruses and Mycoplasma, coding VNTRs), with their unit sequences and
build recipes. None has sinks: each is loaded in whole-file scans, its sinks are learned there and
checked on held-out scans, and only then does its row in `../fetch_menu.tsv` name a sinks file and
become fetchable. The seven files are final for the 1000 Genomes cohort's remaining scans
(5,076,728 k-mers; `candidates/README.md` lists their sha256). `candidates/README.md` gives the
classes, their truth data, what loading them costs a scan, and the steps.
