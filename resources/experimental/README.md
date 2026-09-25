# Experimental resources: dispersed sequence

A positional class (rDNA, the distal junction) has a unit, its reads land in a few places, and a
targeted fetch finds them. The satellite families have neither property: they are spread over
centromere models, decoys and everything else in an alignment, so **only a whole-file scan
measures them** — which is why a cohort that is scanned once should be scanned with these loaded
(NGS-DOSE-1000G's pipeline does; `EXTRA_PANELS` in its `config.sh`). The telomeric repeat is the exception:
the aligner concentrates its reads at the chromosome ends, the bundle's `sinks.bed` carries the
intervals, and a fetch with `-p telomere.k31.panel.tsv.gz` measures it (`FETCH_PANELS`).

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
family's CHM13 arrays, in no other family, and nowhere in CHM13 outside CenSat-annotated
satellite. *Recall* is the share of 150-bp reads drawn from the family's own CHM13 arrays that
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

The first five are measured; `ACRO` nearly; `bSat`, `SST1`, `CER` and `SATR` are relative
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

Two samples say nothing about whether the estimates *track* the assemblies across people, which
is what association work needs. Two hundred samples of the 1000 Genomes cohort have HPRC
assemblies; NGS-DOSE-1000G's `pipeline/04_hprc_satellites.sh` makes the comparison once the cohort is scanned.

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
(`ngsdose sinks --classes TEL`, 10-kb bins holding ≥ 1e-5 of the class and ≥ 25 reads, padded 1 kb)
captured ≥ 99.87% of the class in each of the other 332; the bundle's `sinks.bed` now carries the
set learned from all 372. Two caveats travel with the number: the engine skips unmapped records in
both modes, so a fully unmapped telomeric read is counted by neither (TelSeq counts it); and the
implied length is 4.4 kb per chromosome end counting whole reads or 1.7 kb weighting each read by
its telomeric k-mer share, an under-read that a calibration against TelSeq or NGS-TL on the same
genomes will size.
