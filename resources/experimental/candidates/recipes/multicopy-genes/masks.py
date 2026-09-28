#!/usr/bin/env python3
"""Background masks: every copy of each class in GRCh38 (analysis set: primary, alts, decoys,
HLA, unplaced) and in CHM13, so that k-mers found only at the class's own copies survive and a
k-mer found anywhere else (a paralog, a pseudogene, an interspersed repeat) is dropped.

Copies come from the genome-wide search (find_copies.py; the PAFs of the class units against
1-Mb windows of each assembly). A hit is a copy of the class when its gap-compressed identity and
aligned length pass the class's rule below; the rules were set by reading every hit above 90%
identity (the thresholds sit in the gap between co-copies and paralogs; see the notes). Three
sub-units are searched in their loci only (regional), with the same rule format.

    masks.py GRCh38.paf CHM13.paf GRCh38.fa CHM13.fa OUTDIR
writes OUTDIR/<class>.GRCh38.bed, <class>.CHM13.bed and copies.tsv (every copy with its identity).
"""
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from find_copies import parse  # noqa: E402

# class: (unit searched, min identity, min aligned bp, note)
RULES = {
    "KIV2": ("KIV2", 0.97, 1500, "KIV-2 units 97.8-100%; KIV-3 91%, KIV-1 96% (its first 1.75 kb identical to KIV-2) and KIV-4 84% stay background"),
    "C4": ("C4", 0.99, 2000, "C4A/C4B long and short copies 99.8-100%; HERV-K elements elsewhere 85-89%"),
    "HERVC4": ("C4", 0.99, 2000, "the C4 gene copies (the insertion sits in the long ones); every other HERV-K stays background"),
    "CYP21A2": ("CYP21A2", 0.995, 2500, "CYP21A2 copies 99.6-100%; CYP21A1P 98.1%"),
    "AMY1": ("AMY1", 0.99, 5000, "3 copies GRCh38, 7 CHM13, 99.9%; AMY2A/2B 92%"),
    "AMY2B": ("AMY2B", 0.99, 10000, "one copy each"),
    "AMY2A": ("AMY2A", 0.99, 7000, "the full gene only; the 6.4-kb AMY2Ap partial copies (99.7%) stay background"),
    "SMN": ("SMN", 0.99, 20000, "SMN1, SMN2, three alt copies, two CHM13 copies, all 99.8-100%"),
    "SMN1": ("SMN", 0.99, 20000, "SMN1-like copies only (SMN1 allele at c.840): SMN1, chr5_GL339449v2_alt, chr5_KI270897v1_alt 473-501 kb, CHM13 chr5:71.38 Mb"),
    "RHD": ("RHD", 0.99, 30000, "RHD 99.9-100%; RHCE 98.2% and the upstream Rhesus box stay background"),
    "GSTM1": ("GSTM1", 0.99, 3000, "GRCh38 only (CHM13 carries the deletion); GSTM2/4/5 below 95%"),
    "GSTT1": ("GSTT1", 0.99, 4000, "GRCh38 chr22_KI270879v1_alt and CHM13 chr22"),
    "UGT2B17": ("UGT2B17", 0.99, 5000, "one copy each; UGT2B15 below 97%"),
    "LCE3BC": ("LCE3BC", 0.99, 10000, "one copy each"),
    "APOBEC3B": ("APOBEC3B", 0.99, 5000, "one copy each; APOBEC3A-side homology 96.6%"),
    "HPR": ("HPR", 0.99, 5000, "one copy each; HP 94%"),
    "CCL3L": ("CCLblk", 0.99, 3000, "primary, KI270857 (1), KI270909 (2), CHM13 (2 + a 5.9-kb partial); the CCL3/CCL4 block 94.5%"),
    "ORM1": ("ORM1", 0.99, 3000, "one copy each; ORM2 below 97%"),
    "DEFB": ("DEFB", 0.99, 20000, "REPP, REPD, chr8_KI270813v1_alt, three CHM13 copies, 99.5-100%"),
    "FCGR3": ("FCGR3", 0.975, 5000, "FCGR3B and FCGR3A (98.1%); CHM13 FCGR3A + two FCGR3B"),
    "DEFA1A3": ("DEFA1A3", 0.98, 5000, "2 full units + the 7.7-kb DEFA3 partial unit in each assembly"),
    "NPY4R": ("NPY4R", 0.99, 3000, "both SD copies and the 4-kb 99.0% partial at chr10:47.46 Mb (GRCh38) / 48.35 Mb (CHM13)"),
    "SULT1A": ("SULT1A", 0.99, 3000, "one copy each"),
    "RASA4": ("RASA4", 0.995, 5000, "the two tandem units and the array-contiguous partial third unit; the 44-Mb paralog (98.4%), SPDYE6 (99.4%) and the SPDYE family stay background"),
    "NOTCH2NL": ("NOTCH2NL", 0.99, 50000, "NOTCH2NLA/B/C/R and the NOTCH2 5' region, 99.2-100%"),
    "TBC1D3": ("TBC1D3", 0.99, 8000, "17q12 TBC1D3 paralogs 99.3-99.7% (primary, KI270857, KI270909); 17q23 TBC1D3P1/P2 98% and a 97.3% partial stay background"),
}
# restrict a rule to hits overlapping these loci (SMN1-like copies)
ONLY = {
    "SMN1": {"GRCh38": [("chr5", 70924940, 70953012), ("chr5_GL339449v2_alt", 457778, 485846), ("chr5_KI270897v1_alt", 473375, 501444)],
             "CHM13": [("chr5", 71381728, 71409801)]},
}
# classes whose genome-wide hits are completed by a search of their locus (the genome-wide search
# split the short C4 copy at chr6_GL000252v2_alt:3,262,569-3,276,826 and reported only 2.6 kb of it)
MHC = {"GRCh38": ["chr6:31900001-32100000"] + [f"chr6_GL0002{i}v2_alt" for i in range(50, 57)], "CHM13": ["chr6:31750001-31950000"]}
EXTRA = {"C4": MHC, "HERVC4": MHC}
# sub-units searched only near their locus: class -> {assembly: [regions]}
REGIONAL = {
    "HBA": ({"GRCh38": ["chr16:150001-200000"], "CHM13": ["chr16:145001-195000"]}, 0.97, 1200,
            "HBA2 and HBA1 blocks (98-99%); HBAP1 (the 1.4-kb 95% hit) and HBM/HBQ1 stay background"),
    "HPdup": ({"GRCh38": ["chr16:72040001-72090000"], "CHM13": ["chr16:77860001-77910000"]}, 0.97, 1000,
              "both copies in GRCh38 (Hp2), one in CHM13 (Hp1); HPR stays background"),
}


def regional_hits(unit_fa, fa, regions, tmp):
    q = tmp / "regional.fa"
    with open(q, "w") as fh:
        subprocess.run(["samtools", "faidx", fa] + regions, stdout=fh, check=True)
    paf = tmp / "regional.paf"
    with open(paf, "w") as fh:
        subprocess.run(["minimap2", "-x", "asm20", "-c", "--cs", "-P", str(unit_fa), str(q)], stdout=fh, stderr=subprocess.DEVNULL, check=True)
    return parse(str(paf), 0)


def main():
    g_paf, c_paf, g_fa, c_fa, outdir = sys.argv[1:6]
    here = Path(__file__).resolve().parent
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    rows = {"GRCh38": parse(g_paf, 500), "CHM13": parse(c_paf, 500)}
    fas = {"GRCh38": g_fa, "CHM13": c_fa}
    tab = open(out / "copies.tsv", "w")
    tab.write("class\tassembly\tcontig\tstart\tend\tunit_start\tunit_end\tstrand\tidentity\trule\n")
    for cls in list(RULES) + list(REGIONAL):
        for asm in ("GRCh38", "CHM13"):
            if cls in RULES:
                unit, idy, mlen, note = RULES[cls]
                hits = [r for r in rows[asm] if r["unit"] == unit and r["idy"] >= idy and r["te"] - r["ts"] >= mlen]
                if cls in EXTRA:
                    ex = regional_hits(here / "units" / "C4.fa", fas[asm], EXTRA[cls][asm], out)
                    hits += [r for r in ex if r["idy"] >= idy and r["te"] - r["ts"] >= mlen]
                if cls in ONLY:
                    hits = [r for r in hits if any(r["ctg"] == c and r["gs"] < e and r["ge"] > s for c, s, e in ONLY[cls][asm])]
            else:
                regions, idy, mlen, note = REGIONAL[cls]
                hits = regional_hits(here / "units" / f"{cls}.fa", fas[asm], regions[asm], out)
                hits = [r for r in hits if r["idy"] >= idy and r["te"] - r["ts"] >= mlen]
            ivs = sorted({(r["ctg"], r["gs"], r["ge"]) for r in hits})
            with open(out / f"{cls}.{asm}.bed", "w") as fh:
                for c, s, e in ivs:
                    fh.write(f"{c}\t{s}\t{e}\t{cls}\n")
            for r in sorted(hits, key=lambda r: (r["ctg"], r["gs"])):
                tab.write(f"{cls}\t{asm}\t{r['ctg']}\t{r['gs']}\t{r['ge']}\t{r['ts']}\t{r['te']}\t{r['strand']}\t{r['idy']:.4f}\t{note}\n")
    for f in ("regional.fa", "regional.paf"):
        if (out / f).exists():
            os.remove(out / f)


if __name__ == "__main__":
    main()
