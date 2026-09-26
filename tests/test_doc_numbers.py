"""The numbers the resource READMEs quote about the shipped files are the files' own. Each test
computes a figure from the file and asserts that the README states it as written, so an edit to
sinks.bed, the controls, a panel or the fixture that leaves the text behind fails here instead of
going stale unnoticed. A failure names the text it expected: update the README (and this test,
if the wording changed)."""
import gzip
import json
import os
import struct
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "resources" / "GRCh38"
EXPERIMENTAL = ROOT / "resources" / "experimental"
FIXTURE = ROOT / "tests" / "data" / "NA12878.subsample.bam"


def text(path) -> str:
    """The document with its whitespace collapsed, so a figure may wrap across lines."""
    return " ".join(Path(path).read_text().split())


def states(doc: str, phrase: str, path):
    assert " ".join(phrase.split()) in doc, f"{path.relative_to(ROOT)} should say: {phrase}"


def sinks():
    rows = []
    for line in open(BUNDLE / "sinks.bed"):
        c, s, e, cls = line.rstrip("\n").split("\t")
        rows.append((c, int(s), int(e), cls))
    return rows


def per_class(rows):
    out = defaultdict(lambda: [0, 0])
    for _, s, e, cls in rows:
        out[cls][0] += 1
        out[cls][1] += e - s
    return out


def panel_classes(path):
    """The ##class lines of a panel file, as name -> {key: value}."""
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##class\t"):
                kv = dict(f.split("=", 1) for f in line.rstrip("\n").split("\t")[1:])
                out[kv["name"]] = kv
    return out


def bam(path):
    """A BAM read whole (BGZF reads as gzip): contig names, contig lengths, and each record as (tid, pos, flag, l_seq)."""
    with gzip.open(path, "rb") as fh:
        data = fh.read()
    assert data[:4] == b"BAM\1"
    o = 8 + struct.unpack_from("<i", data, 4)[0]
    n_ref, = struct.unpack_from("<i", data, o)
    o += 4
    names, lengths = [], {}
    for _ in range(n_ref):
        ln, = struct.unpack_from("<i", data, o)
        names.append(data[o + 4:o + 3 + ln].decode())
        lengths[names[-1]], = struct.unpack_from("<i", data, o + 4 + ln)
        o += 8 + ln
    recs = []
    while o < len(data):
        size, tid, pos = struct.unpack_from("<iii", data, o)
        flag, = struct.unpack_from("<H", data, o + 18)
        l_seq, = struct.unpack_from("<i", data, o + 20)
        recs.append((tid, pos, flag, l_seq))
        o += 4 + size
    return names, lengths, recs


def fmt_mb(bp, nd=1):
    return f"{bp / 1e6:.{nd}f} Mb"


def test_bundle_readme_quotes_the_sinks_file():
    path = BUNDLE / "README.md"
    doc, rows = text(path), sinks()
    cls = per_class(rows)
    assert set(cls) == {"rDNA45S", "rDNA5S", "DJ", "TEL"}
    pos = [r for r in rows if r[3] != "TEL"]
    states(doc, f"{len(pos)} intervals for the positional classes (`rDNA45S` {cls['rDNA45S'][0]}, `rDNA5S` {cls['rDNA5S'][0]}, "
                f"`DJ` {cls['DJ'][0]}; {fmt_mb(sum(e - s for _, s, e, _ in pos), 2)} as written", path)
    states(doc, f"{cls['TEL'][0]} intervals ({cls['TEL'][1] / 1e3:.0f} kb) for the telomeric repeat", path)
    # the intervals that run past their contig's end, from the fixture's header (the analysis set's lengths)
    _, lengths, _ = bam(FIXTURE)
    over = [(c, e - lengths[c]) for c, s, e, _ in rows if e > lengths[c]]
    words = {8: "eight"}
    states(doc, f"of which {sum(x for _, x in over) / 1e3:.1f} kb run past the ends of {words.get(len(over), len(over))} short contigs", path)
    decoy = sorted(c for c, _, _, k in rows if k == "DJ" and c.endswith("_decoy"))
    assert {c for c, *_ in rows if c.endswith("_decoy")} == set(decoy), "only DJ sinks are on decoys"
    states(doc, f"Two of the {cls['DJ'][0]} DJ sink intervals are on hs38d1 decoy contigs ({', '.join(decoy)})", path)
    assert len(decoy) == 2


def test_bundle_json_records_the_telomeric_sinks():
    meta = json.loads((BUNDLE / "bundle.json").read_text())
    n, bp = per_class(sinks())["TEL"]
    assert (meta["sinks_learned_from"]["TEL"]["intervals"], meta["sinks_learned_from"]["TEL"]["bp"]) == (n, bp)


def test_experimental_readme_quotes_the_telomeric_sinks():
    path = EXPERIMENTAL / "README.md"
    n, bp = per_class(sinks())["TEL"]
    states(text(path), f"the set learned from all 372 ({n} intervals, {bp / 1e3:.0f} kb)", path)


def test_bundle_readme_quotes_the_controls():
    path = BUNDLE / "README.md"
    doc = text(path)
    bed = [line.rstrip("\n").split("\t") for line in open(BUNDLE / "controls.bed")]
    roles = defaultdict(int)
    for _, _, _, role in bed:
        roles[role] += 1
    ctrl_bp = sum(int(e) - int(s) for _, s, e, r in bed if r == "control")
    flanks = {h.split("flank=")[1].split()[0] for h in gzip.open(BUNDLE / "controls.fa.gz", "rt") if h.startswith(">")}
    assert flanks == {"1000"}
    states(doc, f"{roles['control']} control regions ({fmt_mb(ctrl_bp)}); {roles['test:auto']} held-out autosomal, "
                f"{roles['test:chrX']} chrX and {roles['test:chrY']} chrY known-truth regions; one dosage region each on chrM and chrEBV; "
                "all with 1 kb flanks", path)
    assert roles["dosage:chrM"] == roles["dosage:chrEBV"] == 1 and len(roles) == 6


def test_bundle_readme_quotes_the_panel():
    path = BUNDLE / "README.md"
    doc = text(path)
    k = {n: int(v["kmers_kept"]) for n, v in panel_classes(BUNDLE / "panel.k31.tsv.gz").items()}
    states(doc, f"`rDNA45S` ({k['rDNA45S']:,}), `rDNA5S` ({k['rDNA5S']:,}), `DJ` ({k['DJ']:,})", path)
    states(doc, f"of the {k['DJ']:,} panel k-mers", path)


def test_experimental_readme_quotes_the_panels():
    path = EXPERIMENTAL / "README.md"
    doc = text(path)
    sat = panel_classes(EXPERIMENTAL / "satellites.CHM13v2.k31.panel.tsv.gz")
    assert len(sat) == 10 and all(v["kind"] == "compositional" for v in sat.values())
    total = sum(int(v["kmers_kept"]) for v in sat.values())
    states(doc, f"ten families, {total / 1e6:.2f} M k-mers", path)
    for name, v in sat.items():
        size = int(v["length"]) / 1e6
        # the table's CHM13 (haploid) and k-mers columns
        states(doc, f"| {size:.1f} Mb | {int(v['kmers_kept']):,} |", path)
        assert f"| `{name}` |" in doc, name
    tel = panel_classes(EXPERIMENTAL / "telomere.k31.panel.tsv.gz")
    assert list(tel) == ["TEL"] and int(tel["TEL"]["kmers_kept"]) == 6
    states(doc, "class `TEL`, six k-mers", path)


def test_fixture_readme_quotes_the_fixture():
    path = ROOT / "tests" / "data" / "README.md"
    doc = text(path)
    names, _, recs = bam(FIXTURE)
    states(doc, f"{len(recs):,} reads", path)
    assert all(tid >= 0 for tid, *_ in recs)
    states(doc, "It holds no read without a coordinate", path)

    def intervals(rows):
        d = defaultdict(list)
        for c, s, e in rows:
            d[c].append((s, e))
        out = {}
        for c, v in d.items():
            v.sort()
            m = [list(v[0])]
            for s, e in v[1:]:
                if s <= m[-1][1]:
                    m[-1][1] = max(m[-1][1], e)
                else:
                    m.append([s, e])
            out[c] = (np.array([x[0] for x in m]), np.array([x[1] for x in m]))
        return out

    def depth(recs, iv):
        bases = 0
        for tid, p, flag, l_seq in recs:
            st, en = iv.get(names[tid], (None, None))
            if st is not None:
                i = np.searchsorted(st, p, side="right") - 1
                if i >= 0 and p < en[i]:
                    bases += l_seq
        return bases / sum(int((en - st).sum()) for st, en in iv.values())

    prim = [r for r in recs if not r[2] & 0x904]
    bed = [line.rstrip("\n").split("\t") for line in open(BUNDLE / "controls.bed")]
    ctrl = intervals([(c, int(s), int(e)) for c, s, e, r in bed if r == "control"])
    states(doc, f"~{depth(prim, ctrl):.2f}× over the {sum(len(v[0]) for v in ctrl.values())} control regions "
                f"(duplicates included; ~{depth([r for r in prim if not r[2] & 0x400], ctrl):.2f}× without)", path)
    # the regions make_fixture.sh cut: the controls with 1,600 bp either side, and the sinks of the day (no TEL)
    rows = sinks()
    kept = intervals([(c, max(0, int(s) - 1600), int(e) + 1600) for c, s, e, _ in bed] + [(c, s, e) for c, s, e, k in rows if k != "TEL"])
    states(doc, f"~{depth(prim, kept):.1f}× over all retained regions", path)
    # which TEL sinks hold fixture reads
    tel = [(c, s, e) for c, s, e, k in rows if k == "TEL"]
    starts = defaultdict(list)
    for tid, p, *_ in recs:
        starts[names[tid]].append(p)
    starts = {c: np.sort(v) for c, v in starts.items()}
    held = [(c, s, e) for c, s, e in tel if c in starts and np.searchsorted(starts[c], e) > np.searchsorted(starts[c], s)]
    assert len(held) == 1
    c, s, e = held[0]
    assert (c, s, e, "rDNA45S") in rows
    states(doc, f"sink ({c}:{s:,}-{e:,}); the other {len(tel) - 1} hold no fixture reads", path)


def test_bundle_readme_capture_figures_match_the_cohort_scans():
    """Opt-in (NGSDOSE_COHORT_SCANS = a directory of scan-mode counts files): the capture figures the
    README quotes, recomputed with ngsdose.sinks.capture. Runs only when the directory holds exactly
    the number of scans the README quotes, so it is a check of the figures as of their date."""
    d = os.environ.get("NGSDOSE_COHORT_SCANS")
    if not d:
        pytest.skip("NGSDOSE_COHORT_SCANS not set")
    from ngsdose import io
    from ngsdose import sinks as sk
    files = sorted(Path(d).glob("*.json.gz"))
    path = BUNDLE / "README.md"
    doc = text(path)
    if f"over the {len(files):,} cohort scans" not in doc:
        pytest.skip(f"{d} holds {len(files)} scans, not the number the README quotes")
    bed = sk.read_bed(BUNDLE / "sinks.bed")
    frac = defaultdict(list)
    for f in files:
        for cls, (n, inside) in sk.capture(io.load_counts(f), bed).items():
            frac[cls].append(inside / max(n, 1))
    parts = []
    for cls in ("rDNA45S", "rDNA5S", "DJ", "TEL"):
        v = 100 * np.array(frac[cls])
        parts.append(f"`{cls}` {v.mean():.2f} / {np.percentile(v, 1):.2f} / {v.min():.2f}%")
    states(doc, ", ".join(parts), path)
    tel = 100 * np.array(frac["TEL"])
    states(doc, f"({int((tel < 99.67).sum())} of the {len(files):,} below 99.67%)", path)
    path = EXPERIMENTAL / "README.md"
    states(text(path), f"across the {len(files):,} scans counted by 2026-09-25, a median of {np.median(tel):.2f}%, "
                       f"a 1st percentile of {np.percentile(tel, 1):.2f}% and a minimum of {tel.min():.2f}%", path)
