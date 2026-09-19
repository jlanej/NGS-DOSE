"""Simulate a tiny paired-end WGS with known truth: a single-copy contig and a multi-copy tandem
class, sequenced under a fragment-GC bias. Used by the end-to-end test of engine + estimator."""
from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def random_seq(rng, blocks) -> str:
    out = []
    for gc, n in blocks:
        p = np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])
        out.append("".join(rng.choice(list("ACGT"), size=n, p=p)))
    return "".join(out)


def bias(gc: np.ndarray, strength: float) -> np.ndarray:
    return np.exp(-strength * ((gc - 0.45) / 0.2) ** 2)


def _fragments(rng, template: str, rate: float, frag_mean: int, frag_sd: int, strength: float, lo: int, hi: int):
    """Fragment starts in [lo, hi) of the template, thinned by the GC of the actual fragment."""
    a = np.frombuffer(template.encode(), dtype=np.uint8)
    cg = np.r_[0, np.cumsum((a == 71) | (a == 67))]
    n = rng.poisson(rate * (hi - lo))
    starts = rng.integers(lo, hi, n)
    lens = np.clip(rng.normal(frag_mean, frag_sd, n).round().astype(int), 120, None)
    ok = starts + lens <= len(template)
    starts, lens = starts[ok], lens[ok]
    gc = (cg[starts + lens] - cg[starts]) / lens
    keep = rng.random(len(starts)) < bias(gc, strength)
    return starts[keep], lens[keep]


def _mutate(rng, s: str, err: float) -> str:
    if err <= 0:
        return s
    b = np.frombuffer(s.encode(), dtype=np.uint8).copy()
    hit = np.where(rng.random(len(b)) < err)[0]
    for i in hit:
        b[i] = rng.choice([x for x in b"ACGT" if x != b[i]])
    return b.tobytes().decode()


def satellite_array(rng, monomer_len=171, n_monomers=300, gc=0.38, divergence=0.03) -> str:
    """A tandem array of diverged monomers: a family with no stable unit (a compositional class)."""
    base = random_seq(rng, [(gc, monomer_len)])
    return "".join(_mutate(rng, base, divergence) for _ in range(n_monomers))


def simulate(outdir, seed=3, copies=60.0, depth=30.0, read_len=100, frag_mean=300, frag_sd=25, strength=1.2,
             err=0.002, dup_frac=0.08, ctrl_len=300_000, unit_blocks=((0.42, 900), (0.68, 700), (0.5, 800), (0.33, 600)),
             sat_array_copies=8.0):
    """Returns paths and truth. `copies` is the diploid copy number of the positional class;
    the satellite array is present `sat_array_copies` times per diploid genome."""
    rng = np.random.default_rng(seed)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    ctrl = random_seq(rng, [(gc, ctrl_len // 10) for gc in (0.35, 0.42, 0.55, 0.38, 0.65, 0.45, 0.30, 0.50, 0.60, 0.40)])
    unit = random_seq(rng, unit_blocks)
    U = len(unit)
    n_ref_copies = 3
    sat = satellite_array(rng)
    ref = {"ctrl": ctrl, "sink": unit * n_ref_copies, "satsink": sat}
    (outdir / "ref.fa").write_text("".join(f">{k}\n{v}\n" for k, v in ref.items()))
    (outdir / "unit.fa").write_text(f">unit\n{unit}\n")
    (outdir / "background.fa").write_text(f">ctrl\n{ctrl}\n")
    (outdir / "sat.fa").write_text(f">sat_array\n{sat}\n")
    (outdir / "classes.tsv").write_text("unit\tpositional\tunit.fa\t1\nsat\tcompositional\tsat.fa\t0\t1\n")
    # control regions: five blocks of the single-copy contig, away from the contig ends
    step = ctrl_len // 5
    (outdir / "controls.bed").write_text("".join(f"ctrl\t{i * step + 2000}\t{(i + 1) * step - 2000}\tcontrol\n" for i in range(5)))
    (outdir / "sinks.bed").write_text(f"sink\t0\t{U * n_ref_copies}\tunit\nsatsink\t0\t{len(sat)}\tsat\n")

    # per-strand fragment-start rate that gives `depth` for diploid single-copy DNA before bias:
    # depth = 2 strands * rate * 2 reads * read_len  (fragments counted on both strands)
    rate = depth / (2 * 2 * read_len)
    recs = []

    def emit(name, chrom, pos_f, seq_f, pos_r, seq_r, tlen, mapq, fwd_first, dup):
        f1, f2 = (99, 147) if fwd_first else (163, 83)
        d = 1024 if dup else 0
        recs.append((chrom, pos_f, f"{name}\t{f1 + d}\t{chrom}\t{pos_f + 1}\t{mapq}\t{read_len}M\t=\t{pos_r + 1}\t{tlen}\t{seq_f}\t*"))
        recs.append((chrom, pos_r, f"{name}\t{f2 + d}\t{chrom}\t{pos_r + 1}\t{mapq}\t{read_len}M\t=\t{pos_f + 1}\t{-tlen}\t{seq_r}\t*"))

    # single-copy contig: two haplotypes x two strands are symmetric here, so draw at 2x the rate
    # with a random orientation of the read pair
    starts, lens = _fragments(rng, ctrl, 2 * rate, frag_mean, frag_sd, strength, 0, len(ctrl) - 1000)
    for i, (s, ln) in enumerate(zip(starts, lens)):
        frag = ctrl[s:s + ln]
        emit(f"c{i}", "ctrl", s, _mutate(rng, frag[:read_len], err), s + ln - read_len,
             _mutate(rng, frag[-read_len:], err), ln, 60, rng.random() < 0.5, rng.random() < dup_frac)
    # the class: a long tandem array; every unit position is present `copies` times in the diploid genome
    array = unit * 40
    starts, lens = _fragments(rng, array, 2 * rate * copies / 2 * (1 / 1), frag_mean, frag_sd, strength, U * 10, U * 11)
    for i, (s, ln) in enumerate(zip(starts, lens)):
        frag = array[s:s + ln]
        off = s % U + U * int(rng.integers(0, n_ref_copies - 1))          # the aligner scatters reads over paralogs
        emit(f"u{i}", "sink", off, _mutate(rng, frag[:read_len], err), off + ln - read_len,
             _mutate(rng, frag[-read_len:], err), ln, 0, rng.random() < 0.5, rng.random() < dup_frac / 3)
    n_unit_reads = 2 * len(starts)
    unit_starts, unit_lens = starts, lens
    # the satellite: fragments from the array, placed where they came from (one reference copy)
    starts, lens = _fragments(rng, sat, 2 * rate * sat_array_copies / 2, frag_mean, frag_sd, strength, 0, len(sat) - 600)
    for i, (s, ln) in enumerate(zip(starts, lens)):
        frag = sat[s:s + ln]
        emit(f"s{i}", "satsink", s, _mutate(rng, frag[:read_len], err), s + ln - read_len,
             _mutate(rng, frag[-read_len:], err), ln, 0, rng.random() < 0.5, False)
    # expected ends before GC thinning correspond to this much diploid sequence
    sat_mass = sat_array_copies * (len(sat) - 600)
    starts, lens = unit_starts, unit_lens
    order = {"ctrl": 0, "sink": 1, "satsink": 2}
    recs.sort(key=lambda r: (order[r[0]], r[1]))
    header = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(f"@SQ\tSN:{k}\tLN:{len(v)}\n" for k, v in ref.items()) + "@RG\tID:sim\tSM:simulated\n"
    sam = outdir / "sim.sam"
    with open(sam, "w") as fh:
        fh.write(header)
        for _, _, line in recs:
            fh.write(line + "\tRG:Z:sim\n")
    bam = outdir / "sim.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    subprocess.run(["samtools", "faidx", str(outdir / "ref.fa")], check=True)
    sam.unlink()
    nb = -(-U // 50)
    truth_fwd = np.bincount((starts % U) // 50, minlength=nb)                 # forward read: 5' end = fragment start
    truth_rev = np.bincount(((starts + lens - 1) % U) // 50, minlength=nb)    # reverse read: 5' end = fragment end
    return dict(bam=bam, dir=outdir, unit=unit, copies=copies, n_class_reads=n_unit_reads,
                truth_fwd=truth_fwd, truth_rev=truth_rev, sat_mass_bp=sat_mass)
