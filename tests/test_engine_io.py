"""The engine's inputs and how it fails on them: CRAM (reference from -T or REF_PATH), reads without a
coordinate, a file served over HTTP (with and without range requests, with lost requests, behind a
signed URL; fetched by exact byte ranges and by htslib's own reader), sinks on contigs the input lacks,
and inputs it must refuse. Needs the built binary and samtools; every input is made in the test from the
committed fixture or the simulated genome."""
import hashlib
import json
import os
import re
import shutil
import struct
import subprocess
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import numpy as np
import pytest

from ngsdose import contract, fetchplan, io, resources

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("NGSDOSE_BIN", ROOT / "target" / "release" / "ngs-dose"))
BAM = ROOT / "tests" / "data" / "NA12878.subsample.bam"
BUNDLE = resources.Bundle(ROOT / "resources" / "GRCh38")

pytestmark = pytest.mark.skipif(not BIN.exists() or shutil.which("samtools") is None,
                                reason="needs target/release/ngs-dose (cargo build --release) and samtools")

# what two counts of the same reads differ in: where they were read from, how long it took, and the
# @PG line samtools adds to every copy it writes (the @SQ fingerprint of `pipeline` is compared on its own)
RUN = {"input", "elapsed_sec", "pipeline"}
SCAN_FETCH_DIFFER = {"mode", "sinks", "sinks_sha256", "sinks_skipped", "pad", "records", "primary", "primary_dup_flagged",
                     "unmapped", "contigs", "elapsed_sec"}


def strip(c, drop=RUN):
    return {k: v for k, v in c.items() if k not in drop}


def reads(c):
    return {x["name"]: x["reads"] for x in c["classes"]}


def environ(**kw):
    """The environment without a reference path, plus `kw`; a local server is never reached through a proxy."""
    e = {k: v for k, v in os.environ.items() if k not in ("REF_PATH", "REF_CACHE")}
    e["no_proxy"] = e["NO_PROXY"] = "127.0.0.1,localhost"
    return {**e, **{k: str(v) for k, v in kw.items()}}


def engine(*args, cwd=None, env=None):
    return subprocess.run([str(BIN), *map(str, args)], cwd=cwd, env=env or environ(), capture_output=True, text=True, timeout=300)


def count(inp, mode, out, *extra, threads=2, cwd=None, env=None):
    args = ["count", "-i", inp, "-p", BUNDLE.panel, "-c", BUNDLE.controls, "-m", mode, "-@", threads, "-o", out, *extra]
    if mode == "fetch" and "--sinks" not in map(str, extra):
        args += ["--sinks", BUNDLE.sinks]
    return engine(*args, cwd=cwd, env=env)


def ok(r):
    assert r.returncode == 0, r.stderr
    return r


def samtools(*args, **kw):
    return subprocess.run(["samtools", *map(str, args)], check=True, capture_output=True, **kw)


def own(stderr):
    """The engine's own messages: htslib prints its errors and warnings with the full file name itself."""
    return "\n".join(line for line in stderr.splitlines() if not line.startswith(("[E::", "[W::")))


def bgzf_cut(data: bytes) -> bytes:
    """The file up to the first BGZF block boundary past its middle: a copy that stopped early."""
    o = 0
    while o < len(data) // 2:
        p = o + 12
        while data[p:p + 2] != b"BC":                                  # BSIZE, the block size - 1, is in the BC subfield
            p += 4 + struct.unpack_from("<H", data, p + 2)[0]
        o += struct.unpack_from("<H", data, p + 4)[0] + 1
    return data[:o]


@pytest.fixture(scope="module")
def local(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("local")
    out = {}
    for mode in ("scan", "fetch"):
        ok(count(BAM, mode, tmp / f"{mode}.json"))
        out[mode] = io.load_counts(tmp / f"{mode}.json")
    return out


def test_the_pipeline_record_is_taken_from_the_header(local):
    """`pipeline` can be recomputed from `samtools view -H`: the @PG ids, programs and versions in header
    order, and the sha256 of the sorted "SN<TAB>LN<TAB>M5" lines of the @SQ records."""
    lines = samtools("view", "-H", "--no-PG", BAM, text=True).stdout.splitlines()
    tags = lambda line: dict(f.split(":", 1) for f in line.split("\t")[1:])
    sq = [tags(line) for line in lines if line.startswith("@SQ\t")]
    pg = [{k: t[K] for k, K in (("id", "ID"), ("pn", "PN"), ("vn", "VN")) if K in t} for t in (tags(line) for line in lines if line.startswith("@PG\t"))]
    fp = hashlib.sha256(b"".join(sorted(f"{t['SN']}\t{t['LN']}\t{t.get('M5', '')}\n".encode() for t in sq))).hexdigest()
    want = {"pg": pg, "sq_sha256": fp, "sq_n": len(sq), "sq_m5": sum("M5" in t for t in sq)}
    assert local["scan"]["pipeline"] == local["fetch"]["pipeline"] == want


def test_sinks_on_contigs_the_input_lacks_refuse_the_fetch_unless_told(local, tmp_path):
    """Sinks learned on another reference or pipeline can name contigs the file does not have. A class that
    would lose any interval to them is one `ngsdose estimate` reports as NA (sinks_skipped), so the fetch is
    refused before it reads a byte (engine 0.1.1; 0.1.0 warned and spent the bytes), naming the class and
    what it loses. With --allow-missing-sinks it goes on, leaves those intervals out, says so, and records
    what it left out per class. `--classes` without the class needs no flag, and `plan` reports the drop."""
    bed = tmp_path / "sinks.bed"
    bed.write_text(Path(BUNDLE.sinks).read_text() + "chrAbsent_decoy\t1000\t6000\tDJ\n")
    out = tmp_path / "f.json"
    r = count(BAM, "fetch", out, "--sinks", bed)
    assert r.returncode == 1 and not out.exists(), r.stderr
    assert "class(es) DJ lose sink intervals" in r.stderr and "DJ 1 interval (5.0 kb) on chrAbsent_decoy" in r.stderr, r.stderr
    assert "absent from the alignment header" in r.stderr and "--allow-missing-sinks" in r.stderr and "no interval" not in r.stderr
    r = ok(count(BAM, "fetch", out, "--sinks", bed, "--allow-missing-sinks"))
    assert "WARNING" in r.stderr and "sinks_skipped" in r.stderr
    c = io.load_counts(out)
    assert c["sinks_skipped"] == {"DJ": {"intervals": 1, "bp": 5000}}
    assert "sinks_missing_classes" not in c
    drop = RUN | {"sinks", "sinks_sha256", "sinks_skipped"}
    assert strip(c, drop) == strip(local["fetch"], drop)
    assert set(contract.incomplete_sinks(c)) == {"DJ"} and contract.incomplete_sinks(local["fetch"]) == {}
    r = ok(count(BAM, "fetch", tmp_path / "r.json", "--sinks", bed, "--classes", "rDNA45S,rDNA5S"))
    assert "WARNING" not in r.stderr and "sinks_skipped" not in io.load_counts(tmp_path / "r.json")
    # an interval without a class serves every class (as `ngsdose sinks` writes them): every class loses it
    bed.write_text(Path(BUNDLE.sinks).read_text() + "chrAbsent_decoy\t1000\t6000\n")
    r = count(BAM, "fetch", out, "--sinks", bed)
    assert r.returncode == 1 and "class(es) rDNA45S, rDNA5S, DJ lose sink intervals" in r.stderr, r.stderr
    assert "every class (intervals without a class) 1 interval (5.0 kb) on chrAbsent_decoy" in r.stderr, r.stderr
    ok(count(BAM, "fetch", out, "--sinks", bed, "--allow-missing-sinks"))
    c = io.load_counts(out)
    assert c["sinks_skipped"] == {"": {"intervals": 1, "bp": 5000}} and set(contract.incomplete_sinks(c)) == {"rDNA45S", "rDNA5S", "DJ"}
    r = ok(engine("plan", "-c", BUNDLE.controls, "--sinks", bed, "-i", BAM, "-o", tmp_path / "plan.bed"))
    assert "sink intervals by class: (no class) 1" in r.stderr


def test_a_class_whose_sinks_are_all_on_absent_contigs_is_refused(local, tmp_path):
    """Every DJ sink under a contig name the file does not use (no 'chr' prefix): the fetch would read DJ
    only where other intervals happen to hold it, so it refuses unless told to go on."""
    rows = [line.split("\t") for line in Path(BUNDLE.sinks).read_text().splitlines() if line and not line.startswith("#")]
    dj = [r for r in rows if r[3] == "DJ"]
    bed = tmp_path / "sinks.bed"
    bed.write_text("".join("\t".join([r[0].removeprefix("chr") if r[3] == "DJ" else r[0], *r[1:]]) + "\n" for r in rows))
    out = tmp_path / "f.json"
    r = count(BAM, "fetch", out, "--sinks", bed)
    assert r.returncode == 1 and not out.exists(), r.stderr
    assert "DJ" in r.stderr and "absent from the alignment header" in r.stderr and "--allow-missing-sinks" in r.stderr
    ok(count(BAM, "fetch", out, "--sinks", bed, "--allow-missing-sinks"))
    c = io.load_counts(out)
    assert c["sinks_missing_classes"] == ["DJ"]
    assert c["sinks_skipped"] == {"DJ": {"intervals": len(dj), "bp": sum(int(r[2]) - int(r[1]) for r in dj)}}
    assert reads(c)["DJ"] < 0.1 * reads(local["fetch"])["DJ"]
    assert contract.incomplete_sinks(c)["DJ"] == "no sink intervals in the fetch"


@pytest.fixture(scope="module")
def cram(tmp_path_factory):
    """The fixture as a CRAM that carries its own sequence (no_ref: no reference needed to decode it), in
    small slices so that a fetch's random access stays cheap."""
    d = tmp_path_factory.mktemp("cram")
    samtools("view", "-C", "--output-fmt-option", "no_ref=1", "--output-fmt-option", "seqs_per_slice=1000", "-o", d / "fx.cram", BAM)
    samtools("index", d / "fx.cram")
    (d / "noref").mkdir()
    return d


def test_a_cram_gives_the_counts_of_its_bam(local, cram, tmp_path):
    """Every production input is a CRAM: this runs the fields the engine asks the CRAM decoder for, the .crai
    index and the CRAM end-of-file check. (The engine insists on -T or REF_PATH even for a CRAM that needs no
    reference, rather than let htslib go to the network for one; an empty REF_PATH directory satisfies it.)"""
    env = environ(REF_PATH=f"{cram}/noref/%s", REF_CACHE=f"{cram}/noref/%s")
    for mode in ("scan", "fetch"):
        ok(count(cram / "fx.cram", mode, tmp_path / f"{mode}.json", env=env))
        c = io.load_counts(tmp_path / f"{mode}.json")
        assert strip(c) == strip(local[mode])
        assert c["pipeline"]["sq_sha256"] == local[mode]["pipeline"]["sq_sha256"]


def test_closing_a_cram_reader_touches_no_freed_memory(local, cram, tmp_path):
    """A fetch closes one indexed reader per worker while the other workers still count. A CRAM's index belongs to its
    file handle, and the htslib binding destroys the index after it has closed the file: htslib then reads the freed
    handle. Left alone that memory still reads as it did, and nothing happens; once in a few hundred runs another
    thread has taken it, and the run aborts. With freed memory overwritten (macOS: MallocScribble; glibc:
    MALLOC_PERTURB_) the read fails every time, so this is a test and not a matter of luck: the engine's reader must
    never get there, whatever the number of workers."""
    env = environ(REF_PATH=f"{cram}/noref/%s", REF_CACHE=f"{cram}/noref/%s", MallocScribble="1", MALLOC_PERTURB_="85")
    for threads in (1, 4):
        out = tmp_path / f"t{threads}.json"
        r = count(cram / "fx.cram", "fetch", out, threads=threads, env=env)
        assert r.returncode == 0, (threads, r.returncode, r.stderr[-400:])
        assert strip(io.load_counts(out)) == strip(local["fetch"])


def test_a_truncated_file_is_refused(cram, tmp_path):
    """A copy that stopped early lacks the end-of-file marker: refused, with no counts written, unless
    --allow-truncated, which records the marker as absent."""
    (tmp_path / "cut.bam").write_bytes(bgzf_cut(BAM.read_bytes()))
    data = (cram / "fx.cram").read_bytes()
    (tmp_path / "cut.cram").write_bytes(data[:len(data) // 2])
    env = environ(REF_PATH=f"{cram}/noref/%s", REF_CACHE=f"{cram}/noref/%s")
    for name in ("cut.bam", "cut.cram"):
        out = tmp_path / f"{name}.json"
        r = count(tmp_path / name, "scan", out, env=env)
        assert r.returncode == 1 and "end-of-file marker" in r.stderr and "--allow-truncated" in r.stderr and not out.exists(), r.stderr
    ok(count(tmp_path / "cut.bam", "scan", tmp_path / "cut.json", "--allow-truncated"))
    assert io.load_counts(tmp_path / "cut.json")["eof_marker"] == "absent"


def test_a_reference_cram_is_decoded_with_T_or_REF_PATH(sim, tmp_path):
    """A CRAM stores reads as differences from its reference, found through -T or through REF_PATH (files
    named by the @SQ M5). Both give the BAM's counts. With neither the engine refuses before decoding
    anything, and a -T FASTA that is not the CRAM's reference is refused at once, not retried."""
    d = sim["dir"]
    enc = tmp_path / "enc"
    enc.mkdir()
    shutil.copy(d / "ref.fa", enc)
    samtools("view", "-C", "-T", enc / "ref.fa", "-o", tmp_path / "sim.cram", d / "sim.bam")
    samtools("index", tmp_path / "sim.cram")
    shutil.rmtree(enc)                                    # the @SQ UR tags now lead nowhere: -T or REF_PATH must serve
    cache = tmp_path / "cache"
    cache.mkdir()
    seqs = {}
    for line in (d / "ref.fa").read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip().upper())
    for s in seqs.values():
        s = "".join(s).encode()
        (cache / hashlib.md5(s).hexdigest()).write_bytes(s)
    common = ["-p", d / "panel.tsv.gz", "-c", d / "controls.fa.gz", "--l-grid", "100,200,300,400", "-@", "2"]
    ways = {"-T": (["-T", d / "ref.fa"], environ()),
            "REF_PATH": ([], environ(REF_PATH=f"{cache}/%s", REF_CACHE=f"{cache}/%s"))}
    for how, (extra, env) in ways.items():
        for mode in ("scan", "fetch"):
            out = tmp_path / f"{mode}.json"
            sinks = ["--sinks", d / "sinks.bed"] if mode == "fetch" else []
            ok(engine("count", "-i", tmp_path / "sim.cram", *common, "-m", mode, *sinks, *extra, "-o", out, env=env))
            c = io.load_counts(out)
            assert strip(c) == strip(sim[mode]), how
            assert c["pipeline"]["sq_m5"] == c["pipeline"]["sq_n"] == 3
    out = tmp_path / "none.json"
    # (through a dead proxy, so that an engine which let htslib go to the EBI reference server fails fast)
    r = engine("count", "-i", tmp_path / "sim.cram", *common, "-m", "scan", "-o", out,
               env=environ(http_proxy="http://127.0.0.1:9", https_proxy="http://127.0.0.1:9", REF_CACHE=f"{tmp_path}/enc/%s"))
    assert r.returncode == 1 and "neither -T nor a non-empty REF_PATH" in r.stderr and not out.exists(), r.stderr
    # an empty REF_PATH is unset to htslib (it would go to the server): refused the same way
    r = engine("count", "-i", tmp_path / "sim.cram", *common, "-m", "scan", "-o", out,
               env=environ(REF_PATH="", http_proxy="http://127.0.0.1:9", https_proxy="http://127.0.0.1:9", REF_CACHE=f"{tmp_path}/enc/%s"))
    assert r.returncode == 1 and "neither -T nor a non-empty REF_PATH" in r.stderr and not out.exists(), r.stderr
    lines = (d / "ref.fa").read_text().splitlines()
    i = lines.index(">ctrl") + 1
    lines[i] = lines[i][:150_000] + ("A" if lines[i][150_000] != "A" else "C") + lines[i][150_001:]
    (tmp_path / "mut.fa").write_text("\n".join(lines) + "\n")
    samtools("faidx", tmp_path / "mut.fa")
    env = environ(REF_PATH=f"{tmp_path}/enc/%s", REF_CACHE=f"{tmp_path}/enc/%s")
    for mode in ("scan", "fetch"):
        sinks = ["--sinks", d / "sinks.bed"] if mode == "fetch" else []
        r = engine("count", "-i", tmp_path / "sim.cram", *common, "-m", mode, *sinks, "-T", tmp_path / "mut.fa", "-o", out, env=env)
        assert r.returncode == 1 and "does not match the CRAM on ctrl" in r.stderr and "retry" not in r.stderr, r.stderr
        assert not out.exists()


@pytest.fixture(scope="module")
def unplaced(tmp_path_factory):
    """The fixture with pairs the aligner could not place (flag 77/141, no coordinate), which sit in the
    unmapped section at the end of a file: the fixture has none, so a few hundred carrying rDNA sequence are
    added. Returns the BAM (indexed) and the number of reads added."""
    d = tmp_path_factory.mktemp("unplaced")
    rng = np.random.default_rng(1)
    rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    sam = []
    units = BUNDLE.units()
    for cls, n in (("rDNA45S", 200), ("rDNA5S", 100)):
        u = units[cls].upper()
        for i in range(n):
            s = int(rng.integers(0, len(u) - 450))
            sam += [f"un_{cls}_{i}\t77\t*\t0\t0\t*\t*\t0\t0\t{u[s:s + 150]}\t*\tRG:Z:N",
                    f"un_{cls}_{i}\t141\t*\t0\t0\t*\t*\t0\t0\t{rc(u[s + 300:s + 450])}\t*\tRG:Z:N"]
    whole = samtools("view", "-h", BAM).stdout + ("\n".join(sam) + "\n").encode()
    samtools("sort", "-o", d / "un.bam", "-", input=whole)
    samtools("index", "-c", d / "un.bam")
    return d / "un.bam", len(sam)


def test_reads_without_a_coordinate_are_counted_by_scan_and_by_fetch_unmapped(local, unplaced, tmp_path):
    """A scan classifies the reads without a coordinate; a fetch reads them only with --unmapped, and then
    equals the scan."""
    un, added = unplaced
    s, f, u = (ok(count(un, m, tmp_path / f"{n}.json", *x)) and io.load_counts(tmp_path / f"{n}.json")
               for n, m, x in (("s", "scan", []), ("f", "fetch", []), ("u", "fetch", ["--unmapped"])))
    assert strip(s, SCAN_FETCH_DIFFER | {"unmapped_fetched"}) == strip(u, SCAN_FETCH_DIFFER | {"unmapped_fetched"})
    assert strip(f) == strip(local["fetch"])                                              # plain fetch: as if they were not there
    assert (s["unmapped_fetched"], f["unmapped_fetched"], u["unmapped_fetched"]) == (False, False, True)
    star = lambda c: {p["class"]: p["reads"] for p in c["placements"] if p["contig"] == "*"}
    assert star(s) == star(u) and star(f) == {} and set(star(s)) == {"rDNA45S", "rDNA5S"}
    assert 0 < sum(star(s).values()) <= added
    assert {k: v - reads(local["scan"])[k] for k, v in reads(s).items() if v != reads(local["scan"])[k]} == star(s)
    assert u["unmapped"] - f["unmapped"] == u["primary"] - f["primary"] == added
    assert sum(c["reads"] for c in s["contigs"]) == s["primary"] - s["unmapped"]         # contig tallies leave them out


def test_a_file_without_reads_is_refused(tmp_path):
    """No control read means nothing can be estimated: an empty file, the wrong contig names, or an index
    of another file. The engine refuses rather than write counts of zero."""
    samtools("view", "-H", "-b", "-o", tmp_path / "h.bam", BAM)
    samtools("index", "-c", tmp_path / "h.bam")
    for mode in ("scan", "fetch"):
        out = tmp_path / f"{mode}.json"
        r = count(tmp_path / "h.bam", mode, out)
        assert r.returncode == 1 and "no read of" in r.stderr and "control region" in r.stderr and not out.exists(), r.stderr


@pytest.mark.parametrize("bad", [["--bin", "0"], ["--place-bin", "0"], ["--place-bin=-1"], ["--min-frac", "2"], ["--min-frac=-0.1"],
                                 ["--min-hits", "0"], ["--retries", "0"], ["--l-grid", "100,0"], ["--pad=-5"], ["--pad", "399"]])
def test_out_of_range_parameters_are_refused_by_the_parser(bad, tmp_path):
    out = tmp_path / "o.json"
    r = count(BAM, "fetch", out, *bad)
    assert r.returncode == 2 and not out.exists(), r.stderr


class Server:
    """A threaded HTTP server on 127.0.0.1 over one directory. ranges=False ignores Range headers, as some
    servers do. fail(n, path, range) -> True answers request n (1-based, over all requests) with 503;
    etag(n) -> str sends that entity tag with request n's answer. `log` keeps (path with query, range,
    status) for every request."""

    def __init__(self, root, ranges=True, fail=None, etag=None):
        self.log, lock = [], threading.Lock()
        outer = self

        class Handler(BaseHTTPRequestHandler):
            protocol_version = "HTTP/1.1"

            def log_message(self, *a):
                pass

            def do_HEAD(self):
                self.answer(head=True)

            def do_GET(self):
                self.answer()

            def answer(self, head=False):
                path = self.path.split("?")[0].lstrip("/")
                rng = self.headers.get("Range") if ranges else None
                with lock:
                    n = len(outer.log) + 1
                    bad = bool(fail and fail(n, path, rng))
                    outer.log.append((self.path, rng, 503 if bad else None))
                f = Path(root) / path
                code, body, hdr = 503, b"", []
                if not bad:
                    code = 404
                    if f.is_file():
                        data = f.read_bytes()
                        hdr = [("Accept-Ranges", "bytes")] if ranges else []
                        code, body = 200, data
                        if rng and rng.startswith("bytes="):
                            s, e = rng[6:].split("-")
                            s, e = (len(data) - int(e), len(data) - 1) if s == "" else (int(s), min(int(e or len(data) - 1), len(data) - 1))
                            code, body = (206, data[s:e + 1]) if s < len(data) else (416, b"")
                            hdr.append(("Content-Range", f"bytes {s}-{e}/{len(data)}" if code == 206 else f"bytes */{len(data)}"))
                        if etag:
                            hdr.append(("ETag", etag(n)))
                self.send_response(code)
                for k, v in hdr:
                    self.send_header(k, v)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                if not head:
                    try:
                        self.wfile.write(body)
                    except (BrokenPipeError, ConnectionResetError):          # a reader stops reading when it has enough
                        self.close_connection = True

        class Quiet(ThreadingHTTPServer):
            daemon_threads = True

            def handle_error(self, request, client_address):
                pass

        self.httpd = Quiet(("127.0.0.1", 0), Handler)
        self.url = f"http://127.0.0.1:{self.httpd.server_address[1]}"

    def __enter__(self):
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()
        return self

    def __exit__(self, *a):
        self.httpd.shutdown()
        self.httpd.server_close()

    def ranges(self, name):
        """The byte ranges asked for of one file (not its index), in order: (start, end exclusive), (start, None)
        for an open-ended request, (None, None) for a request of the whole file."""
        out = []
        for path, rng, _ in self.log:
            if path.split("?")[0] == "/" + name:
                if not rng:
                    out.append((None, None))
                    continue
                a, b = rng.removeprefix("bytes=").split("-")
                out.append((int(a), int(b) + 1 if b else None))
        return out


# `count --transport`: the engine's exact byte ranges (the default), and htslib's own reader
TRANSPORTS = ("ranges", "htslib")
# the engine's line on what a fetch by exact byte ranges moved
MOVED = re.compile(r"moved ([\d.]+) MB in .*? (\d+) reads outside the plan")


@pytest.fixture(scope="module")
def www(tmp_path_factory):
    d = tmp_path_factory.mktemp("www")
    shutil.copy(BAM, d / "fx.bam")
    shutil.copy(BAM.with_name(BAM.name + ".csi"), d / "fx.bam.csi")
    (d / "cut.bam").write_bytes(bgzf_cut(BAM.read_bytes()))
    shutil.copy(BAM, d / "noidx.bam")
    # the fixture as a CRAM that needs no reference (REF_PATH at the empty noref/ serves), in small slices
    samtools("view", "-C", "--output-fmt-option", "no_ref=1", "--output-fmt-option", "seqs_per_slice=1000", "-o", d / "fx.cram", BAM)
    samtools("index", d / "fx.cram")
    (d / "noref").mkdir()
    return d


def cram_env(www, **kw):
    return environ(REF_PATH=f"{www}/noref/%s", REF_CACHE=f"{www}/noref/%s", **kw)


@pytest.fixture
def wd(tmp_path):
    """A fresh working directory per run: htslib saves a remote BAM index there and reuses it unchecked."""
    (tmp_path / "wd").mkdir()
    return tmp_path / "wd"


@pytest.mark.parametrize("transport", TRANSPORTS)
def test_counts_over_http_equal_local_counts(local, www, wd, transport):
    with Server(www) as srv:
        for mode in ("scan", "fetch"):
            r = ok(count(srv.url + "/fx.bam", mode, wd / f"{mode}.json", "--transport", transport, cwd=wd))
            c = io.load_counts(wd / f"{mode}.json")
            assert strip(c, {"input", "elapsed_sec"}) == strip(local[mode], {"input", "elapsed_sec"})
            assert c["input"] == srv.url + "/fx.bam" and "retry" not in r.stderr


@pytest.mark.parametrize("name", ["fx.bam", "fx.cram"])
def test_a_remote_fetch_asks_for_exactly_the_bytes_its_queries_read(local, www, wd, name):
    """`--transport ranges`, the default (engine 0.4.0): the engine reads the index, works out which bytes each
    query reads (a CRAM's containers, the BGZF blocks of a BAM's chunks), asks for exactly those with bounded
    range requests, each byte once, and serves them to htslib, which decodes them as it would the file: the
    counts are the local file's, at one and at four workers, with freed memory overwritten (macOS
    MallocScribble, glibc MALLOC_PERTURB_). No request is open-ended, no byte is asked for twice, nothing is
    read outside the plan. htslib's own reader asks for an open-ended range at every query, and counts the same."""
    env = cram_env(www, MallocScribble="1", MALLOC_PERTURB_="85")
    size = (www / name).stat().st_size
    with Server(www) as srv:
        for threads in (1, 4):
            srv.log.clear()
            r = ok(count(srv.url + "/" + name, "fetch", wd / "r.json", threads=threads, cwd=wd, env=env))
            assert strip(io.load_counts(wd / "r.json")) == strip(local["fetch"])
            got = srv.ranges(name)
            assert got and all(a is not None and b is not None for a, b in got), got
            spans = sorted(got)
            assert all(b1 <= a2 for (_, b1), (a2, _) in zip(spans, spans[1:])), f"a byte asked for twice: {spans}"
            assert sum(b - a for a, b in got) <= size
            m = MOVED.search(r.stderr)
            assert m and m[2] == "0", r.stderr
        srv.log.clear()
        ok(count(srv.url + "/" + name, "fetch", wd / "h.json", "--transport", "htslib", cwd=wd, env=env))
        assert strip(io.load_counts(wd / "h.json")) == strip(local["fetch"])
        assert any(b is None for _, b in srv.ranges(name)), "htslib asks for open-ended ranges"


@pytest.mark.parametrize("name", ["fx.bam", "fx.cram"])
def test_a_remote_fetch_of_a_few_regions_moves_only_their_bytes(www, wd, tmp_path, name):
    """A plan that needs a small part of the file - the base control regions of chr1 and the 5S rDNA sinks - moves
    that part: the containers of a CRAM that htslib's iterator visits (and the header of the one after), the BGZF
    blocks of a BAM's chunks (the block a chunk ends inside read to its end, whose length only its header gives),
    and the counts are the local file's."""
    rows = [r.split("\t") for r in (BUNDLE.dir / "controls.base.bed").read_text().splitlines() if r.startswith("chr1\t")]
    fetchplan.subset_fasta(BUNDLE.controls, [f"{r[0]}:{r[1]}-{r[2]}" for r in rows], tmp_path / "chr1.fa.gz")
    env = cram_env(www)
    args = ["-p", BUNDLE.panel, "-c", tmp_path / "chr1.fa.gz", "--sinks", BUNDLE.sinks, "--classes", "rDNA5S", "-m", "fetch", "-@", "2"]
    ok(engine("count", "-i", www / name, *args, "-o", tmp_path / "local.json", env=env))
    size = (www / name).stat().st_size
    with Server(www) as srv:
        r = ok(engine("count", "-i", srv.url + "/" + name, *args, "-o", wd / "r.json", cwd=wd, env=env))
    want = io.load_counts(tmp_path / "local.json")
    assert want["ctrl_reads"] > 0 and reads(want)["rDNA5S"] > 0
    assert strip(io.load_counts(wd / "r.json"), {"input", "elapsed_sec"}) == strip(want, {"input", "elapsed_sec"})
    m = MOVED.search(r.stderr)
    assert m and m[2] == "0", r.stderr
    moved = sum(b - a for a, b in srv.ranges(name))
    assert moved < 0.3 * size, (moved, size, r.stderr)


def test_the_unmapped_bin_is_read_by_exact_ranges_too(unplaced, wd, tmp_path):
    """`--unmapped` over HTTP: the reads without a coordinate, from the first of them to the end of the file, as
    the local file gives them."""
    un, _ = unplaced
    d = tmp_path / "www"
    d.mkdir()
    shutil.copy(un, d / "un.bam")
    shutil.copy(un.with_suffix(".bam.csi"), d / "un.bam.csi")
    ok(count(un, "fetch", tmp_path / "local.json", "--unmapped"))
    with Server(d) as srv:
        r = ok(count(srv.url + "/un.bam", "fetch", wd / "u.json", "--unmapped", cwd=wd))
    c = io.load_counts(wd / "u.json")
    assert c["unmapped_fetched"] and c["unmapped"] > 0
    assert strip(c, {"input", "elapsed_sec"}) == strip(io.load_counts(tmp_path / "local.json"), {"input", "elapsed_sec"})
    assert MOVED.search(r.stderr)[2] == "0", r.stderr


def eof_request(www):
    return f"bytes={(www / 'fx.bam').stat().st_size - 28}-"


def lose(www, what):
    """A rule for Server(fail=...) that loses the requests of one kind of htslib's reader, each once."""
    eof, seen = eof_request(www), {"csi": False, "data": 0, "open": 0, "failed": set()}

    def once(key, cond):
        if cond and key not in seen["failed"]:
            seen["failed"].add(key)
            return True
        return False

    def rule(n, path, rng):
        if path.endswith(".csi"):
            seen["csi"] = True
            return what == "index" and once("index", True)
        is_open = rng is None
        seen["open"] += is_open
        data = seen["csi"] and rng not in (None, eof)
        seen["data"] += data
        if what == "probe":        # the probe's open, then its end-of-file check
            return once("open", n == 1) or once("eof", rng == eof)
        if what == "scan":         # the scan's own open, after the probe's
            return once("open", is_open and seen["open"] == 2)
        if what == "workers":      # a worker's first open; a read of an interval; the reopen after that read failed
            return (once("worker", seen["csi"] and is_open) or once("read", data and seen["data"] == 20)
                    or once("reopen", is_open and "read" in seen["failed"]))
        return False
    return rule


@pytest.mark.parametrize("what", ["probe", "scan", "workers"])
def test_lost_requests_are_retried(local, www, wd, what):
    """A 503 on one request of htslib's reader - the probe's open or its end-of-file check, the scan's open, a
    fetch worker's open, the read of an interval and the reopen after it - costs one attempt, not the run."""
    mode = "scan" if what == "scan" else "fetch"
    with Server(www, fail=lose(www, what)) as srv:
        r = ok(count(srv.url + "/fx.bam", mode, wd / "o.json", "--transport", "htslib", threads=1, cwd=wd))
    lost = [x for x in srv.log if x[2] == 503]
    assert len(lost) == {"probe": 2, "scan": 1, "workers": 3}[what]
    msg = own(r.stderr)
    if what == "workers":
        assert re.search(r"retry 1/5 for interval \S+: cannot open", msg), msg
        m = re.search(r"retry 1/5 for interval (\S+): read error", msg)
        assert m and f"retry 2/5 for interval {m[1]}: cannot open" in msg, msg
    else:
        assert "retry 1/5 to open the input" in msg, msg
        assert what == "scan" or "retry 2/5 to open the input: the end-of-file check" in msg, msg
    c = io.load_counts(wd / "o.json")
    assert strip(c, {"input", "elapsed_sec"}) == strip(local[mode], {"input", "elapsed_sec"})
    assert c["eof_marker"] == "present"


def test_lost_range_requests_are_retried(local, www, wd):
    """By exact ranges, a 503 on any request - the first bytes, the index, the end-of-file marker, a planned range -
    costs one attempt, not the run."""
    size, seen = (www / "fx.bam").stat().st_size, set()

    def rule(n, path, rng):
        kind = ("index" if path.endswith(".csi") else "first" if rng == "bytes=0-65535" else
                "tail" if rng == f"bytes={size - 64}-{size - 1}" else "range" if rng else None)
        if kind and kind not in seen:
            seen.add(kind)
            return True
        return False
    with Server(www, fail=rule) as srv:
        r = ok(count(srv.url + "/fx.bam", "fetch", wd / "o.json", threads=1, cwd=wd))
    assert seen == {"first", "index", "tail", "range"} and len([x for x in srv.log if x[2] == 503]) == 4
    msg = own(r.stderr)
    assert "retry 1/5 to open the input: HTTP 503" in msg and "retry 1/5 to download the index: HTTP 503" in msg, msg
    # the end-of-file marker's range is retried when the input is opened; a planned range by the connection that asked for it,
    # or, when a worker asked for it first, by that worker's interval
    assert len(re.findall(r"retry 1/5 for .*HTTP 503 for bytes \d+-\d+ of ", msg)) == 2, msg
    c = io.load_counts(wd / "o.json")
    assert strip(c, {"input", "elapsed_sec"}) == strip(local["fetch"], {"input", "elapsed_sec"})
    assert MOVED.search(r.stderr)[2] == "0" and c["eof_marker"] == "present"


@pytest.mark.parametrize("transport", TRANSPORTS)
def test_a_lost_index_request_is_retried(local, www, wd, transport):
    """A 503 on the index download is transient: the open is retried and the counts equal the local ones.
    (An index that is really absent is not retried; see test_a_remote_input_without_an_index_is_refused_at_once.)"""
    with Server(www, fail=lose(www, "index")) as srv:
        r = count(srv.url + "/fx.bam", "fetch", wd / "o.json", "--transport", transport, threads=1, cwd=wd)
    assert [x[0] for x in srv.log if x[2] == 503] == ["/fx.bam.csi"]
    assert r.returncode == 0 and "retry 1/5" in r.stderr, r.stderr
    c = io.load_counts(wd / "o.json")
    assert strip(c, {"input", "elapsed_sec"}) == strip(local["fetch"], {"input", "elapsed_sec"})


@pytest.mark.parametrize("transport", TRANSPORTS)
def test_a_remote_input_without_an_index_is_refused_at_once(www, wd, transport):
    """No index next to the file on the server is not a lost request: the fetch exits 1 without retrying
    (not 75, which would have the workflow try the sample again and again)."""
    out = wd / "o.json"
    with Server(www) as srv:
        r = count(srv.url + "/noidx.bam", "fetch", out, "--transport", transport, threads=1, cwd=wd)
    assert r.returncode == 1 and not out.exists(), r.stderr
    assert "no index beside it" in own(r.stderr) and "retry" not in own(r.stderr), r.stderr


@pytest.mark.parametrize("transport", TRANSPORTS)
@pytest.mark.parametrize("start", ["probe", "interval"])
def test_a_remote_input_that_keeps_failing_exits_75(www, wd, start, transport):
    """When every retry fails (a server down from the start, or from the middle of a fetch), the run stops
    with EX_TEMPFAIL, 75, and writes nothing: the workflow retries the sample later."""
    eof, n_data = eof_request(www), [0]

    def down(n, path, rng):
        if start == "probe":
            return True
        if transport == "ranges":             # the first bytes and the end-of-file marker come; every planned range fails
            n_data[0] += bool(rng)
            return n_data[0] > 2
        n_data[0] += rng not in (None, eof)
        return n_data[0] >= 20
    out = wd / "o.json"
    with Server(www, fail=down) as srv:
        r = count(srv.url + "/fx.bam", "fetch", out, "--retries", "2", "--transport", transport, threads=1, cwd=wd)
    assert r.returncode == 75 and "still failing after every retry" in r.stderr and not out.exists(), r.stderr


def test_a_server_without_range_requests(local, www, wd):
    """Without range requests the end-of-file marker cannot be checked up front, so a scan checks it at the
    end of the stream: the whole file passes, a copy cut at a block boundary is refused."""
    with Server(www, ranges=False) as srv:
        ok(count(srv.url + "/fx.bam", "scan", wd / "full.json", cwd=wd))
        c = io.load_counts(wd / "full.json")
        assert strip(c, {"input", "elapsed_sec"}) == strip(local["scan"], {"input", "elapsed_sec"})
        assert c["eof_marker"] == "present"
        out = wd / "cut.json"
        r = count(srv.url + "/cut.bam", "scan", out, cwd=wd)
        assert r.returncode == 1 and "ended without an end-of-file marker" in r.stderr and not out.exists(), r.stderr
        r = ok(count(srv.url + "/cut.bam", "scan", out, "--allow-truncated", cwd=wd))
        assert "WARNING" in r.stderr and io.load_counts(out)["eof_marker"] == "absent"


@pytest.mark.parametrize("transport", TRANSPORTS)
def test_a_fetch_from_a_server_without_range_requests_is_refused_at_once(www, wd, transport):
    """A fetch needs range requests: a server that answers a range with the whole file is refused with exit 1 (no
    retry can change it), and by exact ranges after the first 64 kB of that answer, not the whole file."""
    out = wd / "o.json"
    with Server(www, ranges=False) as srv:
        r = count(srv.url + "/fx.bam", "fetch", out, "--transport", transport, cwd=wd)
    assert r.returncode == 1 and "does not answer range requests" in r.stderr and not out.exists(), r.stderr


def test_a_file_that_changes_during_a_fetch_is_refused(www, wd):
    """The server's entity tag is the file's identity: one that changes between the first answer and a later one (a
    file replaced while it was read) ends the fetch with exit 1 - its bytes would come from two files."""
    out = wd / "o.json"
    with Server(www, etag=lambda n: '"one"' if n <= 3 else '"two"') as srv:
        r = count(srv.url + "/fx.bam", "fetch", out, cwd=wd)
    assert r.returncode == 1 and "changed during the fetch" in r.stderr and not out.exists(), r.stderr


def test_the_downloaded_bytes_need_a_directory_that_can_hold_them(www, wd):
    out = wd / "o.json"
    with Server(www) as srv:
        r = count(srv.url + "/fx.bam", "fetch", out, "--spool-dir", wd / "absent", cwd=wd)
    assert r.returncode == 1 and "--spool-dir" in r.stderr and "retry" not in r.stderr and not out.exists(), r.stderr


SILENCED = "htslib's own messages are silenced for a URL with a query string"


@pytest.mark.parametrize("transport", TRANSPORTS)
def test_a_signed_url_is_kept_out_of_the_output(local, www, wd, transport):
    """A signed URL's query is a credential. It reaches the server (for the file and its index) but not
    the counts or the engine's messages, which carry the URL redacted. htslib's own error lines would print
    it in full on every failed open, index search and retry, so for a URL with a query string the engine
    turns them off before its first open and says so once (engine 0.1.1); a URL without a query keeps them.
    (By exact ranges htslib never sees the URL.)"""
    q = "?X-Amz-Credential=KEYSECRET&X-Amz-Signature=deadbeefSECRET"
    t = ["--transport", transport]
    with Server(www) as srv:
        for mode in ("scan", "fetch"):
            r = ok(count(srv.url + "/fx.bam" + q, mode, wd / f"{mode}.json", *t, cwd=wd))
            text = (wd / f"{mode}.json").read_text()
            c = json.loads(text)
            assert "SECRET" not in text and "SECRET" not in r.stderr, r.stderr
            assert c["input"] == srv.url + "/fx.bam?<redacted>" and r.stderr.count(SILENCED) == 1
            assert strip(c, {"input", "elapsed_sec"}) == strip(local[mode], {"input", "elapsed_sec"})
        assert {p for p, _, _ in srv.log} == {"/fx.bam" + q, "/fx.bam.csi" + q}
        n = len(srv.log)
        r = count(srv.url + "/missing.bam" + q, "scan", wd / "m.json", cwd=wd)
        assert r.returncode == 1 and len(srv.log) == n + 1, r.stderr                 # a 404 is final: no retry, no 75
        assert "SECRET" not in r.stderr and "missing.bam?<redacted>" in own(r.stderr) and SILENCED in r.stderr, r.stderr
        n = len(srv.log)
        r = count(srv.url + "/missing.bam" + q, "fetch", wd / "m.json", *t, cwd=wd)
        assert r.returncode == 1 and len(srv.log) == n + 1, r.stderr                 # so it is for a fetch
        assert "SECRET" not in r.stderr and "missing.bam?<redacted>" in own(r.stderr), r.stderr
        # no index beside the file: htslib's index search names the URL in every message it prints
        r = count(srv.url + "/noidx.bam" + q, "fetch", wd / "n.json", *t, threads=1, cwd=wd)
        assert r.returncode == 1 and "SECRET" not in r.stderr and "noidx.bam?<redacted>" in r.stderr, r.stderr
        assert "no index beside it" in r.stderr and SILENCED in r.stderr, r.stderr
        # the index as a signed URL of its own, the file without one
        r = ok(count(srv.url + "/fx.bam", "fetch", wd / "i.json", "--index", srv.url + "/fx.bam.csi" + q, *t, cwd=wd))
        assert "SECRET" not in r.stderr and SILENCED in r.stderr, r.stderr
        assert strip(io.load_counts(wd / "i.json"), {"input", "elapsed_sec"}) == strip(local["fetch"], {"input", "elapsed_sec"})
    # a server that keeps failing: every attempt's open fails with the URL in htslib's message
    with Server(www, fail=lambda n, path, rng: True) as srv:
        r = count(srv.url + "/fx.bam" + q, "fetch", wd / "d.json", "--retries", "2", *t, cwd=wd)
        assert r.returncode == 75 and "retry 1/2" in r.stderr and "SECRET" not in r.stderr, r.stderr
        assert "fx.bam?<redacted>" in r.stderr and SILENCED in r.stderr, r.stderr
    # and a failed open of a URL without a query string keeps htslib's own line
    with Server(www) as srv:
        r = count(srv.url + "/missing.bam", "scan", wd / "p.json", cwd=wd)
        assert r.returncode == 1 and "[E::" in r.stderr and SILENCED not in r.stderr, r.stderr
