"""The engine's inputs and how it fails on them: CRAM (reference from -T or REF_PATH), reads without a
coordinate, a file served over HTTP (with and without range requests, with lost requests, behind a
signed URL), sinks on contigs the input lacks, and inputs it must refuse. Needs the built binary and
samtools; every input is made in the test from the committed fixture or the simulated genome."""
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

from ngsdose import contract, io, resources

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


def test_sinks_on_contigs_the_input_lacks_are_left_out_and_recorded(local, tmp_path):
    """Sinks learned on another reference or pipeline can name contigs the file does not have. The fetch
    leaves those intervals out, says so, and records what it left out per class."""
    bed = tmp_path / "sinks.bed"
    bed.write_text(Path(BUNDLE.sinks).read_text() + "chrAbsent_decoy\t1000\t6000\tDJ\n")
    r = ok(count(BAM, "fetch", tmp_path / "f.json", "--sinks", bed))
    assert "WARNING" in r.stderr and "sinks_skipped" in r.stderr
    c = io.load_counts(tmp_path / "f.json")
    assert c["sinks_skipped"] == {"DJ": {"intervals": 1, "bp": 5000}}
    assert "sinks_missing_classes" not in c
    drop = RUN | {"sinks", "sinks_sha256", "sinks_skipped"}
    assert strip(c, drop) == strip(local["fetch"], drop)
    assert set(contract.incomplete_sinks(c)) == {"DJ"} and contract.incomplete_sinks(local["fetch"]) == {}
    r = ok(engine("plan", "-c", BUNDLE.controls, "--sinks", bed, "-i", BAM, "-o", tmp_path / "plan.bed"))
    assert "sink intervals by class: DJ 1" in r.stderr


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


def test_reads_without_a_coordinate_are_counted_by_scan_and_by_fetch_unmapped(local, tmp_path):
    """Pairs the aligner could not place (flag 77/141, no coordinate) sit in the unmapped section at the end
    of a file; the fixture has none, so a few hundred carrying rDNA sequence are added. A scan classifies
    them; a fetch reads them only with --unmapped, and then equals the scan."""
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
    samtools("sort", "-o", tmp_path / "un.bam", "-", input=whole)
    samtools("index", "-c", tmp_path / "un.bam")
    s, f, u = (ok(count(tmp_path / "un.bam", m, tmp_path / f"{n}.json", *x)) and io.load_counts(tmp_path / f"{n}.json")
               for n, m, x in (("s", "scan", []), ("f", "fetch", []), ("u", "fetch", ["--unmapped"])))
    assert strip(s, SCAN_FETCH_DIFFER | {"unmapped_fetched"}) == strip(u, SCAN_FETCH_DIFFER | {"unmapped_fetched"})
    assert strip(f) == strip(local["fetch"])                                              # plain fetch: as if they were not there
    assert (s["unmapped_fetched"], f["unmapped_fetched"], u["unmapped_fetched"]) == (False, False, True)
    star = lambda c: {p["class"]: p["reads"] for p in c["placements"] if p["contig"] == "*"}
    assert star(s) == star(u) and star(f) == {} and set(star(s)) == {"rDNA45S", "rDNA5S"}
    assert 0 < sum(star(s).values()) <= len(sam)
    assert {k: v - reads(local["scan"])[k] for k, v in reads(s).items() if v != reads(local["scan"])[k]} == star(s)
    assert u["unmapped"] - f["unmapped"] == u["primary"] - f["primary"] == len(sam)
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
    servers do. fail(n, path, range) -> True answers request n (1-based, over all requests) with 503.
    `log` keeps (path with query, range, status) for every request."""

    def __init__(self, root, ranges=True, fail=None):
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
                    bad = bool(fail and fail(len(outer.log) + 1, path, rng))
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
                self.send_response(code)
                for k, v in hdr:
                    self.send_header(k, v)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                if not head:
                    try:
                        self.wfile.write(body)
                    except (BrokenPipeError, ConnectionResetError):          # htslib stops reading when it has enough
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


@pytest.fixture(scope="module")
def www(tmp_path_factory):
    d = tmp_path_factory.mktemp("www")
    shutil.copy(BAM, d / "fx.bam")
    shutil.copy(BAM.with_name(BAM.name + ".csi"), d / "fx.bam.csi")
    (d / "cut.bam").write_bytes(bgzf_cut(BAM.read_bytes()))
    shutil.copy(BAM, d / "noidx.bam")
    return d


@pytest.fixture
def wd(tmp_path):
    """A fresh working directory per run: htslib saves a remote BAM index there and reuses it unchecked."""
    (tmp_path / "wd").mkdir()
    return tmp_path / "wd"


def test_counts_over_http_equal_local_counts(local, www, wd):
    with Server(www) as srv:
        for mode in ("scan", "fetch"):
            r = ok(count(srv.url + "/fx.bam", mode, wd / f"{mode}.json", cwd=wd))
            c = io.load_counts(wd / f"{mode}.json")
            assert strip(c, {"input", "elapsed_sec"}) == strip(local[mode], {"input", "elapsed_sec"})
            assert c["input"] == srv.url + "/fx.bam" and "retry" not in r.stderr


def eof_request(www):
    return f"bytes={(www / 'fx.bam').stat().st_size - 28}-"


def lose(www, what):
    """A rule for Server(fail=...) that loses the requests of one kind, each once."""
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
    """A 503 on one request of a remote run - the probe's open or its end-of-file check, the scan's open, a
    fetch worker's open, the read of an interval and the reopen after it - costs one attempt, not the run."""
    mode = "scan" if what == "scan" else "fetch"
    with Server(www, fail=lose(www, what)) as srv:
        r = ok(count(srv.url + "/fx.bam", mode, wd / "o.json", threads=1, cwd=wd))
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


def test_a_lost_index_request_is_retried(local, www, wd):
    """A 503 on the index download is transient: the open is retried and the counts equal the local ones.
    (An index that is really absent is not retried; see test_a_remote_input_without_an_index_is_refused_at_once.)"""
    with Server(www, fail=lose(www, "index")) as srv:
        r = count(srv.url + "/fx.bam", "fetch", wd / "o.json", threads=1, cwd=wd)
    assert [x[0] for x in srv.log if x[2] == 503] == ["/fx.bam.csi"]
    assert r.returncode == 0 and "retry 1/5" in r.stderr, r.stderr
    c = io.load_counts(wd / "o.json")
    assert strip(c, {"input", "elapsed_sec"}) == strip(local["fetch"], {"input", "elapsed_sec"})


def test_a_remote_input_without_an_index_is_refused_at_once(www, wd):
    """No index next to the file on the server is not a lost request: the fetch exits 1 without retrying
    (not 75, which would have the workflow try the sample again and again)."""
    out = wd / "o.json"
    with Server(www) as srv:
        r = count(srv.url + "/noidx.bam", "fetch", out, threads=1, cwd=wd)
    assert r.returncode == 1 and not out.exists(), r.stderr
    assert "no index beside it" in own(r.stderr) and "retry" not in own(r.stderr), r.stderr


@pytest.mark.parametrize("start", ["probe", "interval"])
def test_a_remote_input_that_keeps_failing_exits_75(www, wd, start):
    """When every retry fails (a server down from the start, or from the middle of a fetch), the run stops
    with EX_TEMPFAIL, 75, and writes nothing: the workflow retries the sample later."""
    eof, n_data = eof_request(www), [0]

    def down(n, path, rng):
        if start == "probe":
            return True
        n_data[0] += rng not in (None, eof)
        return n_data[0] >= 20
    out = wd / "o.json"
    with Server(www, fail=down) as srv:
        r = count(srv.url + "/fx.bam", "fetch", out, "--retries", "2", threads=1, cwd=wd)
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


def test_a_signed_url_is_kept_out_of_the_output(local, www, wd):
    """A signed URL's query is a credential. It reaches the server (for the file and its index) but not
    the counts or the engine's messages; htslib's own error lines still print it, which the engine cannot stop."""
    q = "?X-Amz-Credential=KEYSECRET&X-Amz-Signature=deadbeefSECRET"
    with Server(www) as srv:
        for mode in ("scan", "fetch"):
            r = ok(count(srv.url + "/fx.bam" + q, mode, wd / f"{mode}.json", cwd=wd))
            text = (wd / f"{mode}.json").read_text()
            c = json.loads(text)
            assert "SECRET" not in text and "SECRET" not in r.stderr, r.stderr
            assert c["input"] == srv.url + "/fx.bam?<redacted>"
            assert strip(c, {"input", "elapsed_sec"}) == strip(local[mode], {"input", "elapsed_sec"})
        assert {p for p, _, _ in srv.log} == {"/fx.bam" + q, "/fx.bam.csi" + q}
        n = len(srv.log)
        r = count(srv.url + "/missing.bam" + q, "scan", wd / "m.json", cwd=wd)
        assert r.returncode == 1 and len(srv.log) == n + 1, r.stderr                 # a 404 is final: no retry, no 75
        assert "SECRET" not in own(r.stderr) and "missing.bam?<redacted>" in own(r.stderr), r.stderr
