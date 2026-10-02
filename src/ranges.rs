//! A remote fetch by exact byte ranges: over HTTP(S) it moves the bytes its index prices.
//!
//! htslib reads a URL with open-ended range requests (hfile_libcurl): a seek starts a request from the new
//! offset to the end of the file and abandons the one before, and what was in flight on it is thrown away.
//! A fetch seeks at every query, so a remote fetch moved several times what its index prices (HG02300's
//! fetch of the karyotype windows: 1,517 to 2,184 MB on the wire for 422 MB of containers). Here the engine reads the index
//! itself, works out which bytes each query will read - the CRAM containers htslib's iterator visits, the
//! BGZF blocks of a BAM query's chunks - and downloads exactly those, each once, with bounded range requests
//! over kept-alive connections, several at a time, into a spool file. htslib reads the input through an
//! hFILE backend of this module's (`ngsdose:` names), which serves those bytes at their offsets: htslib
//! decodes what it would have decoded from the server, byte for byte, so the counts are those of htslib's
//! own transport. A read the plan did not foresee is fetched when it happens, and counted ("outside the
//! plan"): a plan that misses costs bytes, never counts.

use crate::count::{backoff, redact, scrub_urls, set_errno, tick, Indexed, Permanent, TempFail};
use anyhow::{bail, Context, Result};
use rust_htslib::htslib;
use rustc_hash::FxHashMap;
use std::collections::BTreeMap;
use std::ffi::CStr;
use std::os::raw::{c_char, c_int, c_void};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering::Relaxed};
use std::sync::{Arc, Condvar, Mutex, MutexGuard, Weak};
use std::time::Duration;

/// The first request: a CRAM's file definition and the start of its header, a BAM's header, the format.
const HEAD: u64 = 64 << 10;
/// The last request: the end-of-file marker (38 bytes in a CRAM, 28 in a BAM).
const TAIL: u64 = 64;
/// What htslib reads of a container it does not enter (one past the query's end, or one wholly before
/// its start): the container's header, a few dozen bytes and one landmark per slice.
const PEEK: u64 = 1 << 10;
/// CRAM ranges closer than this are asked for in one request.
const JOIN: u64 = 16 << 10;
/// No request is larger: the data reach the readers in pieces, and a failed request costs one piece.
const MAX_REQUEST: u64 = 4 << 20;
/// A read outside the plan fetches at least this much.
const OUTSIDE: u64 = 64 << 10;
/// A BGZF block's header: gzip's ten bytes, XLEN, and the BC subfield that holds BSIZE (the block's length - 1).
const BGZF_HEADER: u64 = 18;
/// An index larger than this is not one.
const MAX_INDEX: usize = 1 << 30;

const EIO: i32 = 5;
const ENOENT: i32 = 2;
const EINVAL: i32 = 22;
const EROFS: i32 = 30;
const EBADF: i32 = 9;

// ------------------------------------------------------------------------------------------------ HTTP

/// Why a request failed, and whether asking again can help.
#[derive(Debug, Clone, PartialEq)]
enum Failure {
    /// a lost connection, a timeout, a 5xx: another attempt may pass
    Transient(String),
    /// an answer that will not change for this request: a 404, a refused range (the interval's retries may still ask again)
    Final(String),
    /// what no attempt can cure, which ends the run: a file that changed, a server without range requests, a
    /// spool that cannot be written
    Fatal(String),
}

impl Failure {
    fn msg(&self) -> &str {
        match self {
            Failure::Transient(m) | Failure::Final(m) | Failure::Fatal(m) => m,
        }
    }
}

/// A Content-Range: (first, last) or None for "*", and the total or None for "*".
type ContentRange = (Option<(u64, u64)>, Option<u64>);

/// One response: status, the headers this module reads, and the body.
#[derive(Default, Debug)]
struct Answer {
    status: u32,
    range: Option<ContentRange>,
    etag: Option<String>,
    body: Vec<u8>,
    /// more than this and the transfer is stopped: a server that ignores Range sends the whole file
    limit: usize,
    over: bool,
}

struct Collect {
    a: Answer,
    stop: Arc<AtomicBool>,
    since_tick: usize,
}

impl curl::easy::Handler for Collect {
    fn write(&mut self, data: &[u8]) -> Result<usize, curl::easy::WriteError> {
        if self.a.body.len() + data.len() > self.a.limit {
            self.a.over = true;
            return Ok(0); // stops the transfer
        }
        self.a.body.extend_from_slice(data);
        self.since_tick += data.len();
        if self.since_tick >= 1 << 18 {
            self.since_tick = 0;
            tick();
        }
        Ok(data.len())
    }
    fn header(&mut self, data: &[u8]) -> bool {
        let line = String::from_utf8_lossy(data);
        let line = line.trim_end();
        if line.starts_with("HTTP/") {
            // a response of its own (after a redirect, the next one): its headers replace the last one's
            self.a.status = line.split_whitespace().nth(1).and_then(|s| s.parse().ok()).unwrap_or(0);
            self.a.range = None;
            self.a.etag = None;
        } else if let Some((k, v)) = line.split_once(':') {
            if k.trim().eq_ignore_ascii_case("content-range") {
                self.a.range = content_range(v);
            } else if k.trim().eq_ignore_ascii_case("etag") {
                self.a.etag = Some(v.trim().to_string());
            }
        }
        true
    }
    fn progress(&mut self, _: f64, _: f64, _: f64, _: f64) -> bool {
        !self.stop.load(Relaxed)
    }
}

/// "bytes 0-99/1000" -> (Some((0, 99)), Some(1000)); "bytes */1000" -> (None, Some(1000)).
fn content_range(v: &str) -> Option<ContentRange> {
    let rest = v.trim().strip_prefix("bytes")?.trim_start();
    let (span, total) = rest.split_once('/')?;
    let total = match total.trim() {
        "*" => None,
        t => Some(t.parse().ok()?),
    };
    let span = match span.trim() {
        "*" => None,
        s => {
            let (a, b) = s.split_once('-')?;
            Some((a.trim().parse().ok()?, b.trim().parse().ok()?))
        }
    };
    Some((span, total))
}

/// A connection to the server: one handle, whose connection is kept alive from request to request.
struct Client {
    easy: curl::easy::Easy2<Collect>,
}

fn curl_failure(e: &curl::Error) -> Failure {
    let msg = scrub_urls(&match e.extra_description() {
        Some(x) if !x.is_empty() => format!("{} ({})", e.description(), x.trim()),
        _ => e.description().to_string(),
    });
    if e.is_url_malformed()
        || e.is_unsupported_protocol()
        || e.is_ssl_cacert()
        || e.is_ssl_cacert_badfile()
        || e.is_peer_failed_verification()
        || e.is_ssl_certproblem()
        || e.is_too_many_redirects()
    {
        Failure::Final(msg)
    } else {
        Failure::Transient(msg)
    }
}

impl Client {
    fn new(stop: Arc<AtomicBool>) -> Result<Client, Failure> {
        let mut e = curl::easy::Easy2::new(Collect { a: Answer::default(), stop, since_tick: 0 });
        let set = |e: &mut curl::easy::Easy2<Collect>| -> Result<(), curl::Error> {
            e.get(true)?;
            e.follow_location(true)?;
            e.max_redirections(10)?;
            e.connect_timeout(Duration::from_secs(30))?;
            // a connection that moves less than 1 kB/s for a minute is dead: fail it, and ask again
            e.low_speed_limit(1024)?;
            e.low_speed_time(Duration::from_secs(60))?;
            e.useragent(concat!("ngs-dose/", env!("CARGO_PKG_VERSION")))?;
            // the CA store, as htslib takes it (main points CURL_CA_BUNDLE at the system's)
            if let Some(ca) = std::env::var_os("CURL_CA_BUNDLE") {
                e.cainfo(ca)?;
            }
            e.tcp_keepalive(true)?;
            e.progress(true)?;
            Ok(())
        };
        set(&mut e).map_err(|x| curl_failure(&x))?;
        Ok(Client { easy: e })
    }

    /// GET `url`, bytes [a, b) or the whole resource; a body larger than `limit` is cut off (`over`).
    fn get(&mut self, url: &str, range: Option<(u64, u64)>, limit: usize) -> Result<Answer, Failure> {
        // a range's body has room for all of it from the start
        let body = Vec::with_capacity(if range.is_some() { limit } else { 0 });
        self.easy.get_mut().a = Answer { limit, body, ..Default::default() };
        let set = |e: &mut curl::easy::Easy2<Collect>| -> Result<(), curl::Error> {
            e.url(url)?;
            if let Some((a, b)) = range {
                e.range(&format!("{}-{}", a, b - 1))?;
            }
            Ok(())
        };
        set(&mut self.easy).map_err(|x| curl_failure(&x))?;
        let res = self.easy.perform();
        let a = std::mem::take(&mut self.easy.get_mut().a);
        match res {
            Ok(()) => Ok(a),
            Err(_) if a.over => Ok(a),
            Err(e) => Err(curl_failure(&e)),
        }
    }
}

/// A server that answers a request for bytes [a, b) with the whole file: nothing a fetch can use.
fn no_ranges(shown: &str, a: u64, b: u64) -> Failure {
    Failure::Fatal(format!(
        "{} does not answer range requests (HTTP 200 to a request for bytes {}-{}): a fetch needs them. Scan the file instead, or copy it (or a \
         samtools cut along `ngs-dose plan`) and fetch from the copy",
        shown,
        a,
        b - 1
    ))
}

/// What a status that is not a success says.
fn status_failure(status: u32, what: &str) -> Failure {
    match status {
        408 | 425 | 429 | 500..=599 => Failure::Transient(format!("HTTP {} for {}", status, what)),
        0 => Failure::Transient(format!("no answer for {}", what)),
        _ => Failure::Final(format!("HTTP {} for {}", status, what)),
    }
}

// ------------------------------------------------------------------------------------------------ the index

/// One line of a .crai: a slice, or a reference's run within a slice of several references.
#[derive(Debug, Clone, Copy, PartialEq)]
struct CraiLine {
    refid: i32,
    /// 1-based start and inclusive end of its alignment span, as htslib reads them (start + span - 1)
    start: i64,
    end: i64,
    container: u64,
}

/// A .crai (gzip-compressed or not): every line six integers.
fn parse_crai(raw: &[u8]) -> Result<Vec<CraiLine>> {
    let text = if raw.starts_with(&[0x1f, 0x8b]) {
        let mut s = Vec::new();
        std::io::Read::read_to_end(&mut flate2::read::MultiGzDecoder::new(raw), &mut s).context("the CRAM index does not decompress")?;
        s
    } else {
        raw.to_vec()
    };
    let text = std::str::from_utf8(&text).context("the CRAM index is not text")?;
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate().filter(|(_, l)| !l.trim().is_empty()) {
        let f: Vec<i64> = line
            .split_ascii_whitespace()
            .map(|x| x.parse::<i64>())
            .collect::<Result<_, _>>()
            .ok()
            .filter(|v: &Vec<i64>| v.len() == 6)
            .with_context(|| format!("line {} of the CRAM index is not six integers", i + 1))?;
        if f[0] < -1 || f[3] < 0 || f[0] > i32::MAX as i64 {
            bail!("line {} of the CRAM index has reference {} and offset {}", i + 1, f[0], f[3]);
        }
        out.push(CraiLine { refid: f[0] as i32, start: f[1], end: f[1] + f[2] - 1, container: f[3] as u64 });
    }
    Ok(out)
}

#[derive(Debug, Clone)]
struct Container {
    off: u64,
    /// the next container's offset; the last one's, the end of the file
    end: u64,
    /// the references its index lines name (several for a container of several references)
    refs: Vec<i32>,
    /// its alignment span over its lines (1-based, inclusive), for a container of one reference
    start: i64,
    last: i64,
}

/// An entry of the list htslib searches for a reference: the entries no earlier entry of the reference
/// contains (cram_index_load nests the others under the entry that contains them).
#[derive(Debug, Clone, Copy)]
struct Top {
    start: i64,
    end: i64,
    c: usize,
}

/// A CRAM index as htslib loads it, with the containers in file order.
#[derive(Debug, Default)]
pub struct Crai {
    containers: Vec<Container>,
    tops: FxHashMap<i32, Vec<Top>>,
}

impl Crai {
    fn new(lines: &[CraiLine], size: u64) -> Result<Crai> {
        let mut offs: Vec<u64> = lines.iter().map(|l| l.container).collect();
        offs.sort_unstable();
        offs.dedup();
        if let Some(&o) = offs.last() {
            if o >= size {
                bail!("the CRAM index names a container at byte {}, but the file has {} bytes: it is the index of another file", o, size);
            }
        }
        let at: FxHashMap<u64, usize> = offs.iter().enumerate().map(|(i, &o)| (o, i)).collect();
        let mut containers: Vec<Container> = offs
            .iter()
            .enumerate()
            .map(|(i, &o)| Container { off: o, end: offs.get(i + 1).copied().unwrap_or(size), refs: Vec::new(), start: i64::MAX, last: i64::MIN })
            .collect();
        for l in lines {
            let c = &mut containers[at[&l.container]];
            if !c.refs.contains(&l.refid) {
                c.refs.push(l.refid);
            }
            c.start = c.start.min(l.start);
            c.last = c.last.max(l.end);
        }
        // cram_index_load: lines in file order; each reference's list starts afresh when its lines begin, and an
        // entry is nested under the last entry that contains it (an unmapped entry under none)
        let (root_start, root_end) = (i32::MIN as i64, i32::MAX as i64);
        let mut tops: FxHashMap<i32, Vec<Top>> = FxHashMap::default();
        let mut cur: Option<i32> = None;
        let mut stack: Vec<(i64, i64)> = Vec::new();
        for l in lines {
            if cur != Some(l.refid) {
                cur = Some(l.refid);
                tops.insert(l.refid, Vec::new());
                stack.clear();
                stack.push((root_start, root_end));
            }
            while stack.len() > 1 && {
                let t = stack[stack.len() - 1];
                !(l.start >= t.0 && l.end <= t.1) || (t.0 == 0 && l.refid == -1)
            } {
                stack.pop();
            }
            if stack.len() == 1 {
                tops.get_mut(&l.refid).expect("inserted above").push(Top { start: l.start, end: l.end, c: at[&l.container] });
            }
            stack.push((l.start, l.end));
        }
        Ok(Crai { containers, tops })
    }

    /// The container htslib starts a query at (cram_index_query): `pos` is the query's 1-based start (0
    /// with refid -1, the reads without a coordinate). None: htslib reads nothing.
    fn start_of(&self, refid: i32, pos: i64) -> Option<usize> {
        let top = self.tops.get(&refid).filter(|t| !t.is_empty())?;
        let n = top.len() as i64;
        let (mut i, mut j) = (0i64, n - 1);
        let mut k = j / 2;
        while k != i {
            if top[k as usize].start >= pos {
                j = k;
            } else {
                i = k;
            }
            k = (j - i) / 2 + i;
        }
        if j >= 0 && top[j as usize].start < pos {
            i = j;
        }
        while i > 0 && top[(i - 1) as usize].end >= pos {
            i -= 1;
        }
        while i + 1 < n && top[i as usize].end < pos {
            i += 1;
        }
        Some(top[i as usize].c)
    }

    fn whole(&self, c: usize) -> Need {
        let ct = &self.containers[c];
        Need { start: ct.off, end: ct.end, extend: None }
    }

    fn peek(&self, c: usize) -> Need {
        let ct = &self.containers[c];
        Need { start: ct.off, end: ct.end.min(ct.off + PEEK), extend: None }
    }

    /// The bytes one query reads, as htslib's CRAM iterator (cram_next_slice and cram_get_seq) walks the
    /// containers from its start: a container of the query's reference that overlaps [beg1, end1] whole; one
    /// wholly before it, or the first past it or of another reference, its header; a container of several
    /// references whole, and the walk ends in it when it holds a later reference. `rref` -1: the reads without
    /// a coordinate, every container from the first of them to the end. A single-reference container in range
    /// may end the walk inside (a read past the end): the next container's header is read for nothing.
    fn walk(&self, rref: i32, beg1: i64, end1: i64, out: &mut Vec<Need>) {
        let Some(mut c) = self.start_of(rref, beg1) else {
            return;
        };
        while let Some(ct) = self.containers.get(c) {
            if let [r] = ct.refs[..] {
                if r != rref || (rref != -1 && ct.start > end1) {
                    out.push(self.peek(c));
                    return;
                }
                if rref != -1 && ct.last < beg1 {
                    out.push(self.peek(c));
                    c += 1;
                    continue;
                }
                out.push(self.whole(c));
            } else {
                out.push(self.whole(c));
                if rref != -1 && ct.refs.iter().any(|&r| r > rref || r == -1) {
                    return;
                }
            }
            c += 1;
        }
    }

    /// The first container's offset: where the header ends.
    fn first(&self) -> Option<u64> {
        self.containers.first().map(|c| c.off)
    }
}

// ------------------------------------------------------------------------------------------------ the plan

/// Bytes a query reads: [start, end), and with `extend` Some(t) the BGZF block at t, whose header these are,
/// must be read whole (its length is in the header).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Need {
    start: u64,
    end: u64,
    extend: Option<u64>,
}

/// A BAM query's chunk [u, v) in virtual offsets: the blocks it covers whole, and the block v points into,
/// whose length only its header tells.
fn bam_needs(u: u64, v: u64, out: &mut Vec<Need>) {
    let (cu, cv) = (u >> 16, v >> 16);
    if cv > cu {
        out.push(Need { start: cu, end: cv, extend: None });
    }
    if v & 0xffff != 0 {
        out.push(Need { start: cv, end: cv + BGZF_HEADER, extend: Some(cv) });
    }
}

/// The requests for a set of needs: ranges merged where they overlap or lie less than `join` apart, cut at
/// MAX_REQUEST bytes. A BGZF block to be read whole that a range of whole blocks already covers is dropped
/// (a range of whole blocks ends on a block boundary, so one that starts at or before the block and ends past
/// its start holds it); the others are requests of their own.
fn requests(mut needs: Vec<Need>, join: u64) -> Vec<Need> {
    needs.sort_unstable();
    let mut ranges: Vec<Need> = Vec::new();
    for n in needs.iter().filter(|n| n.extend.is_none() && n.end > n.start) {
        match ranges.last_mut() {
            Some(r) if n.start <= r.end + join => r.end = r.end.max(n.end),
            _ => ranges.push(*n),
        }
    }
    let covered = |t: u64| {
        let i = ranges.partition_point(|r| r.start <= t);
        i > 0 && ranges[i - 1].end > t
    };
    let mut tails: Vec<Need> = needs.iter().filter(|n| n.extend.is_some_and(|t| !covered(t))).copied().collect();
    tails.dedup();
    let mut out = Vec::new();
    for r in ranges {
        let mut a = r.start;
        while a < r.end {
            let b = r.end.min(a + MAX_REQUEST);
            out.push(Need { start: a, end: b, extend: None });
            a = b;
        }
    }
    out.extend(tails);
    out.sort_unstable();
    out
}

/// The BSIZE of a BGZF block header (the block's length - 1), checked: gzip, deflate, FEXTRA, the BC subfield.
fn bgzf_bsize(h: &[u8]) -> Option<u64> {
    if h.len() < BGZF_HEADER as usize || h[0] != 0x1f || h[1] != 0x8b || h[2] != 8 || h[3] & 4 == 0 {
        return None;
    }
    let xlen = u16::from_le_bytes([h[10], h[11]]);
    (xlen >= 6 && h[12] == b'B' && h[13] == b'C' && h[14] == 2 && h[15] == 0).then(|| u16::from_le_bytes([h[16], h[17]]) as u64)
}

// ------------------------------------------------------------------------------------------------ the store

#[derive(Debug, Clone, Copy, PartialEq)]
enum State {
    Queued,
    Running,
    Done,
    Failed,
}

#[derive(Debug)]
struct Job {
    end: u64,
    extend: Option<u64>,
    state: State,
    /// where its bytes are in the spool file, once done
    at: u64,
}

#[derive(Default)]
struct Book {
    /// by start; jobs do not overlap
    jobs: BTreeMap<u64, Job>,
    /// the planned jobs' starts, in the order the downloaders take them
    order: Vec<u64>,
    cursor: usize,
    spool_end: u64,
    /// the last failure, for messages
    failure: Option<String>,
    /// a failure that ends the run (the spool cannot be written, the file changed)
    fatal: Option<String>,
}

#[derive(Default)]
struct Stats {
    bytes: AtomicU64,
    requests: AtomicU64,
    outside_bytes: AtomicU64,
    outside_requests: AtomicU64,
    planned_bytes: AtomicU64,
    planned_requests: AtomicU64,
    retries: AtomicU64,
}

/// What every reader and downloader of one input shares: the file's identity, the jobs, the spool.
struct Shared {
    /// the URL as curl gets it (signed URLs keep their query)
    url: String,
    /// the URL for messages
    shown: String,
    size: u64,
    etag: Option<String>,
    retries: usize,
    spool: std::fs::File,
    spool_dir: PathBuf,
    book: Mutex<Book>,
    cv: Condvar,
    stop: Arc<AtomicBool>,
    clients: Mutex<Vec<Client>>,
    stats: Stats,
}

impl Shared {
    fn lock(&self) -> MutexGuard<'_, Book> {
        self.book.lock().unwrap_or_else(|p| p.into_inner())
    }

    fn client(&self) -> Result<Client, Failure> {
        let c = self.clients.lock().unwrap_or_else(|p| p.into_inner()).pop();
        match c {
            Some(c) => Ok(c),
            None => Client::new(self.stop.clone()),
        }
    }

    fn give_back(&self, c: Client) {
        self.clients.lock().unwrap_or_else(|p| p.into_inner()).push(c);
    }

    /// Bytes [a, b) of the file, checked against what the first answer said of it.
    fn get(&self, client: &mut Client, a: u64, b: u64) -> Result<Vec<u8>, Failure> {
        let what = format!("bytes {}-{} of {}", a, b - 1, self.shown);
        let ans = client.get(&self.url, Some((a, b)), (b - a) as usize)?;
        self.stats.requests.fetch_add(1, Relaxed);
        self.stats.bytes.fetch_add(ans.body.len() as u64, Relaxed);
        match ans.status {
            206 => match ans.range {
                Some((_, Some(t))) if t != self.size => {
                    return Err(Failure::Fatal(format!(
                        "{}: the server now gives the file as {} bytes, not {}: it changed during the fetch",
                        what, t, self.size
                    )))
                }
                Some((Some((x, y)), _)) if x == a && y + 1 == b => {}
                None => return Err(Failure::Fatal(format!("{}: the server answered without saying which bytes it sent (no Content-Range)", what))),
                Some(r) => return Err(Failure::Transient(format!("{}: the server answered with another range ({:?})", what, r))),
            },
            200 if a == 0 && b == self.size => {}
            200 => return Err(no_ranges(&self.shown, a, b)),
            416 => return Err(Failure::Final(format!("{}: the server says these bytes are past the end of the file (HTTP 416)", what))),
            s => return Err(status_failure(s, &what)),
        }
        if ans.over || ans.body.len() as u64 != b - a {
            return Err(Failure::Transient(format!("{}: {} bytes came", what, ans.body.len())));
        }
        if let (Some(want), Some(got)) = (&self.etag, &ans.etag) {
            if want != got {
                return Err(Failure::Fatal(format!("{}: the file changed during the fetch (ETag {} at the start, {} now)", what, want, got)));
            }
        }
        Ok(ans.body)
    }

    /// Run the job that starts at `s` (Running, claimed by the caller), with up to `attempts` attempts.
    fn run(&self, s: u64, attempts: usize) {
        let (mut end, extend) = {
            let b = self.lock();
            let j = &b.jobs[&s];
            (j.end, j.extend)
        };
        let mut client = match self.client() {
            Ok(c) => c,
            Err(f) => return self.settle(s, Err(f)),
        };
        let mut tries = 0;
        let res = loop {
            tries += 1;
            let got = self.get(&mut client, s, end).and_then(|mut data| {
                // a BGZF block whose header ends the range: read the rest of the block too, up to the next job
                if let Some(t) = extend {
                    let bsize = bgzf_bsize(&data[(t - s) as usize..]).ok_or_else(|| {
                        Failure::Final(format!(
                            "bytes {}-{} of {} are not a BGZF block header: the index is not this file's?",
                            t,
                            t + BGZF_HEADER - 1,
                            self.shown
                        ))
                    })?;
                    let mut b = self.lock();
                    let next = b.jobs.range(end..).next().map_or(self.size, |(&n, _)| n);
                    let stop = (t + bsize + 1).min(next).min(self.size);
                    if stop > end {
                        let j = b.jobs.get_mut(&s).expect("the job is ours");
                        j.end = stop;
                        j.extend = None;
                        drop(b);
                        data.extend(self.get(&mut client, end, stop)?);
                        end = stop;
                    }
                }
                Ok(data)
            });
            match got {
                Ok(d) => break Ok(d),
                Err(Failure::Transient(m)) if tries < attempts && !self.stop.load(Relaxed) => {
                    self.stats.retries.fetch_add(1, Relaxed);
                    eprintln!("[count] retry {}/{} for {}", tries, attempts, m);
                    // the job may have grown to the end of its block: ask for all of it next time
                    end = self.lock().jobs[&s].end;
                    self.pause(tries);
                }
                Err(f) => break Err(f),
            }
        };
        self.give_back(client);
        self.settle(s, res);
    }

    /// Wait before attempt `attempt + 1` as `backoff` does (1, 2, 4, 8, 16, then 16 s), but no longer than the run
    /// lasts: a stop ends the wait.
    fn pause(&self, attempt: usize) {
        let until = std::time::Instant::now() + Duration::from_millis(500 << attempt.min(5));
        while std::time::Instant::now() < until && !self.stop.load(Relaxed) {
            std::thread::sleep(Duration::from_millis(100));
            tick();
        }
    }

    /// Record a job's outcome: its bytes in the spool, or its failure; wake the readers.
    fn settle(&self, s: u64, res: Result<Vec<u8>, Failure>) {
        let res = res.and_then(|data| {
            let at = {
                let mut b = self.lock();
                let at = b.spool_end;
                b.spool_end += data.len() as u64;
                at
            };
            self.spool.write_all_at(&data, at).map_err(|e| {
                Failure::Fatal(format!(
                    "cannot keep the downloaded bytes in a file in {} ({}): give --spool-dir a directory with room",
                    self.spool_dir.display(),
                    e
                ))
            })?;
            Ok((at, data.len() as u64))
        });
        let mut b = self.lock();
        match res {
            Ok((at, len)) => {
                let j = b.jobs.get_mut(&s).expect("the job is ours");
                j.state = State::Done;
                j.at = at;
                j.end = s + len;
            }
            Err(f) => {
                if let Failure::Fatal(m) = &f {
                    b.fatal.get_or_insert_with(|| m.clone());
                }
                b.failure = Some(f.msg().to_string());
                b.jobs.get_mut(&s).expect("the job is ours").state = State::Failed;
            }
        }
        drop(b);
        self.cv.notify_all();
    }

    /// Up to `buf.len()` bytes at `pos`: from the spool, after waiting for or running the job that holds them,
    /// or fetched now if no job does. Ok(0) only at the end of the file.
    fn read_at(&self, pos: u64, buf: &mut [u8]) -> std::io::Result<usize> {
        if pos >= self.size || buf.is_empty() {
            return Ok(0);
        }
        let mut b = self.lock();
        loop {
            if let Some(f) = &b.fatal {
                b.failure = Some(f.clone());
                return Err(std::io::Error::from_raw_os_error(EIO));
            }
            let found = b.jobs.range(..=pos).next_back().filter(|(_, j)| j.end > pos).map(|(&s, j)| (s, j.state, j.end, j.at));
            match found {
                Some((s, State::Done, end, at)) => {
                    drop(b);
                    let n = (buf.len() as u64).min(end - pos) as usize;
                    self.spool.read_exact_at(&mut buf[..n], at + (pos - s))?;
                    return Ok(n);
                }
                Some((_, State::Running, _, _)) => {
                    if self.stop.load(Relaxed) {
                        return Err(std::io::Error::from_raw_os_error(EIO));
                    }
                    b = self.cv.wait_timeout(b, Duration::from_secs(1)).unwrap_or_else(|p| p.into_inner()).0;
                }
                Some((s, State::Queued, _, _)) => {
                    // a planned job no downloader has reached: the reader runs it now, once (the interval's retries retry it)
                    b.jobs.get_mut(&s).expect("found above").state = State::Running;
                    drop(b);
                    self.run(s, 1);
                    b = self.lock();
                }
                Some((s, State::Failed, _, _)) => {
                    // the read fails; the next attempt (the interval's retry) asks again
                    b.jobs.get_mut(&s).expect("found above").state = State::Queued;
                    return Err(std::io::Error::from_raw_os_error(EIO));
                }
                None => {
                    // a read the plan did not foresee: fetch it now, up to the next job
                    let next = b.jobs.range(pos..).next().map_or(self.size, |(&n, _)| n);
                    let end = next.min(pos + (buf.len() as u64).max(OUTSIDE)).min(self.size);
                    b.jobs.insert(pos, Job { end, extend: None, state: State::Running, at: 0 });
                    drop(b);
                    self.stats.outside_requests.fetch_add(1, Relaxed);
                    self.stats.outside_bytes.fetch_add(end - pos, Relaxed);
                    self.run(pos, 1);
                    b = self.lock();
                }
            }
        }
    }

    /// Add jobs for the parts of [a, b) no job covers yet (an `extend` job only where none covers any of it).
    fn add(b: &mut Book, n: Need) -> u64 {
        let mut added = 0;
        if n.extend.is_some() {
            let free = b.jobs.range(..n.end).next_back().is_none_or(|(_, j)| j.end <= n.start);
            if free {
                // a queued job that ends where the block starts asks for its header too: one request, not two
                let before =
                    b.jobs.range_mut(..n.start).next_back().filter(|(_, j)| j.end == n.start && j.state == State::Queued && j.extend.is_none());
                match before {
                    Some((_, j)) => {
                        j.end = n.end;
                        j.extend = n.extend;
                    }
                    None => {
                        b.jobs.insert(n.start, Job { end: n.end, extend: n.extend, state: State::Queued, at: 0 });
                        b.order.push(n.start);
                    }
                }
                added += n.end - n.start;
            }
            return added;
        }
        let mut a = n.start;
        if let Some((_, j)) = b.jobs.range(..=a).next_back() {
            a = a.max(j.end);
        }
        while a < n.end {
            let next = b.jobs.range(a..).next().map(|(&s, j)| (s, j.end));
            let stop = next.map_or(n.end, |(s, _)| s.min(n.end));
            if stop > a {
                b.jobs.insert(a, Job { end: stop, extend: None, state: State::Queued, at: 0 });
                b.order.push(a);
                added += stop - a;
            }
            match next {
                Some((s, e)) if s < n.end => a = e.max(stop),
                _ => break,
            }
        }
        added
    }

    /// A downloader: planned jobs in order, each with every attempt it is allowed.
    fn download(self: Arc<Self>) {
        loop {
            let s = {
                let mut b = self.lock();
                let mut found = None;
                while b.cursor < b.order.len() && !self.stop.load(Relaxed) {
                    let s = b.order[b.cursor];
                    b.cursor += 1;
                    if let Some(j) = b.jobs.get_mut(&s).filter(|j| j.state == State::Queued) {
                        j.state = State::Running;
                        found = Some(s);
                        break;
                    }
                }
                match found {
                    Some(s) => s,
                    None => return,
                }
            };
            self.run(s, self.retries);
        }
    }
}

fn spool_file(dir: &Path) -> Result<std::fs::File> {
    let nanos = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.subsec_nanos()).unwrap_or(0);
    let name = dir.join(format!(".ngs-dose.{}.{}.spool", std::process::id(), nanos));
    let f = std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create_new(true)
        .open(&name)
        .with_context(|| format!("cannot create a file in {} for the downloaded bytes (--spool-dir, or TMPDIR)", dir.display()))?;
    // the file lives while it is open, and goes when the run ends, however it ends
    let _ = std::fs::remove_file(&name);
    Ok(f)
}

// ------------------------------------------------------------------------------------------------ htslib's view

/// htslib's hFILE backend vector (hfile_internal.h): read, write, seek, flush, close.
#[repr(C)]
struct Backend {
    read: unsafe extern "C" fn(*mut htslib::hFILE, *mut c_void, usize) -> isize,
    write: unsafe extern "C" fn(*mut htslib::hFILE, *const c_void, usize) -> isize,
    seek: unsafe extern "C" fn(*mut htslib::hFILE, htslib::off_t, c_int) -> htslib::off_t,
    flush: unsafe extern "C" fn(*mut htslib::hFILE) -> c_int,
    close: unsafe extern "C" fn(*mut htslib::hFILE) -> c_int,
}

/// A URL scheme's handler (hfile_internal.h): open, isremote, provider, priority (below 2000: no vopen).
#[repr(C)]
struct SchemeHandler {
    open: unsafe extern "C" fn(*const c_char, *const c_char) -> *mut htslib::hFILE,
    isremote: unsafe extern "C" fn(*const c_char) -> c_int,
    provider: *const c_char,
    priority: c_int,
}

struct SyncHandler(SchemeHandler);
// SAFETY: the handler is never written; `provider` points at a static C string
unsafe impl Sync for SyncHandler {}

extern "C" {
    fn hfile_init(struct_size: usize, mode: *const c_char, capacity: usize) -> *mut htslib::hFILE;
    fn hfile_add_scheme_handler(scheme: *const c_char, handler: *const SchemeHandler);
}

static BACKEND: Backend = Backend { read: hf_read, write: hf_write, seek: hf_seek, flush: hf_flush, close: hf_close };
static HANDLER: SyncHandler = SyncHandler(SchemeHandler { open: hf_open, isremote: hf_isremote, provider: c"ngs-dose".as_ptr(), priority: 50 });
/// htslib keeps the scheme's name by pointer: it must live for ever
const SCHEME: &CStr = c"ngsdose";

/// An open hFILE: htslib's part first, as a backend's struct must begin.
#[repr(C)]
struct HFile {
    base: htslib::hFILE,
    handle: *mut Handle,
}

enum Source {
    Remote(Arc<Shared>),
    Bytes(Arc<Vec<u8>>),
}

struct Handle {
    src: Source,
    pos: u64,
}

impl Handle {
    fn size(&self) -> u64 {
        match &self.src {
            Source::Remote(s) => s.size,
            Source::Bytes(b) => b.len() as u64,
        }
    }
}

/// An input htslib can open: `ngsdose:ID` (the file) and `ngsdose:ID.idx` (its index).
struct Served {
    id: u64,
    file: Weak<Shared>,
    index: Arc<Vec<u8>>,
}

static REGISTRY: Mutex<Vec<Served>> = Mutex::new(Vec::new());
static NEXT_ID: AtomicU64 = AtomicU64::new(1);

fn lookup(name: &str) -> Option<Source> {
    let rest = name.strip_prefix("ngsdose:")?;
    let (id, index) = match rest.strip_suffix(".idx") {
        Some(id) => (id, true),
        None => (rest, false),
    };
    let id: u64 = id.parse().ok()?;
    let reg = REGISTRY.lock().unwrap_or_else(|p| p.into_inner());
    let s = reg.iter().find(|s| s.id == id)?;
    if index {
        Some(Source::Bytes(s.index.clone()))
    } else {
        s.file.upgrade().map(Source::Remote)
    }
}

/// Register the `ngsdose:` scheme with htslib, once, before any thread opens a file: htslib's table of schemes
/// is not locked against a registration.
pub fn register() {
    static ONCE: std::sync::Once = std::sync::Once::new();
    ONCE.call_once(|| unsafe {
        // the table exists once htslib has loaded its plugins
        htslib::hfile_has_plugin(c"libcurl".as_ptr());
        hfile_add_scheme_handler(SCHEME.as_ptr(), &HANDLER.0);
    });
}

unsafe fn handle<'a>(fp: *mut htslib::hFILE) -> &'a mut Handle {
    &mut *(*(fp as *mut HFile)).handle
}

unsafe extern "C" fn hf_open(name: *const c_char, mode: *const c_char) -> *mut htslib::hFILE {
    let r = std::panic::catch_unwind(|| {
        let name = CStr::from_ptr(name).to_str().map_err(|_| ENOENT)?;
        if CStr::from_ptr(mode).to_bytes().iter().any(|c| matches!(c, b'w' | b'a' | b'+')) {
            return Err(EROFS);
        }
        let src = lookup(name).ok_or(ENOENT)?;
        let fp = hfile_init(std::mem::size_of::<HFile>(), mode, 0);
        if fp.is_null() {
            return Err(EIO);
        }
        (*fp).backend = &BACKEND as *const Backend as *const htslib::hFILE_backend;
        (*(fp as *mut HFile)).handle = Box::into_raw(Box::new(Handle { src, pos: 0 }));
        Ok(fp)
    });
    match r {
        Ok(Ok(fp)) => fp,
        Ok(Err(e)) => {
            set_errno(e);
            std::ptr::null_mut()
        }
        Err(_) => {
            set_errno(EIO);
            std::ptr::null_mut()
        }
    }
}

unsafe extern "C" fn hf_isremote(_: *const c_char) -> c_int {
    // what it serves is local to htslib: no index is downloaded or cached for it
    0
}

unsafe extern "C" fn hf_read(fp: *mut htslib::hFILE, buf: *mut c_void, n: usize) -> isize {
    let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let h = handle(fp);
        let out = std::slice::from_raw_parts_mut(buf as *mut u8, n);
        let got = match &h.src {
            Source::Remote(s) => s.read_at(h.pos, out).map_err(|e| e.raw_os_error().unwrap_or(EIO))?,
            Source::Bytes(b) => {
                let p = (h.pos as usize).min(b.len());
                let k = n.min(b.len() - p);
                out[..k].copy_from_slice(&b[p..p + k]);
                k
            }
        };
        h.pos += got as u64;
        Ok::<isize, i32>(got as isize)
    }));
    match r {
        Ok(Ok(k)) => k,
        Ok(Err(e)) => {
            set_errno(e);
            -1
        }
        Err(_) => {
            set_errno(EIO);
            -1
        }
    }
}

unsafe extern "C" fn hf_write(_: *mut htslib::hFILE, _: *const c_void, _: usize) -> isize {
    set_errno(EBADF);
    -1
}

unsafe extern "C" fn hf_seek(fp: *mut htslib::hFILE, offset: htslib::off_t, whence: c_int) -> htslib::off_t {
    let h = handle(fp);
    let base = match whence {
        0 => 0,
        1 => h.pos as i64,
        2 => h.size() as i64,
        _ => {
            set_errno(EINVAL);
            return -1;
        }
    };
    match base.checked_add(offset) {
        Some(p) if p >= 0 => {
            h.pos = p as u64;
            p
        }
        _ => {
            set_errno(EINVAL);
            -1
        }
    }
}

unsafe extern "C" fn hf_flush(_: *mut htslib::hFILE) -> c_int {
    0
}

unsafe extern "C" fn hf_close(fp: *mut htslib::hFILE) -> c_int {
    let p = (*(fp as *mut HFile)).handle;
    if !p.is_null() {
        drop(Box::from_raw(p));
        (*(fp as *mut HFile)).handle = std::ptr::null_mut();
    }
    0
}

// ------------------------------------------------------------------------------------------------ the input

/// Whether a URL is read by exact byte ranges: http and https.
pub fn supports(url: &str) -> bool {
    let scheme = url.split("://").next().unwrap_or("").to_ascii_lowercase();
    url.contains("://") && (scheme == "http" || scheme == "https")
}

/// A remote BAM or CRAM read by exact byte ranges: htslib opens it as `htslib_url()`.
pub struct Remote {
    id: u64,
    shared: Arc<Shared>,
    crai: Option<Crai>,
    index_len: u64,
    index_downloaded: bool,
    downloaders: Mutex<Vec<std::thread::JoinHandle<()>>>,
    started: std::time::Instant,
}

/// Run `f` up to `retries` times while it fails transiently, backing off in between; a final failure stops
/// at once. Messages say `what`.
fn with_retries<T>(retries: usize, what: &str, mut f: impl FnMut() -> Result<T, Failure>) -> Result<T> {
    let tries = retries.max(1);
    let mut attempt = 0;
    loop {
        attempt += 1;
        match f() {
            Ok(v) => return Ok(v),
            Err(Failure::Final(m) | Failure::Fatal(m)) => return Err(anyhow::Error::new(Permanent(m))),
            Err(Failure::Transient(m)) if attempt >= tries => return Err(anyhow::anyhow!("{}: {}", what, m).context(TempFail)),
            Err(Failure::Transient(m)) => {
                eprintln!("[count] retry {}/{} to {}: {}", attempt, tries, what, m);
                backoff(attempt);
            }
        }
    }
}

impl Remote {
    /// Open `url` for a fetch: its first bytes (the format, the size, the header), its index (`index`, a path or
    /// URL, or beside the file as htslib would look for it), and its end-of-file marker.
    pub fn open(url: &str, index: Option<&str>, retries: usize, spool_dir: Option<&Path>) -> Result<Arc<Remote>> {
        let shown = redact(url);
        let curl_url = url::Url::parse(url).map_err(|e| anyhow::Error::new(Permanent(format!("bad URL {}: {}", shown, e))))?.to_string();
        let stop = Arc::new(AtomicBool::new(false));
        // the first bytes: whether the server answers ranges, how large the file is, what it is
        let (head, size, etag) = with_retries(retries, "open the input", || {
            let mut c = Client::new(stop.clone())?;
            let a = c.get(&curl_url, Some((0, HEAD)), HEAD as usize)?;
            match a.status {
                206 => match a.range {
                    Some((Some((0, last)), Some(total))) if a.body.len() as u64 == last + 1 && !a.over => Ok((a.body, total, a.etag)),
                    None | Some((_, None)) => Err(Failure::Final(format!(
                        "the first bytes of {}: the server answered without saying which bytes of how many it sent (Content-Range)",
                        shown
                    ))),
                    r => Err(Failure::Transient(format!("the first bytes of {}: the server answered with another range ({:?})", shown, r))),
                },
                // a file no longer than the request comes whole
                200 if !a.over => {
                    let n = a.body.len() as u64;
                    Ok((a.body, n, a.etag))
                }
                200 => Err(no_ranges(&shown, 0, HEAD)),
                s => Err(status_failure(s, &shown)),
            }
        })?;
        let cram = head.starts_with(b"CRAM");
        if head.is_empty() {
            return Err(anyhow::Error::new(Permanent(format!("{} is empty", shown))));
        }
        if !cram && !head.starts_with(&[0x1f, 0x8b]) {
            return Err(anyhow::Error::new(Permanent(format!("{} is neither a CRAM nor a BAM (its first bytes)", shown))));
        }
        let dir = spool_dir.map(Path::to_path_buf).unwrap_or_else(std::env::temp_dir);
        let shared = Arc::new(Shared {
            url: curl_url.clone(),
            shown: shown.clone(),
            size,
            etag,
            retries,
            spool: spool_file(&dir).map_err(|e| anyhow::Error::new(Permanent(format!("{:#}", e))))?,
            spool_dir: dir,
            book: Mutex::new(Book::default()),
            cv: Condvar::new(),
            stop: stop.clone(),
            clients: Mutex::new(Vec::new()),
            stats: Stats::default(),
        });
        shared.stats.requests.fetch_add(1, Relaxed);
        shared.stats.bytes.fetch_add(head.len() as u64, Relaxed);
        {
            let mut b = shared.lock();
            Shared::add(&mut b, Need { start: 0, end: head.len() as u64, extend: None });
            b.order.clear();
            b.jobs.get_mut(&0).expect("added").state = State::Running;
        }
        shared.settle(0, Ok(head));
        // the index
        let (raw, downloaded) = load_index(url, index, cram, retries, &shown, &stop)?;
        let crai = if cram {
            let lines = parse_crai(&raw).map_err(|e| anyhow::Error::new(Permanent(format!("the index of {}: {:#}", shown, e))))?;
            Some(Crai::new(&lines, size).map_err(|e| anyhow::Error::new(Permanent(format!("the index of {}: {:#}", shown, e))))?)
        } else {
            None
        };
        let id = NEXT_ID.fetch_add(1, Relaxed);
        REGISTRY.lock().unwrap_or_else(|p| p.into_inner()).push(Served { id, file: Arc::downgrade(&shared), index: Arc::new(raw.clone()) });
        let r = Remote {
            id,
            shared,
            crai,
            index_len: raw.len() as u64,
            index_downloaded: downloaded,
            downloaders: Mutex::new(Vec::new()),
            started: std::time::Instant::now(),
        };
        // the rest of the header, and the end-of-file marker, which every open reads
        let mut first: Vec<Need> = vec![Need { start: size.saturating_sub(TAIL), end: size, extend: None }];
        if let Some(h) = r.crai.as_ref().and_then(Crai::first).filter(|&h| h > HEAD) {
            first.push(Need { start: HEAD, end: h, extend: None });
        }
        let jobs = {
            let mut b = r.shared.lock();
            let from = b.order.len();
            for n in requests(first, 0) {
                Shared::add(&mut b, n);
            }
            let jobs: Vec<u64> = b.order.drain(from..).collect();
            for s in &jobs {
                b.jobs.get_mut(s).expect("added").state = State::Running;
            }
            jobs
        };
        for s in jobs {
            r.shared.run(s, retries);
            let b = r.shared.lock();
            if b.jobs[&s].state == State::Failed {
                let why = b.failure.clone().unwrap_or_default();
                return Err(match &b.fatal {
                    Some(f) => anyhow::Error::new(Permanent(f.clone())),
                    None => anyhow::anyhow!("reading the header and the end of {}: {}", shown, why).context(TempFail),
                });
            }
        }
        Ok(Arc::new(r))
    }

    /// What htslib opens: the file and its index, both served by this module.
    pub fn htslib_url(&self) -> url::Url {
        url::Url::parse(&format!("ngsdose:{}##idx##ngsdose:{}.idx", self.id, self.id)).expect("a valid URL")
    }

    /// Plan the queries of a fetch - (reference id, 0-based start, end) each, and the reads without a coordinate
    /// with `unmapped` - and start `downloaders` connections on the bytes they read. A BAM is planned with
    /// htslib's own iterator, on a reader `open` gives.
    pub fn prefetch(&self, queries: &[(i32, i64, i64)], unmapped: bool, open: impl FnOnce() -> Result<Indexed>, downloaders: usize) -> Result<()> {
        let mut needs = Vec::new();
        let join = match &self.crai {
            Some(crai) => {
                for &(tid, beg, end) in queries {
                    crai.walk(tid, beg + 1, end, &mut needs);
                }
                if unmapped {
                    crai.walk(-1, 0, 0, &mut needs);
                }
                JOIN
            }
            None => {
                let rd = open()?;
                let idx = rd.index().inner_ptr();
                let mut chunks = |itr: *mut htslib::hts_itr_t| unsafe {
                    if itr.is_null() {
                        return;
                    }
                    let it = &*itr;
                    if it.read_rest() != 0 {
                        // the reads without a coordinate: from the first of them to the end
                        needs.push(Need { start: it.curr_off >> 16, end: self.shared.size, extend: None });
                    } else if !it.off.is_null() {
                        for k in 0..it.n_off.max(0) as usize {
                            let o = *it.off.add(k);
                            bam_needs(o.u, o.v, &mut needs);
                        }
                    }
                    htslib::hts_itr_destroy(itr);
                };
                for &(tid, beg, end) in queries {
                    chunks(unsafe { htslib::sam_itr_queryi(idx, tid, beg, end) });
                }
                if unmapped {
                    chunks(unsafe { htslib::sam_itr_queryi(idx, htslib::HTS_IDX_NOCOOR, 0, 0) });
                }
                0
            }
        };
        let (n, bytes) = {
            let mut b = self.shared.lock();
            let from = b.order.len();
            let mut bytes = 0;
            for n in requests(needs, join) {
                bytes += Shared::add(&mut b, n);
            }
            (b.order.len() - from, bytes)
        };
        self.shared.stats.planned_bytes.fetch_add(bytes, Relaxed);
        self.shared.stats.planned_requests.fetch_add(n as u64, Relaxed);
        eprintln!(
            "[count] {} queries read {:.1} MB of {} ({:.2}% of {:.2} GB) in {} byte ranges, by exact range requests",
            queries.len() + unmapped as usize,
            bytes as f64 / 1e6,
            self.shared.shown,
            100.0 * bytes as f64 / self.shared.size.max(1) as f64,
            self.shared.size as f64 / 1e9,
            n
        );
        let mut d = self.downloaders.lock().unwrap_or_else(|p| p.into_inner());
        for _ in 0..downloaders.max(1) {
            let s = self.shared.clone();
            d.push(std::thread::spawn(move || s.download()));
        }
        Ok(())
    }

    /// Stop the downloaders (a request in flight is abandoned) and wait for them.
    pub fn finish(&self) {
        self.shared.stop.store(true, Relaxed);
        self.shared.cv.notify_all();
        for h in self.downloaders.lock().unwrap_or_else(|p| p.into_inner()).drain(..) {
            let _ = h.join();
        }
    }

    /// An error of a read through this input, with what the downloads said: a failure no attempt can cure
    /// (the file changed, the spool cannot be written) is Permanent.
    pub fn explain(&self, e: anyhow::Error) -> anyhow::Error {
        let mut b = self.shared.lock();
        match (b.fatal.clone(), b.failure.take()) {
            (Some(f), _) => e.context(Permanent(f)),
            (None, Some(f)) => anyhow::anyhow!("{} ({})", e, f),
            _ => e,
        }
    }

    /// One line on what the fetch moved.
    pub fn summary(&self) -> String {
        let st = &self.shared.stats;
        let mb = |n: u64| n as f64 / 1e6;
        let moved = st.bytes.load(Relaxed) + if self.index_downloaded { self.index_len } else { 0 };
        format!(
            "[count] moved {:.1} MB in {} requests ({:.2}% of the {:.2}-GB file, {:.0} s): {:.1} MB planned in {} byte ranges, {:.1} MB in {} reads \
             outside the plan, the header and end-of-file marker{}; {} requests retried",
            mb(moved),
            st.requests.load(Relaxed) + self.index_downloaded as u64,
            100.0 * moved as f64 / self.shared.size.max(1) as f64,
            self.shared.size as f64 / 1e9,
            self.started.elapsed().as_secs_f64(),
            mb(st.planned_bytes.load(Relaxed)),
            st.planned_requests.load(Relaxed),
            mb(st.outside_bytes.load(Relaxed)),
            st.outside_requests.load(Relaxed),
            if self.index_downloaded { format!(", and the index ({:.1} MB)", mb(self.index_len)) } else { String::new() },
            st.retries.load(Relaxed)
        )
    }
}

impl Drop for Remote {
    fn drop(&mut self) {
        self.finish();
        REGISTRY.lock().unwrap_or_else(|p| p.into_inner()).retain(|s| s.id != self.id);
    }
}

/// The index: read from disk, downloaded from `index`, or found beside the file as htslib looks for it (X.crai
/// or the stem's .crai for a CRAM; X.csi, X.bai, then the stem's, for a BAM; a signed URL's query kept).
/// Returns the bytes and whether they came over the network.
fn load_index(url: &str, index: Option<&str>, cram: bool, retries: usize, shown: &str, stop: &Arc<AtomicBool>) -> Result<(Vec<u8>, bool)> {
    let get = |u: &str| -> Result<Vec<u8>, Failure> {
        let mut c = Client::new(stop.clone())?;
        let curl_u = url::Url::parse(u).map_err(|e| Failure::Final(format!("bad URL {}: {}", redact(u), e)))?.to_string();
        let a = c.get(&curl_u, None, MAX_INDEX)?;
        match a.status {
            200 if !a.over => Ok(a.body),
            200 => Err(Failure::Final(format!("{} is larger than an index can be", redact(u)))),
            s => Err(status_failure(s, &redact(u))),
        }
    };
    match index {
        Some(i) if i.contains("://") => Ok((with_retries(retries, "download the index", || get(i))?, true)),
        Some(i) => {
            let raw = std::fs::read(i).map_err(|e| anyhow::Error::new(Permanent(format!("cannot read the index {}: {}", i, e))))?;
            Ok((raw, false))
        }
        None => {
            let (base, query) = match url.split_once('?') {
                Some((b, q)) => (b, format!("?{}", q)),
                None => (url, String::new()),
            };
            let file = base.rsplit('/').next().unwrap_or("");
            let stem = file.rfind('.').map(|d| &base[..base.len() - file.len() + d]);
            let exts: &[&str] = if cram { &[".crai"] } else { &[".csi", ".bai"] };
            let tried: Vec<String> =
                std::iter::once(base).chain(stem).flat_map(|b| exts.iter().map(|e| format!("{}{}{}", b, e, query)).collect::<Vec<_>>()).collect();
            for u in &tried {
                let r = with_retries(retries, "download the index", || match get(u) {
                    // not there: the next name
                    Err(Failure::Final(m)) if ["HTTP 404", "HTTP 403", "HTTP 410", "HTTP 401"].iter().any(|s| m.starts_with(s)) => Ok(None),
                    r => r.map(Some),
                })?;
                if let Some(raw) = r {
                    return Ok((raw, true));
                }
            }
            Err(anyhow::Error::new(Permanent(format!(
                "cannot open {} with its index (no index beside it: {}, after the file name or in place of its extension, are not there, or refused)",
                shown,
                if cram { ".crai" } else { ".csi and .bai" }
            ))))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn crai(text: &str, size: u64) -> Crai {
        Crai::new(&parse_crai(text.as_bytes()).unwrap(), size).unwrap()
    }

    #[test]
    fn content_ranges_are_read() {
        assert_eq!(content_range(" bytes 0-99/1000"), Some((Some((0, 99)), Some(1000))));
        assert_eq!(content_range("bytes */1000"), Some((None, Some(1000))));
        assert_eq!(content_range("bytes 5-9/*"), Some((Some((5, 9)), None)));
        assert_eq!(content_range("items 0-1/2"), None);
        assert_eq!(content_range("bytes 0-x/2"), None);
    }

    #[test]
    fn a_crai_is_read_as_htslib_nests_it() {
        // ref 0: A [1,100] at 100, B [50,80] nested in A (200), C [150,300] (300); ref 1: D (400), a container
        // of refs 1 and 2 (500); unmapped (600)
        let c = crai("0\t1\t100\t100\t10\t5\n0\t50\t31\t200\t10\t5\n0\t150\t151\t300\t10\t5\n1\t1\t50\t400\t10\t5\n1\t60\t10\t500\t10\t5\n2\t1\t10\t500\t10\t5\n-1\t0\t150\t600\t10\t5\n", 700);
        assert_eq!(
            c.containers.iter().map(|x| (x.off, x.end)).collect::<Vec<_>>(),
            vec![(100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700)]
        );
        assert_eq!(c.containers[4].refs, vec![1, 2]);
        let tops: Vec<(i64, i64)> = c.tops[&0].iter().map(|t| (t.start, t.end)).collect();
        assert_eq!(tops, vec![(1, 100), (150, 300)], "B is nested in A");
        // htslib's start for a query: the first top entry whose end reaches the position
        assert_eq!(c.start_of(0, 1), Some(0));
        assert_eq!(c.start_of(0, 90), Some(0));
        assert_eq!(c.start_of(0, 101), Some(2));
        assert_eq!(c.start_of(0, 5000), Some(2), "past every entry: the last");
        assert_eq!(c.start_of(3, 1), None, "a reference without entries reads nothing");
        assert_eq!(c.start_of(-1, 0), Some(5));
        // gzip-compressed, as the .crai files are
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        std::io::Write::write_all(&mut gz, b"0\t1\t100\t100\t10\t5\n").unwrap();
        assert_eq!(parse_crai(&gz.finish().unwrap()).unwrap()[0], CraiLine { refid: 0, start: 1, end: 100, container: 100 });
        assert!(parse_crai(b"0\t1\t100\t100\t10\n").is_err());
        assert!(Crai::new(&parse_crai(b"0\t1\t100\t900\t10\t5\n").unwrap(), 700).is_err(), "a container past the end of the file");
    }

    #[test]
    fn a_query_reads_the_containers_htslib_walks() {
        // containers of 10 kB: five along ref 0 (10 kb of it each), one of ref 1, one of refs 1-3, one of ref 3,
        // one of the reads without a coordinate
        let mut t = String::new();
        for i in 0..5u64 {
            t += &format!("0\t{}\t10000\t{}\t10\t5\n", 1 + i * 10_000, 10_000 + i * 10_000);
        }
        t += "1\t1\t10000\t60000\t10\t5\n1\t10001\t500\t70000\t10\t5\n2\t1\t300\t70000\t10\t5\n3\t1\t200\t70000\t10\t5\n";
        t += "3\t201\t10000\t80000\t10\t5\n-1\t0\t150\t90000\t10\t5\n";
        let c = crai(&t, 100_000);
        let walk = |r: i32, beg: i64, end: i64| {
            let mut v = Vec::new();
            c.walk(r, beg + 1, end, &mut v);
            v.into_iter().map(|n| (n.start, n.end)).collect::<Vec<_>>()
        };
        let peek = |o: u64| (o, o + PEEK);
        // inside the second container: it whole, then the next one's header
        assert_eq!(walk(0, 12_000, 13_000), vec![(20_000, 30_000), peek(30_000)]);
        // across two containers; to the end of the reference, where the next reference's first header ends it
        assert_eq!(walk(0, 19_000, 21_000), vec![(20_000, 30_000), (30_000, 40_000), peek(40_000)]);
        assert_eq!(walk(0, 45_000, 60_000), vec![(50_000, 60_000), peek(60_000)]);
        // past the reference's last read: htslib starts at its last entry, a container wholly before the query
        assert_eq!(walk(0, 70_000, 71_000), vec![peek(50_000), peek(60_000)]);
        // ref 1 ends in a container of several references: read whole, and its records of ref 2 end the walk
        assert_eq!(walk(1, 10_200, 10_300), vec![(70_000, 80_000)]);
        // ref 3 starts in that container (its last reference): whole; then ref 3's own container, its header only
        // when the query ends before it, whole when the query reaches into it
        assert_eq!(walk(3, 0, 150), vec![(70_000, 80_000), peek(80_000)]);
        assert_eq!(walk(3, 0, 300), vec![(70_000, 80_000), (80_000, 90_000), peek(90_000)]);
        // the reads without a coordinate: every container from the first of them to the end
        assert_eq!(walk(-1, -1, 0), vec![(90_000, 100_000)]);
    }

    #[test]
    fn needs_become_requests() {
        let n = |s: u64, e: u64| Need { start: s, end: e, extend: None };
        // overlapping and close ranges join; far ones do not; a long range is cut
        let r = requests(vec![n(0, 10), n(5, 20), n(20 + JOIN, 40_000), n(200_000, 200_000 + 2 * MAX_REQUEST + 7)], JOIN);
        assert_eq!(
            r.iter().map(|x| (x.start, x.end)).collect::<Vec<_>>(),
            vec![
                (0, 40_000),
                (200_000, 200_000 + MAX_REQUEST),
                (200_000 + MAX_REQUEST, 200_000 + 2 * MAX_REQUEST),
                (200_000 + 2 * MAX_REQUEST, 200_000 + 2 * MAX_REQUEST + 7)
            ]
        );
        // BAM: a chunk [u, v) covers its whole blocks, and the block v points into is read to its end
        let mut v = Vec::new();
        bam_needs(1000 << 16 | 5, 5000 << 16 | 7, &mut v);
        bam_needs(5000 << 16 | 9, 5000 << 16 | 300, &mut v); // inside one block
        bam_needs(9000 << 16, 12000 << 16, &mut v); // ends on a block boundary
        bam_needs(11000 << 16 | 3, 11000 << 16 | 40, &mut v); // a block a range of whole blocks holds
        let r = requests(v, 0);
        assert_eq!(
            r,
            vec![n(1000, 5000), Need { start: 5000, end: 5000 + BGZF_HEADER, extend: Some(5000) }, n(9000, 12000)],
            "the tail block at 5000 once; the one at 11000 is inside [9000, 12000)"
        );
    }

    #[test]
    fn a_bgzf_header_gives_its_block_size() {
        let mut h = vec![0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0x34, 0x12];
        assert_eq!(bgzf_bsize(&h), Some(0x1234));
        h[12] = b'X';
        assert_eq!(bgzf_bsize(&h), None);
        assert_eq!(bgzf_bsize(&h[..10]), None);
    }

    /// A store over bytes in memory, as if a server held them: jobs, outside reads, failures, without a network.
    fn shared_over(data: &[u8], dir: &Path) -> Arc<Shared> {
        Arc::new(Shared {
            url: "http://127.0.0.1:9/none".into(),
            shown: "test".into(),
            size: data.len() as u64,
            etag: None,
            retries: 1,
            spool: spool_file(dir).unwrap(),
            spool_dir: dir.to_path_buf(),
            book: Mutex::new(Book::default()),
            cv: Condvar::new(),
            stop: Arc::new(AtomicBool::new(false)),
            clients: Mutex::new(Vec::new()),
            stats: Stats::default(),
        })
    }

    #[test]
    fn jobs_cover_only_what_no_job_covers() {
        let dir = std::env::temp_dir();
        let s = shared_over(&[0u8; 1000], &dir);
        let mut b = s.lock();
        let n = |a: u64, e: u64| Need { start: a, end: e, extend: None };
        assert_eq!(Shared::add(&mut b, n(100, 200)), 100);
        assert_eq!(Shared::add(&mut b, n(50, 300)), 150, "the parts before and after the job already there");
        assert_eq!(Shared::add(&mut b, n(120, 180)), 0);
        assert_eq!(b.jobs.iter().map(|(&a, j)| (a, j.end)).collect::<Vec<_>>(), vec![(50, 100), (100, 200), (200, 300)]);
        let t = Need { start: 250, end: 268, extend: Some(250) };
        assert_eq!(Shared::add(&mut b, t), 0, "a block header inside a job is not asked for again");
        let t = Need { start: 400, end: 418, extend: Some(400) };
        assert_eq!(Shared::add(&mut b, t), 18);
        assert_eq!(b.jobs[&400].extend, Some(400));
        // a block header right after a queued job is asked for with it: the job grows by the header
        assert_eq!(Shared::add(&mut b, n(500, 600)), 100);
        assert_eq!(Shared::add(&mut b, Need { start: 600, end: 618, extend: Some(600) }), 18);
        assert_eq!((b.jobs[&500].end, b.jobs[&500].extend, b.jobs.contains_key(&600)), (618, Some(600), false));
    }

    /// The store and htslib's view of it, without a network: a BAM in memory served as `ngsdose:` to htslib's
    /// reader, every job filled from the bytes, gives the same records as the file read from disk.
    #[test]
    fn htslib_reads_a_bam_through_the_backend() {
        use rust_htslib::bam::Read;
        register();
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let bam = root.join("tests/data/NA12878.subsample.bam");
        let data = std::fs::read(&bam).unwrap();
        let csi = std::fs::read(root.join("tests/data/NA12878.subsample.bam.csi")).unwrap();
        let dir = std::env::temp_dir();
        let s = shared_over(&data, &dir);
        // every byte as one done job, written to the spool as a download would be
        {
            let mut b = s.lock();
            Shared::add(&mut b, Need { start: 0, end: data.len() as u64, extend: None });
            b.jobs.get_mut(&0).unwrap().state = State::Running;
        }
        s.settle(0, Ok(data.clone()));
        let id = NEXT_ID.fetch_add(1, Relaxed);
        REGISTRY.lock().unwrap().push(Served { id, file: Arc::downgrade(&s), index: Arc::new(csi) });
        let url = url::Url::parse(&format!("ngsdose:{}##idx##ngsdose:{}.idx", id, id)).unwrap();
        let mut via = rust_htslib::bam::IndexedReader::from_url(&url).unwrap();
        let mut local = rust_htslib::bam::IndexedReader::from_path(&bam).unwrap();
        let tid = local.header().tid(b"chr1").unwrap();
        for (a, b) in [(1_000_000i64, 1_100_000i64), (50_000_000, 50_500_000), (0, 1 << 29)] {
            let mut keys = Vec::new();
            for rd in [&mut via, &mut local] {
                rd.fetch((tid, a, b)).unwrap();
                let mut k = Vec::new();
                for r in rd.records() {
                    let r = r.unwrap();
                    k.push((r.pos(), r.qname().to_vec(), r.flags()));
                }
                keys.push(k);
            }
            assert!(!keys[0].is_empty() || b - a < 1 << 20);
            assert_eq!(keys[0], keys[1]);
        }
        unsafe { assert_eq!(htslib::hts_check_EOF(via.htsfile()), 1, "the end-of-file marker, through a seek from the end") };
        drop(via);
        REGISTRY.lock().unwrap().retain(|s| s.id != id);
    }
}
