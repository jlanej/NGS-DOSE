//! 2-bit k-mer encoding and canonical k-mer iteration (k <= 32).

/// Encode an ASCII base to 2 bits; anything that is not A/C/G/T returns 4.
#[inline]
pub fn enc(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// BAM 4-bit nibble (=ACMGRSVTWYHKDBN) to 2 bits; ambiguity codes return 4.
pub const NIB2BIT: [u8; 16] = [4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4];

/// Decode a BAM packed sequence into 2-bit codes (4 = not ACGT). Returns the GC count.
pub fn unpack_bam_seq(encoded: &[u8], len: usize, out: &mut Vec<u8>) -> u32 {
    out.clear();
    out.reserve(len);
    let mut gc = 0u32;
    for i in 0..len {
        let byte = encoded[i >> 1];
        let nib = if i & 1 == 0 { byte >> 4 } else { byte & 0xf };
        let c = NIB2BIT[nib as usize];
        gc += (c == 1 || c == 2) as u32;
        out.push(c);
    }
    gc
}

pub fn encode_ascii(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| enc(b)).collect()
}

/// One k-mer occurrence: start offset in the sequence, canonical value, and whether the
/// forward-strand k-mer is the canonical one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Kmer {
    pub offset: usize,
    pub canon: u64,
    pub fwd_is_canon: bool,
}

/// Rolling canonical k-mer iterator over 2-bit codes; k-mers spanning a non-ACGT base are skipped.
pub struct KmerIter<'a> {
    codes: &'a [u8],
    k: usize,
    mask: u64,
    rc_shift: u32,
    i: usize,
    valid: usize,
    fwd: u64,
    rc: u64,
}

impl<'a> KmerIter<'a> {
    pub fn new(codes: &'a [u8], k: usize) -> Self {
        assert!((1..=32).contains(&k), "k must be in 1..=32");
        let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
        KmerIter { codes, k, mask, rc_shift: (2 * (k - 1)) as u32, i: 0, valid: 0, fwd: 0, rc: 0 }
    }
}

impl<'a> Iterator for KmerIter<'a> {
    type Item = Kmer;
    #[inline]
    fn next(&mut self) -> Option<Kmer> {
        while self.i < self.codes.len() {
            let c = self.codes[self.i];
            self.i += 1;
            if c > 3 {
                self.valid = 0;
                self.fwd = 0;
                self.rc = 0;
                continue;
            }
            self.fwd = ((self.fwd << 2) | c as u64) & self.mask;
            self.rc = (self.rc >> 2) | ((3 - c as u64) << self.rc_shift);
            self.valid += 1;
            if self.valid >= self.k {
                let fwd_is_canon = self.fwd <= self.rc;
                return Some(Kmer { offset: self.i - self.k, canon: if fwd_is_canon { self.fwd } else { self.rc }, fwd_is_canon });
            }
        }
        None
    }
}

pub fn kmer_to_string(mut v: u64, k: usize) -> String {
    let mut s = vec![b'A'; k];
    for i in (0..k).rev() {
        s[i] = b"ACGT"[(v & 3) as usize];
        v >>= 2;
    }
    String::from_utf8(s).unwrap()
}

pub fn string_to_kmer(s: &str) -> Option<u64> {
    let mut v = 0u64;
    for b in s.bytes() {
        let c = enc(b);
        if c > 3 {
            return None;
        }
        v = (v << 2) | c as u64;
    }
    Some(v)
}

/// Reverse complement of a 2-bit packed k-mer.
pub fn revcomp_kmer(v: u64, k: usize) -> u64 {
    debug_assert!((1..=32).contains(&k));
    // complement, reverse the 2-bit groups within the u64, then drop the unused low bits
    let mut x = !v;
    x = ((x >> 2) & 0x3333_3333_3333_3333) | ((x & 0x3333_3333_3333_3333) << 2);
    x = ((x >> 4) & 0x0F0F_0F0F_0F0F_0F0F) | ((x & 0x0F0F_0F0F_0F0F_0F0F) << 4);
    x = x.swap_bytes();
    x >> (64 - 2 * k)
}

/// The canonical form of a packed k-mer: the smaller of it and its reverse complement.
pub fn canonical(v: u64, k: usize) -> u64 {
    v.min(revcomp_kmer(v, k))
}

#[cfg(test)]
pub fn revcomp_ascii(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_is_strand_symmetric() {
        let s = b"ACGTTGCATGCCGATAGCTAGCTAGGATCCGATCGATTTACG";
        let r = revcomp_ascii(s);
        let k = 11;
        let a: Vec<u64> = KmerIter::new(&encode_ascii(s), k).map(|x| x.canon).collect();
        let mut b: Vec<u64> = KmerIter::new(&encode_ascii(&r), k).map(|x| x.canon).collect();
        b.reverse();
        assert_eq!(a, b);
        assert_eq!(a.len(), s.len() - k + 1);
    }

    #[test]
    fn n_resets_window() {
        let s = b"ACGTACGTNACGTACGT";
        let v: Vec<usize> = KmerIter::new(&encode_ascii(s), 4).map(|x| x.offset).collect();
        assert_eq!(v, vec![0, 1, 2, 3, 4, 9, 10, 11, 12, 13]);
    }

    #[test]
    fn roundtrip_and_k32() {
        let s = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        assert_eq!(kmer_to_string(string_to_kmer(s).unwrap(), 32), s);
        let it: Vec<Kmer> = KmerIter::new(&encode_ascii(s.as_bytes()), 32).collect();
        assert_eq!(it.len(), 1);
    }

    #[test]
    fn packed_revcomp_matches_the_ascii_one() {
        let s = b"ACGTTGCATGCCGATAGCTAGCTAGGATCCGAT";
        for k in [1, 5, 11, 31, 32] {
            let f = std::str::from_utf8(&s[..k]).unwrap();
            let r = String::from_utf8(revcomp_ascii(&s[..k])).unwrap();
            let (fv, rv) = (string_to_kmer(f).unwrap(), string_to_kmer(&r).unwrap());
            assert_eq!(revcomp_kmer(fv, k), rv, "k={}", k);
            assert_eq!(canonical(fv, k), canonical(rv, k));
            let it = KmerIter::new(&encode_ascii(&s[..k]), k).next().unwrap();
            assert_eq!(it.canon, canonical(fv, k));
        }
    }

    #[test]
    fn unpack_matches_ascii() {
        // "ACGTN" packed: A=1 C=2 G=4 T=8 N=15 -> 0x12 0x48 0xF0
        let enc = [0x12u8, 0x48, 0xF0];
        let mut out = Vec::new();
        let gc = unpack_bam_seq(&enc, 5, &mut out);
        assert_eq!(out, vec![0, 1, 2, 3, 4]);
        assert_eq!(gc, 2);
    }
}
