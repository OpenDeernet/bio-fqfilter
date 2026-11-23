//! Phred score helpers for FASTQ qualities.
//!
//! FASTQ Phred+33 encodes Q0..Q93 as `!`..`~`. Bytes outside that range are
//! clamped rather than rejected so that slightly non-conforming
//! (Phred+64, padded) data still gets a conservative score.

pub const OFFSET: u8 = 33;
pub const MAX_SCORE: u8 = 93;

/// Converts one encoded byte to a Phred score.
#[inline]
pub fn phred(byte: u8) -> u8 {
    byte.saturating_sub(OFFSET).min(MAX_SCORE)
}

/// Arithmetic mean Phred score of `qual`; empty input scores 0.
#[inline]
pub fn mean_phred(qual: &[u8]) -> f64 {
    if qual.is_empty() {
        return 0.0;
    }
    let sum: u32 = qual.iter().map(|&b| u32::from(phred(b))).sum();
    f64::from(sum) / qual.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn maps_ascii_range() {
        assert_eq!(phred(b'!'), 0);
        assert_eq!(phred(b'I'), 40);
        assert_eq!(phred(b'~'), MAX_SCORE);
        assert_eq!(phred(0), 0);
    }

    #[test]
    fn averages_scores() {
        assert_eq!(mean_phred(b"I!I!"), 20.0);
        assert_eq!(mean_phred(b""), 0.0);
    }
}
