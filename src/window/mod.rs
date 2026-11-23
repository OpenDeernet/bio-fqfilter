//! Sliding-window trimming of low-quality read ends.
//!
//! A base is trimmed from a 5' end while the fixed-size window starting there
//! scores below the threshold, from the 3' end equivalently. Windows are
//! evaluated with the same mean Phred score used by the filter, so a window
//! is only kept when the whole window clears the threshold.

use crate::fastq::Record;
use crate::quality::mean_phred;

/// Trims low-quality ends of `record` in place and returns the number of bases
/// removed (0 when nothing moved).
///
/// A window wider than the read trims nothing, so callers should validate the
/// requested width against the shortest expected read.
pub fn trim_ends(record: &mut Record, window: usize, min_mean: f64) -> usize {
    let len = record.qual.len();
    if window == 0 || window > len {
        return 0;
    }

    let mut start = 0;
    while start + window <= len && mean_phred(&record.qual[start..start + window]) < min_mean {
        start += 1;
    }

    // Fewer bases left than one window: everything remaining sat in windows
    // that scored below the threshold, so there is no core to keep.
    if start + window > len {
        record.seq.clear();
        record.qual.clear();
        return len;
    }

    let mut end = len;
    while end >= start + window && mean_phred(&record.qual[end - window..end]) < min_mean {
        end -= 1;
    }

    record.seq.truncate(end);
    record.seq.drain(..start);
    record.qual.truncate(end);
    record.qual.drain(..start);
    len - (end - start)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(seq: &str, qual: &str) -> Record {
        Record {
            header: "r".into(),
            seq: seq.as_bytes().to_vec(),
            plus: "+".into(),
            qual: qual.as_bytes().to_vec(),
        }
    }

    #[test]
    fn trims_both_ends_but_keeps_core() {
        let mut rec = record("AAAACCCCGGGGTTTT", "!!!!IIIIIIII!!!!");
        // Boundary windows mix junk and strong bases, so one base per side
        // stays: the window mean there is exactly the threshold.
        let removed = trim_ends(&mut rec, 4, 20.0);
        assert_eq!(removed, 4);
        assert_eq!(rec.seq, b"AACCCCGGGGTT");
        assert_eq!(rec.qual, b"!!IIIIIIII!!");
    }

    #[test]
    fn leaves_clean_read_untouched() {
        let mut rec = record("ACGT", "IIII");
        assert_eq!(trim_ends(&mut rec, 2, 20.0), 0);
        assert_eq!(rec.seq, b"ACGT");
    }

    #[test]
    fn oversize_window_is_a_noop() {
        let mut rec = record("AC", "II");
        assert_eq!(trim_ends(&mut rec, 5, 20.0), 0);
        assert_eq!(rec.seq, b"AC");
    }

    #[test]
    fn all_low_quality_collapses_to_empty() {
        let mut rec = record("ACGT", "!!!!");
        assert_eq!(trim_ends(&mut rec, 2, 20.0), 4);
        assert!(rec.seq.is_empty());
        assert!(rec.qual.is_empty());
    }
}
