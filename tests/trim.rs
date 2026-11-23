//! Integration test: trimming feeds the keep/drop decision.
//!
//! A read with a strong core but junk ends must survive `-q 20` once a window
//! is applied, and must be dropped without one.

use std::io::Cursor;

use fqfilter::window::trim_ends;
use fqfilter::{classify, fastq, Config, Decision, Reader, Reason};

fn first_record(input: &str) -> fastq::Record {
    Reader::new(Cursor::new(input))
        .next_record()
        .unwrap()
        .unwrap()
}

const JUNK_ENDS: &str = "@r\nAAAACCCCGGGGTTTT\n+\n!!!!IIIIIIII!!!!\n";

#[test]
fn trimming_lifts_a_borderline_read_over_threshold() {
    let mut record = first_record(JUNK_ENDS);
    let cfg = Config { min_mean: 22.0, min_length: 4 };

    assert_eq!(classify(&record, &cfg), Decision::Drop(Reason::LowQuality));

    let removed = trim_ends(&mut record, 4, 22.0);
    assert_eq!(removed, 6);
    assert_eq!(classify(&record, &cfg), Decision::Keep);
}

#[test]
fn trimmed_read_shorter_than_min_length_is_still_dropped() {
    let mut record = first_record(JUNK_ENDS);
    let cfg = Config { min_mean: 22.0, min_length: 11 };
    trim_ends(&mut record, 4, 22.0);
    assert_eq!(record.seq.len(), 10);
    assert_eq!(classify(&record, &cfg), Decision::Drop(Reason::TooShort));
}
