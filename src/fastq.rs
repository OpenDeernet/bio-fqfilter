//! Minimal streaming FASTQ reader.
//!
//! Only the record being parsed is allocated; callers can reuse a `Vec`-free
//! API by handing the record onward or letting it drop.

use std::fmt;
use std::io::{self, BufRead};

/// A single FASTQ record, stored as raw bytes to avoid UTF-8 validation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// Header line without the leading `@` and the line terminator.
    pub header: String,
    pub seq: Vec<u8>,
    /// Separator line, including any text after `+`.
    pub plus: String,
    /// Phred+33 qualities, same length as `seq`.
    pub qual: Vec<u8>,
}

impl Record {
    /// Writes the record back out in FASTQ form.
    pub fn write_to<W: io::Write>(&self, out: &mut W) -> io::Result<()> {
        out.write_all(b"@")?;
        out.write_all(self.header.as_bytes())?;
        out.write_all(b"\n")?;
        out.write_all(&self.seq)?;
        out.write_all(b"\n")?;
        out.write_all(self.plus.as_bytes())?;
        out.write_all(b"\n")?;
        out.write_all(&self.qual)?;
        out.write_all(b"\n")
    }
}

/// Parse failure, tagged with the 1-based line where it was detected.
#[derive(Debug)]
pub enum ParseError {
    Io(io::Error),
    Malformed { line: u64, reason: &'static str },
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::Io(e) => write!(f, "{e}"),
            ParseError::Malformed { line, reason } => {
                write!(f, "line {line}: {reason}")
            }
        }
    }
}

impl std::error::Error for ParseError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ParseError::Io(e) => Some(e),
            ParseError::Malformed { .. } => None,
        }
    }
}

impl From<io::Error> for ParseError {
    fn from(e: io::Error) -> Self {
        ParseError::Io(e)
    }
}

/// Streaming reader; `next_record` borrows no state beyond the current line.
pub struct Reader<R: BufRead> {
    inner: R,
    line: u64,
}

impl<R: BufRead> Reader<R> {
    pub fn new(inner: R) -> Self {
        Self { inner, line: 0 }
    }

    fn read_line(&mut self, buf: &mut Vec<u8>) -> Result<usize, ParseError> {
        buf.clear();
        self.line += 1;
        let n = self.inner.read_until(b'\n', buf)?;
        while matches!(buf.last(), Some(b'\n') | Some(b'\r')) {
            buf.pop();
        }
        Ok(n)
    }

    /// Returns the next record, or `None` at a clean end of input.
    pub fn next_record(&mut self) -> Result<Option<Record>, ParseError> {
        let mut header = Vec::new();
        if self.read_line(&mut header)? == 0 {
            return Ok(None);
        }
        if header.first() != Some(&b'@') {
            return Err(self.malformed("header does not start with '@'"));
        }

        let mut seq = Vec::new();
        let mut plus = Vec::new();
        let mut qual = Vec::new();
        self.read_line(&mut seq)?;
        self.read_line(&mut plus)?;
        self.read_line(&mut qual)?;

        if plus.first() != Some(&b'+') {
            return Err(self.malformed("third line does not start with '+'"));
        }
        if seq.len() != qual.len() {
            return Err(self.malformed("sequence and quality lengths differ"));
        }
        if seq.is_empty() {
            return Err(self.malformed("empty sequence"));
        }

        header.remove(0);
        Ok(Some(Record {
            header: String::from_utf8_lossy(&header).into_owned(),
            seq,
            plus: String::from_utf8_lossy(&plus).into_owned(),
            qual,
        }))
    }

    fn malformed(&self, reason: &'static str) -> ParseError {
        ParseError::Malformed {
            line: self.line,
            reason,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reader(input: &str) -> Reader<&[u8]> {
        Reader::new(input.as_bytes())
    }

    #[test]
    fn parses_record_and_trims_crlf() {
        let mut r = reader("@read1 note\r\nACGT\r\n+\r\nIIII\r\n");
        let rec = r.next_record().unwrap().unwrap();
        assert_eq!(rec.header, "read1 note");
        assert_eq!(rec.seq, b"ACGT");
        assert_eq!(rec.qual, b"IIII");
        assert!(r.next_record().unwrap().is_none());
    }

    #[test]
    fn reports_mismatched_quality_length() {
        let mut r = reader("@r\nACGT\n+\nIII\n");
        let err = r.next_record().unwrap_err();
        assert_eq!(err.to_string(), "line 4: sequence and quality lengths differ");
    }

    #[test]
    fn reports_bad_header() {
        let mut r = reader("r\nACGT\n+\nIIII\n");
        assert!(matches!(
            r.next_record(),
            Err(ParseError::Malformed { line: 1, .. })
        ));
    }
}
