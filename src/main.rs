use std::env;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::process::ExitCode;

use fqfilter::{classify, fastq, Config, Decision, Reason, Reader, Stats};

const USAGE: &str = "\
usage: fqfilter [-q MIN_MEAN] [-l MIN_LENGTH] [INPUT]

Streaming FASTQ quality filter. Reads INPUT (or stdin) and writes passing
records to stdout. Summary counters go to stderr.

  -q MIN_MEAN    minimum mean Phred score to keep a read (default 20)
  -l MIN_LENGTH  minimum sequence length to keep a read (default 1)
  -h, --help     show this help

Exit status: 0 on success, 1 on usage or parse error.
";

fn parse_args() -> Result<(Config, Option<String>), String> {
    let mut cfg = Config::default();
    let mut input = None;
    let mut args = env::args().skip(1);

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => return Err(USAGE.to_string()),
            "-q" => {
                let value = args.next().ok_or("-q needs a value")?;
                cfg.min_mean = value
                    .parse()
                    .map_err(|_| format!("invalid -q value: {value}"))?;
            }
            "-l" => {
                let value = args.next().ok_or("-l needs a value")?;
                cfg.min_length = value
                    .parse()
                    .map_err(|_| format!("invalid -l value: {value}"))?;
            }
            other if other.starts_with('-') => return Err(format!("unknown flag: {other}")),
            other => {
                if input.replace(other.to_string()).is_some() {
                    return Err("more than one input file given".into());
                }
            }
        }
    }

    Ok((cfg, input))
}

fn open_input(path: Option<&str>) -> io::Result<Box<dyn Read>> {
    match path {
        Some(p) => Ok(Box::new(File::open(p)?)),
        None => Ok(Box::new(io::stdin())),
    }
}

fn run() -> Result<(), String> {
    let (cfg, input) = parse_args()?;
    let source = |e: io::Error| format!("cannot read input: {e}");
    let mut reader = Reader::new(BufReader::with_capacity(
        1 << 16,
        open_input(input.as_deref()).map_err(source)?,
    ));
    let stdout = io::stdout();
    let mut out = BufWriter::with_capacity(1 << 16, stdout.lock());
    let mut stats = Stats::default();

    loop {
        let record = match reader.next_record() {
            Ok(Some(record)) => record,
            Ok(None) => break,
            Err(fastq::ParseError::Io(e)) => return Err(format!("read error: {e}")),
            Err(e) => return Err(e.to_string()),
        };
        let decision = classify(&record, &cfg);
        stats.record(decision);
        if decision == Decision::Keep {
            record
                .write_to(&mut out)
                .map_err(|e| format!("write error: {e}"))?;
        }
    }

    out.flush().map_err(|e| format!("write error: {e}"))?;
    eprintln!(
        "seen={} kept={} low_quality={} too_short={}",
        stats.seen, stats.kept, stats.low_quality, stats.too_short
    );
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(msg) => {
            let help_only = msg.starts_with("usage:");
            if help_only {
                print!("{msg}");
                ExitCode::SUCCESS
            } else {
                eprintln!("fqfilter: {msg}");
                eprintln!("try 'fqfilter --help'");
                let code = u8::try_from(Reason::LowQuality as u8 + 1).unwrap_or(1);
                ExitCode::from(code)
            }
        }
    }
}
