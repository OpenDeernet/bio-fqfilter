# bio-fqfilter

Streaming FASTQ quality filter for the front of a read-QC pipeline.

Long-read and single-cell runs routinely arrive with a tail of low-quality
reads that are cheap to drop before assembly or quantification. `fqfilter`
does exactly that one job: it reads FASTQ (plain or gzipped through a pipe),
keeps records whose mean Phred score clears a threshold, and writes the rest
unchanged to stdout. Nothing is buffered beyond the current record, so a
multi-gigabyte file costs the same memory as the first read.

## Build

```sh
cargo build --release
# binary: target/release/fqfilter
```

Requires Rust 1.70 or newer; there are no crate dependencies.

## Usage

```sh
fqfilter -q 20 reads.fastq > reads.q20.fastq
zcat reads.fastq.gz | fqfilter -q 25 -l 500 | gzip > reads.q25.fastq.gz
fqfilter -q 22 -w 4 reads.fastq > reads.trimmed.fastq
```

Options:

| Flag | Meaning | Default |
| --- | --- | --- |
| `-q MIN_MEAN` | keep reads with mean Phred at least this value | `20` |
| `-l MIN_LENGTH` | drop reads shorter than this before scoring | `1` |
| `-w WINDOW` | trim low-quality ends with a sliding window of this width | off |
| `-h`, `--help` | usage | |

Passing records are byte-identical to the input (line endings normalised to
LF). The per-run counters `seen`, `kept`, `low_quality` and `too_short` go to
stderr so stdout stays pipeable; with `-w` the count also reports
`trimmed_bases`. Malformed input aborts with a line number and
exit status 1.

### Score semantics

Qualities are decoded as Phred+33 and clamped to the valid Q0-Q93 range. The
score is the arithmetic mean over the whole read; a single low base in an
otherwise strong read does not drop it. For adapter or poly-G trimming, run a
trimmer first and filter the trimmed file.

With `-w`, ends are trimmed before scoring: while the window at either end
scores below `MIN_MEAN`, one base is dropped from that end. The same threshold
is used for window and read, so trimming and filtering stay consistent. A read
whose remaining length is shorter than one window is emptied and therefore
dropped. Boundary windows can leave a low base or two next to the core; that
window already met the threshold as a whole.

### Memory bounds

One record plus two 64 KiB buffers (`BufReader`, `BufWriter`). Peak RSS scales
with the longest single read, not with the file size.

## Layout

```
src/lib.rs      config, keep/drop decision, run statistics
src/fastq.rs    streaming FASTQ parser and writer
src/quality.rs  Phred decoding and mean score
src/window/     sliding-window end trimming
src/main.rs     CLI: argument parsing, I/O wiring, summary
```

The library is usable without the binary: `fqfilter::classify` scores a
parsed record against a `Config`.

## Testing

```sh
cargo test
```

Tests cover parser edge cases (CRLF, truncated quality lines, bad headers), the
decision order in `classify`, and trimming end-to-end in `tests/trim.rs`.

## Citation

If this filter contributes to published work, cite the repository:

> Nair, R. (2026). bio-fqfilter: streaming FASTQ quality filtering
> (version 0.1.0). https://github.com/d3orr56hoq/bio-fqfilter

## License

MIT, see [LICENSE](LICENSE).
