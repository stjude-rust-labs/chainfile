//! Input and output helpers.

use std::fs::File;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Cursor;
use std::io::Read;
use std::io::Write;
use std::path::Path;

use anyhow::Context as _;
use anyhow::Result;
use chainfile::liftover::Machine;
use chainfile::liftover::machine::builder::Builder;
use flate2::read::MultiGzDecoder;
use tracing::debug;
use tracing::info;

/// The path that refers to standard input or standard output.
const STDIO: &str = "-";

/// The magic bytes at the start of every gzip (and BGZF) stream.
const GZIP_MAGIC: [u8; 2] = [0x1F, 0x8B];

/// The capacity of input and output buffers.
const BUFFER_SIZE: usize = 128 * 1024;

/// Returns whether a path refers to standard input or standard output.
fn is_stdio(path: &Path) -> bool {
    path.as_os_str() == STDIO
}

/// Opens an input for reading.
///
/// A path of `-` reads from standard input. Gzip-compressed (including BGZF)
/// inputs are detected from their contents and decompressed transparently.
pub fn open_input(path: &Path) -> Result<Box<dyn BufRead>> {
    debug!(path = %path.display(), "opening input");

    let inner: Box<dyn Read> = if is_stdio(path) {
        Box::new(io::stdin().lock())
    } else {
        Box::new(File::open(path).with_context(|| format!("opening `{}`", path.display()))?)
    };

    decompress(inner).with_context(|| format!("reading `{}`", path.display()))
}

/// Wraps a reader so that gzip-compressed data is decompressed.
///
/// The bytes consumed while sniffing the format are replayed.
fn decompress(mut inner: Box<dyn Read>) -> io::Result<Box<dyn BufRead>> {
    let mut magic = Vec::with_capacity(GZIP_MAGIC.len());
    (&mut inner)
        .take(GZIP_MAGIC.len() as u64)
        .read_to_end(&mut magic)?;

    let replayed = Cursor::new(magic.clone()).chain(inner);

    if magic == GZIP_MAGIC {
        debug!("detected gzip-compressed input");
        Ok(Box::new(BufReader::with_capacity(
            BUFFER_SIZE,
            MultiGzDecoder::new(replayed),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(BUFFER_SIZE, replayed)))
    }
}

/// Opens an output for writing.
///
/// [`None`] or a path of `-` writes to standard output.
pub fn open_output(path: Option<&Path>) -> Result<Box<dyn Write>> {
    debug!(
        path = %path.map(|p| p.display().to_string()).unwrap_or_else(|| String::from(STDIO)),
        "opening output"
    );

    match path {
        None => Ok(Box::new(BufWriter::with_capacity(
            BUFFER_SIZE,
            io::stdout().lock(),
        ))),
        Some(path) if is_stdio(path) => Ok(Box::new(BufWriter::with_capacity(
            BUFFER_SIZE,
            io::stdout().lock(),
        ))),
        Some(path) => File::create(path)
            .map(|file| Box::new(BufWriter::with_capacity(BUFFER_SIZE, file)) as Box<dyn Write>)
            .with_context(|| format!("creating `{}`", path.display())),
    }
}

/// Opens a chain file for reading.
pub fn open_chain(path: &Path) -> Result<chainfile::Reader<Box<dyn BufRead>>> {
    open_input(path).map(chainfile::Reader::new)
}

/// Builds a liftover machine from a chain file.
pub fn load_machine(path: &Path) -> Result<Machine> {
    info!(path = %path.display(), "loading chain file");

    let reader = open_chain(path)?;
    let machine = Builder
        .try_build_from(reader)
        .with_context(|| format!("building liftover machine from `{}`", path.display()))?;

    info!(
        reference_contigs = machine.reference_chromosomes().len(),
        "loaded chain file"
    );

    Ok(machine)
}

#[cfg(test)]
mod tests {
    use flate2::Compression;
    use flate2::write::GzEncoder;

    use super::*;

    fn read_all(reader: Box<dyn Read>) -> String {
        let mut result = String::new();
        decompress(reader)
            .unwrap()
            .read_to_string(&mut result)
            .unwrap();
        result
    }

    fn gzip(data: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(data).unwrap();
        encoder.finish().unwrap()
    }

    /// A reader that returns at most one byte per read.
    struct Trickle(Cursor<Vec<u8>>);

    impl Read for Trickle {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = buf.len().min(1);
            self.0.read(&mut buf[..n])
        }
    }

    #[test]
    fn plain_input_is_replayed() {
        assert_eq!(read_all(Box::new(&b"hello"[..])), "hello");
        assert_eq!(read_all(Box::new(&b"h"[..])), "h");
        assert_eq!(read_all(Box::new(&b""[..])), "");
    }

    #[test]
    fn gzip_input_is_decompressed() {
        let data = gzip(b"hello");
        assert_eq!(read_all(Box::new(Cursor::new(data))), "hello");
    }

    #[test]
    fn concatenated_gzip_members_are_decompressed() {
        let mut data = gzip(b"hello, ");
        data.extend(gzip(b"world"));
        assert_eq!(read_all(Box::new(Cursor::new(data))), "hello, world");
    }

    #[test]
    fn short_reads_are_handled() {
        let data = gzip(b"hello");
        assert_eq!(read_all(Box::new(Trickle(Cursor::new(data)))), "hello");
        assert_eq!(
            read_all(Box::new(Trickle(Cursor::new(b"hello".to_vec())))),
            "hello"
        );
    }

    #[test]
    fn truncated_gzip_input_is_an_error() {
        let data = gzip(b"hello, world");
        let truncated = data[..data.len() / 2].to_vec();

        let mut result = String::new();
        assert!(
            decompress(Box::new(Cursor::new(truncated)))
                .unwrap()
                .read_to_string(&mut result)
                .is_err()
        );
    }
}
