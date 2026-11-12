//! Input helpers.

use std::fs::File;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Cursor;
use std::io::Read;
use std::path::Path;

use anyhow::Context as _;
use anyhow::Result;
use flate2::read::MultiGzDecoder;

/// The magic bytes at the start of every gzip (and BGZF) stream.
const GZIP_MAGIC: [u8; 2] = [0x1F, 0x8B];

/// Opens a file for reading.
///
/// Gzip-compressed (including BGZF) files are detected from their contents
/// and decompressed transparently.
pub fn open_input(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening `{}`", path.display()))?;
    decompress(Box::new(file)).with_context(|| format!("reading `{}`", path.display()))
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
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(replayed))))
    } else {
        Ok(Box::new(BufReader::new(replayed)))
    }
}

/// Reads an entire (optionally gzip-compressed) file into memory.
pub fn read_to_end(path: &Path) -> Result<Vec<u8>> {
    let mut data = Vec::new();
    open_input(path)?
        .read_to_end(&mut data)
        .with_context(|| format!("reading `{}`", path.display()))?;
    Ok(data)
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use flate2::Compression;
    use flate2::write::GzEncoder;

    use super::*;

    #[test]
    fn plain_and_gzip_inputs_are_read() {
        let directory = tempfile::tempdir().unwrap();

        let plain = directory.path().join("plain.txt");
        std::fs::write(&plain, b"hello").unwrap();
        assert_eq!(read_to_end(&plain).unwrap(), b"hello");

        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b"hello").unwrap();
        let compressed = directory.path().join("compressed.txt.gz");
        std::fs::write(&compressed, encoder.finish().unwrap()).unwrap();
        assert_eq!(read_to_end(&compressed).unwrap(), b"hello");

        let empty = directory.path().join("empty.txt");
        std::fs::write(&empty, b"").unwrap();
        assert_eq!(read_to_end(&empty).unwrap(), b"");
    }
}
