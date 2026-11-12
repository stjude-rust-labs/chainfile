//! Ordered, multithreaded processing of line-based records.
//!
//! The input is split into chunks of whole lines by a reader thread. Chunks
//! are processed by a pool of worker threads, and their outputs are written
//! by the calling thread in input order, so the output is identical to
//! processing the records one at a time.

use std::collections::VecDeque;
use std::io;
use std::io::BufRead as _;
use std::mem;
use std::num::NonZeroUsize;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::PoisonError;
use std::sync::mpsc;
use std::sync::mpsc::Receiver;
use std::sync::mpsc::SyncSender;
use std::sync::mpsc::TryRecvError;
use std::thread;

use anyhow::Context as _;
use anyhow::Result;
use anyhow::anyhow;
use tracing::debug;

use crate::io::Input;

/// The maximum number of lines in a chunk.
const CHUNK_LINES: usize = 4096;

/// The size at which a chunk is considered full.
const CHUNK_BYTES: usize = 1024 * 1024;

/// The number of chunks that may be in flight per worker.
const CHUNKS_PER_WORKER: usize = 4;

/// A chunk of whole lines (each of which is one record).
pub struct Chunk {
    /// The 1-based number of the first record in the chunk.
    pub first: u64,

    /// The number of lines in the chunk.
    pub lines: usize,

    /// The lines, including their line endings.
    pub data: Vec<u8>,
}

/// The result of processing a chunk: its output and whether it succeeded.
type Processed<O> = (O, Result<()>);

/// Gets the number of worker threads to use.
///
/// When no number is provided, the available parallelism is used.
pub fn threads(threads: Option<NonZeroUsize>) -> NonZeroUsize {
    threads.unwrap_or_else(|| thread::available_parallelism().unwrap_or(NonZeroUsize::MIN))
}

/// Processes the records of an input in parallel.
///
/// One worker thread is started for each state. Each chunk is processed by
/// `process` into a default-initialized output, and the outputs are passed to
/// `write` in input order. The output of a chunk that fails to process is
/// written before the error is returned, so the records before the failing
/// record are written as they would be when processing serially.
///
/// `kind` names the records (for example, `BED`) in error messages. Lines that
/// start with `comment` (when provided) are skipped and are not counted as
/// records.
pub fn process<S, O, P, W>(
    input: Input,
    kind: &'static str,
    comment: Option<u8>,
    states: Vec<S>,
    process: P,
    mut write: W,
) -> Result<()>
where
    S: Send,
    O: Default + Send,
    P: Fn(&mut S, &Chunk, &mut O) -> Result<()> + Sync,
    W: FnMut(O) -> Result<()>,
{
    let in_flight = states.len() * CHUNKS_PER_WORKER;
    debug!(workers = states.len(), "processing records in parallel");

    // NOTE: the reader thread is not scoped, so returning early (for example,
    // after an invalid record) does not wait for a read that may block
    // indefinitely (such as from standard input).
    let (chunk_tx, chunk_rx) = mpsc::sync_channel(in_flight);
    thread::spawn(move || read(input, kind, comment, &chunk_tx));

    thread::scope(|scope| {
        let (job_tx, job_rx) = mpsc::channel::<(Chunk, SyncSender<Processed<O>>)>();
        let job_rx = Arc::new(Mutex::new(job_rx));

        for mut state in states {
            let job_rx = Arc::clone(&job_rx);
            let process = &process;

            scope.spawn(move || {
                loop {
                    let job = job_rx.lock().unwrap_or_else(PoisonError::into_inner).recv();

                    let Ok((chunk, tx)) = job else {
                        break;
                    };

                    let mut output = O::default();
                    let result = process(&mut state, &chunk, &mut output);

                    // NOTE: the writer stops receiving after an error, so a
                    // failed send is not an error.
                    let _ = tx.send((output, result));
                }
            });
        }

        let mut pending = VecDeque::<Receiver<Processed<O>>>::new();
        let mut read_error = None;
        let mut reading = true;

        loop {
            while reading && pending.len() < in_flight {
                // NOTE: the next chunk is only waited for when no results are
                // pending, as the reader may block indefinitely (for example,
                // on standard input) while pending results are ready to write.
                let next = if pending.is_empty() {
                    chunk_rx.recv().map_err(|_| TryRecvError::Disconnected)
                } else {
                    chunk_rx.try_recv()
                };

                match next {
                    Ok(Ok(chunk)) => {
                        let (tx, rx) = mpsc::sync_channel(1);
                        job_tx
                            .send((chunk, tx))
                            .map_err(|_| anyhow!("the worker threads exited unexpectedly"))?;
                        pending.push_back(rx);
                    }
                    Ok(Err(error)) => {
                        read_error = Some(error);
                        reading = false;
                    }
                    Err(TryRecvError::Empty) => break,
                    Err(TryRecvError::Disconnected) => reading = false,
                }
            }

            let Some(rx) = pending.pop_front() else {
                break;
            };

            let (output, result) = rx
                .recv()
                .map_err(|_| anyhow!("a worker thread exited unexpectedly"))?;
            write(output)?;
            result?;
        }

        read_error.map_or(Ok(()), Err)
    })
}

/// Reads an input into chunks of whole lines.
///
/// Before each read that may block (for example, on standard input), the
/// complete lines read so far are sent, so records are processed as soon as
/// they are available. Comment lines are discarded as they are read.
///
/// A read error is sent after the chunk of lines read before it.
fn read(mut input: Input, kind: &str, comment: Option<u8>, tx: &SyncSender<Result<Chunk>>) {
    let mut chunk = Chunk::new(1);

    // The length of the complete lines in the chunk; any bytes after it are
    // the start of the next line.
    let mut complete = 0;

    // Whether the input is at the start of a line.
    let mut line_start = true;

    // Whether the current line is a comment that is being discarded.
    let mut skipping = false;

    loop {
        if input.buffer().is_empty() && chunk.lines > 0 {
            let next = chunk.split_off(complete);
            complete = 0;

            if tx.send(Ok(mem::replace(&mut chunk, next))).is_err() {
                return;
            }
        }

        let buffer = match input.fill_buf() {
            Ok(buffer) => buffer,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => {
                // NOTE: the bytes of a partially read line are discarded.
                chunk.data.truncate(complete);
                let record = chunk.first + chunk.lines as u64;

                if chunk.lines == 0 || tx.send(Ok(chunk)).is_ok() {
                    let _ =
                        tx.send(Err(e).with_context(|| format!("reading {kind} record {record}")));
                }

                return;
            }
        };

        if buffer.is_empty() {
            // NOTE: a final line without a line ending is still a line.
            if chunk.data.len() > complete {
                chunk.lines += 1;
            }

            if chunk.lines > 0 {
                let _ = tx.send(Ok(chunk));
            }

            return;
        }

        if line_start {
            skipping = comment.is_some() && buffer.first() == comment.as_ref();
        }

        let (length, line_end) = match memchr::memchr(b'\n', buffer) {
            Some(i) => (i + 1, true),
            None => (buffer.len(), false),
        };

        if !skipping {
            chunk.data.extend_from_slice(&buffer[..length]);
        }

        input.consume(length);
        line_start = line_end;

        if line_end && !skipping {
            chunk.lines += 1;
            complete = chunk.data.len();

            if chunk.lines >= CHUNK_LINES || complete >= CHUNK_BYTES {
                let first = chunk.first + chunk.lines as u64;

                if tx
                    .send(Ok(mem::replace(&mut chunk, Chunk::new(first))))
                    .is_err()
                {
                    return;
                }

                complete = 0;
            }
        }
    }
}

impl Chunk {
    /// Creates an empty chunk whose first record has the given number.
    fn new(first: u64) -> Self {
        Self {
            first,
            lines: 0,
            data: Vec::with_capacity(CHUNK_BYTES + CHUNK_BYTES / 4),
        }
    }

    /// Splits off the bytes after the complete lines into a new chunk for the
    /// following records.
    fn split_off(&mut self, complete: usize) -> Self {
        let mut next = Self::new(self.first + self.lines as u64);
        next.data.extend_from_slice(&self.data[complete..]);
        self.data.truncate(complete);
        next
    }
}

#[cfg(test)]
mod tests {
    use std::io;
    use std::io::Cursor;
    use std::io::Read;

    use super::*;

    /// Collects the record numbers and lines seen by the workers, in output
    /// order.
    fn run(input: impl Read + Send + 'static, workers: usize) -> (Vec<(u64, String)>, Result<()>) {
        let mut seen = Vec::new();

        let result = process(
            io::BufReader::new(Box::new(input)),
            "TEST",
            Some(b'#'),
            vec![(); workers],
            |_, chunk, output: &mut Vec<(u64, String)>| {
                for (i, line) in chunk.data.split_inclusive(|&b| b == b'\n').enumerate() {
                    let line = String::from_utf8(line.to_vec())?;

                    if line.starts_with("bad") {
                        anyhow::bail!("bad record {}", chunk.first + i as u64);
                    }

                    output.push((chunk.first + i as u64, line));
                }

                assert_eq!(output.len(), chunk.lines);
                Ok(())
            },
            |output| {
                seen.extend(output);
                Ok(())
            },
        );

        (seen, result)
    }

    fn lines(n: usize) -> String {
        (1..=n).map(|i| format!("{i}\n")).collect()
    }

    #[test]
    fn records_are_processed_in_order() {
        let n = CHUNK_LINES * 5 + 7;
        let (seen, result) = run(Cursor::new(lines(n)), 4);
        result.unwrap();

        assert_eq!(seen.len(), n);
        for (i, (record, line)) in seen.iter().enumerate() {
            assert_eq!(*record, i as u64 + 1);
            assert_eq!(*line, format!("{}\n", i + 1));
        }
    }

    #[test]
    fn a_final_line_without_a_newline_is_processed() {
        let (seen, result) = run(Cursor::new("1\n2"), 2);
        result.unwrap();
        assert_eq!(seen, [(1, String::from("1\n")), (2, String::from("2"))]);
    }

    #[test]
    fn comment_lines_are_skipped() {
        let (seen, result) = run(Cursor::new("#a\n1\n#b\n#c\n2\n#d"), 2);
        result.unwrap();
        assert_eq!(seen, [(1, String::from("1\n")), (2, String::from("2\n"))]);
    }

    #[test]
    fn empty_input_is_processed() {
        let (seen, result) = run(Cursor::new(""), 2);
        result.unwrap();
        assert!(seen.is_empty());
    }

    #[test]
    fn records_before_an_error_are_written() {
        let mut input = lines(CHUNK_LINES * 3);
        input.push_str("bad\n");
        input.push_str(&lines(CHUNK_LINES * 3));

        let (seen, result) = run(Cursor::new(input), 4);

        assert_eq!(
            result.unwrap_err().to_string(),
            format!("bad record {}", CHUNK_LINES * 3 + 1)
        );
        assert_eq!(seen.len(), CHUNK_LINES * 3);
    }

    /// A reader that fails after its data.
    struct Failing(Cursor<Vec<u8>>);

    impl Read for Failing {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            match self.0.read(buf)? {
                0 => Err(io::Error::other("failed")),
                n => Ok(n),
            }
        }
    }

    /// A reader that blocks after its data until its sender is dropped.
    struct Blocking(Cursor<Vec<u8>>, mpsc::Receiver<()>);

    impl Read for Blocking {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            match self.0.read(buf)? {
                0 => {
                    let _ = self.1.recv();
                    Ok(0)
                }
                n => Ok(n),
            }
        }
    }

    #[test]
    fn errors_are_reported_while_the_input_blocks() {
        let (unblock, blocked) = mpsc::channel();
        let input = Blocking(Cursor::new(b"1\nbad\npar".to_vec()), blocked);

        let (tx, rx) = mpsc::channel();
        thread::spawn(move || tx.send(run(input, 2)));

        let (seen, result) = rx
            .recv_timeout(std::time::Duration::from_secs(10))
            .expect("processing should not wait for the blocked input");
        drop(unblock);

        assert_eq!(result.unwrap_err().to_string(), "bad record 2");
        assert_eq!(seen, [(1, String::from("1\n"))]);
    }

    #[test]
    fn read_errors_name_the_record() {
        let (seen, result) = run(Failing(Cursor::new(b"1\n2\n3".to_vec())), 2);

        assert_eq!(result.unwrap_err().to_string(), "reading TEST record 3");
        assert_eq!(seen, [(1, String::from("1\n")), (2, String::from("2\n"))]);
    }
}
