//! `chainfile` is a crate for reading a processing genomic chain files.
//!
//! The crate provides two main points of entry:
//!
//! - Parsing and reading chain files directly.
//! - Providing a machine for lifting over intervals given a chain file.
//!
//! Since the main purpose of a chain file is to lift over intervals from one
//! genome build to another, we expect that most users will be interested in the
//! latter functionality. However, we have exposed the former functionality in
//! the event that it is needed for some other purpose.
//!
//! ## Parsing and reading chain files
//!
//! If you're interested in parsing and reading chain files directly, you can
//! use the [`reader::Reader`] facility to accomplish that. Most users will want
//! to read the parsed [`reader::AlignmentDataSections`] using
//! [`reader::Reader::sections()`]. For each data section, you can access the
//! [`record::HeaderRecord`] (via [`reader::AlignmentDataSection::header()`])
//! and the subsequent [`record::AlignmentDataRecord`]s for that section (via
//! [`reader::AlignmentDataSection::alignment_data_records()`]). However, most
//! users will not be interested in working with the raw alignment data records.
//!
//! Generally, what one _actually_ wants is the mapping between contiguous
//! regions of the reference and query genomes that are defined by the alignment
//! data section. The translation between a raw alignment data records and this
//! mapping can be tricky, especially considering gotchas such as coordinates on
//! the reverse strand being stored as the reverse complement of the sequence.
//! Instead of computing these yourself, you should use the
//! [`liftover::StepThrough`] facility that can be obtained from each alignment
//! data section via [`reader::AlignmentDataSection::stepthrough()`].
//!
//! Iterating over this stepthrough provides a series of
//! [`liftover::stepthrough::interval_pair::ContiguousIntervalPair`]s that
//! represent contiguous mappings between the two genomes. This struct includes
//! the ever-important
//! [`liftover::stepthrough::interval_pair::ContiguousIntervalPair::liftover()`]
//! method to translate a [`core::Coordinate`] from the reference
//! [`core::Interval`] to the query [`core::Interval`].
//!
//! Below is a representative example of how you might want to access and
//! explore a chain file with the facilities discussed above.
//!
//! ```
//! use chainfile as chain;
//!
//! let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
//! let mut reader = chain::Reader::new(&data[..]);
//!
//! for result in reader.sections() {
//!     let section = result?;
//!     println!("{}", section.header());
//!
//!     for result in section.stepthrough()? {
//!         let pair = result?;
//!         println!("{} -> {}", pair.reference(), pair.query());
//!     }
//! }
//!
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! ## Liftover Machine
//!
//! Most often, users will be interested in the [`liftover::Machine`] facility.
//! This struct provides the capability to attempt a lift over of a
//! [`core::Interval`] from the reference genome to the query genome via
//! [`liftover::Machine::liftover()`]. Perhaps importantly (and different from
//! most other liftover tools that the author is aware of), this method provides
//! the complete list of mapped contiguous interval pairs that are encompassed
//! by the provided interval rather than providing an inexact mapping and/or
//! lifting over a single position. Note that, if you're interested in lifting
//! over a single position, you can still construct a 1-sized [`core::Interval`]
//! that contains your position and lift that over instead.
//!
//! A [`liftover::Machine`] cannot be instantiated directly. Instead, you should
//! use [`liftover::machine::Builder`] and the associated
//! [`liftover::machine::Builder::try_build_from()`] method to construct a
//! liftover machine.
//!
//! Below is a representative example of how one might read in a chain file,
//! construct a liftover machine, parse an interval of interest, and then lift
//! over that interval of interest from the reference genome to the query
//! genome.
//!
//! ```
//! use chainfile as chain;
//!
//! let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
//! let mut reader = chain::Reader::new(&data[..]);
//! let machine = chain::liftover::machine::Builder::default().try_build_from(reader)?;
//!
//! let interval = "seq0:3-4".parse::<chain::core::Interval>()?;
//! for result in machine.liftover(&interval).unwrap() {
//!     println!("{} -> {}", result.reference(), result.query());
//! }
//!
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

#![warn(missing_docs)]
#![warn(rust_2018_idioms)]
#![warn(rust_2021_compatibility)]

pub mod core;
pub mod liftover;
pub mod line;
pub mod reader;
pub mod record;

pub use self::reader::Reader;
