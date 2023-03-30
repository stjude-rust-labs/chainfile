//! `chainfile` is a crate for reading a processing genomic chain files.

#![warn(missing_docs)]
#![warn(rust_2018_idioms)]
#![warn(rust_2021_compatibility)]

pub mod line;
pub mod reader;
pub mod record;

pub use self::reader::Reader;
