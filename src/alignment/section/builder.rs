//! A builder for an alignment section.

use nonempty::NonEmpty;

use crate::alignment::Section;
use crate::alignment::section::data;
use crate::alignment::section::header;

/// An error that occurs when a required field was never provided to the
/// [`Builder`].
#[derive(Debug)]
pub enum MissingError {
    /// No data was provided to the [`Builder`].
    Data,

    /// No header was provided to the [`Builder`].
    Header,
}

impl std::fmt::Display for MissingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MissingError::Data => write!(f, "data"),
            MissingError::Header => write!(f, "header"),
        }
    }
}

impl std::error::Error for MissingError {}

/// An error that occurs when a singular field was provided multiple times to
/// the [`Builder`].
#[derive(Debug)]
pub enum MultipleError {
    /// The header field was provided multiple times to the [`Builder`].
    Header,
}

impl std::fmt::Display for MultipleError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MultipleError::Header => write!(f, "header"),
        }
    }
}

impl std::error::Error for MultipleError {}

/// An error related to a [`Builder`].
#[derive(Debug)]
pub enum Error {
    /// An error where a required field was never provided to the [`Builder`].
    Missing(MissingError),

    /// An error where a singular field was provided to the [`Builder`] more
    /// than once.
    Multiple(MultipleError),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Missing(err) => write!(f, "missing required field: {err}"),
            Error::Multiple(err) => write!(f, "singular field set multiple times: {err}"),
        }
    }
}

impl std::error::Error for Error {}

/// A [`Result`](std::result::Result) with an [`Error`].
type Result<T> = std::result::Result<T, Error>;

/// A builder for a [`Section`].
#[derive(Debug, Default)]
pub struct Builder {
    /// The header record.
    header: Option<header::Record>,

    /// The data records.
    data: Option<NonEmpty<data::Record>>,
}

impl Builder {
    /// Sets the header record for the [`Builder`].
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::Builder;
    /// use chainfile::alignment::section::header::Record;
    ///
    /// let header = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse::<Record>()?;
    /// let builder = Builder::default().header(header)?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn header(mut self, record: header::Record) -> Result<Self> {
        if self.header.is_some() {
            return Err(Error::Multiple(MultipleError::Header));
        }

        self.header = Some(record);
        Ok(self)
    }

    /// Pushes a [data record](crate::alignment::section::data::Record) into the
    /// [`Builder`].
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::Builder;
    /// use chainfile::alignment::section::data::Record;
    ///
    /// let data = "10".parse::<Record>()?;
    /// let builder = Builder::default().push_data(data);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn push_data(mut self, record: data::Record) -> Self {
        let data = match self.data {
            Some(mut data) => {
                data.push(record);
                data
            }
            None => NonEmpty::new(record),
        };

        self.data = Some(data);
        self
    }

    /// Consumes `self` to attempt to build a [`Section`].
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::Builder;
    /// use chainfile::alignment::section::data;
    /// use chainfile::alignment::section::header;
    ///
    /// let mut section = Builder::default()
    ///     .header("chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse()?)?
    ///     .push_data("3\t0\t1".parse()?)
    ///     .push_data("1".parse()?)
    ///     .try_build()?;
    ///
    /// assert_eq!(section.header().score(), 0);
    /// assert_eq!(section.header().id(), 1);
    /// assert_eq!(section.data().len(), 2);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_build(self) -> Result<Section> {
        let header = self.header.ok_or(Error::Missing(MissingError::Header))?;

        let alignment_data = self.data.ok_or(Error::Missing(MissingError::Data))?;

        Ok(Section {
            header,
            data: alignment_data,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_fails_to_produce_a_section_when_no_header_is_provided()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = Builder::default()
            .push_data("1".parse()?)
            .try_build()
            .unwrap_err();

        assert_eq!(err.to_string(), "missing required field: header");

        Ok(())
    }

    #[test]

    fn it_fails_to_produce_a_section_when_the_header_field_is_provided_more_than_once()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = Builder::default()
            .header("chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse()?)?
            .header("chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse()?)
            .unwrap_err();

        assert_eq!(err.to_string(), "singular field set multiple times: header");

        Ok(())
    }

    #[test]
    fn it_fails_to_produce_a_section_when_no_data_is_provided()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = Builder::default()
            .header("chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse()?)?
            .try_build()
            .unwrap_err();

        assert_eq!(err.to_string(), "missing required field: data");

        Ok(())
    }
}
