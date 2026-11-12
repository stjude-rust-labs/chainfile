//! Swapping the reference and alternate alleles of a biallelic VCF record.
//!
//! When a record's alternate allele matches the query genome, the alleles are
//! swapped so that the record's reference allele matches the query genome.
//! Fields that refer to alleles by their index are rewritten to match:
//!
//! * `GT` allele indices `0` and `1` are exchanged.
//! * `Number=R` fields (such as `AD`) have their two values exchanged.
//! * FORMAT `Number=G` fields (such as `PL`) are reversed. For a biallelic
//!   record, genotypes are ordered by their number of alternate alleles, so
//!   reversing the values matches the swapped alleles for any ploidy.
//! * The INFO `AF` field becomes `1 - AF`, and the INFO `AC` field becomes `AN
//!   - AC`.
//!
//! Other `Number=A` fields, local allele fields (`LAA`, `LA`, `LR`, and `LG`),
//! INFO `Number=G` fields, and fields known to describe specific alleles (such
//! as `ANN`, `CSQ`, `SB`, and `PGT`) cannot be rewritten in general, so they
//! are dropped. All other fields (including `QUAL` and `FILTER`) are kept
//! unchanged.

use noodles::vcf::Header;
use noodles::vcf::header::record::value::Map;
use noodles::vcf::header::record::value::map::Format;
use noodles::vcf::header::record::value::map::Info;
use noodles::vcf::header::record::value::map::format::Number as FormatNumber;
use noodles::vcf::header::record::value::map::info::Number as InfoNumber;
use noodles::vcf::variant::record_buf::Info as InfoFields;
use noodles::vcf::variant::record_buf::RecordBuf;
use noodles::vcf::variant::record_buf::Samples;
use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
use noodles::vcf::variant::record_buf::info::field::value::Array as InfoArray;
use noodles::vcf::variant::record_buf::samples::Keys;
use noodles::vcf::variant::record_buf::samples::sample::Value as SampleValue;
use noodles::vcf::variant::record_buf::samples::sample::value::Array as SampleArray;
use noodles::vcf::variant::record_buf::samples::sample::value::Genotype;

/// The genotype FORMAT field.
const GENOTYPE_KEY: &str = "GT";

/// The INFO field for allele frequencies.
const ALLELE_FREQUENCY_KEY: &str = "AF";

/// The INFO field for allele counts.
const ALLELE_COUNT_KEY: &str = "AC";

/// The INFO field for the total number of alleles.
const TOTAL_ALLELE_COUNT_KEY: &str = "AN";

/// INFO fields that describe specific alleles without being declared as
/// per-allele fields.
const ALLELE_SPECIFIC_INFO_KEYS: &[&str] = &["ANN", "BCSQ", "CSQ"];

/// FORMAT fields that describe specific alleles without being declared as
/// per-allele fields.
const ALLELE_SPECIFIC_FORMAT_KEYS: &[&str] = &["LAA", "PGT", "SB"];

/// How a FORMAT field is rewritten.
#[derive(Clone, Copy, Debug)]
enum FormatAction {
    /// The field is kept unchanged.
    Keep,

    /// The field is dropped.
    Drop,

    /// The field is the genotype.
    Genotype,

    /// The field has one value per allele.
    PerAllele,

    /// The field has one value per genotype.
    PerGenotype,
}

/// Swaps the reference and alternate alleles in a biallelic record's fields.
///
/// The record's alleles themselves are not changed. Fields are looked up in
/// `header`, falling back to the definitions reserved by the VCF
/// specification.
///
/// Returns `false` without changing the record when its fields are not valid
/// for a biallelic record (for example, a `GT` that refers to a third allele).
pub fn swap_alleles(header: &Header, record: &mut RecordBuf) -> bool {
    let Some(info) = swap_info(header, record.info()) else {
        return false;
    };

    let Some(samples) = swap_samples(header, record.samples()) else {
        return false;
    };

    *record.info_mut() = info;
    *record.samples_mut() = samples;

    true
}

/// Gets the number of values of an INFO field (if it is defined).
fn info_number(header: &Header, key: &str) -> Option<InfoNumber> {
    match header.infos().get(key) {
        Some(info) => Some(info.number()),
        None => {
            let info = Map::<Info>::from((header.file_format(), key));
            (!info.description().is_empty()).then(|| info.number())
        }
    }
}

/// Gets the number of values of a FORMAT field (if it is defined).
fn format_number(header: &Header, key: &str) -> Option<FormatNumber> {
    match header.formats().get(key) {
        Some(format) => Some(format.number()),
        None => {
            let format = Map::<Format>::from((header.file_format(), key));
            (!format.description().is_empty()).then(|| format.number())
        }
    }
}

/// Swaps the alleles in the INFO fields.
///
/// Returns `None` when a per-allele field does not have two values.
fn swap_info(header: &Header, info: &InfoFields) -> Option<InfoFields> {
    let total_allele_count = match info.get(TOTAL_ALLELE_COUNT_KEY) {
        Some(Some(InfoValue::Integer(n))) if *n >= 0 => Some(*n),
        _ => None,
    };

    let mut swapped = InfoFields::default();

    for (key, value) in info.as_ref() {
        if ALLELE_SPECIFIC_INFO_KEYS.contains(&key.as_str()) {
            continue;
        }

        let value = match info_number(header, key) {
            None | Some(InfoNumber::Count(_) | InfoNumber::Unknown) => value.clone(),
            Some(InfoNumber::ReferenceAlternateBases) => match value {
                Some(InfoValue::Array(array)) if info_array_len(array) == 2 => {
                    let mut array = array.clone();
                    reverse_info_array(&mut array);
                    Some(InfoValue::Array(array))
                }
                Some(_) => return None,
                None => None,
            },
            Some(InfoNumber::AlternateBases) => {
                let value = match key.as_str() {
                    ALLELE_FREQUENCY_KEY => swap_allele_frequency(value.as_ref()),
                    ALLELE_COUNT_KEY => swap_allele_count(value.as_ref(), total_allele_count),
                    _ => None,
                };

                match value {
                    Some(value) => value,
                    None => continue,
                }
            }
            Some(InfoNumber::Samples) => continue,
        };

        swapped.insert(key.clone(), value);
    }

    Some(swapped)
}

/// Swaps an INFO `AF` value.
///
/// Returns `None` when the value cannot be swapped, in which case the field is
/// dropped.
fn swap_allele_frequency(value: Option<&InfoValue>) -> Option<Option<InfoValue>> {
    let Some(value) = value else {
        return Some(None);
    };

    let InfoValue::Array(InfoArray::Float(values)) = value else {
        return None;
    };

    let [frequency] = values.as_slice() else {
        return None;
    };

    let frequency = match frequency {
        Some(frequency) if frequency.is_finite() && (0.0..=1.0).contains(frequency) => {
            // NOTE: the complement is computed with double precision so that
            // it rounds to the nearest single-precision value (e.g., `1 - 0.3`
            // is written as `0.7` rather than `0.70000005`).
            Some((1.0 - f64::from(*frequency)) as f32)
        }
        Some(_) => return None,
        None => None,
    };

    Some(Some(InfoValue::Array(InfoArray::Float(vec![frequency]))))
}

/// Swaps an INFO `AC` value using the total number of alleles (`AN`).
///
/// Returns `None` when the value cannot be swapped, in which case the field is
/// dropped.
fn swap_allele_count(
    value: Option<&InfoValue>,
    total_allele_count: Option<i32>,
) -> Option<Option<InfoValue>> {
    let Some(value) = value else {
        return Some(None);
    };

    let InfoValue::Array(InfoArray::Integer(values)) = value else {
        return None;
    };

    let [count] = values.as_slice() else {
        return None;
    };

    let count = match count {
        Some(count) => {
            let total = total_allele_count?;

            if !(0..=total).contains(count) {
                return None;
            }

            Some(total - count)
        }
        None => None,
    };

    Some(Some(InfoValue::Array(InfoArray::Integer(vec![count]))))
}

/// Swaps the alleles in the FORMAT fields of each sample.
///
/// Returns `None` when a genotype refers to an allele other than the first two
/// or a per-allele or per-genotype field has the wrong number of values.
fn swap_samples(header: &Header, samples: &Samples) -> Option<Samples> {
    let (keys, values): (Keys, Vec<Vec<Option<SampleValue>>>) = samples.clone().into();

    let actions = keys
        .as_ref()
        .iter()
        .map(|key| format_action(header, key))
        .collect::<Vec<_>>();

    let swapped_keys = keys
        .as_ref()
        .iter()
        .zip(&actions)
        .filter(|(_, action)| !matches!(action, FormatAction::Drop))
        .map(|(key, _)| key.clone())
        .collect::<Keys>();

    let genotype_index = keys.as_ref().get_index_of(GENOTYPE_KEY);

    let mut swapped_values = Vec::with_capacity(values.len());

    for sample in values {
        let ploidy = genotype_index
            .and_then(|i| sample.get(i))
            .and_then(|value| match value {
                Some(SampleValue::Genotype(genotype)) => Some(genotype.as_ref().len()),
                _ => None,
            });

        let mut swapped = Vec::with_capacity(sample.len());

        for (value, action) in sample.into_iter().zip(&actions) {
            let value = match (action, value) {
                (FormatAction::Drop, _) => continue,
                (FormatAction::Keep, value) | (_, value @ None) => value,
                (FormatAction::Genotype, Some(SampleValue::Genotype(genotype))) => {
                    Some(SampleValue::Genotype(swap_genotype(genotype)?))
                }
                (FormatAction::PerAllele, Some(SampleValue::Array(mut array)))
                    if sample_array_len(&array) == 2 =>
                {
                    reverse_sample_array(&mut array);
                    Some(SampleValue::Array(array))
                }
                (FormatAction::PerGenotype, Some(SampleValue::Array(mut array)))
                    if ploidy.is_none_or(|ploidy| sample_array_len(&array) == ploidy + 1) =>
                {
                    reverse_sample_array(&mut array);
                    Some(SampleValue::Array(array))
                }
                (_, Some(_)) => return None,
            };

            swapped.push(value);
        }

        swapped_values.push(swapped);
    }

    Some(Samples::new(swapped_keys, swapped_values))
}

/// Gets how a FORMAT field is rewritten.
fn format_action(header: &Header, key: &str) -> FormatAction {
    if key == GENOTYPE_KEY {
        return FormatAction::Genotype;
    }

    if ALLELE_SPECIFIC_FORMAT_KEYS.contains(&key) {
        return FormatAction::Drop;
    }

    match format_number(header, key) {
        None
        | Some(
            FormatNumber::Count(_)
            | FormatNumber::Unknown
            | FormatNumber::Ploidy
            | FormatNumber::BaseModifications,
        ) => FormatAction::Keep,
        Some(FormatNumber::ReferenceAlternateBases) => FormatAction::PerAllele,
        Some(FormatNumber::Samples) => FormatAction::PerGenotype,
        Some(
            FormatNumber::AlternateBases
            | FormatNumber::LocalAlternateBases
            | FormatNumber::LocalReferenceAlternateBases
            | FormatNumber::LocalSamples,
        ) => FormatAction::Drop,
    }
}

/// Exchanges allele indices `0` and `1` in a genotype.
///
/// Returns `None` when the genotype refers to any other allele.
fn swap_genotype(mut genotype: Genotype) -> Option<Genotype> {
    for allele in genotype.as_mut() {
        let position = match allele.position() {
            Some(0) => Some(1),
            Some(1) => Some(0),
            Some(_) => return None,
            None => None,
        };

        *allele.position_mut() = position;
    }

    Some(genotype)
}

/// Gets the number of values in an INFO array.
fn info_array_len(array: &InfoArray) -> usize {
    match array {
        InfoArray::Integer(values) => values.len(),
        InfoArray::Float(values) => values.len(),
        InfoArray::Character(values) => values.len(),
        InfoArray::String(values) => values.len(),
    }
}

/// Reverses the values in an INFO array.
fn reverse_info_array(array: &mut InfoArray) {
    match array {
        InfoArray::Integer(values) => values.reverse(),
        InfoArray::Float(values) => values.reverse(),
        InfoArray::Character(values) => values.reverse(),
        InfoArray::String(values) => values.reverse(),
    }
}

/// Gets the number of values in a FORMAT array.
fn sample_array_len(array: &SampleArray) -> usize {
    match array {
        SampleArray::Integer(values) => values.len(),
        SampleArray::Float(values) => values.len(),
        SampleArray::Character(values) => values.len(),
        SampleArray::String(values) => values.len(),
    }
}

/// Reverses the values in a FORMAT array.
fn reverse_sample_array(array: &mut SampleArray) {
    match array {
        SampleArray::Integer(values) => values.reverse(),
        SampleArray::Float(values) => values.reverse(),
        SampleArray::Character(values) => values.reverse(),
        SampleArray::String(values) => values.reverse(),
    }
}

#[cfg(test)]
mod tests {
    use noodles::vcf;
    use noodles::vcf::variant::io::Write as _;

    use super::*;

    /// The header lines of the test VCF. `AD` is intentionally undeclared so
    /// that its reserved definition is used.
    const HEADER: &[&str] = &[
        "##fileformat=VCFv4.3",
        r#"##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">"#,
        r#"##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">"#,
        r#"##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">"#,
        r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">"#,
        r#"##INFO=<ID=XA,Number=A,Type=Integer,Description="Per-ALT">"#,
        r#"##INFO=<ID=XR,Number=R,Type=Integer,Description="Per-allele">"#,
        r#"##INFO=<ID=ANN,Number=.,Type=String,Description="Annotations">"#,
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        r#"##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">"#,
        r#"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">"#,
        r#"##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Strand bias">"#,
        r#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">"#,
    ];

    /// Swaps a biallelic `A`/`C` record with the given INFO field and sample
    /// columns (the FORMAT keys followed by one column per sample).
    ///
    /// Returns the swapped INFO and sample columns or `None` if the record
    /// could not be swapped.
    fn swap(info: &str, samples: &[&str]) -> Option<Vec<String>> {
        let mut names = [
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        ]
        .map(String::from)
        .to_vec();
        let mut record = vec!["chr1", "1", ".", "A", "C", ".", ".", info];

        if let Some((_, sample_columns)) = samples.split_first() {
            names.push(String::from("FORMAT"));
            names.extend((1..=sample_columns.len()).map(|i| format!("sample{i}")));
            record.extend(samples);
        }

        let mut lines = HEADER.to_vec();
        let (names, record) = (names.join("\t"), record.join("\t"));
        lines.extend([names.as_str(), record.as_str()]);
        let text = lines.join("\n") + "\n";

        let mut reader = vcf::io::Reader::new(text.as_bytes());
        let header = reader.read_header().unwrap();
        let mut record = reader.record_bufs(&header).next().unwrap().unwrap();

        if !swap_alleles(&header, &mut record) {
            return None;
        }

        let mut writer = vcf::io::Writer::new(Vec::new());
        writer.write_variant_record(&header, &record).unwrap();

        let line = String::from_utf8(writer.into_inner()).unwrap();
        Some(
            line.trim_end()
                .split('\t')
                .skip(7)
                .map(String::from)
                .collect(),
        )
    }

    #[test]
    fn genotypes() {
        assert_eq!(
            swap(".", &["GT", "0/0", "0|1", "1", "./.", "0/.", "1|1"]).unwrap(),
            [".", "GT", "1/1", "1|0", "0", "./.", "1/.", "0|0"]
        );

        assert_eq!(swap(".", &["GT", "0/2"]), None);
    }

    #[test]
    fn per_allele_and_per_genotype_fields() {
        assert_eq!(
            swap(
                ".",
                &[
                    "GT:AD:PL:GQ",
                    "0/1:3,5:20,0,100:20",
                    "0/0/1:7,2:0,10,20,30:9",
                    "./.:.:.:.",
                ]
            )
            .unwrap(),
            [
                ".",
                "GT:AD:PL:GQ",
                "1/0:5,3:100,0,20:20",
                "1/1/0:2,7:30,20,10,0:9",
                "./.:.:.:.",
            ]
        );

        assert_eq!(swap(".", &["GT:AD", "0/1:1,2,3"]), None);
        assert_eq!(swap(".", &["GT:PL", "0/1:0,10"]), None);
    }

    #[test]
    fn dropped_format_fields() {
        assert_eq!(
            swap(".", &["GT:AF:SB:GQ", "0/1:0.4:1,2,3,4:30"]).unwrap(),
            [".", "GT:GQ", "1/0:30"]
        );
    }

    #[test]
    fn info_fields() {
        assert_eq!(
            swap("AF=0.3;AC=1;AN=4;DP=10;XA=7;XR=1,2;ANN=C|missense", &[]).unwrap(),
            ["AF=0.7;AC=3;AN=4;DP=10;XR=2,1"]
        );

        assert_eq!(swap("AF=1.5;AC=1", &[]).unwrap(), ["."]);
        assert_eq!(swap("AC=5;AN=4", &[]).unwrap(), ["AN=4"]);
        assert_eq!(swap("XR=1,2,3", &[]), None);
    }
}
