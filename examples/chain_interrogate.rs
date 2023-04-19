use std::env;
use std::fs::File;
use std::io::BufReader;

use chain::core::Coordinate;
use chain::core::Interval;
use chainfile as chain;
use flate2::read::GzDecoder;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let interval = env::args().nth(1).expect("missing interval");
    let src = env::args().nth(2).expect("missing src");

    let interval = interval.parse::<Interval>()?;

    let mut reader = File::open(src)
        .map(GzDecoder::new)
        .map(BufReader::new)
        .map(chain::Reader::new)?;

    for result in reader.sections() {
        let section = result?;
        let section_interval = Interval::try_new(
            Coordinate::try_new(
                section
                    .header()
                    .reference_sequence()
                    .chromosome_name()
                    .clone(),
                section.header().reference_sequence().alignment_start(),
                section.header().reference_sequence().strand().clone(),
            )?,
            Coordinate::try_new(
                section
                    .header()
                    .reference_sequence()
                    .chromosome_name()
                    .clone(),
                section.header().reference_sequence().alignment_end(),
                section.header().reference_sequence().strand().clone(),
            )?,
        )?;

        if section_interval.contains(interval.start()) {
            println!("{}", section.header());
        }
    }

    Ok(())
}
