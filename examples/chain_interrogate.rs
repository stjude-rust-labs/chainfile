use std::env;
use std::fs::File;
use std::io::BufReader;

use chainfile as chain;
use flate2::read::GzDecoder;
use omics::coordinate::interval::interbase::Interval;

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
        let section_interval = section.reference_sequence().interval()?;

        if section_interval.contains_coordinate(interval.start())
            || section_interval.contains_coordinate(interval.end())
        {
            println!("{}", section.header());
        }
    }

    Ok(())
}
