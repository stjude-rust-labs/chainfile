use std::env;
use std::fs::File;
use std::io::BufReader;

use chain::liftover;
use chainfile as chain;
use flate2::read::GzDecoder;
use omics::coordinate::interval::interbase::Interval;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let interval = env::args().nth(1).expect("missing interval");
    let src = env::args().nth(2).expect("missing src");

    let interval = interval.parse::<Interval>()?;

    let reader = File::open(src)
        .map(GzDecoder::new)
        .map(BufReader::new)
        .map(chain::Reader::new)?;

    let machine = liftover::machine::builder::Builder.try_build_from(reader)?;
    let results = machine.liftover(interval);

    match results {
        Some(results) => {
            for result in results {
                println!("{}", result);
            }
        }
        None => println!("Does not exist in new genome"),
    }

    Ok(())
}
