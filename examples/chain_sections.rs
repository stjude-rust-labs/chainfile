use std::env;
use std::fs::File;
use std::io::BufReader;

use chainfile as chain;
use flate2::read::GzDecoder;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(GzDecoder::new)
        .map(BufReader::new)
        .map(chain::Reader::new)?;

    for result in reader.sections() {
        let section = result?;

        // Print all alignment data headers where the chromosome has changed
        if section.header().reference_sequence().chromosome_name()
            != section.header().query_sequence().chromosome_name()
        {
            println!("{:?}", section.header());
        }
    }

    Ok(())
}
