use std::env;
use std::fs::File;
use std::io::BufReader;

use chain::Line;
use chainfile as chain;
use flate2::read::GzDecoder;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(GzDecoder::new)
        .map(BufReader::new)
        .map(chain::Reader::new)?;

    let mut sections = 0;
    let mut alignment_data_entries = 0;

    for result in reader.lines() {
        let line = result?;

        match line {
            Line::Empty => {}
            Line::Header(_) => sections += 1,
            Line::AlignmentData(_) => alignment_data_entries += 1,
        }
    }

    println!(
        "Counted {sections} alignment sections and {alignment_data_entries} alignment data entries."
    );

    Ok(())
}
