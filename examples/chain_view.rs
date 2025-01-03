use std::env;
use std::fs::File;
use std::io::BufReader;

use chainfile as chain;
use flate2::read::GzDecoder;
use tabled::builder::Builder;
use tabled::settings::Alignment;
use tabled::settings::Style;
use tabled::settings::object::Rows;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let chain_ids: Vec<usize> = env::args()
        .skip(2)
        .map(|s| {
            s.parse::<usize>()
                .unwrap_or_else(|_| panic!("could not parse chain id: {s}"))
        })
        .collect();

    let chain_ids = if chain_ids.is_empty() {
        None
    } else {
        Some(chain_ids)
    };

    let mut reader = File::open(src)
        .map(GzDecoder::new)
        .map(BufReader::new)
        .map(chain::Reader::new)?;

    let mut builder = Builder::default();
    builder.push_record([
        "Reference",
        "--",
        "--",
        "-->",
        "Query",
        "--",
        "--",
        "-->",
        "Data",
        "--",
        "-->",
    ]);
    builder.push_record([
        "Contig", "Strand", "Start", "End", "Contig", "Strand", "Start", "End", "Size", "Dt", "Dq",
    ]);

    for result in reader.sections() {
        let section = result?;

        if let Some(ref chain_ids) = chain_ids {
            if !chain_ids.contains(&section.header().id()) {
                continue;
            }
        }

        for result in section.stepthrough_with_data()? {
            let (pair, data) = result?;

            builder.push_record([
                pair.reference().contig().as_str(),
                &pair.reference().strand().to_string(),
                &pair.reference().start().to_string(),
                &pair.reference().end().to_string(),
                pair.query().contig().as_str(),
                &pair.query().strand().to_string(),
                &pair.query().start().to_string(),
                &pair.query().end().to_string(),
                &data.size().to_string(),
                &data
                    .dt()
                    .map(|v| v.to_string())
                    .unwrap_or(String::from("<None>")),
                &data
                    .dq()
                    .map(|v| v.to_string())
                    .unwrap_or(String::from("<None>")),
            ]);
        }
    }

    let table = builder
        .build()
        .with(Style::rounded())
        .modify(Rows::new(1..), Alignment::left())
        .to_string();

    println!("{}", table);

    Ok(())
}
