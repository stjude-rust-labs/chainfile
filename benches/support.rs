use std::fmt::Write;

/// Generates synthetic chain file data in memory.
///
/// Each contig gets `blocks_per_contig` aligned blocks of `block_size` bases,
/// separated by gaps of `gap_size` in both reference and query.
pub fn generate_chain_data(
    num_contigs: usize,
    blocks_per_contig: usize,
    block_size: u64,
    gap_size: u64,
) -> Vec<u8> {
    let mut buf = String::new();

    for contig_idx in 0..num_contigs {
        let ref_name = format!("chr{}", contig_idx + 1);
        let query_name = format!("qchr{}", contig_idx + 1);

        let total_ref_size =
            blocks_per_contig as u64 * block_size + (blocks_per_contig as u64 - 1) * gap_size;
        let total_query_size = total_ref_size;

        writeln!(
            buf,
            "chain 100 {ref_name} {total_ref_size} + 0 {total_ref_size} {query_name} \
             {total_query_size} + 0 {total_query_size} {contig_idx}"
        )
        .unwrap();

        for block_idx in 0..blocks_per_contig {
            if block_idx < blocks_per_contig - 1 {
                writeln!(buf, "{block_size}\t{gap_size}\t{gap_size}").unwrap();
            } else {
                writeln!(buf, "{block_size}").unwrap();
            }
        }

        if contig_idx < num_contigs - 1 {
            writeln!(buf).unwrap();
        }
    }

    buf.into_bytes()
}
