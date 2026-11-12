#!/usr/bin/env bash
#
# Downloads the UCSC chain files used for comparisons into a directory,
# skipping any that are already present.

set -euo pipefail

directory="${1:?usage: $0 <directory>}"
mkdir -p "$directory"

for chain in hg19/liftOver/hg19ToHg38 hg38/liftOver/hg38ToHs1; do
  path="$directory/$(basename "$chain").over.chain.gz"

  if [[ ! -f "$path" ]]; then
    curl -fsSL -o "$path" "https://hgdownload.soe.ucsc.edu/goldenPath/$chain.over.chain.gz"
  fi
done
