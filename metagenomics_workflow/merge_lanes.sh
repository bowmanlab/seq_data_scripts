#!/usr/bin/env bash
# merge_lanes.sh
# Concatenates L001 and L002 FASTQ.gz files for each sample/read direction.
# Usage: bash merge_lanes.sh <input_dir> <output_dir>

set -euo pipefail

INPUT_DIR="./"
OUTPUT_DIR="./merged"

mkdir -p "$OUTPUT_DIR"

# Find all L001 files and pair with corresponding L002
while IFS= read -r l001; do
    l002="${l001//_L001_/_L002_}"

    if [[ ! -f "$l002" ]]; then
        echo "WARNING: No L002 counterpart found for $(basename "$l001"), skipping." >&2
        continue
    fi

    # Build output filename by removing the lane tag entirely
    # e.g. CHABO_DNA_13_S35_L001_R1_001.fastq.gz -> CHABO_DNA_13_S35_R1_001.fastq.gz
    basename_out=$(basename "$l001" | sed 's/_L001//')
    out="$OUTPUT_DIR/$basename_out"

    echo "Merging:"
    echo "  $(basename "$l001")"
    echo "  $(basename "$l002")"
    echo "  -> $basename_out"

    cat "$l001" "$l002" > "$out"

done < <(find "$INPUT_DIR" -maxdepth 1 -name "*_L001_*.fastq.gz" | sort)

echo ""
echo "Done. Merged files written to: $OUTPUT_DIR"