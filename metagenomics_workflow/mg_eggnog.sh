#!/bin/bash

export EGGNOG_DATA_DIR=/home/jsbowman/eggnog-mapper-data

INPUT_DIR="/data_store/shared_work/botswana_chabo/high_qual_draft"
OUTPUT_DIR="/data_store/shared_work/botswana_chabo/high_qual_draft_eggnog"

rm -rf "$OUTPUT_DIR"

mkdir -p "$OUTPUT_DIR"

cat "$INPUT_DIR"/*.faa > "$OUTPUT_DIR"/combined.faa

CPUS=48

emapper.py \
    -i "$OUTPUT_DIR"/combined.faa \
    -o "$OUTPUT_DIR"/combined \
    --itype proteins \
	-m diamond \
	--evalue 0.001 \
	--score 60 \
	--seed_ortholog_score 60 \
	--dbmem \
	--cpu $CPUS \
	--output_dir "$OUTPUT_DIR"
