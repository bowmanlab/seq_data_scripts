#!/bin/bash

## Fastq files from ANL may be labeled as .gz, but will not be gzipped.  Rename before running this script.
## Modify parse_anl_map_file.py as appropriate for file names.

#python parse_anl_map_file.py

## Make sure that the oligos file has no spaces in sample names (python script should do this in the future)

mkdir demultiplexed

MAPFILE=230326_Bowman_16s_DG_230317.oligos

iu-demultiplex -s $MAPFILE --r1 Undetermined_S0_L001_R1_001.fastq --r2 Undetermined_S0_L001_R2_001.fastq --index Undetermined_S0_L001_I1_001.fastq -o demultiplexed/

## usually dada2.r gets run manually
#./dada2.r

seqtable2fasta.py seqtable_16S.csv
