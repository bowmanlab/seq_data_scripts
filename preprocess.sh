#!/bin/bash

## Fastq files from ANL may be labeled as .gz, but will not be gzipped.  Rename before running this script.
## Modify parse_anl_map_file.py as appropriate for file names.

#python parse_anl_map_file.py

## Make sure that the oligos file has no spaces in sample names (python script should do this in the future)

#mkdir demultiplexed

#iu-demultiplex -s 191021_Mapping_File_Bowman_18S.oligos --r1 Undetermined_S0_L001_R1_001.fastq --r2 Undetermined_S0_L001_R2_001.fastq --index Undetermined_S0_L001_I1_001.fastq -o demultiplexed/

## usually dada2.r gets run manually
#./dada2.r

ls *csv|parallel ./deunique_dada2.py {}
