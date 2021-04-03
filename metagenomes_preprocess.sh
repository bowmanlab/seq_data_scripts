#!/bin/bash

## This script is for demultiplexing and merging PE reads from the Illumina HiSeq platform for genomes or metagenomes.
## It assumes that your barcode file has two columns: <sample_name>/t<barcode>.  If your barcode has 10 characters but
## the sequencing center only provides 8 this method still works.  This method also considers only the forward read index
## primer.

## Demultiplex.

R1=bowman_S1_L001_R1_001.fastq.gz
R2=bowman_S1_L001_R2_001.fastq.gz
I1=bowman_S1_L001_I1_001.fastq.gz
MAP=ucsd_core_2021331_mapping_F.oligos

mkdir demultiplexed

iu-demultiplex -s $MAP --r1 $R1 --r2 $R2 --index $I1 -o demultiplexed/

mkdir merged

## Get a list of all sample names.

for f in merged/*-R1.fastq; do
	printf '%s\n' "${f%-R1.fastq}" >> samples.txt
done

## Merge reads.

while read f;do
	pear -f demultiplexed/${f}-R1.fastq -r demultiplexed/${f}-R2.fastq -o merged/$f -j 72 -n 100 > demultiplexed/${f}.pear.txt
done < samples.txt

