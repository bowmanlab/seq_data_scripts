#!/bin/bash

#### run this script inside the conda env "imagine" ####

### testing ###

#name=CHABO_DNA_91

#TOMAP="CHABO_DNA_1 CHABO_DNA_13 CHABO_DNA_36 CHABO_DNA_45 CHABO_DNA_51 CHABO_DNA_58 CHABO_DNA_63 CHABO_DNA_69 CHABO_DNA_84 CHABO_DNA_89 CHABO_DNA_90 CHABO_DNA_91"

### running ###

name="$1"

TOMAP=""

while read f;do
	TOMAP="${TOMAP} ${f}"
done < names4coverage.txt &&

cd binning_${name} || exit 1

bwa index contigs.fasta &&

for fmap in ${TOMAP}; do
    bwa mem -t 32 contigs.fasta ../${fmap}_filt_R1.fq.gz ../${fmap}_filt_R2.fq.gz | \
        samtools sort -@ 16 -m 2G -o ${fmap}_map_sorted.bam -
    
    TOMAPSORT="${TOMAPSORT} ${fmap}_map_sorted.bam"
done

jgi_summarize_bam_contig_depths --outputDepth ${name}_depth.txt ${TOMAPSORT}
rm -f *.bam
