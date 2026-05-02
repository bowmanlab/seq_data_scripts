#!/bin/bash

DIR=merged
ODIR=repaired

mkdir -p "$ODIR"

for R1 in ${DIR}/*R1*fastq.gz; do
    SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
    R2=${DIR}/${SAMPLE}_R2_001.fastq.gz
    repair.sh \
        in1="$R1" \
        in2="$R2" \
        out1=${ODIR}/${SAMPLE}_repaired_R1.fastq.gz \
        out2=${ODIR}/${SAMPLE}_repaired_R2.fastq.gz \
        outs=${ODIR}/${SAMPLE}_singletons.fastq.gz \
        repair=t
done
