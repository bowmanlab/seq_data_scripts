#!/bin/bash

#### run this script inside the conda env "imagine" ####

### testing ###

#name=CHABO_DNA_91

### running ###

name="$1"
MEMORY=512
CPUS=64

##reformat.sh in1=${name}_filt_R1.fq.gz in2=${name}_filt_R2.fq.gz out1=${name}_filt_sub_R1.fq.gz out2=${name}_filt_sub_R2.fq.gz sample=50000000 

### --only-assembler flag used because error correction seems to hang for many samples ###

if [ ! -d "metaspades_output_${name}" ]; then
    metaspades.py -1 ${name}_filt_R1.fq.gz -2 ${name}_filt_R2.fq.gz -k 21,33,55 -o metaspades_output_${name} -t ${CPUS} -m ${MEMORY} --only-assembler
else
    metaspades.py -o metaspades_output_${name} --restart-from last -m ${MEMORY}
fi &&

quast.py metaspades_output_${name}/contigs.fasta -o quast_output_${name} &&

mkdir binning_${name} &&

seqtk seq -L 1500 metaspades_output_${name}/contigs.fasta > binning_${name}/contigs.fasta


