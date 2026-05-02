#!/bin/bash

#### run this script inside the conda env "imagine" ####

### testing ###

#name=CHABO_DNA_91

### running ###

name="$1"

echo "metaspades.py -1 ${name}_filt_R1.fq.gz -2 ${name}_filt_R2.fq.gz -k 21,33,55 -o metaspades_output_${name} -t ${CPUS} -m ${MEMORY} --only-assembler"