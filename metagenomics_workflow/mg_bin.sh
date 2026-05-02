#!/bin/bash

#### run this script inside the conda env "imagine" ####

### testing ###

#name=CHABO_DNA_91

### running ###

name=$3

metabat2 -i contigs.fasta -a ${name}_depth.txt -o bins_dir/${name}_bin --seed 1234 -m 1500 -v &&

checkm data setRoot /home/jsbowman/checkm_files &&

checkm lineage_wf bins_dir/ checkm/ -x .fa -t ${CPUS} &&

mag_extract.py

echo $name "binning complete"