#!/bin/bash

export GTDBTK_DATA_PATH=/data_store/gtdbtk_files/release232
PWD=`pwd`

CPUS=48

## note the need for full links or symlinks created by gtdbk will fail

gtdbtk identify --genes --genome_dir  ${PWD}/high_qual_draft/ --out_dir ${PWD}/high_qual_identify_output --cpus ${CPUS} --extension faa
gtdbtk align --identify_dir ${PWD}/high_qual_identify_output --out_dir ${PWD}/high_qual_align_output --cpus ${CPUS}
gtdbtk classify --genome_dir  ${PWD}/high_qual_draft/ --align_dir ${PWD}/high_qual_align_output/ --out_dir ${PWD}/high_qual_classify_output --cpus ${CPUS} --extension fa
gtdbtk ani_rep --genome_dir  ${PWD}/high_qual_draft/ --out_dir high_qual_ani_rep/ --cpus ${CPUS} --extension fa