#!/bin/bash

./drep_combined.sh

python3 select_genomes.py

## gtdbk

rm -rf high_qual_identify_output
rm -rf high_qual_classify_output
rm -rf high_qual_align_output
rm -rf high_qual_ani_rep

./gtdbk_combined.sh

python3 aggregate_genome_info.py

## eggnog mapper

./eggnog_combined.sh

python3 make_annotation_table.py 

## map 

python3 make_combined_bins_fasta.py

bwa index high_qual_draft.fasta

REF=high_qual_draft.fasta

rm -rf mag_map
mkdir mag_map

while read f;do
	f_NAME=`basename $f`
	R1=${f}_R1_001.fastq.gz
	R2=${f}_R2_001.fastq.gz
	
	bwa mem -t 64 ${REF} ${R1} ${R2} | samtools sort -l 0 | samtools view -F 260 -q 20 | gzip > mag_map/${f_NAME}_map.q20.sam.gz
	
done < raw_reads_combined.txt

## count reads in original sample files

python3 count_reads.py

## return reads per million and breadth

python3 tally_sam_RPM.py

rm -rf map_dets
ls mag_map/*gz|parallel -j 24 --ungroup --joblog tally_sam_FC.log python3 tally_sam_FC.py {}
parallel -j 6 --ungroup --retry-failed --joblog tally_sam_FC.log python3 tally_sam_FC.py {}

python_3 tally_sam_FC_combine.py

#### cluster ####

conda activate imagine

rm -rf clustering

mkdir clustering

cat high_qual_draft_prodigal/*faa > clustering/high_qual_draft.faa

mmseqs createdb clustering/high_qual_draft.faa clustering/high_qual_draft.db

mmseqs cluster clustering/high_qual_draft.db clustering/high_qual_draft.clu tmp

mmseqs createsubdb clustering/high_qual_draft.clu clustering/high_qual_draft.db clustering/high_qual_draft.clu.rep

mmseqs createtsv clustering/high_qual_draft.db clustering/high_qual_draft.db clustering/high_qual_draft.clu clustering/high_qual_draft.clu.tsv

mmseqs convert2fasta clustering/high_qual_draft.clu.rep clustering/high_qual_draft.clu.rep.fasta

#### structural prediction - nimrod ####

python3 predict.py --fasta high_qual_draft.clu.rep.fasta --output-dir 2024_WATL_combined_clu_rep

## back on fram/gjoa - requires predicted_topologies.3line from nimrod ##

python3 parse_deeptmhmm.py


