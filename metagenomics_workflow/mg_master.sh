#!/bin/bash

## consider whether samples need to be merged or repaired prior to QC

export DIR=/data_store/seq_data/igm_20260118/repaired

## samples.txt is a tab delimited files of read 1, read 2, sample name

## The QC and assembly steps can and should be split up by creating separate samples.txt files for fram and gjoa, and executing the workflows simulataneously in the working directory. 
## Subsampling is also handled here.

export NREADS=50000000

parallel -u --colsep '\t' \
    'fastp --reads_to_process '"$NREADS"' \
    -i '"$DIR"'/{1} \
    -I '"$DIR"'/{2} \
    -o {3}_filt_R1.fq.gz \
    -O {3}_filt_R2.fq.gz \
    -y -e 30 -w 1 \
    -j fastp_{3}.json \
    -h fastp_{3}.html' :::: samples.txt

## mg_check_complete.py looks for the final directory created by mg_assemble.sh.  If it doesn't find it, assumes the
## script needs to be run and adds to samples_resume.txt.

python3 mg_check_complete.py samples.txt

parallel -u ./mg_assemble.sh {} :::: samples_resume.txt

## mg_check_complete.py should be run manually until all assemblies complete without error (or investigate individual cases)

python3 mg_check_complete.py samples_resume.txt

### mapping, make sure that names4coverage.txt is present with the samples that you want to use for mapping ###

while read -r col1 col2 col3; do
    ./mg_map.sh "$col3"
done < samples.txt

### binning ###

while read -r col1 col2 col3;do
	name="$col3"
	metabat2 -i binning_${name}/contigs.fasta -a binning_${name}/${name}_depth.txt -o bins_dir/${name}_bin --seed 1234 -m 1500 -v 
done < samples.txt

## checkm2 is difficult to install and requires mamba.  If installed differently, change. 

mamba run -n checkm2 checkm2 predict -i bins_dir/ -o checkm2_output/ -x fa --force --threads 48

awk 'BEGIN{FS="\t"; OFS=","} NR==1{print "genome","completeness","contamination","strain_heterogeneity","length","N50"} NR>1{print $1".fa",$2,$3,"0",$9,$8}' checkm2_output/quality_report.tsv > checkm2_output/drep_compatible_genomeInfo.csv

### dreplicate, high quality bins only ###

rm -rf drep_output

dRep dereplicate drep_output \
  -comp 50 \
  -con 10 \
  -sa 0.95 \
  -g bins_dir/*.fa \
  -p 64 \
  --genomeInfo checkm2_output/drep_compatible_genomeInfo.csv
  
python3 mg_aggregate_files.py
  
### annotate ###

./mg_eggnog.sh

python3 mg_disaggregate_eggnog.py

## classify

rm -rf high_qual_identify_output
rm -rf high_qual_classify_output
rm -rf high_qual_align_output
rm -rf high_qual_ani_rep

conda run -n gtdbtk ./mg_gtdbk.sh

## decorate drep_compatible_genomeInfo.csv

python3 mg_collect_genome_info.py

### map reads and calculate RPMs (make sure that global $DIR variable is set!) ###

./mg_make_mag_index.sh

## to make workflow more consistent, mg_map_mags.sh loop should be moved outside to master script
## this is a major bottleneck, would be good place to use fram and gjoa.  Note that this does not
## overwrite summary_stats.txt, as the last server running will create that fille for all MAGs.

./mg_map_mags.sh samples.txt

## count reads - probably not necessary

(echo "sample,read_count"; ls *_filt_R1.fq.gz | awk 'BEGIN{OFS=","} 
{
    sample_name = $1
    gsub(/_filt_R1\.fq\.gz$/, "", sample_name)
    
    cmd = "seqtk seq " $1 " | wc -l"
    cmd | getline line_count
    close(cmd)
    read_count = line_count / 4
    print sample_name, read_count
}') > read_counts.csv

## sadly samtools has faster tools for this

samtools depth "${sample}_vs_$(basename $bin .fa).bam" | \
            awk '{sum+=$3; count++} END{print sum/count}' > coverage.txt
			
# Compare positions with coverage to genome length
samtools depth file.bam | wc -l    # Positions covered
samtools faidx reference.fa        # Shows genome length

ls mag_map *sam.gz|python3 tally_sam_FC.py {}

python3 tally_sam_FC_combine.py

## next steps, convert sample mapping files to bam and assess above commands.
## If they work for bam files, restructure workflow

bwa mem -t 16 ${REF} ${R1} ${R2} | \
    samtools view -F 260 -q 20 -b - | \
    tee >(samtools view -c - > mapped_reads.count) | \
    samtools sort - | \
    samtools depth -a - | \
    awk '{sum+=$3; count++} END{print sum/count}' > mean_coverage.txt

# Calculate RPM
MAPPED_READS=$(cat mapped_reads.count)
TOTAL_READS=$(zcat ${R1} | wc -l | awk '{print $1/4}')  # If you know total reads
RPM=$(awk "BEGIN {printf \"%.2f\", $MAPPED_READS / $TOTAL_READS * 1000000}")
REF_LENGTH=$(grep -v '^>' ${REF} | wc -c) 
