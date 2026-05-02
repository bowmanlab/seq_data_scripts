#!/bin/bash

CPUS=16
SAMPLES_FILE="$1"

mkdir -p mag_map

map_bins(){
    local REF="$1"
    local R1="$2"
    local R2="$3"
    local SAMPLE_NAME="$4"
    local CPUS="$5"
    local TOTAL_READS="$6"
    
    local ref_name=$(basename "$REF" .fa)
    local base_name="mag_map/${SAMPLE_NAME}__${ref_name}"

    # Check if input files exist
    if [[ ! -f "$R1" ]]; then
        echo "Warning: $R1 not found, skipping ${SAMPLE_NAME} vs ${ref_name}"
        return 1
    fi
    
    if [[ ! -f "$R2" ]]; then
        echo "Warning: $R2 not found, skipping ${SAMPLE_NAME} vs ${ref_name}"
        return 1
    fi
    
    # Check if BWA index exists
    if [[ ! -f "${REF}.bwt" ]]; then
        echo "Warning: BWA index not found for $REF, skipping"
        return 1
    fi
    
    echo "Mapping $SAMPLE_NAME to $ref_name..."
    
    # Index the reference if not already done
    if [[ ! -f "${REF}.fai" ]]; then
        samtools faidx ${REF}
    fi
	
	# Map reads, count mapped reads, calculate coverage stats
	
    REF_LENGTH=$(awk '{sum+=$2} END{print sum}' ${REF}.fai)

    local mapped_count_file="${base_name}_mapped.count"
    local stats_file="${base_name}_stats.txt"

    bwa mem -t "$CPUS" "$REF" "$R1" "$R2" | \
        samtools view -F 260 -q 20 -b - | \
        tee >(samtools view -c - > "$mapped_count_file") | \
        samtools sort - | \
        samtools depth -a - | \
        awk -v ref_len="$REF_LENGTH" \
            -v sample="$SAMPLE_NAME" \
            -v ref="$ref_name" \
            -v total_reads="$TOTAL_READS" \
            -v count_file="$mapped_count_file" '
        {
            sum += $3
            if ($3 > 0) covered++
        }
        END {
            mean_cov = sum / ref_len
            frac_covered = covered / ref_len
            
            # Read mapped count from file
            getline mapped_reads < count_file
            
            # Calculate RPM
            rpm = mapped_reads / total_reads * 1000000

            print "sample\t" sample
            print "reference\t" ref
            print "ref_length\t" ref_len
            print "mapped_reads\t" mapped_reads
            print "total_reads\t" total_reads
            print "rpm\t" rpm
            print "mean_coverage\t" mean_cov
            print "fraction_covered\t" frac_covered
        }' > "$stats_file"
    
    if [[ $? -eq 0 ]]; then
        echo "Completed: $stats_file"
        return 0
    else
        echo "Error: Failed mapping $SAMPLE_NAME to $ref_name"
        return 1
    fi
}
 
export -f map_bins

while read -r col1 col2 col3; do
    SAMPLE_NAME="$col3"
    R1="${SAMPLE_NAME}_filt_R1.fq.gz"
    R2="${SAMPLE_NAME}_filt_R2.fq.gz"

    # Calculate total reads from R1 file prior to mapping
    echo "Counting reads in $R1..."
    TOTAL_READS=$(seqtk seq "$R1" | wc -l | awk '{print $1/4}')
    echo "$SAMPLE_NAME: $TOTAL_READS total reads"
    
    ls bwa_mag_index/*.fa | parallel -j 4 map_bins {} "$R1" "$R2" "$SAMPLE_NAME" "$CPUS" "$TOTAL_READS"
    
done < "$SAMPLES_FILE"

# Aggregate all stats into a single summary file
echo "Aggregating stats..."
echo -e "sample\treference\tref_length\tmapped_reads\ttotal_reads\trpm\tmean_coverage\tfraction_covered" > mag_map/summary_stats.txt

for stats_file in mag_map/*_stats.txt; do
    awk 'NR==1{sample=$2} NR==2{ref=$2} NR==3{len=$2} NR==4{mapped=$2} NR==5{total=$2} NR==6{rpm=$2} NR==7{cov=$2} NR==8{frac=$2}
    END{print sample"\t"ref"\t"len"\t"mapped"\t"total"\t"rpm"\t"cov"\t"frac}' "$stats_file"
done >> mag_map/summary_stats.txt

echo "=== Mapping complete ==="
echo "Summary stats written to mag_map/summary_stats.txt"