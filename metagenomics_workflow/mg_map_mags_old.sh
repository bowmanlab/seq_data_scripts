#!/bin/bash

CPUS=16

mkdir -p mag_map

map_bins(){
    local REF="$1"
    local R1="$2"
    local R2="$3"
    local SAMPLE_NAME="$4"
    local CPUS="$5"
    
    local ref_name=$(basename "$REF" .fa)
    local output_file="mag_map/${SAMPLE_NAME}__${ref_name}_map.q20.sam.gz"
    
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
    
    # bwa mem -t "$CPUS" "$REF" "$R1" "$R2" | \
        # samtools view -F 260 -q 20 | \
        # gzip > "$output_file"
		
	# First index the reference if not already done
	samtools faidx ${REF}
	REF_LENGTH=$(awk '{sum+=$2} END{print sum}' ${REF}.fai)

	bwa mem -t 16 ${REF} ${R1} ${R2} | \
		samtools view -F 260 -q 20 -b - | \
		tee >(samtools view -c - > mapped_reads.count) | \
		samtools sort - | \
		samtools depth -a - | \
		awk -v ref_len="$REF_LENGTH" '
		{
			sum += $3
			if ($3 > 0) covered++
		}
		END {
			mean_cov = sum / ref_len
			frac_covered = covered / ref_len
			print "mean_coverage\t" mean_cov
			print "fraction_covered\t" frac_covered
		}' > coverage_stats.txt
    
    if [[ $? -eq 0 ]]; then
        echo "Completed: $output_file"
        return 0
    else
        echo "Error: Failed mapping $SAMPLE_NAME to $ref_name"
        return 1
    fi
}
 
export -f map_bins

while read -r col1 col2 col3; do
	SAMPLE_NAME="$col3"
	R1=${SAMPLE_NAME}_filt_R1.fq.gz
	R2=${SAMPLE_NAME}_filt_R2.fq.gz
	SAMPLE_NAME="$col3"
	
	ls bwa_mag_index/*fa|parallel -j 4 map_bins {} "$R1" "$R2" "$SAMPLE_NAME" "$CPUS"
	
done < "$1"