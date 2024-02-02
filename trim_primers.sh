#!/bin/bash
# Author: Benjamin Klempay

# This script trims primers and adapters from the 3 prime end of amplicon
# sequence reads (fastq) using cutadapt. It is optimized for our standard
# ANL amplicon sequencing workflow, but it could easilty be adapted for use
# with different naming schemes or sequencing archetecture.

# Most of this functionality is present in the dada2.r script, but this script
# is included for standalone functionality.

#### Default parameters: ####
# can be changed manually or by running the script with
# the flags -i -o -t -f -r & -l (respectively)

indir=.
outdir=trimmed
threads=1
FWD=GTACACACCGCCCGTC		# 1391f and EukBr eukaryotic 18S primers
REV=TGATCCTTCTGCAGGTTCACCTAC	# used by ANL Environmental Sequencing Facility
minLen=50 # must be >= 1

# parse options from command line
while getopts 'i:o:t:f:r:l:' opt; do
        case $opt in
                i) indir=$OPTARG ;;
                o) outdir=$OPTARG ;;
		t) threads=$OPTARG ;;
		f) FWD=$OPTARG ;;
		r) REV=$OPTARG ;;
		l) minLen=$OPTARG ;;
        esac
done

FWD_RC=$(echo $FWD | tr "ATCGWSMKRYBDHVN" "TAGCWSKMYRVHDBN" | rev)
REV_RC=$(echo $REV | tr "ATCGWSMKRYBDHVN" "TAGCWSKMYRVHDBN" | rev)

mkdir -p $outdir

for i in $indir/*-R1.fastq; do
	name=$(basename $i "-R1.fastq")
	echo -e "\n\n#### Removing primers/adapters from $name ####"

	cutadapt -j $threads -a $REV_RC -A $FWD_RC --minimum-length $minLen \
	-o $outdir/$name"-R1.fastq" -p $outdir/$name"-R2.fastq" \
	$indir/$name"-R1.fastq" $indir/$name"-R2.fastq"
	# Add the flags '-g $FWD -G $REV' to remove primers/adapters from 5 prime end
	# (not required for our standard ANL amplicon sequencing workflow)
done
