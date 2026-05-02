#!/bin/bash



mkdir bwa_mag_index
cp high_qual_draft/*fa bwa_mag_index/

for f in `ls bwa_mag_index/*fa`;do
	bwa index $f
done