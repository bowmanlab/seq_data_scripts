# seq_data_scripts
Basic scripts for QC and analysis of Illumina MiSeq data

preprocess.sh - basic steps for running all the scripts necessary to take Illumina MiSeq data from ANL through all preprocessing and QC steps.  I don't recommend running the script, just use it as a guide.

dada2.r - the script to execute the dada2 package in R, this will QC, denoise, and merge PE reads.  The trim parameters should be used as is unless there is a clear need to change them.

deunique_dada2.py - this script will take a csv file as output by dada2.r and expand it to create a redundant fasta file.  This redundant fasta file is suitable for analysis with paprica or another tool.
