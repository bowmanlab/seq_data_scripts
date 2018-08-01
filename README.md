# seq_data_scripts
Basic scripts for QC and analysis of Illumina MiSeq data

dada2.r - the script to execute the dada2 package in R, this will QC, denoise, and merge PE reads

deunique_dada2.py - this script will take a csv file as output by dada2.r and expand it to create a redundant fasta file
