#!/usr/bin/Rscript

# See tutorial at https://benjjneb.github.io/dada2/tutorial.html.  This script will
# QC and identify valid unique reads, then assemble.  It will execute on all files in
# the directory "multiplexed", and create the directories "filtered", and "merged".
# The Python script deunique_dada2.py should be used to inflate all of the output tables
# into redundant fasta files for paprica.

library(dada2)

path <- 'demultiplexed_18S'
gene <- '18S'

#fnFs <- sort(list.files(path, pattern = 'Mock.*-R1.fastq', full.names = T))
#fnRs <- sort(list.files(path, pattern = 'Mock.*-R2.fastq', full.names = T))

fnFs <- sort(list.files(path, pattern = '-R1.fastq', full.names = T))
fnRs <- sort(list.files(path, pattern = '-R2.fastq', full.names = T))

sample.names <- sapply(strsplit(basename(fnFs), "-R"), `[`, 1)

pdf(paste0(path, '/', 'quality_profiles.pdf'), width = 6, height = 6)

for(i in 1:length(fnFs)){
	print(plotQualityProfile(fnFs[i]))
	print(plotQualityProfile(fnRs[i]))
}
	
dev.off()

file_path <- file.path(paste0(path, '/','filtered')) # Place filtered files in filtered/ subdirectory
filtFs <- file.path(file_path, paste0(sample.names, "_R1.filt.fastq.gz"))
filtRs <- file.path(file_path, paste0(sample.names, "_R2.filt.fastq.gz"))

## multithreading only useful if multiple fastq files

out <- filterAndTrim(fnFs,
	filtFs,
	fnRs,
	filtRs,
	multithread = T,
	trimLeft = 15,
	truncLen = 150,
	verbose = T)

plotQualityProfile(filtFs[50])

## need distribution of lengths for filtFs and filtRs
	
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

pdf(paste0(path, '/', 'error_rates.pdf'), width = 6, height = 6)
print(plotErrors(errF, nominalQ = T))
print(plotErrors(errR, nominalQ = T))
dev.off()

## redefine names in case some files didn't pass QC

filtFs <- sort(list.files(paste0(path, '/filtered/'), pattern = '_R1.filt.fastq.gz', full.names = T))
filtRs <- sort(list.files(paste0(path, '/filtered/'), pattern = '_R2.filt.fastq.gz', full.names = T))

sample.names <- sapply(strsplit(basename(filtFs), "_R1"), `[`, 1)

derepFs <- derepFastq(filtFs, verbose = T)
derepRs <- derepFastq(filtRs, verbose = T)

## wouldn't it make more sense to assemble here?

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs,
                      derepFs,
                      dadaRs,
                      derepRs,
                      maxMismatch = 0,
                      verbose=TRUE)
					  
seqtab <- makeSequenceTable(mergers)

## remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    method = "consensus",
                                    multithread = TRUE,
                                    verbose = TRUE)
					  
## remove ASVs associated with mitos and chloros for 16S only
#!!! You should save discarded seqs in a different seqtable in case they're needed

if(gene == '16S'){
  
  ## classify
  
  taxa <- assignTaxonomy(seqtab.nochim,
                         '/data_store/silva_databases/silva_nr99_v138.1_train_set.fa.gz',
                         multithread = T)

  chloros <- row.names(taxa)[grep('Chloroplast', taxa[,'Order'])]
  mitos <- row.names(taxa)[grep('Mitochondria', taxa[,'Family'])]
  
  seqtab.discard <- seqtab.nochim[,colnames(seqtab.nochim) %in% chloros|
				   colnames(seqtab.nochim) %in% mitos]

  seqtab.reduced <- seqtab.nochim[,!colnames(seqtab.nochim) %in% chloros]
  seqtab.reduced <- seqtab.reduced[,!colnames(seqtab.reduced) %in% mitos]
  
  write.csv(t(seqtab.reduced), paste0('seqtab.reduced_', gene, '.csv'), quote = F)
  write.csv(taxa, paste0('taxa_', gene, '.csv'), quote = F)
  write.csv(t(seqtab.discard), paste0('seqtab.discard_', gene, '.csv'), quote = F)
}

## ZymoBIOMICS Microbial Community Standard should have 8 bacterial strains, 2 fungal strains,
## so if you have more than 10 strains you probably QC'd insufficiently.

write.csv(t(seqtab.nochim), paste0('seqtable_', gene, '.csv'), quote = F)

