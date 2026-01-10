


library(dada2)
packageVersion('dada2')
path <- '/mnt'
setwd('/mnt')
list.files(path)
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

pdf("quality_profile_forward.pdf", width = 10, height = 6)
plotQualityProfile(fnFs[1:2])
dev.off()

pdf("quality_profile_backward.pdf", width = 10, height = 6)
plotQualityProfile(fnRs[1:2])
dev.off()

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,260),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)  
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
write.table(seqtab.nochim,'df.txt',sep='\t',col.names=NA)

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE)
write.table(taxa,'df_taxa.txt',sep='\t',col.names=NA)
