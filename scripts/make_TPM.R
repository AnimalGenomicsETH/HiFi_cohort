library(tximport)
library(tidyverse)
library(preprocessCore)
library(RNOmni)
library(argparser)
library(PCAForQTL)

args <- commandArgs(trailingOnly = TRUE)

name <- args[1]
gene_locations <- args[2]
dir <- args[3]
TPM_unfiltered_transcript <- args[4]
counts_unfiltered_transcript <- args[5]
TPM_unfiltered_gene <- args[6]
counts_unfiltered_gene <- args[7]
TPM_filtered_transcript <- args[8]
counts_filtered_transcript <- args[9]
TPM_filtered_gene <- args[10]
counts_filtered_gene <- args[11]
TPM_normalized_transcript <- args[12]
TPM_normalized_gene <- args[13]
TPM_filter <- as.numeric(args[14])
read_filter <- as.numeric(args[15])
samp_filter <- as.numeric(args[16])
PC_all_transcript <- args[17]
PC_cov_transcript <- args[18]
var_explain_transcript <-args[19]
PC_all_gene <- args[20]
PC_cov_gene <- args[21]
var_explain_gene <-args[22]

# Gene/transcript name list (from kallisto index)
name <- read.table(name, header = T)
names(name) <- c("pid", "gid", "gname")

# Get gene location from NCBI location table, format TSS for QTLtools
gene_locos <- read.table(gene_locations, header =T, sep = "\t")
gene_locos <- gene_locos[-c(3)]
names(gene_locos) <- c("gid", "#chr", "start", "end", "strand")
gene_locos <- gene_locos[!(gene_locos$`#chr` %in% c('Un', '')),]
gene_locos <- gene_locos %>% mutate(`#chr` = str_replace(`#chr`,"X,Y", "Y"))
gene_locos <- gene_locos %>% mutate(strand = str_replace(strand,"plus", "+"))
gene_locos <- gene_locos %>% mutate(strand = str_replace(strand, "minus", "-"))
neg <- gene_locos %>% filter(strand == "-")					# Separate by strand
pos <- gene_locos %>% filter(strand == "+")
neg$start <- neg$end - 1
pos$end <- pos$start
pos$start <- pos$start - 1
gene_locos <- rbind(neg, pos)						#Final table for merging 

# Make TPM matrix from kallisto output (transcripts)
dir <- dir
samples <- system(paste("ls -d ", dir, "/*/ | xargs -n 1 basename", sep=""), intern=TRUE)
files <- file.path(dir, c(samples), "abundance.tsv")
names(files) <- c(samples)
trans2gene <- name[, c(1,2)]
txi.kallisto_transcripts <- tximport(files, type = "kallisto", txOut = TRUE)
tpm_transcript <- as.data.frame(txi.kallisto_transcripts$abundance)
tpm_transcript <- cbind(rownames(tpm_transcript), tpm_transcript)
colnames(tpm_transcript)[1] <- "pid"
tpm_transcript <- merge(trans2gene, tpm_transcript, by = "pid", all = F)
tpm_transcript <- merge(gene_locos, tpm_transcript, by = "gid", all = F)
tpm_transcript <- tpm_transcript[c(2,3,4,6,1,5,7:as.numeric(length(tpm_transcript)))]
tpm_transcript <- tpm_transcript %>% distinct(pid, .keep_all = T)

# Make TPM matrix from kallisto output (gene)
txi.kallisto_gene <- tximport(files, type = "kallisto", txOut = FALSE, tx2gene=trans2gene)
tpm_gene <- as.data.frame(txi.kallisto_gene$abundance)
tpm_gene <- cbind(rownames(tpm_gene), tpm_gene)
colnames(tpm_gene)[1] <- "gid"
name_gene <- as.data.frame(unique(name[-c(1)]))
tpm_gene <- merge(gene_locos, tpm_gene, by = "gid", all = F)
tpm_gene <- tpm_gene[c(2,3,4,1,1,5,6,7:as.numeric(length(tpm_gene)))]
colnames(tpm_gene)[5] <- "gid"
colnames(tpm_gene)[4] <- "pid"
tpm_gene <- tpm_gene %>% distinct(gid, .keep_all = T)

# Make count matrix (same process)
counts_transcript <- as.data.frame(txi.kallisto_transcripts$counts)
counts_transcript <- cbind(rownames(counts_transcript), counts_transcript)
colnames(counts_transcript)[1] <- "pid"
counts_transcript <- merge(trans2gene, counts_transcript, by = "pid", all = F)
counts_transcript <- merge(gene_locos, counts_transcript, by = "gid", all = F)
counts_transcript <- counts_transcript[c(2,3,4,6,1,5,7:as.numeric(length(counts_transcript)))]
counts_transcript <- counts_transcript %>% distinct(pid, .keep_all = T)

# Make count matrix gene
counts_gene <- as.data.frame(txi.kallisto_gene$counts)
counts_gene <- cbind(rownames(counts_gene), counts_gene)
colnames(counts_gene)[1] <- "gid"
counts_gene <- merge(gene_locos, counts_gene, by = "gid", all = F)
counts_gene <- counts_gene[c(2,3,4,1,1,5,6,7:as.numeric(length(counts_gene)))]
colnames(counts_gene)[5] <- "gid"
colnames(counts_gene)[4] <- "pid"
counts_gene <- counts_gene %>% distinct(gid, .keep_all = T)

# Write unfiltered matrices 
write.table(tpm_transcript, TPM_unfiltered_transcript, col.names = T, row.names = F, quote = F, sep = "\t")
write.table(counts_transcript, counts_unfiltered_transcript, col.names = T, row.names = F, quote = F, sep = "\t")
write.table(tpm_gene, TPM_unfiltered_gene, col.names = T, row.names = F, quote = F, sep = "\t")
write.table(counts_gene, counts_unfiltered_gene, col.names = T, row.names = F, quote = F, sep = "\t")

# Filtering TPM matrix and count matrix
# Keep transcripts with > 0.1 TPM in > 10% of samples
tpm_transcript_filtered <- tpm_transcript[rowSums(tpm_transcript[,c(7:ncol(tpm_transcript))] >= TPM_filter) / (ncol(tpm_transcript)-6) >= samp_filter,]
tpm_gene_filtered <- tpm_gene[rowSums(tpm_gene[,c(7:ncol(tpm_gene))] >= TPM_filter) / (ncol(tpm_gene)-6) >= samp_filter,]

# Keep transcripts with > 6 reads in > 10% of samples
counts_transcript_good <- counts_transcript[rowSums(counts_transcript[,c(7:ncol(counts_transcript))] >= read_filter) / (ncol(counts_transcript)-6) >= samp_filter,]
tpm_transcript_filtered = tpm_transcript_filtered[tpm_transcript_filtered$pid %in% counts_transcript_good$pid,]
tpm_transcript_filtered <- tpm_transcript_filtered %>% arrange(`#chr`, start)
counts_gene_good <- counts_gene[rowSums(counts_gene[,c(7:ncol(counts_gene))] >= read_filter) / (ncol(counts_gene)-6) >= samp_filter,]
tpm_gene_filtered = tpm_gene_filtered[tpm_gene_filtered$gid %in% counts_gene_good$gid,]
tpm_gene_filtered <- tpm_gene_filtered %>% arrange(`#chr`, start)

# Subset counts matrix to match TPM matrix
counts_transcript_filtered = counts_transcript[counts_transcript$pid %in% tpm_transcript_filtered$pid,]
counts_transcript_filtered <- counts_transcript_filtered %>% arrange(`#chr`, start) 
counts_gene_filtered = tpm_gene[tpm_gene$gid %in% tpm_gene_filtered$gid,]
counts_gene_filtered <- counts_gene_filtered %>% arrange(`#chr`, start) 

# Write filtered matrices
write.table(tpm_transcript_filtered, TPM_filtered_transcript, col.names = T, row.names = F, quote = F, sep = "\t")
write.table(counts_transcript_filtered, counts_filtered_transcript, col.names = T, row.names = F, quote = F, sep = "\t") 
write.table(tpm_gene_filtered, TPM_filtered_gene, col.names = T, row.names = F, quote = F, sep = "\t")
write.table(counts_gene_filtered, counts_filtered_gene, col.names = T, row.names = F, quote = F, sep = "\t") 

# Transform and normalize TPM values
tpm_transcript_filtered_norm <- as.matrix(tpm_transcript_filtered[7:ncol(tpm_transcript_filtered)])  
tpm_transcript_filtered_norm <- normalize.quantiles(tpm_transcript_filtered_norm, copy = TRUE)
tpm_transcript_filtered_norm  <- as.data.frame(t(apply(tpm_transcript_filtered_norm, 1, RNOmni::RankNorm)))
colnames(tpm_transcript_filtered_norm) <- colnames(tpm_transcript_filtered[7:ncol(tpm_transcript_filtered)])
tpm_transcript_filtered_norm <- cbind(tpm_transcript_filtered[c(1:6)], tpm_transcript_filtered_norm)

# Transform and normalize gene matrix
tpm_gene_filtered_norm <- as.matrix(tpm_gene_filtered[7:ncol(tpm_gene_filtered)])  
tpm_gene_filtered_norm <- normalize.quantiles(tpm_gene_filtered_norm, copy = TRUE)
tpm_gene_filtered_norm  <- as.data.frame(t(apply(tpm_gene_filtered_norm, 1, RNOmni::RankNorm)))
colnames(tpm_gene_filtered_norm) <- colnames(tpm_gene_filtered[7:ncol(tpm_gene_filtered)])
tpm_gene_filtered_norm <- cbind(tpm_gene_filtered[c(1:6)], tpm_gene_filtered_norm)

# Write filtered and normalized matrix
write.table(tpm_transcript_filtered_norm, TPM_normalized_transcript, col.names = T, row.names = F, quote = F, sep = "\t") 
write.table(tpm_gene_filtered_norm, TPM_normalized_gene, col.names = T, row.names = F, quote = F, sep = "\t") 

# Make PCA
expr_t <- tpm_transcript_filtered_norm[-c(1:6)]
rownames(expr_t) <- tpm_transcript_filtered_norm[,4]
expr_pc_t <- prcomp(t(expr_t))
var_explained_t <- expr_pc_t$sdev^2/sum(expr_pc_t$sdev^2)
resultRunElbow_t <- as.numeric(PCAForQTL::runElbow(prcompResult=expr_pc_t))
PCs_all_t <- as.data.frame(expr_pc_t$x)
PCs_cov_t <- PCs_all_t[c(1:resultRunElbow_t)]
var_explained_t <- as.data.frame(var_explained_t)
row.names(var_explained_t) <- colnames(PCs_all_t)
write.table(PCs_all_t, PC_all_transcript, row.names = T, col.names = T, quote = F, sep = "\t")
write.table(PCs_cov_t, PC_cov_transcript, row.names = T, col.names = T, quote = F, sep = "\t")
write.table(var_explained_t, var_explain_transcript, row.names = T, col.names = F, quote = F, sep = "\t")

# Make PCA (Gene)
expr_g <- tpm_gene_filtered_norm[-c(1:6)]
rownames(expr_g) <- tpm_gene_filtered_norm[,4]
expr_pc_g <- prcomp(t(expr_g))
var_explained_g <- expr_pc_g$sdev^2/sum(expr_pc_g$sdev^2)
resultRunElbow_g <- as.numeric(PCAForQTL::runElbow(prcompResult=expr_pc_g))
PCs_all_g <- as.data.frame(expr_pc_g$x)
PCs_cov_g <- PCs_all_g[c(1:resultRunElbow_g)]
var_explained_g <- as.data.frame(var_explained_g)
row.names(var_explained_g) <- colnames(PCs_all_g)
write.table(PCs_all_g, PC_all_gene, row.names = T, col.names = T, quote = F, sep = "\t")
write.table(PCs_cov_g, PC_cov_gene, row.names = T, col.names = T, quote = F, sep = "\t")
write.table(var_explained_g, var_explain_gene, row.names = T, col.names = F, quote = F, sep = "\t")