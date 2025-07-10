library(tidyverse)
library(argparser, quietly=TRUE)
args <- commandArgs(trailingOnly = TRUE)

rin_age_file <- args[1]
rna_pca <- args[2]
geno_pcs <-args[3]
output <- args[4]

rin_age <- read.table(rin_age_file, header = T, row.names = 1)
rin_age <- rin_age[c("rin_epi_h", "age")]
rna_pcs <- read.table(rna_pca, header = T, row.names = 1)
geno_pcs  <- read.table(geno_pcs, header = T, comment.char = "", row.names = 1)
#geno_pcs <- geno_pcs[c(1:5)]
colnames(geno_pcs)[1:10]<-paste0("genotypePC",1:10)

all_covs <- merge(rin_age, geno_pcs, by = 0, all = F)
row.names(all_covs) <- all_covs$Row.names
all_covs <- all_covs[-c(1)]
keep <- row.names(rna_pcs)
all_covs <- all_covs[(row.names(all_covs) %in% keep),]

knownCovariatesFiltered<-PCAForQTL::filterKnownCovariates(all_covs,rna_pcs,unadjustedR2_cutoff=0.5)

all_covs <- merge(rna_pcs, knownCovariatesFiltered, by = 0, all = F)
colnames(all_covs)[1] <- "ID"
all_covs <- as.data.frame(t(all_covs))
write.table(all_covs, output, col.names = F, row.names = T, quote = F, sep = "\t")
