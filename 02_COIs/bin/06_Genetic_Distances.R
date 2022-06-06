#!/usr/bin/Rscript
# 06_Genetic_Distancess.R
# Made by Luis Rodrigo Arce Valdés, to estimate genetic distances between COIs of hybridising species
rm(list = ls())

# Calling up libraries
library(ape)
library(dplyr)
library(tidyr)

# Genetic Distance Model (help: dist.dna):
model <- "raw"

# Reading FASTAs
COIs <- list()
consensus <- list()
# Listing all .fasta files per group
temp <- list.files(path = "../data/07_Paired_Muscle/", pattern = "*.fasta")
# Editing names
temp <- gsub(".fasta","",temp)
# Reading all .fasta files and assigning them to their name
for (i in temp) {
  print(i)
  COIs[[i]] <- read.FASTA(paste0("../data/07_Paired_Muscle/",i,".fasta"), type = "DNA")
  consensus[[i]] <- dist.dna(COIs[[i]], model = model, variance = F, as.matrix = T)
  consensus[[i]] <- consensus[[i]][1,2]
}
consensus <- as.data.frame(t.data.frame(as.data.frame(consensus)))
consensus$Cross <- row.names(consensus)
row.names(consensus) <- 1:nrow(consensus)
consensus <- consensus[,c(2,1)]
colnames(consensus)[2] <- "Distance"

# Writing results
write.table(consensus, "../results/01_Genetic_Distances.tsv", sep = "\t", quote = F, row.names = F)
