# 02_COIs_Table.R
# With this script we will create a table with each BOLD and GeneBank record key for each species
rm(list = ls())
library(tidyr)
library(stringr)

# COIs meta information ####
# BOLD
BOLD <- read.table("../data/08_Meta/Bold.txt", header = F, sep = ">")
BOLD <- data.frame(BOLD[,-1])

# Removing duplicates
BOLD <- unique(BOLD)

# Splitting columns
BOLD <- separate(BOLD,BOLD....1.,c("BOLD","Sp","Gen","GenBank"), sep = "\\|")
BOLD <- BOLD[,c(2,1,4)]
BOLD$Sp <- gsub("_"," ",BOLD$Sp)

# GenBank
GenBank <- read.table("../data/08_Meta/GB.txt", header = F, sep = ">")
GenBank <- data.frame(GenBank[,-1])

# Removing duplicates
GenBank <- unique(GenBank)

# Splitting columns
GenBank <- separate(GenBank,GenBank....1.,as.character(1:5), sep = "\\|")
GenBank <- GenBank[,c(5,4)]
GenBank$`5` <- gsub("^ ","", GenBank$`5`)
GenBank <- separate(GenBank,`5`, c("Gen","sp"), sep = " ", extra = "drop")
GenBank$sp <- paste(GenBank$Gen, GenBank$sp)
GenBank$BOLD <- "NA"
GenBank <- GenBank[,c(2,4,3)]
colnames(GenBank) <- colnames(BOLD)

# Binding data.frames
BOLD <- rbind(BOLD, GenBank)
rm(GenBank)

# Ordering
BOLD <- BOLD[order(BOLD$Sp),]

# Exporting
write.table(BOLD, "../results/COIs_Meta_Info.tsv", quote = F, row.names = F, sep = "\t")
