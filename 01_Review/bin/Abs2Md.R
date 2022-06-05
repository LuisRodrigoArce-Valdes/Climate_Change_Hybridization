# Script to noquote out title, authors, years and abstract of all found manuscripts
rm(list = ls())
library(dplyr)

# Add search information
engine <- "Web of Science"
string <- paste0("See attached Excel file")
date <- "30/05/2022"

# Reading files
files <- list()
for (i in 1:3) {
files[[i]] <- read.table(paste0("../data/",i,".txt"), header = T, sep = "\t", fill = T, quote = "", stringsAsFactors = F, colClasses = "character")}
files <- bind_rows(files)

# Removing duplicates
files <- unique(files)

# Adding bold to key terms
files$AB <- gsub("hybrid","**hybrid**", files$AB)
files$AB <- gsub("Hybrid","**Hybrid**", files$AB)

files$AB <- gsub("introgression","**introgression**", files$AB)
files$AB <- gsub("Introgression","**Introgression**", files$AB)

# printing abstracts
sink("../results/Abstracts.md")

# Additional information
cat(noquote(paste0("###", engine, "\n\n")))
cat(noquote(paste0("*", string, "*\n\n")))
cat(noquote(paste0(date, "\n\n")))

# Printing abstracts
for (i in 1:nrow(files)) {
  cat(noquote(paste0("**", files[i,"TI"],"**\n\n")))
  cat(noquote(paste0(files[i, "AU"], ", ")))
  cat(noquote(paste0(files[i, "PY"], "\n\n")))
  cat(noquote(paste0(files[i, "AB"], "\n\n")))
  cat(noquote("---\n\n"))
}
sink()

# Convert to pdf using https://www.markdowntopdf.com/