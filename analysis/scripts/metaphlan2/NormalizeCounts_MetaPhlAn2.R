#Author: Alicia Sorgen
#Date: 10-12-21
#Description: Normalizes taxa count tables from Kraken2 output.

R <- sessionInfo()
message(R$R.version$version.string)

## Libraries
library(stringr); message("stringr:", packageVersion("stringr"))
library(data.table); message("data.table:", packageVersion("data.table"))

rm(list=ls())

##### Prep #####
params <- commandArgs(trailingOnly = TRUE)
inputDir = paste0(params[1], "/"); message("inputDir = ", inputDir)
outputDir <- paste0(params[2], "/"); message("outputDir = ", outputDir)

funcScript <- params[3]
source(funcScript)

levels <- c("phylum", "class", "order", "family", "genus", "species")

for (level in levels) {
  
  countFile <- paste0(level, "_rawCounts.tsv")
  counts <- read.delim(paste0(inputDir, countFile), sep="\t",header = TRUE)
  rownames(counts) <- counts$SampleID
  counts <- counts[,-1]
  
  SampleID <- rownames(counts)
  
  normTable <- getNormTable(counts)
  
  normTable <- cbind(SampleID, normTable)
  
  write.table(normTable, paste0(outputDir,level,"_LogNormalizedCounts.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
  

}
