#Author: Alicia Sorgen
#Date: 10-12-21
#Description: Normalizes taxa count tables from Kraken2 output.

## Libraries
library(stringr)
library(data.table)

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
  relabTable <- getRelAbun(counts)
  
  normTable <- cbind(SampleID, normTable)
  relabTable <- cbind(SampleID, relabTable)
  
  write.table(normTable, paste0(outputDir,level,"_LogNormalizedCounts.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
  write.table(relabTable, paste0(outputDir,level,"_RelativeAbundanceCounts.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
  

}