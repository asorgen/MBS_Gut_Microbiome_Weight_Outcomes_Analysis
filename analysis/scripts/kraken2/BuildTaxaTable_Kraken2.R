#Author: Alicia Sorgen
#Date: 10-11-21
#Description: Builds raw taxa count tables from Kraken2 output.

R <- sessionInfo()
message(R$R.version$version.string)

## Libraries
library(data.table); message("data.table:", packageVersion("data.table"))
library(stringr); message("stringr:", packageVersion("stringr"))

rm(list=ls())

##### Prep #####
params <- commandArgs(trailingOnly = TRUE)
inputDir = paste0(params[1], "/"); message("inputDir = ", inputDir)
outputDir <- paste0(params[2], "/"); message("outputDir = ", outputDir)
dir.create(outputDir, showWarnings = FALSE)

funcScript <- params[3]
source(funcScript)

# files <- list.files(inputDir, pattern = "1_AB")
files <- list.files(inputDir, pattern = "k2_REPORT.txt")

index <- 1

for (file in files) {
 
  file.path <- paste0(inputDir, file) 
  sampleID <- sapply(strsplit(file, "_k2"), `[`, 1)
  table <- read.delim(file.path, sep = "\t", header = FALSE)
  
  phylum <- getphylumTableKraken2(table)
  class <- getclassTableKraken2(table)
  order <- getorderTableKraken2(table)
  family <- getfamilyTableKraken2(table)
  genus <- getgenusTableKraken2(table)
  species <- getspeciesTableKraken2(table)
  
  if (index == 1) {
    phylum.df <- phylum
    class.df <- class
    order.df <- order
    family.df <- family
    genus.df <- genus
    species.df <- species
  } else {
    phylum.df <- merge(phylum.df, phylum, by = "V1", all = TRUE)
    class.df <- merge(class.df, class, by = "V1", all = TRUE)
    order.df <- merge(order.df, order, by = "V1", all = TRUE)
    family.df <- merge(family.df, family, by = "V1", all = TRUE)
    genus.df <- merge(genus.df, genus, by = "V1", all = TRUE)
    species.df <- merge(species.df, species, by = "V1", all = TRUE)
  }
  
  index <- index + 1
}

df <- parseKraken2Taxa(phylum.df)
write.table(df, paste0(outputDir, "phylum_rawCounts.tsv"),sep="\t",quote = FALSE, row.names = FALSE)

df <- parseKraken2Taxa(class.df)
write.table(df, paste0(outputDir, "class_rawCounts.tsv"),sep="\t",quote = FALSE, row.names = FALSE)

df <- parseKraken2Taxa(order.df)
write.table(df, paste0(outputDir, "order_rawCounts.tsv"),sep="\t",quote = FALSE, row.names = FALSE)

df <- parseKraken2Taxa(family.df)
write.table(df, paste0(outputDir, "family_rawCounts.tsv"),sep="\t",quote = FALSE, row.names = FALSE)

df <- parseKraken2Taxa(genus.df)
write.table(df, paste0(outputDir, "genus_rawCounts.tsv"),sep="\t",quote = FALSE, row.names = FALSE)

df <- parseKraken2Taxa(species.df)
write.table(df, paste0(outputDir, "species_rawCounts.tsv"),sep="\t",quote = FALSE, row.names = FALSE)

