#Author: Alicia Sorgen
#Date: 10-12-21
#Description: Builds raw taxa count tables from MetaPhlAn2 output.

## Libraries
library(data.table); message("data.table:", packageVersion("data.table"))
library(stringr); message("stringr:", packageVersion("stringr"))

rm(list=ls())

##### Prep #####
params <- commandArgs(trailingOnly = TRUE)
inputDir = paste0(params[1], "/"); #message("inputDir = ", inputDir)
outputDir <- paste0(params[2], "/");# message("outputDir = ", outputDir)
dir.create(outputDir, showWarnings = FALSE)

funcScript <- params[3]
source(funcScript)


##### Relative abundance #####
file <- "merged_abundance_table.txt"

file.path <- paste0(inputDir, file) ; #message("file.path = ", file.path)
table <- read.table(file.path, header = T, sep = "\t")


phylum <- getphylumTableMetaPhlAn2(table)
class <- getclassTableMetaPhlAn2(table)
order <- getorderTableMetaPhlAn2(table)
family <- getfamilyTableMetaPhlAn2(table)
genus <- getgenusTableMetaPhlAn2(table)
species <- getspeciesTableMetaPhlAn2(table)

taxonomy <- species$species
taxonomy <- gsub("[|]", "_", taxonomy)
taxonomy <- str_split(taxonomy, pattern = "_[pcofgst]__")
taxonomy=do.call(rbind,taxonomy)
taxonomy <- as.data.frame(taxonomy)
colnames(taxonomy)[colnames(taxonomy)=="V1"] <- "Domain"
colnames(taxonomy)[colnames(taxonomy)=="V2"] <- "Phylum"
colnames(taxonomy)[colnames(taxonomy)=="V3"] <- "Class"
colnames(taxonomy)[colnames(taxonomy)=="V4"] <- "Order"
colnames(taxonomy)[colnames(taxonomy)=="V5"] <- "Family"
colnames(taxonomy)[colnames(taxonomy)=="V6"] <- "Genus"
colnames(taxonomy)[colnames(taxonomy)=="V7"] <- "Species"
colnames(taxonomy)[colnames(taxonomy)=="V8"] <- "Strain"
write.table(taxonomy, paste0(outputDir, "Full_Taxonomy.tsv"),sep="\t", quote = FALSE, row.names = F)


tableList <- list(phylum, class, order, family, genus, species)
count <- 1

for (i in tableList) {

 df <- parseMetaPhlAn2Taxa(i)
 df$SampleID <- gsub(pattern = "_profile", replacement = "", df$SampleID)
 df$SampleID <- gsub(pattern = "[.]", replacement = "-", df$SampleID)
 taxaLevel <- names(tableList[[count]])[1]

 write.table(df, paste0(outputDir, taxaLevel, "_RelativeAbundanceCounts.tsv"),sep="\t", quote = FALSE, row.names = F)

 count <- count + 1
}
message("Relative abundance tables COMPLETE!\n")






##### Raw Counts #####
profileDir <- paste0(inputDir, "profileOutput/"); #message("profileDir = ", profileDir)
files <- list.files(profileDir)

index <- 1

for (file in files) {
  
  file.path <- paste0(profileDir, file); #message("file.path = ", file.path)
  sampleID <- sapply(strsplit(file, "_profile"), `[`, 1); #message("sampleID = ",sampleID)
  table <- read.delim(file.path, sep = "\t", header = FALSE)
  table <- table[-(1:2),]
  table <- table[,-(2:4)]
  names(table)[names(table) == "V5"] <- "V2"
  

  phylum.df <- getphylumTableKraken2(table)
  class.df <- getclassTableKraken2(table)
  order.df <- getorderTableKraken2(table)
  family.df <- getfamilyTableKraken2(table)
  genus.df <- getgenusTableKraken2(table)
  species.df <- getspeciesTableKraken2(table)
  # message("***** NO ERROR *****\n")
  
  if (index == 1) {
    phylum <- phylum.df
    class <- class.df
    order <- order.df
    family <- family.df
    genus <- genus.df
    species <- species.df
  } else {
    phylum <- merge(phylum, phylum.df, by = "V1", all = TRUE)
    class <- merge(class, class.df, by = "V1", all = TRUE)
    order <- merge(order, order.df, by = "V1", all = TRUE)
    family <- merge(family, family.df, by = "V1", all = TRUE)
    genus <- merge(genus, genus.df, by = "V1", all = TRUE)
    species <- merge(species, species.df, by = "V1", all = TRUE)
  }
  
  index <- index + 1
}


tableList <- list(phylum, class, order, family, genus, species)
tableName <- list("phylum", "class", "order", "family", "genus", "species")
count <- 1

for (i in tableList) {
  
  df <- transposeTaxaTable(i)
  df[is.na(df)] <- 0
  df <- parseMetaPhlAn2Taxa(i)
  taxaLevel <- tableName[count]
  
  write.table(df, paste0(outputDir, taxaLevel, "_rawCounts.tsv"),sep="\t", quote = FALSE, row.names = F)
  
  count <- count + 1
}
message("Raw count tables COMPLETE!\n")

