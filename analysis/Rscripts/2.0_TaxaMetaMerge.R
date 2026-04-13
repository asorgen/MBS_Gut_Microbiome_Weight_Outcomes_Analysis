#Author: Alicia Sorgen
#Date: 10-12-21
#Description: Merge taxa count tables with metadata

##### Libraries #####
R <- sessionInfo()
message(R$R.version$version.string)

library(stringr); message("stringr: Version ", packageVersion("stringr"))
# library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
# library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))
# library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))

##### Edits for script #####
rm(list=ls())
params <- vector()

levels <- vector()
levels <- c(levels, "phylum")
levels <- c(levels, "class")
levels <- c(levels, "order")
levels <- c(levels, "family")
levels <- c(levels, "genus")
levels <- c(levels, "species")

ANALYSIS <- "microbiome_n124"
params <- c(params, "~/UNCC/Projects/Bariatric_Surgery/Git_Repositories/gut-microbiome-bariatric-weight-outcomes")
params <- c(params, "MetaPhlAn2")
# params <- c(params, "Kraken2")

moduleRoot <- "2.0_TaxaMetaMerge"

##### Set up working environment #####
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- params
}

if (args[1] == "BLJ") {
  message("\n************* Running in BioLockJ *************")
} else {
  message("\n************* Running locally *************")
  gitRoot    <- args[1]
  inputEnv <- Sys.getenv("INPUT_ROOT"); gitInput <- if (nchar(inputEnv) > 0) inputEnv else file.path(gitRoot, "..", "Data")
  gitScripts <- file.path(gitRoot, "analysis", "Rscripts")

  resultsEnv <- Sys.getenv("RESULTS_ROOT"); root <- if (nchar(resultsEnv) > 0) resultsEnv else file.path(gitRoot, "Results")
  dir.create(root, showWarnings = FALSE, recursive = TRUE)

  if (length(list.files(file.path(root, "input"), recursive = TRUE)) == 0) {
    dir.create(file.path(root, "input"), showWarnings = FALSE, recursive = TRUE)
    invisible(file.copy(list.files(gitInput, full.names = TRUE, include.dirs = TRUE),
                        file.path(root, "input"), recursive = TRUE))
  }

  module <- moduleRoot
  moduleDir <- file.path(root, module)
  dir.create(moduleDir, showWarnings = FALSE)

  outputDir <- file.path(moduleDir, "output")
  dir.create(outputDir, showWarnings = FALSE)
  unlink(list.files(outputDir, full.names = TRUE, recursive = TRUE))

  pipeRoot <- root
}

if (args[1] == "BLJ") {
  pipeRoot  <- dirname(dirname(getwd()))
  moduleDir <- dirname(getwd())
}

##### Set up functions file #####
funcScript <- if (args[1] == "BLJ") file.path(moduleDir, "resources", "functions.R") else file.path(gitScripts, "functions.R")
source(funcScript)


##### Set up input #####
classifier <- args[2]

inputDir = file.path(pipeRoot, "input")
rawDir   = file.path(inputDir, "microbiome", "taxa_raw")
normDir  = file.path(inputDir, "microbiome", "taxa_normalized")
# sampleID_corrections.R defines sample ID corrections specific to this dataset.
# This file is not included in the public repository; see .gitignore.
source(file.path(inputDir, "metadata", "sampleID_corrections.R"))
rawFile <- "_rawCounts"
logFile <- "_LogNormalizedCounts"
relFile <- "_RelativeAbundanceCounts"

metaDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "WeightMetaMerge"),"/output/")
metaFile <- "metadata.tsv"

##### Set up output #####
output = file.path(moduleDir,"output/")

##### Merge data #####
metaTable <- read.delim(paste0(metaDir, metaFile), sep="\t",header = TRUE)


for (level in levels) {
  
  for (File in c(rawFile, logFile, relFile)) {

    if (File == logFile && level == "species") next  # no species-level log-normalized file in Data

    if (File == rawFile)      Table <- read.delim(file.path(rawDir,  paste0("bs_taxa_", level, "_raw.tsv")),         sep="\t", header = TRUE, row.names = 1)
    else if (File == logFile) Table <- read.delim(file.path(normDir, paste0("bs_taxa_", level, "_norm_log10.tsv")), sep="\t", header = TRUE, row.names = 1)
    else                      Table <- getRelAbun(read.delim(file.path(rawDir, paste0("bs_taxa_", level, "_raw.tsv")), sep="\t", header = TRUE, row.names = 1))
    
    SampleID <- sapply(strsplit(rownames(Table), "_"), "[", 1)
    # Strip aliquot suffixes (e.g. "1-241-00-Aliq1" -> "1-241-00") before format conversion
    SampleID <- sub("-Aliq[0-9]+$", "", SampleID)
    # Convert taxa SampleID format "site-patient-timepoint" -> "BIO-site-patient-ZZtimepoint"
    SampleID <- sapply(strsplit(SampleID, "-"), function(x) {
      if (length(x) == 3) sprintf("BIO-%s-%s-%02d", x[1], x[2], as.integer(x[3]))
      else paste(x, collapse="-")
    })
    SampleID <- gsub(pattern = sampleID_mislabeled, replacement = sampleID_corrected, SampleID)
    
    Table <- cbind(SampleID, Table)
    
    if (File == rawFile) {
      dup <- SampleID[which(duplicated(SampleID) == TRUE)]
      dupRows <- vector()
      
      for (i in dup) {
        
        df <- Table[Table$SampleID == i,]
        rowSum <- rowSums(df[2:ncol(df)])
        dupRows <- c(dupRows, which(rownames(Table) %in% rownames(df)[which(rowSum != max(rowSum))]))
        
      } # for (i in dup)
      
    }    
    Table <- Table[-(dupRows),]
    
    df <- merge(metaTable, Table, by = "SampleID", all = TRUE)
    
    for (i in 1:nrow(df)) {
      
      if (is.na(df$PatientID[i]) && !is.na(df$SampleID[i])) {
        string <- strsplit(as.character(df$SampleID[i]),split = "-")
        temp_string=do.call(rbind,string)
        if (ncol(temp_string) < 4) next  # skip IDs that don't match BIO-site-patient-tp format
        df$PatientID[i] <- paste0(temp_string[,1], "-", temp_string[,2], "-", temp_string[,3])
        df$Timepoint[i] <- temp_string[,4]
        df$Site[i] <- temp_string[,2]
      } # if (is.na(df$PatientID[i]))
      
    } # for (i in 1:nrow(df))
    
    df$Timepoint <- gsub(pattern = "24", replacement = "TWENTY_FOUR", df$Timepoint)
    df$Timepoint <- gsub(pattern = "18", replacement = "EIGHTEEN", df$Timepoint)
    df$Timepoint <- gsub(pattern = "12", replacement = "TWELVE", df$Timepoint)
    df$Timepoint <- gsub(pattern = "06", replacement = "SIX", df$Timepoint)
    df$Timepoint <- gsub(pattern = "01", replacement = "ONE", df$Timepoint)
    df$Timepoint <- gsub(pattern = "00", replacement = "BL", df$Timepoint)
    
    df$Site <- gsub(pattern = "1", replacement = "Fargo", df$Site)
    df$Site <- gsub(pattern = "2", replacement = "Cleveland", df$Site)
    
    for (i in 1:nrow(df)) {
      
      if (is.na(df$Sex[i])) {
        
        ID <- df$PatientID[i]
        df1 <- df[df$PatientID == ID,]
        df1 <- na.omit(df1)
        
        if (nrow(df1) > 0) {
          df$Sex[i] <- df1$Sex
          df$TypeofSurgery[i] <- df1$TypeofSurgery
          df$Ethnicity[i] <- df1$Ethnicity
          
        } # if (nrow(df1) > 0)
        
      } # if (is.na(df$Sex[i]))
      
    } # for (i in 1:nrow(df))
    
    fileOutput <- paste0(output, level, File, "_", classifier, ".tsv"); fileOutput
    write.table(df, fileOutput, sep="\t",quote = FALSE, row.names = FALSE)
    
  } # for (File in c(rawFile, logFile, relFile))  
  
} # for (level in levels)
