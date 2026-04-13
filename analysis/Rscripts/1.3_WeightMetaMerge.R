#Author: Alicia Sorgen
#Date: 10-12-21
#Description: Merge taxa count tables with metadata

##### Edits for script #####
rm(list=ls())
params <- vector()

ANALYSIS <- "microbiome_n124"
params <- c(params, "~/UNCC/Projects/Bariatric_Surgery/Git_Repositories/gut-microbiome-bariatric-weight-outcomes")

moduleRoot <- "1.3_WeightMetaMerge"

##### Libraries #####
R <- sessionInfo()
message(R$R.version$version.string)

library(stringr); message("stringr: Version ", packageVersion("stringr"))
# library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
# library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
# library(scales); message("scales: Version ", packageVersion("scales"))
# library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
# library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
# library(nlme); message("nlme: Version ", packageVersion("nlme"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))

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


##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Set input  #####
prevModule <- "ExcessWeightLoss"
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
weightFile <- "updated_weight2.tsv"
weightTable <- read.delim(paste0(inputDir, weightFile), sep="\t",header = TRUE)
metaFile <- args[2]

if (ANALYSIS == "microbiome_n124") {
  
  mergeTable <- weightTable
  
} else {
  
  inputDir = paste0(pipeRoot, "/input/clinical/")
  metaTable <- read.delim(paste0(inputDir, metaFile), sep="\t",header = TRUE)
  
  metaIDs <- unique(metaTable$SampleID)
  weightIDs <- unique(weightTable$SampleID)
  combinedIDs <- c(metaIDs, weightIDs)
  
  duplicateIDs <- combinedIDs[duplicated(combinedIDs)]
  
  not_in_meta <- weightIDs[! weightIDs %in% duplicateIDs]
  not_in_weight <- metaIDs[! metaIDs %in% duplicateIDs]
  
  
  mergeTable <- merge(metaTable, weightTable, by = "SampleID", all = TRUE)
  row.names(mergeTable) = mergeTable$SampleID
  
  
  if (ANALYSIS == "ASA24") {
    ##### Edit HEI data for merging #####
    HEI.fileName <- args[3]
    file.path <- paste0(inputDir, HEI.fileName)
    HEI.df <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
    HEI.df <- na.omit(HEI.df)
    
    HEI.df$Timepoint <- gsub(pattern = "0", replacement = "00", HEI.df$Timepoint)
    HEI.df$Timepoint <- gsub(pattern = "1$", replacement = "01", HEI.df$Timepoint)
    HEI.df$Timepoint <- gsub(pattern = "6", replacement = "06", HEI.df$Timepoint)
    
    SampleID <- paste0(HEI.df$ParticipantID_, "-", HEI.df$Timepoint)
    HEI <- HEI.df$HEI2015_TOTAL_SCORE
    
    df <- data.frame(SampleID, HEI)
    
    
    mergeTable <- merge(df, mergeTable, by = "SampleID", all = TRUE)
    
  }
  
  mergeTable <- mergeTable[!(mergeTable$SampleID %in% not_in_weight),]
  
}

write.table(mergeTable, paste0(outputDir,"metadata.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
