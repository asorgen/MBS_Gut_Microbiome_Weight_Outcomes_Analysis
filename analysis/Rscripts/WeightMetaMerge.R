#Author: Alicia Sorgen
#Date: 10-12-21
#Description: Merge taxa count tables with metadata

##### Edits for script #####
rm(list=ls())
params <- vector()

ANALYSIS <- "microbiome_n124"
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")

moduleRoot <- paste0("WeightMetaMerge")

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
# library(data.table); message("data.table: Version ", packageVersion("data.table"))

##### Set up working environment #####
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- params
}

if (args[1] == "BLJ") {
  message("\n************* Running in BioLockJ *************")
} else {
  message("\n************* Running locally *************")
  gitRoot <- args[1]
  gitInput <- file.path(gitRoot, "analysis", "input")
  gitScripts <- file.path(gitRoot, "analysis", "Rscripts")
  # message("gitRoot = ", gitRoot)
  
  root <- paste0("~/BioLockJ_pipelines/")
  pipeline <- paste0(ANALYSIS,"_analysis_")
  
  if (any(dir(root) %like% pipeline) == TRUE) {
    root <- paste0(root,"/",str_subset(dir(root), pipeline), "/")
  } else {
    today <- as.character(format(Sys.Date(), "%Y%b%d"))
    root <- paste0(root,ANALYSIS,"_analysis_", today, "/")
    dir.create(root, showWarnings = FALSE)
  }
  
  if (any(dir(root) == "input") == FALSE) {
    # rootInput <- paste0(root, "input/")
    # dir.create(rootInput, showWarnings = FALSE)
    
    file.copy(gitInput,
              root,
              recursive = TRUE)
    
  }
  
  module <- moduleRoot
  
  if (args[2] == "MetaPhlAn2" | args[2] == "Kraken2") {
    module <- paste0(args[2], module)
  } 
  
  if (any(dir(root) %like% module) == TRUE) {
    moduleDir <- paste0(root,str_subset(dir(root), module), "/")
  } else {
    moduleDir <- paste0(root, module, "/")
    dir.create(moduleDir, showWarnings = FALSE)
  }
  
  scriptDir <- paste0(moduleDir, "script/")
  if (any(dir(moduleDir) == "script") == FALSE) {
    dir.create(scriptDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/", moduleRoot, ".R"); script
    file.copy(script,
              scriptDir,
              recursive = TRUE)
  } 
  
  outputDir <- paste0(moduleDir, "output/")
  if (any(dir(moduleDir) == "output") == FALSE) {
    dir.create(outputDir, showWarnings = FALSE)
  }
  
  files <- list.files(outputDir, recursive = TRUE, full.names = TRUE)
  file.remove(files)
  
  resourcesDir <- paste0(moduleDir, "resources/")
  if (any(dir(moduleDir) == "resources") == FALSE) {
    dir.create(resourcesDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/functions.R")
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
  }
  setwd(paste0(moduleDir, "script/"))
  # rm(files, gitInput, gitRoot, gitScripts, module, moduleDir, outputDir, params, pipeline, resourcesDir, root, rootInput, script, today)
  
}

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())

##### Set up functions file #####
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]
rm(str)

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
  
  inputDir = paste0(pipeRoot, "/input/metadataTables/")
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
