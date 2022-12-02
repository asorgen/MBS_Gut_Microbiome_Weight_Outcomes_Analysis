#Author: Alicia Sorgen
#Date: 10-12-21
#Description: Merge taxa count tables with metadata

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
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
params <- c(params, "MetaPhlAn2")
# params <- c(params, "Kraken2")

moduleRoot <- paste0("TaxaMetaMerge")

##### Libraries #####
R <- sessionInfo()
message(R$R.version$version.string)

library(stringr); message("stringr: Version ", packageVersion("stringr"))
# library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
# library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
# library(data.table); message("data.table: Version ", packageVersion("data.table"))
# library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))

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
    module <- paste0(args[2], "_", module)
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

##### Set up input #####
classifier <- args[2]

inputDir = file.path(pipeRoot, "input")
count.InputDir = file.path(inputDir, paste0(classifier, "TaxaTables/"))
rawFile <- "_rawCounts.tsv"
logFile <- "_LogNormalizedCounts.tsv"
relFile <- "_RelativeAbundanceCounts.tsv"

metaDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "WeightMetaMerge"),"/output/")
metaFile <- "metadata.tsv"

##### Set up output #####
output = file.path(moduleDir,"output/")
rawFileOutput <- paste0("_rawCounts_", classifier, ".tsv")
logFileOutput <- paste0("_LogNormalizedCounts_", classifier, ".tsv")
relFileOutput <- paste0("_RelativeAbundanceCounts_", classifier, ".tsv")

##### Merge data #####
metaTable <- read.delim(paste0(metaDir, metaFile), sep="\t",header = TRUE)


for (level in levels) {
  rawTable <- read.delim(paste0(count.InputDir, level, rawFile), sep="\t",header = TRUE, row.names = 1)
  logTable <- read.delim(paste0(count.InputDir, level, logFile), sep="\t",header = TRUE, row.names = 1)
  relTable <- read.delim(paste0(count.InputDir, level, relFile), sep="\t",header = TRUE, row.names = 1)
  
  dup <- which(rownames(rawTable) == "BIO-2-013-01_ACTGAGGT-CTGTTGAC_S051_L005")

  rawTable <- rawTable[-dup,]
  logTable <- logTable[-dup,]
  relTable <- relTable[-dup,]

  dup <- which(rownames(rawTable) == "BIO-1-270-12_CAATGATG-CACGGACG_S56_L006")
  
  rawTable <- rawTable[-dup,]
  logTable <- logTable[-dup,]
  relTable <- relTable[-dup,]

  dup <- which(rownames(rawTable) == "BIO-2-035-24_ATCGTGGA-GGCCGTTG_S50_L005")
  
  rawTable <- rawTable[-dup,]
  logTable <- logTable[-dup,]
  relTable <- relTable[-dup,]
  
  ## extract sample IDs from knead file names
  SampleID <- row.names(rawTable)
  SampleID <- str_split(SampleID, "_[ACGT]")
  SampleID <- do.call(rbind, SampleID)
  SampleID <- SampleID[,1]
  SampleID <- gsub(pattern = "BIO-2-033A-06", replacement = "BIO-2-033-06", SampleID)
  SampleID <- gsub(pattern = "-Aliq2_B6772", replacement = "", SampleID)
  
  ## Duplicate sample 1-241-00 - Aliq1 has lower read depth
  Aliq1 <- which(SampleID == "BIO-1-241-00-Aliq1_B6761")
  SampleID <- SampleID[-Aliq1]
  rawTable <- rawTable[-Aliq1,]
  logTable <- logTable[-Aliq1,]
  relTable <- relTable[-Aliq1,]
  
  rawTable <- cbind(SampleID, rawTable)
  logTable <- cbind(SampleID, logTable)
  relTable <- cbind(SampleID, relTable)
  
  df <- merge(metaTable, rawTable, by = "SampleID", all = TRUE)
  for (i in 1:nrow(df)) {
    
    if (is.na(df$PatientID[i])) {
      string <- strsplit(as.character(df$SampleID[i]),split = "-")
      temp_string=do.call(rbind,string)
      df$PatientID[i] <- paste0(temp_string[,1], "-", temp_string[,2], "-", temp_string[,3])
      df$Timepoint[i] <- temp_string[,4]
      df$Site[i] <- temp_string[,2]
    }
    
  }
  
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
        
        
      }
    }
    
  }
  
  write.table(df, paste0(output,level,rawFileOutput),sep="\t",quote = FALSE, row.names = FALSE)

  df <- merge(metaTable, logTable, by = "SampleID", all = TRUE)
  for (i in 1:nrow(df)) {
    
    if (is.na(df$PatientID[i])) {
      string <- strsplit(as.character(df$SampleID[i]),split = "-")
      temp_string=do.call(rbind,string)
      df$PatientID[i] <- paste0(temp_string[,1], "-", temp_string[,2], "-", temp_string[,3])
      df$Timepoint[i] <- temp_string[,4]
      df$Site[i] <- temp_string[,2]
    }
    
  }
  
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
        
        
      }
    }
    
  }
  write.table(df, paste0(output,level,logFileOutput),sep="\t",quote = FALSE, row.names = FALSE)

  df <- merge(metaTable, relTable, by = "SampleID", all = TRUE)
  for (i in 1:nrow(df)) {
    
    if (is.na(df$PatientID[i])) {
      string <- strsplit(as.character(df$SampleID[i]),split = "-")
      temp_string=do.call(rbind,string)
      df$PatientID[i] <- paste0(temp_string[,1], "-", temp_string[,2], "-", temp_string[,3])
      df$Timepoint[i] <- temp_string[,4]
      df$Site[i] <- temp_string[,2]
    }
    
  }
  
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
        
        
      }
    }
    
  }
  write.table(df, paste0(output,level,relFileOutput),sep="\t",quote = FALSE, row.names = FALSE)
  
}
