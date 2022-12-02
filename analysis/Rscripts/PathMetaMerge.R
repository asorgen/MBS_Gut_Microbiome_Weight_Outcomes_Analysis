#Author: Alicia Sorgen
#Date: 2022 Nov 15
#Description: Merge HumanN2 pathway abundance table with metadata.

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
params <- c(params, 0)
params <- c(params, 1)
params <- c(params, 6)
params <- c(params, 12)
params <- c(params, 18)
params <- c(params, 24)
# end_month <- params[length(params)]

moduleRoot <- "PathMetaMerge"

##### Libraries #####
R <- sessionInfo()
message(R$R.version$version.string)

library(stringr); message("stringr: Version ", packageVersion("stringr"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
library(ggrepel); message("ggrepel: Version ", packageVersion("ggrepel"))

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
  rm(gitRoot)
  
  root <- paste0("~/BioLockJ_pipelines/")
  pipeline <- paste0(ANALYSIS,"_analysis_")
  
  if (any(dir(root) %like% pipeline) == TRUE) {
    root <- paste0(root,"/",str_subset(dir(root), pipeline), "/")
  } else {
    today <- as.character(format(Sys.Date(), "%Y%b%d"))
    root <- paste0(root,ANALYSIS,"_analysis_", today, "/")
    dir.create(root, showWarnings = FALSE)
  }
  rm(pipeline)
  
  if (any(dir(root) == "input") == FALSE) {
    file.copy(gitInput,
              root,
              recursive = TRUE)
    
  }
  rm(gitInput)
  module <- moduleRoot
  
  if (args[2] == "MetaPhlAn2" | args[2] == "Kraken2") {
    module <- paste0(args[2], "_", module)
  } 
  
  if (args[3] %in% c("Quintile", "Quartile", "Thirds", "Half")) {
    module <- paste0(module, "_", args[3])
  } 
  
  if (exists("end_month") == TRUE) {
    module <- paste0(module, "_BLto", end_month, "M")
  }
  
  if (any(dir(root) %like% module) == TRUE) {
    moduleDir <- paste0(root,str_subset(dir(root), module), "/")
  } else {
    moduleDir <- paste0(root, module, "/")
    dir.create(moduleDir, showWarnings = FALSE)
  }
  rm(module, root)
  
  scriptDir <- paste0(moduleDir, "script/")
  if (any(dir(moduleDir) == "script") == FALSE) {
    dir.create(scriptDir, showWarnings = FALSE)
  }
  rm(scriptDir)
  
  outputDir <- paste0(moduleDir, "output/")
  if (any(dir(moduleDir) == "output") == FALSE) {
    dir.create(outputDir, showWarnings = FALSE)
  }
  
  files <- list.files(outputDir, recursive = TRUE, full.names = TRUE)
  file.remove(files)
  rm(outputDir, files)
  
  resourcesDir <- paste0(moduleDir, "resources/")
  if (any(dir(moduleDir) == "resources") == FALSE) {
    dir.create(resourcesDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/functions.R")
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
    
    script = paste0(gitScripts,"/", moduleRoot, ".R"); script
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
    rm(script)
    
  }
  rm(resourcesDir, gitScripts)
  
  setwd(paste0(moduleDir, "script/"))
  rm(moduleDir)
  
}
rm(params, moduleRoot, ANALYSIS, end_month, R)

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
included <- args[2:length(args)]

inputDir = paste0(pipeRoot,"/input/")
humanN2File <- "humanN2_pathabundance_cpm.tsv"

metaDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "WeightMetaMerge"),"/output/")
metaFile <- "metadata.tsv"

month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")




##### Read in metadata file #####
metaTable <- read.delim(paste0(metaDir, metaFile), sep="\t",header = TRUE)

##### Read in HumanN2 file #####
humanN <- read.delim(paste0(inputDir, humanN2File), sep="\t",header = TRUE, row.names = 1)
humanN2 <- as.data.frame(t(humanN))

##### Filter out unclassified/unintegrated reads #####
message(">> ", ncol(humanN2), " pathways")
colsToRemove <- vector()

UNMAPPED <- colnames(humanN2)[which(colnames(humanN2) %like% "UNMAPPED")]
message(">> ", length(UNMAPPED), " UNMAPPED pathway")
colsToRemove <- c(colsToRemove, UNMAPPED)

UNINTEGRATED <- colnames(humanN2)[which(colnames(humanN2) %like% "UNINTEGRATED")]
message(">> ", length(UNINTEGRATED), " UNINTEGRATED pathways")
colsToRemove <- c(colsToRemove, UNINTEGRATED)

unclassified <- colnames(humanN2)[which(colnames(humanN2) %like% "unclassified")]
message(">> ", length(unclassified), " unclassified species")
colsToRemove <- c(colsToRemove, unclassified)

noTaxa <- colnames(humanN2)[which(!(colnames(humanN2) %like% ".s__"))]
message(">> ", length(noTaxa), " pathways with no taxa")
colsToRemove <- c(colsToRemove, noTaxa)

colsToRemove <- unique(colsToRemove)
colsToRemove <- which(colnames(humanN2) %in% colsToRemove)
message(">> ", length(colsToRemove), " unique pathways removed")

humanN3 <- humanN2[, -colsToRemove]
message(">> ", ncol(humanN3), " pathways remaining")

##### extract sample IDs from knead file names #####
SampleID <- row.names(humanN2)
SampleID <- str_split(SampleID, "_[ACGT]")
SampleID <- do.call(rbind, SampleID)
SampleID <- SampleID[,1]
SampleID <- gsub("X", "", SampleID)

## Duplicate sample 1-241-00 - Aliq1 has lower read depth
Aliq1 <- which(SampleID %like% "Aliq1_B6761")
SampleID <- SampleID[-Aliq1]
humanN2 <- humanN2[-Aliq1,]

tempStr <- str_split(SampleID, "[.]")
tempStr <- do.call(rbind, tempStr)
tempStr[,3] <- str_pad(tempStr[,3], 3, pad = "0") # add leading zeros to subject ID
tempStr[,4] <- str_pad(tempStr[,4], 2, pad = "0") # add leading zeros to time point


SampleID <- paste0(tempStr[,1], "-", tempStr[,2], "-", tempStr[,3], "-", tempStr[,4])
# PatientID <- paste0(tempStr[,1], "-", tempStr[,2], "-", tempStr[,3])
humanN2 <- cbind(SampleID, humanN2)
# humanN2 <- cbind(PatientID, humanN2)

df <- merge(metaTable, humanN2, by = "SampleID")
write.table(df, paste0(outputDir,"metaHumanN.tsv"),sep="\t",quote = FALSE, row.names = FALSE)
