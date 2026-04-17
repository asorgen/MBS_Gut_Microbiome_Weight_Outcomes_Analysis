#Author: Alicia Sorgen
#Date: 10-12-21
#Description: Merge taxa count tables with metadata

# Set up ------------------------------------------------------------------
rm(list=ls())
set.seed(1989)
H1 = function(comment) {
   delim = "#"
   side = "#"
   fill = " "
   maxLen <- 90
   lineLen <- 60
   message("")
   message(paste0(strrep(delim, maxLen)))
   if (nchar(comment) > lineLen) {
      cut <- strsplit(comment, " ")
      line <- ""
      for (word in cut[[1]]) {
         if ((nchar(line) + 1 + nchar(word)) > lineLen) {
            edge1 <- round((maxLen - nchar(line))/2, 0) - 2
            edge2 <- maxLen - edge1 - nchar(line) - 4
            line_print <- paste0(side, strrep(fill, edge1), " ", line, " ", strrep(fill, edge2), side)
            message(line_print)
            line <- word
         } else {
            line <- paste(line, word)
         }
      }
      edge1 <- round((maxLen - nchar(line))/2, 0) - 2
      edge2 <- maxLen - edge1 - nchar(line) - 4
      line_print <- paste0(side, strrep(fill, edge1), " ", line, " ", strrep(fill, edge2), side)
      message(line_print)
   } else {
      edge1 <- round((maxLen - nchar(comment))/2, 0) - 2
      edge2 <- maxLen - edge1 - nchar(comment) - 4
      line_print <- paste0(side, strrep(fill, edge1), " ", comment, " ", strrep(fill, edge2), side)
      message(line_print)
   }
   message(paste0(strrep(delim, maxLen)))
   message("")
}


# Libraries ----------------------------------------------------------------
H1("Libraries")
R <- sessionInfo()
message(R$R.version$version.string); rm(R)

library(stringr); message("stringr: Version ", packageVersion("stringr"))
# library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
# library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
# library(scales); message("scales: Version ", packageVersion("scales"))
# library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
# library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
# library(nlme); message("nlme: Version ", packageVersion("nlme"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))


# Script-specific edits -------------------------------------------------
ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, getOption("mbs.pipe_root", default = stop("Set options(mbs.pipe_root = ...) in your local .Rprofile (gitignored)")))

module <- "1.3_WeightMetaMerge"


# Set working environment ---------------------------------------------------------
H1("Working Environment")
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
   args <- params
}

message("\n************* Running locally *************")
proj_root    <- args[1]
message("Project root directory: ", proj_root, "\n")
inputEnv <- Sys.getenv("INPUT_ROOT"); input_root <- if (nchar(inputEnv) > 0) inputEnv else file.path(dirname(proj_root), "Data")
script_root <- file.path(proj_root, basename(proj_root), "analysis", "Rscripts")

resultsEnv <- Sys.getenv("RESULTS_ROOT"); pipeRoot <- if (nchar(resultsEnv) > 0) resultsEnv else file.path(proj_root, gsub('Analysis', 'Results', basename(proj_root)))
# dir.create(pipeRoot, showWarnings = FALSE, recursive = TRUE)

# if (length(list.files(file.path(pipeRoot, "input"), recursive = TRUE)) == 0) {
#    dir.create(file.path(pipeRoot, "input"), showWarnings = FALSE, recursive = TRUE)
#    invisible(file.copy(list.files(input_root, full.names = TRUE, include.dirs = TRUE),
#                        file.path(pipeRoot, "input"), recursive = TRUE))
# }

moduleDir <- file.path(pipeRoot, module)
# dir.create(moduleDir, showWarnings = FALSE)

outputDir <- file.path(moduleDir, "output/")
# dir.create(outputDir, showWarnings = FALSE)
unlink(list.files(outputDir, full.names = TRUE, recursive = TRUE))

rm(proj_root, inputEnv, resultsEnv, module, params)


# Set functions ------------------------------------------------
H1("Functions")
message("Loading functions from: ", script_root)
funcScript <- file.path(script_root, "functions.R")
source(funcScript); rm(funcScript)


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
