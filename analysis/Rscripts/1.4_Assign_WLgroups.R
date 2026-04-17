#Author: Alicia Sorgen
#Date: 2022 October 21
#Description: Generate bar plots


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
library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
library(scales); message("scales: Version ", packageVersion("scales"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))


# Script-specific edits -------------------------------------------------
params <- vector()
params <- c(params, getOption("mbs.pipe_root", default = stop("Set options(mbs.pipe_root = ...) in your local .Rprofile (gitignored)")))
# params <- c(params, "Quintile")
# params <- c(params, "Quartile")
params <- c(params, "Tertile")
# params <- c(params, "Half")

module <- "1.4_Assign_WLgroups"


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


##### Set up input #####
prevModule <- "WeightMetaMerge"
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
inputFile <- "metadata.tsv"

month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")
divisions <- c("Quintile", "Quartile", "Tertile", "Half")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Script variables #####
divisions <- args[2:length(args)]

##### Assign weight loss groups #####

# Read in table
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

# Order Timepoint as factor
myTable$Timepoint <- as.factor(myTable$Timepoint)
myTable$Timepoint <- factor(myTable$Timepoint, levels = c("BL", "ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR"))

# Assign SampleID to as dataframe row names
row.names(myTable) = myTable$SampleID

# Convert the timepoint data to numeric values.
myTable$time = as.numeric(myTable$time)

months <- c("12", "18", "24")
included <- c("0", "1", "6")

# Weight loss group assignment
myTable2 <- myTable[ !is.na( myTable$Percent_Loss_kg ), ] # remove rows with percent change NAs
months <- c(1, 6, 12, 18, 24)

for (m in months) {
  
  df.end <- myTable2[myTable2$time == m,]
  divAssignments_endM <- vector()
  
  for (division in divisions) {
    message(paste0("\n************* ", division, " assignment *************\n"))
    
    if (division == "Quintile") {
      div_endM <- quantile(df.end$Percent_Loss_kg , probs=c(0,.2,.4,.6,.8,1), na.rm = TRUE)
      
      for( i in 1:nrow(df.end)) {
        df.end$divAssignments_endM[i] = getQuintileGroup(  df.end$Percent_Loss_kg[i], div_endM )
      }
      
    }
    if (division == "Quartile") {
      div_endM <- quantile(df.end$Percent_Loss_kg , probs=c(0,.25,.5,.75,1), na.rm = TRUE)
      
      for( i in 1:nrow(df.end)) {
        df.end$divAssignments_endM[i] = getQuartileGroup(  df.end$Percent_Loss_kg[i], div_endM )
      }
      
    }
    if (division == "Tertile") {
      div_endM <- quantile(df.end$Percent_Loss_kg , probs=c(0, (1/3), (2/3), 1), na.rm = TRUE)
      
      for( i in 1:nrow(df.end)) {
        df.end$divAssignments_endM[i] = getTertileGroup(  df.end$Percent_Loss_kg[i], div_endM )
      }
      
    }
    if (division == "Half") {
      div_endM <- quantile(df.end$Percent_Loss_kg , probs=c(0,.5,1), na.rm = TRUE)
      
      for( i in 1:nrow(df.end)) {
        df.end$divAssignments_endM[i] = getHalfGroup(  df.end$Percent_Loss_kg[i], div_endM )
      }
      
    }
    
    PatientID <- df.end$PatientID
    Division <- df.end$divAssignments_endM
    PWL_endM <- df.end$Percent_Loss_kg
    df.div <- data.frame(PatientID, Division)
    names(df.div)[names(df.div) == "Division"] <- paste0(division, "Assignment_", m, "M")
    myTable <- merge(df.div, myTable, by = "PatientID", all = TRUE)
    
    
  }
  
  
} # for (m in months)


file.path <- paste0(outputDir, inputFile)
write.table(myTable, file.path, sep="\t",row.names = FALSE,quote = FALSE)
