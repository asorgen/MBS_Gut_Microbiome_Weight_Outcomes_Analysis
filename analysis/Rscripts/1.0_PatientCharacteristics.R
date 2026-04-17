#Author: Alicia Sorgen
#Date: 2022 June 13
#Description: Summarize patient demographics

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
library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))


# Script-specific edits -------------------------------------------------
params <- vector()
params <- c(params, getOption("mbs.pipe_root", default = stop("Set options(mbs.pipe_root = ...) in your local .Rprofile (gitignored)")))

module <- "1.0_PatientCharacteristics"


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
inputDir = paste0(pipeRoot, "/input/clinical/")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")
sumTableName <- "SummaryTable.tsv"

##### Read in tables #####
weightTable <- read.table(file.path(inputDir, "bs_weight_bmi.tsv"), sep="\t", header=TRUE, check.names=FALSE)
demoTable   <- read.table(file.path(inputDir, "bs_patient_demographics.tsv"), sep="\t", header=TRUE, check.names=FALSE)
myTable <- merge(weightTable, demoTable, by="PatientID")
myTable$time = as.numeric(gsub("BIO-.-...-","",myTable$SampleID))
myTable$Surgery <- myTable$TypeofSurgery
myTable$Surgery <- ifelse(myTable$Surgery == "Gastric Bypass", "RYGB",
                          ifelse(myTable$Surgery == "Sleeve Gastrectomy", "SG",
                                 NA))

row.names(myTable) = myTable$SampleID
BL.df <- myTable[myTable$Timepoint == "BL",]

##### Start summary #####
BL.df$Sex[is.na(BL.df$Sex)] <- "Not recorded"
var <- table(BL.df$Sex)
len <- length(var)

Variable <- vector()
Type <- vector()
Count <- vector()

for (i in 1:len) {
  
  Variable[i] <- "Sex"
  Type[i] <- names(var[i])
  Count[i] <- var[[i]]
  
}

varTable <- data.frame(Variable, Type, Count)
varTable$Percent <- (varTable$Count / sum(varTable$Count)) * 100
final <- varTable

BL.df$Ethnicity[is.na(BL.df$Ethnicity)] <- "Not recorded"
var <- table(BL.df$Ethnicity)
len <- length(var)

Variable <- vector()
Type <- vector()
Count <- vector()

for (i in 1:len) {
  
  Variable[i] <- "Ethnicity"
  Type[i] <- names(var[i])
  Count[i] <- var[[i]]
  
}

varTable <- data.frame(Variable, Type, Count)
varTable$Percent <- (varTable$Count / sum(varTable$Count)) * 100
final <- rbind(final, varTable)

BL.df$TypeofSurgery[is.na(BL.df$TypeofSurgery)] <- "Not recorded"
var <- table(BL.df$TypeofSurgery)
len <- length(var)

Variable <- vector()
Type <- vector()
Count <- vector()

for (i in 1:len) {
  
  Variable[i] <- "TypeofSurgery"
  Type[i] <- names(var[i])
  Count[i] <- var[[i]]
  
}

varTable <- data.frame(Variable, Type, Count)
varTable$Percent <- (varTable$Count / sum(varTable$Count)) * 100
final <- rbind(final, varTable)

BL.df$Site[is.na(BL.df$Site)] <- "Not recorded"
var <- table(BL.df$Site)
len <- length(var)

Variable <- vector()
Type <- vector()
Count <- vector()

for (i in 1:len) {
  
  Variable[i] <- "Site"
  Type[i] <- names(var[i])
  Count[i] <- var[[i]]
  
}

varTable <- data.frame(Variable, Type, Count)
varTable$Percent <- (varTable$Count / sum(varTable$Count)) * 100
final <- rbind(final, varTable)

ages <- BL.df$Age
Variable <- "Age"
Type <- "Avg & St Dev"
Count <- mean(ages, na.rm = TRUE)
Percent <- sd(ages, na.rm = TRUE)

varTable <- data.frame(Variable, Type, Count, Percent)
final <- rbind(final, varTable)

ages <- BL.df$Weight_kg
Variable <- "Weight (kg)"
Type <- "Avg & St Dev"
Count <- mean(ages, na.rm = TRUE)
Percent <- sd(ages, na.rm = TRUE)

varTable <- data.frame(Variable, Type, Count, Percent)
final <- rbind(final, varTable)

ages <- BL.df$BMI_kgm2
Variable <- "BMI (kg/m2)"
Type <- "Avg & St Dev"
Count <- mean(ages, na.rm = TRUE)
Percent <- sd(ages, na.rm = TRUE)

varTable <- data.frame(Variable, Type, Count, Percent)
final <- rbind(final, varTable)

BL.df$Ideal_BL_kg <- BL.df$Baseline_m^2 * 25
BL.df$Excess_kg <- BL.df$Weight_kg - BL.df$Ideal_BL_kg

wgt <- BL.df$Excess_kg
Variable <- "Excess weight (kg)"
Type <- "Avg & St Dev"
Count <- mean(wgt, na.rm = TRUE)
Percent <- sd(wgt, na.rm = TRUE)

varTable <- data.frame(Variable, Type, Count, Percent)
final <- rbind(final, varTable)

file.path <- paste0(outputDir, sumTableName)
write.table(final, file.path, sep="\t",quote = FALSE, row.names = FALSE)
