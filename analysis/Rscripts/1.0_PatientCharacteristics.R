#Author: Alicia Sorgen
#Date: 2022 June 13
#Description: Summarize patient demographics

##### Edits for script #####
rm(list=ls())
params <- vector()

ANALYSIS <- "microbiome_n124"
params <- c(params, "~/UNCC/Projects/Bariatric_Surgery/Git_Repositories/gut-microbiome-bariatric-weight-outcomes")

moduleRoot <- "1.0_PatientCharacteristics"

##### Libraries #####
R <- sessionInfo()
message(R$R.version$version.string)

library(stringr); message("stringr: Version ", packageVersion("stringr"))
library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
library(scales); message("scales: Version ", packageVersion("scales"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
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
  gitInput   <- file.path(gitRoot, "..", "Data")
  gitScripts <- file.path(gitRoot, "analysis", "Rscripts")

  root <- file.path(gitRoot, "Results")
  dir.create(root, showWarnings = FALSE, recursive = TRUE)

  if (!dir.exists(file.path(root, "input"))) {
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
