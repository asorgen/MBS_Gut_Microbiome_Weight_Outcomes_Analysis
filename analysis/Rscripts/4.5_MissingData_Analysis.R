#' ---
#' title: "Missing data"
#' author: "Alicia Sorgen"
#' date: December 12, 2023
#' ---

## ----Setup, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo=FALSE)
rm(list=ls())
set.seed(1989)


## ----Library, include=FALSE-----------------------------------
R <- sessionInfo()
message(R$R.version$version.string)

library(MASS); message("MASS: Version ", packageVersion("MASS"))
library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
library(lcmm); message("lcmm: Version ", packageVersion("lcmm"))
library(stringr); message("stringr: Version ", packageVersion("stringr"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
# library(scales); message("scales: Version ", packageVersion("scales"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))
library(dplyr); message("dplyr: Version ", packageVersion("dplyr"))


## ----Script-Edits---------------------------------------------
ANALYSIS <- "MetaPhlAn2_microbiome"
date = "2023Dec14"
moduleRoot <- "4.5_MissingData_Analysis"
params <- vector()
params <- c(params, "~/UNCC/Projects/Bariatric_Surgery/Git_Repositories/gut-microbiome-bariatric-weight-outcomes")
params <- c(params, "MetaPhlAn2")
# params <- c(params, "Kraken2")
# params <- c(params, 0)
# params <- c(params, 1)
# params <- c(params, 6)
# params <- c(params, 12)
# params <- c(params, 18)
params <- c(params, 24)

end_month <- params[length(params)]


level <- "genus"


## ----Working-Environment, include=FALSE-----------------------
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


#' 
## ----Functions------------------------------------------------
funcScript <- if (args[1] == "BLJ") file.path(moduleDir, "resources", "functions.R") else file.path(gitScripts, "functions.R")
source(funcScript)


## ----Input----------------------------------------------------
prevModule <- str_subset(dir(pipeRoot), "WeightMetaMerge")
inputDir = paste0(pipeRoot,"/",prevModule,"/output/")
inputFile <- "metadata.tsv"

# classifier <- args[2]
# endM <- args[3]

# prevModule2 <- paste0("TaxaMetaMerge")
# inputDir2 = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule2),"/output/")
# logCountFile <- paste0(level, "_LogNormalizedCounts_", classifier, ".tsv")


## ----Output----------------------------------------------------
outputDir = file.path(moduleDir,"output/")

## ----Script-Variables-----------------------------------------
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

outcomes <- c("Percent_Loss_kg", "Weight_kg")
outcomeLabel <- c("Weight Loss from BL (%)", "Weight_kg")
repeatedMeasure <- c("time")
covariates <- c("Surgery")
subjects <- c("SampleID", "PatientID")
Palette <- c("orange", "blue", "green3", "black", "red")

## ----Read-Table, echo=FALSE-----------------------------------

# Read in table
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

# file.path <- paste0(inputDir2, logCountFile)
# taxa.df <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

# Trim variables
df.1 <- myTable[ , which(colnames(myTable) %in% c(outcomes, repeatedMeasure, covariates, subjects))]
months <- unique(df.1$time)

## ----Data-Summary---------------------------------------------
n <- nrow(df.1)
n.subjects <- length(unique(df.1$PatientID))

# Percent weight
sum.stats <- df.1 %>%
  group_by(time) %>%
  get_summary_stats(Percent_Loss_kg, type = "mean_sd")

knitr::kable(sum.stats[, c(1,3)], "pipe", caption = "RYGB Patients at each timepoint", booktabs = TRUE)

outputName <- paste0("Percent_Loss_kg_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)


# Weight
sum.stats <- df.1 %>%
  group_by(time, Surgery) %>%
  get_summary_stats(Weight_kg, type = "mean_sd")

outputName <- paste0("Weight_kg_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
weight.table <- sum.stats[,c("time", "Surgery", "n")]

## ----Age-Summary-1, include=FALSE-----------------------------
# Trim variables
df.1a <- myTable[ , which(colnames(myTable) %in% c("Age", "Weight_kg", repeatedMeasure, covariates, subjects))]


df.2 <- df.1a[,-which(colnames(df.1a) == "SampleID")] # remove SampleID column
df.3 <- spread(df.2, time, Weight_kg)
df.3$Missing <- ifelse(rowSums(is.na(df.3[,which(colnames(df.3) %in% months)])) == 0, "Full", "Missing")


# RYGB
surg <- "RYGB"
df.4 <- df.3[df.3$Surgery == surg,]
stat.test <- df.4 %>%
  wilcox_test(Age ~ Missing)
p_interpret <- ifelse(stat.test$p < 0.05, "did", "did NOT")

sum.stats <- df.4 %>%
  group_by(Missing) %>%
  get_summary_stats(Age, type = "full")

outputName <- paste0(surg, "_Age_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)


line <- (paste0(surg, " participants (n = ", nrow(df.4), ") with a weight at every time point (n = ", 
                sum.stats$n[1], " (", round((sum.stats$n[1] / nrow(df.4)) * 100, 1), 
                "%); median (IQR) age ", sum.stats$median[1], " (", sum.stats$q1[1], " - ", 
                sum.stats$q3[1], ") years) ", p_interpret, " significantly differ in age from those missing weight records (n = ", 
                sum.stats$n[2], " (", round((sum.stats$n[2] / nrow(df.4)) * 100, 1), 
                "%); median age ", sum.stats$median[2], " (", sum.stats$q1[2],  "- ", 
                sum.stats$q3[2], ") years; p = ", stat.test$p, ").\n"))

write(line,file=file.path(outputDir, "result_interpretations.txt"),append=TRUE)

# SG
surg <- "SG"
df.4 <- df.3[df.3$Surgery == surg,]

stat.test <- df.4 %>%
  wilcox_test(Age ~ Missing)
p_interpret <- ifelse(stat.test$p < 0.05, "did", "did NOT")

sum.stats <- df.4 %>%
  group_by(Missing) %>%
  get_summary_stats(Age, type = "full")
outputName <- paste0(surg, "_Age_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

line <- (paste0(surg, " participants (n = ", nrow(df.4), ") with a weight at every time point (n = ", 
                sum.stats$n[1], " (", round((sum.stats$n[1] / nrow(df.4)) * 100, 1), 
                "%); median (IQR) age ", sum.stats$median[1], " (", sum.stats$q1[1], " - ", 
                sum.stats$q3[1], ") years) ", p_interpret, " significantly differ in age from those missing weight records (n = ", 
                sum.stats$n[2], " (", round((sum.stats$n[2] / nrow(df.4)) * 100, 1), 
                "%); median age ", sum.stats$median[2], " (", sum.stats$q1[2],  "- ", 
                sum.stats$q3[2], ") years; p = ", stat.test$p, ").\n"))
write(line,file=file.path(outputDir, "result_interpretations.txt"),append=TRUE)

## ----Baseline_BMI-Summary-1, include=FALSE-----------------------------
# Trim variables
df.1a <- myTable[ , which(colnames(myTable) %in% c("Baseline_BMI", "Weight_kg", repeatedMeasure, covariates, subjects))]


df.2 <- df.1a[,-which(colnames(df.1a) == "SampleID")]
df.3 <- spread(df.2, time, Weight_kg)
df.3$Missing <- ifelse(rowSums(is.na(df.3[,which(colnames(df.3) %in% months)])) == 0, "Full", "Missing")


# RYGB
surg <- "RYGB"
df.4 <- df.3[df.3$Surgery == surg,]
stat.test <- df.4 %>%
  wilcox_test(Baseline_BMI ~ Missing)
p_interpret <- ifelse(stat.test$p < 0.05, "did", "did NOT")

sum.stats <- df.4 %>%
  group_by(Missing) %>%
  get_summary_stats(Baseline_BMI, type = "full")

outputName <- paste0(surg, "_Baseline_BMI_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

line <- (paste0(surg, " participants (n = ", nrow(df.4), ") with a weight at every time point (n = ", 
                sum.stats$n[1], " (", round((sum.stats$n[1] / nrow(df.4)) * 100, 1), 
                "%); median (IQR) baseline BMI ", sum.stats$median[1], " (", sum.stats$q1[1], " - ", 
                sum.stats$q3[1], ") kg/m^2) ", p_interpret, " significantly differ in baseline BMI from those missing weight records (n = ", 
                sum.stats$n[2], " (", round((sum.stats$n[2] / nrow(df.4)) * 100, 1), 
                "%); median baseline BMI ", sum.stats$median[2], " (", sum.stats$q1[2],  "- ", 
                sum.stats$q3[2], ") kg/m^2; p = ", stat.test$p, ").\n"))
write(line,file=file.path(outputDir, "result_interpretations.txt"),append=TRUE)


# SG
surg <- "SG"
df.4 <- df.3[df.3$Surgery == surg,]

stat.test <- df.4 %>%
  wilcox_test(Baseline_BMI ~ Missing)
p_interpret <- ifelse(stat.test$p < 0.05, "did", "did NOT")

sum.stats <- df.4 %>%
  group_by(Missing) %>%
  get_summary_stats(Baseline_BMI, type = "full")
outputName <- paste0(surg, "_Baseline_BMI_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

line <- (paste0(surg, " participants (n = ", nrow(df.4), ") with a weight at every time point (n = ", 
                sum.stats$n[1], " (", round((sum.stats$n[1] / nrow(df.4)) * 100, 1), 
                "%); median (IQR) baseline BMI ", sum.stats$median[1], " (", sum.stats$q1[1], " - ", 
                sum.stats$q3[1], ") kg/m^2) ", p_interpret, " significantly differ in baseline BMI from those missing weight records (n = ", 
                sum.stats$n[2], " (", round((sum.stats$n[2] / nrow(df.4)) * 100, 1), 
                "%); median baseline BMI ", sum.stats$median[2], " (", sum.stats$q1[2],  "- ", 
                sum.stats$q3[2], ") kg/m^2; p = ", stat.test$p, ").\n"))
write(line,file=file.path(outputDir, "result_interpretations.txt"),append=TRUE)

## ----Site-Summary-1, include=FALSE-----------------------------
# Trim variables
df.1a <- myTable[ , which(colnames(myTable) %in% c("Site", "Weight_kg", repeatedMeasure, covariates, subjects))]


df.2 <- df.1a[,-which(colnames(df.1a) == "SampleID")]
df.3 <- spread(df.2, time, Weight_kg)
df.3$Missing <- ifelse(rowSums(is.na(df.3[,which(colnames(df.3) %in% months)])) == 0, "Full", "Missing")


# RYGB
surg <- "RYGB"
df.4 <- df.3[df.3$Surgery == surg,]
contingency_table <- table(df.4$Site, df.4$Missing)
chi_square_result <- chisq.test(contingency_table)
p_interpret <- ifelse(chi_square_result$p.value < 0.05, "did", "did NOT")

name <- vector()
n_full <- vector()
percent_full <- vector()
n_missing <- vector()
percent_missing <- vector()

for (i in 1:nrow(contingency_table)) {
  
  name[i] <- rownames(contingency_table)[i]
  n_full[i] <- contingency_table[i,1]
  percent_full[i] <- contingency_table[i,1] / sum(contingency_table[i,])*100
  n_missing[i] <- contingency_table[i,2]
  percent_missing[i] <- contingency_table[i,2] / sum(contingency_table[i,])*100
  
}
sum.stats <- data.frame(name, n_full, percent_full, n_missing, percent_missing)
outputName <- paste0(surg, "_Site_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

line <- (paste0("The percentage of ", surg, " particpants with complete weight data ", p_interpret,
                " vary between study sites (", 
                name[1], " ", round(percent_full[1], 0), "%; ", 
                name[2], " ", round(percent_full[2], 0), 
                "%; p = ", round(chi_square_result$p.value, 3), 
                "; chi-square = ", round(chi_square_result$statistic, 3), ").\n"))
write(line,file=file.path(outputDir, "result_interpretations.txt"),append=TRUE)


# SG
surg <- "SG"
df.4 <- df.3[df.3$Surgery == surg,]
contingency_table <- table(df.4$Site, df.4$Missing)
chi_square_result <- chisq.test(contingency_table)
p_interpret <- ifelse(chi_square_result$p.value < 0.05, "did", "did NOT")

name <- vector()
n_full <- vector()
percent_full <- vector()
n_missing <- vector()
percent_missing <- vector()

for (i in 1:nrow(contingency_table)) {
  
  name[i] <- rownames(contingency_table)[i]
  n_full[i] <- contingency_table[i,1]
  percent_full[i] <- contingency_table[i,1] / sum(contingency_table[i,])*100
  n_missing[i] <- contingency_table[i,2]
  percent_missing[i] <- contingency_table[i,2] / sum(contingency_table[i,])*100
  
}
sum.stats <- data.frame(name, n_full, percent_full, n_missing, percent_missing)
outputName <- paste0(surg, "_Site_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

line <- (paste0("The percentage of ", surg, " particpants with complete weight data ", p_interpret,
                " vary between study sites (", 
                name[1], " ", round(percent_full[1], 0), "%; ", 
                name[2], " ", round(percent_full[2], 0), 
                "%; p = ", round(chi_square_result$p.value, 3), 
                "; chi-square = ", round(chi_square_result$statistic, 3), ").\n"))
write(line,file=file.path(outputDir, "result_interpretations.txt"),append=TRUE)

## ----Sex-Summary-1, include=FALSE-----------------------------
# Trim variables
df.1a <- myTable[ , which(colnames(myTable) %in% c("Sex", "Weight_kg", repeatedMeasure, covariates, subjects))]


df.2 <- df.1a[,-which(colnames(df.1a) == "SampleID")]
df.3 <- spread(df.2, time, Weight_kg)
df.3$Missing <- ifelse(rowSums(is.na(df.3[,which(colnames(df.3) %in% months)])) == 0, "Full", "Missing")


# RYGB
surg <- "RYGB"
df.4 <- df.3[df.3$Surgery == surg,]
contingency_table <- table(df.4$Sex, df.4$Missing)
chi_square_result <- chisq.test(contingency_table)
p_interpret <- ifelse(chi_square_result$p.value < 0.05, "did", "did NOT")

name <- vector()
n_full <- vector()
percent_full <- vector()
n_missing <- vector()
percent_missing <- vector()

for (i in 1:nrow(contingency_table)) {
  
  name[i] <- rownames(contingency_table)[i]
  n_full[i] <- contingency_table[i,1]
  percent_full[i] <- contingency_table[i,1] / sum(contingency_table[i,])*100
  n_missing[i] <- contingency_table[i,2]
  percent_missing[i] <- contingency_table[i,2] / sum(contingency_table[i,])*100
  
}
sum.stats <- data.frame(name, n_full, percent_full, n_missing, percent_missing)
outputName <- paste0(surg, "_Sex_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

line <- (paste0("The percentage of ", surg, " particpants with complete weight data ", p_interpret,
                " vary between sexes (", 
                name[1], " ", round(percent_full[1], 0), "%; ", 
                name[2], " ", round(percent_full[2], 0), 
                "%; p = ", round(chi_square_result$p.value, 3), 
                "; chi-square = ", round(chi_square_result$statistic, 3), ").\n"))
write(line,file=file.path(outputDir, "result_interpretations.txt"),append=TRUE)


# SG
surg <- "SG"
df.4 <- df.3[df.3$Surgery == surg,]
contingency_table <- table(df.4$Sex, df.4$Missing)
chi_square_result <- chisq.test(contingency_table)
p_interpret <- ifelse(chi_square_result$p.value < 0.05, "did", "did NOT")

name <- vector()
n_full <- vector()
percent_full <- vector()
n_missing <- vector()
percent_missing <- vector()

for (i in 1:nrow(contingency_table)) {
  
  name[i] <- rownames(contingency_table)[i]
  n_full[i] <- contingency_table[i,1]
  percent_full[i] <- contingency_table[i,1] / sum(contingency_table[i,])*100
  n_missing[i] <- contingency_table[i,2]
  percent_missing[i] <- contingency_table[i,2] / sum(contingency_table[i,])*100
  
}
sum.stats <- data.frame(name, n_full, percent_full, n_missing, percent_missing)
outputName <- paste0(surg, "_Sex_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

line <- (paste0("The percentage of ", surg, " particpants with complete weight data ", p_interpret,
                " vary between sexes (", 
                name[1], " ", round(percent_full[1], 0), "%; ", 
                name[2], " ", round(percent_full[2], 0), 
                "%; p = ", round(chi_square_result$p.value, 3), 
                "; chi-square = ", round(chi_square_result$statistic, 3), ").\n"))
write(line,file=file.path(outputDir, "result_interpretations.txt"),append=TRUE)

## ----Ethnicity-Summary-1, include=FALSE-----------------------------
# Trim variables
df.1a <- myTable[ , which(colnames(myTable) %in% c("Ethnicity", "Weight_kg", repeatedMeasure, covariates, subjects))]


df.2 <- df.1a[,-which(colnames(df.1a) == "SampleID")]
df.3 <- spread(df.2, time, Weight_kg)
df.3$Missing <- ifelse(rowSums(is.na(df.3[,which(colnames(df.3) %in% months)])) == 0, "Full", "Missing")


# RYGB
surg <- "RYGB"
df.4 <- df.3[df.3$Surgery == surg,]
contingency_table <- table(df.4$Ethnicity, df.4$Missing)
chi_square_result <- chisq.test(contingency_table)
p_interpret <- ifelse(chi_square_result$p.value < 0.05, "did", "did NOT")

name <- vector()
n_full <- vector()
percent_full <- vector()
n_missing <- vector()
percent_missing <- vector()

for (i in 1:nrow(contingency_table)) {
  
  name[i] <- rownames(contingency_table)[i]
  n_full[i] <- contingency_table[i,1]
  percent_full[i] <- contingency_table[i,1] / sum(contingency_table[i,])*100
  n_missing[i] <- contingency_table[i,2]
  percent_missing[i] <- contingency_table[i,2] / sum(contingency_table[i,])*100
  
}
sum.stats <- data.frame(name, n_full, percent_full, n_missing, percent_missing)
outputName <- paste0(surg, "_Ethnicity_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

line <- (paste0("The percentage of ", surg, " particpants with complete weight data ", p_interpret,
                " vary between ethnicities (p = ", round(chi_square_result$p.value, 3), 
                "; chi-square = ", round(chi_square_result$statistic, 3), ").\n"))
write(line,file=file.path(outputDir, "result_interpretations.txt"),append=TRUE)



# SG
surg <- "SG"
df.4 <- df.3[df.3$Surgery == surg,]
contingency_table <- table(df.4$Ethnicity, df.4$Missing)
chi_square_result <- chisq.test(contingency_table)
p_interpret <- ifelse(chi_square_result$p.value < 0.05, "did", "did NOT")

name <- vector()
n_full <- vector()
percent_full <- vector()
n_missing <- vector()
percent_missing <- vector()

for (i in 1:nrow(contingency_table)) {
  
  name[i] <- rownames(contingency_table)[i]
  n_full[i] <- contingency_table[i,1]
  percent_full[i] <- contingency_table[i,1] / sum(contingency_table[i,])*100
  n_missing[i] <- contingency_table[i,2]
  percent_missing[i] <- contingency_table[i,2] / sum(contingency_table[i,])*100
  
}
sum.stats <- data.frame(name, n_full, percent_full, n_missing, percent_missing)
outputName <- paste0(surg, "_Ethnicity_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

line <- (paste0("The percentage of ", surg, " particpants with complete weight data ", p_interpret,
                " vary between ethnicities (p = ", round(chi_square_result$p.value, 3), 
                "; chi-square = ", round(chi_square_result$statistic, 3), ").\n"))
write(line,file=file.path(outputDir, "result_interpretations.txt"),append=TRUE)


