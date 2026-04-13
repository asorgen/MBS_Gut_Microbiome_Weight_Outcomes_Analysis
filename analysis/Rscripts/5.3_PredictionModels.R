#' ---
#' title: "Prediciton Models"
#' author: "Alicia Sorgen"
#' date: December 7, 2023
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
# library(MuMIn); message("MuMIn: Version ", packageVersion("MuMIn"))


## ----Script-Edits---------------------------------------------
ANALYSIS <- "MetaPhlAn2_microbiome"
date = "2023Dec14"
moduleRoot <- "5.3_PredictionModels"
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
params <- c(params, "LCGA")
# params <- c(params, "GMM-1")
# params <- c(params, "GMM-2")

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
# modelType <- args[4]

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

Formula <- vector()
Model <- vector()
df <- vector()
AIC <- vector()
BIC <- vector()
logLik <- vector()
# R2m <- vector()
# R2c <- vector()
Residual_SE <- vector()
index <- 0
## ----Multivariate_ML, echo=FALSE, fig.height=5, fig.width=10, warning=FALSE----
index <- index +1

df.1a <- df.1
model <- lme( Weight_kg ~ time + Surgery, method = "ML", random = ~1 | PatientID, data = df.1a, na.action = na.omit )
sm <- summary(model)
# Save the summary results to a text file
writeLines(capture.output(summary(model)), 
           con = file.path(outputDir, 
                           paste0(gsub(" ", "", model$call[2]), "_", model$call[5], "_", paste0(model$call[3]), ".txt")))


Formula[index] <- paste0(model$call[2])
Model[index] <- paste0(model$call[5])
df[index] <- paste0(model$call[3])
AIC[index] <- sm$AIC
BIC[index] <- sm$BIC
logLik[index] <- sm$logLik
# R2m[index] <- r.squaredGLMM(model)[1]
# R2c[index] <- r.squaredGLMM(model)[2]
Residual_SE[index] <- model$sigma

## ----Multivariate_REML, echo=FALSE, fig.height=5, fig.width=10, warning=FALSE----
index <- index + 1

df.1a <- df.1
model <- lme( Weight_kg ~ time + Surgery, method = "REML", random = ~1 | PatientID, data = df.1a, na.action = na.omit )
sm <- summary(model)

# Save the summary results to a text file
writeLines(capture.output(summary(model)), 
           con = file.path(outputDir, 
                           paste0(gsub(" ", "", model$call[2]), "_", model$call[5], "_", paste0(model$call[3]), ".txt")))


Formula[index] <- paste0(model$call[2])
Model[index] <- paste0(model$call[5])
df[index] <- paste0(model$call[3])
AIC[index] <- sm$AIC
BIC[index] <- sm$BIC
logLik[index] <- sm$logLik
# R2m[index] <- r.squaredGLMM(model)[1]
# R2c[index] <- r.squaredGLMM(model)[2]
Residual_SE[index] <- model$sigma



## ----Univariate_ML, echo=FALSE, fig.height=5, fig.width=10, warning=FALSE----
index <- index +1

df.1a <- df.1
model <- lme( Weight_kg ~ time, method = "ML", random = ~1 | PatientID, data = df.1a, na.action = na.omit )
sm <- summary(model)
# Save the summary results to a text file
writeLines(capture.output(summary(model)), 
           con = file.path(outputDir, 
                           paste0(gsub(" ", "", model$call[2]), "_", model$call[5], "_", paste0(model$call[3]), ".txt")))


Formula[index] <- paste0(model$call[2])
Model[index] <- paste0(model$call[5])
df[index] <- paste0(model$call[3])
AIC[index] <- sm$AIC
BIC[index] <- sm$BIC
logLik[index] <- sm$logLik
# R2m[index] <- r.squaredGLMM(model)[1]
# R2c[index] <- r.squaredGLMM(model)[2]
Residual_SE[index] <- model$sigma


## ----Univariate_REML, echo=FALSE, fig.height=5, fig.width=10, warning=FALSE----
index <- index + 1

df.1a <- df.1
model <- lme( Weight_kg ~ time, method = "REML", random = ~1 | PatientID, data = df.1a, na.action = na.omit )
sm <- summary(model)

# Save the summary results to a text file
writeLines(capture.output(summary(model)), 
           con = file.path(outputDir, 
                           paste0(gsub(" ", "", model$call[2]), "_", model$call[5], "_", paste0(model$call[3]), ".txt")))


Formula[index] <- paste0(model$call[2])
Model[index] <- paste0(model$call[5])
df[index] <- paste0(model$call[3])
AIC[index] <- sm$AIC
BIC[index] <- sm$BIC
logLik[index] <- sm$logLik
# R2m[index] <- r.squaredGLMM(model)[1]
# R2c[index] <- r.squaredGLMM(model)[2]
Residual_SE[index] <- model$sigma



## ----Multivariate_ML, echo=FALSE, fig.height=5, fig.width=10, warning=FALSE----
index <- index +1

df.1a <- df.1
model <- lme( Weight_kg ~ Surgery + time, method = "ML", random = ~1 | PatientID, data = df.1a, na.action = na.omit )
sm <- summary(model)
# Save the summary results to a text file
writeLines(capture.output(summary(model)), 
           con = file.path(outputDir, 
                           paste0(gsub(" ", "", model$call[2]), "_", model$call[5], "_", paste0(model$call[3]), ".txt")))


Formula[index] <- paste0(model$call[2])
Model[index] <- paste0(model$call[5])
df[index] <- paste0(model$call[3])
AIC[index] <- sm$AIC
BIC[index] <- sm$BIC
logLik[index] <- sm$logLik
# R2m[index] <- r.squaredGLMM(model)[1]
# R2c[index] <- r.squaredGLMM(model)[2]
Residual_SE[index] <- model$sigma


## ----Multivariate_REML, echo=FALSE, fig.height=5, fig.width=10, warning=FALSE----
index <- index + 1

df.1a <- df.1
model <- lme( Weight_kg ~ Surgery + time, method = "REML", random = ~1 | PatientID, data = df.1a, na.action = na.omit )
sm <- summary(model)

# Save the summary results to a text file
writeLines(capture.output(summary(model)), 
           con = file.path(outputDir, 
                           paste0(gsub(" ", "", model$call[2]), "_", model$call[5], "_", paste0(model$call[3]), ".txt")))


Formula[index] <- paste0(model$call[2])
Model[index] <- paste0(model$call[5])
df[index] <- paste0(model$call[3])
AIC[index] <- sm$AIC
BIC[index] <- sm$BIC
logLik[index] <- sm$logLik
# R2m[index] <- r.squaredGLMM(model)[1]
# R2c[index] <- r.squaredGLMM(model)[2]
Residual_SE[index] <- model$sigma

## ----Univariate_ML_Surgery, echo=FALSE, fig.height=5, fig.width=10, warning=FALSE----

index <- index + 1
df.RYGB <- df.1[df.1$Surgery == "RYGB",]
model <- lme( Weight_kg ~ time, method = "ML", random = ~1 | PatientID, data = df.RYGB, na.action = na.omit )
sm <- summary(model)

# Save the summary results to a text file
writeLines(capture.output(summary(model)), 
           con = file.path(outputDir, 
                           paste0(gsub(" ", "", model$call[2]), "_", model$call[5], "_", paste0(model$call[3]), ".txt")))

Formula[index] <- paste0(model$call[2])
Model[index] <- paste0(model$call[5])
df[index] <- paste0(model$call[3])
AIC[index] <- sm$AIC
BIC[index] <- sm$BIC
logLik[index] <- sm$logLik
# R2m[index] <- r.squaredGLMM(model)[1]
# R2c[index] <- r.squaredGLMM(model)[2]
Residual_SE[index] <- model$sigma




index <- index + 1
df.SG <- df.1[df.1$Surgery == "SG",]
model <- lme( Weight_kg ~ time, method = "ML", random = ~1 | PatientID, data = df.SG, na.action = na.omit )
sm <- summary(model)

# Save the summary results to a text file
writeLines(capture.output(summary(model)), 
           con = file.path(outputDir, 
                           paste0(gsub(" ", "", model$call[2]), "_", model$call[5], "_", paste0(model$call[3]), ".txt")))

Formula[index] <- paste0(model$call[2])
Model[index] <- paste0(model$call[5])
df[index] <- paste0(model$call[3])
AIC[index] <- sm$AIC
BIC[index] <- sm$BIC
logLik[index] <- sm$logLik
# R2m[index] <- r.squaredGLMM(model)[1]
# R2c[index] <- r.squaredGLMM(model)[2]
Residual_SE[index] <- model$sigma



## ----Univariate_REML_Surgery, echo=FALSE, fig.height=5, fig.width=10, warning=FALSE----

index <- index + 1
df.RYGB <- df.1[df.1$Surgery == "RYGB",]
model <- lme( Weight_kg ~ time, method = "REML", random = ~1 | PatientID, data = df.RYGB, na.action = na.omit )
sm <- summary(model)

# Save the summary results to a text file
writeLines(capture.output(summary(model)), 
           con = file.path(outputDir, 
                           paste0(gsub(" ", "", model$call[2]), "_", model$call[5], "_", paste0(model$call[3]), ".txt")))

Formula[index] <- paste0(model$call[2])
Model[index] <- paste0(model$call[5])
df[index] <- paste0(model$call[3])
AIC[index] <- sm$AIC
BIC[index] <- sm$BIC
logLik[index] <- sm$logLik
# R2m[index] <- r.squaredGLMM(model)[1]
# R2c[index] <- r.squaredGLMM(model)[2]
Residual_SE[index] <- model$sigma




index <- index + 1
df.SG <- df.1[df.1$Surgery == "SG",]
model <- lme( Weight_kg ~ time, method = "REML", random = ~1 | PatientID, data = df.SG, na.action = na.omit )
sm <- summary(model)

# Save the summary results to a text file
writeLines(capture.output(summary(model)), 
           con = file.path(outputDir, 
                           paste0(gsub(" ", "", model$call[2]), "_", model$call[5], "_", paste0(model$call[3]), ".txt")))

Formula[index] <- paste0(model$call[2])
Model[index] <- paste0(model$call[5])
df[index] <- paste0(model$call[3])
AIC[index] <- sm$AIC
BIC[index] <- sm$BIC
logLik[index] <- sm$logLik
# R2m[index] <- r.squaredGLMM(model)[1]
# R2c[index] <- r.squaredGLMM(model)[2]
Residual_SE[index] <- model$sigma





## ----Summary--------------------------------------------------
results <- data.frame(Formula, Model, df, AIC, BIC, logLik, 
                      # R2m, R2c, 
                      Residual_SE)
# results <- results[order(results$R2c, decreasing = TRUE),]

outputName <- paste0("model_summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(results, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

