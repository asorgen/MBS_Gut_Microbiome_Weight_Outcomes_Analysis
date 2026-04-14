#' Author: Alicia Sorgen
#' Date: December 7, 2023
#' Description: Prediction Models

# Setup -----
rm(list=ls())
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
set.seed(1989)


# Libraries -----
H1("Libraries"); start <- Sys.time()
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


# Script parameters -----
H1("Script Parameters"); start <- Sys.time()
params <- vector()
params <- c(params, getOption("mbs.pipe_root", default = stop("Set options(mbs.pipe_root = ...) in your local .Rprofile (gitignored)")))
params <- c(params, "MetaPhlAn2")
# params <- c(params, "Kraken2")
params <- c(params, 24)
params <- c(params, "LCGA")
# params <- c(params, "GMM-1")
# params <- c(params, "GMM-2")

module <- "5.3_PredictionModels"
level <- "genus"


# Working environment -----
H1("Working Environment"); start <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- params
}

message("\n************* Running locally *************")
proj_root    <- args[1]
inputEnv <- Sys.getenv("INPUT_ROOT"); input_root <- if (nchar(inputEnv) > 0) inputEnv else file.path(dirname(proj_root), "Data")
script_root <- file.path(proj_root, basename(proj_root), "analysis", "Rscripts")

resultsEnv <- Sys.getenv("RESULTS_ROOT"); pipeRoot <- if (nchar(resultsEnv) > 0) resultsEnv else file.path(proj_root, gsub('Analysis', 'Results', basename(proj_root)))
# dir.create(pipeRoot, showWarnings = FALSE, recursive = TRUE)

# if (length(list.files(file.path(pipeRoot, "input"), recursive = TRUE)) == 0) {
#   dir.create(file.path(pipeRoot, "input"), showWarnings = FALSE, recursive = TRUE)
#   invisible(file.copy(list.files(input_root, full.names = TRUE, include.dirs = TRUE),
#                       file.path(pipeRoot, "input"), recursive = TRUE))
# }

moduleDir <- file.path(pipeRoot, module)
# dir.create(moduleDir, showWarnings = FALSE)

outputDir <- file.path(moduleDir, "output")
# dir.create(outputDir, showWarnings = FALSE)
unlink(list.files(outputDir, full.names = TRUE, recursive = TRUE))

rm(proj_root, inputEnv, resultsEnv, module, params)


# Functions -----
H1("Functions"); start <- Sys.time()
source(file.path(script_root, "functions.R"))


# Input -----
H1("Input"); start <- Sys.time()
prevModule <- str_subset(dir(pipeRoot), "WeightMetaMerge")
inputDir = paste0(pipeRoot,"/",prevModule,"/output/")
inputFile <- "metadata.tsv"

# classifier <- args[2]
# endM <- args[3]
# modelType <- args[4]

# prevModule2 <- paste0("TaxaMetaMerge")
# inputDir2 = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule2),"/output/")
# logCountFile <- paste0(level, "_LogNormalizedCounts_", classifier, ".tsv")


# Output -----
H1("Output"); start <- Sys.time()
outputDir = file.path(moduleDir,"output/")

# Script variables -----
H1("Script Variables"); start <- Sys.time()
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

outcomes <- c("Percent_Loss_kg", "Weight_kg")
outcomeLabel <- c("Weight Loss from BL (%)", "Weight_kg")
repeatedMeasure <- c("time")
covariates <- c("Surgery")
subjects <- c("SampleID", "PatientID")
Palette <- c("orange", "blue", "green3", "black", "red")

# Read data -----
H1("Read Data"); start <- Sys.time()
# Read in table
file.path <- paste0(inputDir, inputFile)
myTable <- import_file(file.path)

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

# Multivariate ML: time + Surgery -----
H1("Multivariate ML: time + Surgery"); start <- Sys.time()
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

# Multivariate REML: time + Surgery -----
H1("Multivariate REML: time + Surgery"); start <- Sys.time()
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


# Univariate ML -----
H1("Univariate ML: time"); start <- Sys.time()
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


# Univariate REML -----
H1("Univariate REML: time"); start <- Sys.time()
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


# Multivariate ML: Surgery + time -----
H1("Multivariate ML: Surgery + time"); start <- Sys.time()
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


# Multivariate REML: Surgery + time -----
H1("Multivariate REML: Surgery + time"); start <- Sys.time()
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

# Univariate ML by surgery -----
H1("Univariate ML by surgery"); start <- Sys.time()
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


# Univariate REML by surgery -----
H1("Univariate REML by surgery"); start <- Sys.time()
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


# Summary -----
H1("Summary"); start <- Sys.time()
results <- data.frame(Formula, Model, df, AIC, BIC, logLik,
                      # R2m, R2c,
                      Residual_SE)
# results <- results[order(results$R2c, decreasing = TRUE),]

outputName <- paste0("model_summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(results, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
