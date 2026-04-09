#Author: Alicia Sorgen
#Date: 2022 July 18
#Description: Plot patient rankings at 12- and 24-months.

##### Edits for script #####
rm(list=ls())
params <- vector()

ANALYSIS <- "microbiome_n124"
params <- c(params, "~/UNCC/Projects/Bariatric_Surgery/Git_Repositories/gut-microbiome-bariatric-weight-outcomes")

moduleRoot <- "1.6_Tertile_Ranks"

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
prevModule <- "WeightMetaMerge"
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
inputFile <- "metadata.tsv"

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Formatting variables #####
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")
tags <- c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l.", "m.", "n.", "o.", "p.", "q.", "r.", "s.", "t.", "u.", "v.", "w.", "x.", "y.", "z.")
statusColors <- c("tomato", "steelblue")
surgeryColors <- c("#3171BC", "#E8C241")

##### Prep data table for analysis #####
# Read in table
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

# Data columns
PatientID <- myTable$PatientID
Timepoint <- myTable$time
Percent_Loss_kg <- myTable$Percent_Loss_kg

df <- data.frame(PatientID, Timepoint, Percent_Loss_kg)

# Filter only 12-month weights
df.12M <- df[df$Timepoint == 12,]
df.12M <- na.omit(df.12M)
df.12M <- df.12M[order(df.12M$Percent_Loss_kg, decreasing = FALSE),]
df.12M$Rank_12M <- 1:nrow(df.12M)
# quintiles_12M <- quantile(df.12M$Percent_Loss_kg , probs=c(0,.2,.4,.6,.8,1), na.rm = TRUE)
Q_12M <- nrow(df.12M) / 3

# Filter only 18-month weights
df.18M <- df[df$Timepoint == 18,]
df.18M <- na.omit(df.18M)
df.18M <- df.18M[order(df.18M$Percent_Loss_kg, decreasing = FALSE),]
df.18M$Rank_18M <- 1:nrow(df.18M)
# quintiles_18M <- quantile(df.18M$Percent_Loss_kg , probs=c(0,.2,.4,.6,.8,1), na.rm = TRUE)
Q_18M <- nrow(df.18M) / 3

# Filter only 24-month weights
df.24M <- df[df$Timepoint == 24,]
df.24M <- na.omit(df.24M)
df.24M <- df.24M[order(df.24M$Percent_Loss_kg, decreasing = FALSE),]
df.24M$Rank_24M <- 1:nrow(df.24M)
# quintiles_24M <- quantile(df.24M$Percent_Loss_kg , probs=c(0,.2,.4,.6,.8,1), na.rm = TRUE)
Q_24M <- nrow(df.24M) / 3

plotList <- list()

df.merge <- merge(df.12M, df.24M, by = "PatientID")

plot <- ggscatter(df.merge, x = "Rank_12M", y = "Rank_24M",
                  shape = 21, size = 2.5, # Points shape and size
                  # add = "reg.line",  # Add regression line
                  # conf.int = TRUE, # Add confidence interval
                  # cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  # cor.coeff.args = list(method = "kendall"),
                  add.params = list(fill = "lightgray")
); plot

x.lab <- "12-month weight loss rank"
y.lab <- "24-month weight loss rank"
plot <- plot + labs(x=x.lab, y = y.lab); plot

plot <- plot + geom_hline(yintercept=Q_24M, linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_hline(yintercept=(Q_24M*2), linetype="dashed",
                          color = "red", size=0.5); plot
# plot <- plot + geom_hline(yintercept=(Q_24M*3), linetype="dashed", 
#                           color = "red", size=0.5); plot
# plot <- plot + geom_hline(yintercept=(Q_24M*4), linetype="dashed", 
#                           color = "red", size=0.5); plot

plot <- plot + geom_vline(xintercept=Q_12M, linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_vline(xintercept=(Q_12M*2), linetype="dashed",
                          color = "red", size=0.5); plot
# plot <- plot + geom_vline(xintercept=(Q_12M*3), linetype="dashed", 
#                           color = "red", size=0.5); plot
# plot <- plot + geom_vline(xintercept=(Q_12M*4), linetype="dashed", 
#                           color = "red", size=0.5); plot

plot <- plot + annotate(geom="text", x=10, y=(Q_24M/2), label="Bottom 3rd",
                          color="red"); plot
plot <- plot + annotate(geom="text", x=nrow(df.12M)-10, y=nrow(df.24M), label="Top 3rd",
                        color="red"); plot
plotList[[1]] <- plot


df.merge <- merge(df.18M, df.24M, by = "PatientID")

plot <- ggscatter(df.merge, x = "Rank_18M", y = "Rank_24M",
                  shape = 21, size = 2.5, # Points shape and size
                  # add = "reg.line",  # Add regression line
                  # conf.int = TRUE, # Add confidence interval
                  # cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  # cor.coeff.args = list(method = "kendall"),
                  add.params = list(fill = "lightgray")
); plot

x.lab <- "18-month weight loss rank"
y.lab <- "24-month weight loss rank"
plot <- plot + labs(x=x.lab, y = y.lab); plot

plot <- plot + geom_hline(yintercept=Q_24M, linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_hline(yintercept=(Q_24M*2), linetype="dashed",
                          color = "red", size=0.5); plot
# plot <- plot + geom_hline(yintercept=(Q_24M*3), linetype="dashed", 
#                           color = "red", size=0.5); plot
# plot <- plot + geom_hline(yintercept=(Q_24M*4), linetype="dashed", 
#                           color = "red", size=0.5); plot

plot <- plot + geom_vline(xintercept=Q_18M, linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_vline(xintercept=(Q_18M*2), linetype="dashed",
                          color = "red", size=0.5); plot
# plot <- plot + geom_vline(xintercept=(Q_18M*3), linetype="dashed", 
#                           color = "red", size=0.5); plot
# plot <- plot + geom_vline(xintercept=(Q_18M*4), linetype="dashed", 
#                           color = "red", size=0.5); plot

plot <- plot + annotate(geom="text", x=10, y=(Q_24M/2), label="Bottom 3rd",
                        color="red"); plot
plot <- plot + annotate(geom="text", x=nrow(df.18M)-10, y=nrow(df.24M), label="Top 3rd",
                        color="red"); plot
plotList[[2]] <- plot


df.merge <- merge(df.12M, df.18M, by = "PatientID")

plot <- ggscatter(df.merge, x = "Rank_12M", y = "Rank_18M",
                  shape = 21, size = 2.5, # Points shape and size
                  # add = "reg.line",  # Add regression line
                  # conf.int = TRUE, # Add confidence interval
                  # cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  # cor.coeff.args = list(method = "kendall"),
                  add.params = list(fill = "lightgray")
); plot

x.lab <- "12-month weight loss rank"
y.lab <- "18-month weight loss rank"
plot <- plot + labs(x=x.lab, y = y.lab); plot

plot <- plot + geom_hline(yintercept=Q_18M, linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_hline(yintercept=(Q_18M*2), linetype="dashed",
                          color = "red", size=0.5); plot
# plot <- plot + geom_hline(yintercept=(Q_18M*3), linetype="dashed", 
#                           color = "red", size=0.5); plot
# plot <- plot + geom_hline(yintercept=(Q_18M*4), linetype="dashed", 
#                           color = "red", size=0.5); plot

plot <- plot + geom_vline(xintercept=Q_12M, linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_vline(xintercept=(Q_12M*2), linetype="dashed",
                          color = "red", size=0.5); plot
# plot <- plot + geom_vline(xintercept=(Q_12M*3), linetype="dashed", 
#                           color = "red", size=0.5); plot
# plot <- plot + geom_vline(xintercept=(Q_12M*4), linetype="dashed", 
#                           color = "red", size=0.5); plot

plot <- plot + annotate(geom="text", x=10, y=(Q_18M/2), label="Bottom 3rd",
                        color="red"); plot
plot <- plot + annotate(geom="text", x=nrow(df.12M)-10, y=nrow(df.18M), label="Top 3rd",
                        color="red"); plot
plotList[[3]] <- plot


plotFileName <- "weight_loss_ranking_plots_third.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 5, height = 5)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]],
               ncol = 1, nrow = 1)
}
dev.off()
