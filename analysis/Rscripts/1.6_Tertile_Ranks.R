#Author: Alicia Sorgen
#Date: 2022 July 18
#Description: Plot patient rankings at 12- and 24-months.

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
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))


# Script-specific edits -------------------------------------------------
params <- vector()
params <- c(params, getOption("mbs.pipe_root", default = stop("Set options(mbs.pipe_root = ...) in your local .Rprofile (gitignored)")))

module <- "1.6_Tertile_Ranks"


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
