#Author: Alicia Sorgen
#Date: 2022 July 18
#Description: Plot patient rankings at 12- and 24-months.

##### Edits for script #####
rm(list=ls())
params <- vector()

ANALYSIS <- "microbiome_n124"
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")

moduleRoot <- paste0("Quintile_Ranks")

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
  gitRoot <- args[1]
  gitInput <- file.path(gitRoot, "analysis", "input")
  gitScripts <- file.path(gitRoot, "analysis", "Rscripts")
  # message("gitRoot = ", gitRoot)
  
  root <- paste0("~/BioLockJ_pipelines/")
  pipeline <- paste0(ANALYSIS,"_analysis_")
  
  if (any(dir(root) %like% pipeline) == TRUE) {
    root <- paste0(root,"/",str_subset(dir(root), pipeline), "/")
  } else {
    today <- as.character(format(Sys.Date(), "%Y%b%d"))
    root <- paste0(root,ANALYSIS,"_analysis_", today, "/")
    dir.create(root, showWarnings = FALSE)
  }
  
  if (any(dir(root) == "input") == FALSE) {
    # rootInput <- paste0(root, "input/")
    # dir.create(rootInput, showWarnings = FALSE)
    
    file.copy(gitInput,
              root,
              recursive = TRUE)
    
  }
  
  module <- moduleRoot
  
  if (args[2] == "MetaPhlAn2" | args[2] == "Kraken2") {
    module <- paste0(args[2], module)
  } 
  
  if (any(dir(root) %like% module) == TRUE) {
    moduleDir <- paste0(root,str_subset(dir(root), module), "/")
  } else {
    moduleDir <- paste0(root, module, "/")
    dir.create(moduleDir, showWarnings = FALSE)
  }
  
  scriptDir <- paste0(moduleDir, "script/")
  if (any(dir(moduleDir) == "script") == FALSE) {
    dir.create(scriptDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/", moduleRoot, ".R"); script
    file.copy(script,
              scriptDir,
              recursive = TRUE)
  } 
  
  outputDir <- paste0(moduleDir, "output/")
  if (any(dir(moduleDir) == "output") == FALSE) {
    dir.create(outputDir, showWarnings = FALSE)
  }
  
  files <- list.files(outputDir, recursive = TRUE, full.names = TRUE)
  file.remove(files)
  
  resourcesDir <- paste0(moduleDir, "resources/")
  if (any(dir(moduleDir) == "resources") == FALSE) {
    dir.create(resourcesDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/functions.R")
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
  }
  setwd(paste0(moduleDir, "script/"))
  # rm(files, gitInput, gitRoot, gitScripts, module, moduleDir, outputDir, params, pipeline, resourcesDir, root, rootInput, script, today)
  
}

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
Q_12M <- nrow(df.12M) / 5

# Filter only 18-month weights
df.18M <- df[df$Timepoint == 18,]
df.18M <- na.omit(df.18M)
df.18M <- df.18M[order(df.18M$Percent_Loss_kg, decreasing = FALSE),]
df.18M$Rank_18M <- 1:nrow(df.18M)
# quintiles_18M <- quantile(df.18M$Percent_Loss_kg , probs=c(0,.2,.4,.6,.8,1), na.rm = TRUE)
Q_18M <- nrow(df.18M) / 5

# Filter only 24-month weights
df.24M <- df[df$Timepoint == 24,]
df.24M <- na.omit(df.24M)
df.24M <- df.24M[order(df.24M$Percent_Loss_kg, decreasing = FALSE),]
df.24M$Rank_24M <- 1:nrow(df.24M)
# quintiles_24M <- quantile(df.24M$Percent_Loss_kg , probs=c(0,.2,.4,.6,.8,1), na.rm = TRUE)
Q_24M <- nrow(df.24M) / 5

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
plot <- plot + geom_hline(yintercept=(Q_24M*3), linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_hline(yintercept=(Q_24M*4), linetype="dashed", 
                          color = "red", size=0.5); plot

plot <- plot + geom_vline(xintercept=Q_12M, linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_vline(xintercept=(Q_12M*2), linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_vline(xintercept=(Q_12M*3), linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_vline(xintercept=(Q_12M*4), linetype="dashed", 
                          color = "red", size=0.5); plot

plot <- plot + annotate(geom="text", x=0, y=(Q_24M/2), label="Q1",
                          color="red"); plot
plot <- plot + annotate(geom="text", x=nrow(df.12M)-10, y=nrow(df.24M), label="Q5",
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
plot <- plot + geom_hline(yintercept=(Q_24M*3), linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_hline(yintercept=(Q_24M*4), linetype="dashed", 
                          color = "red", size=0.5); plot

plot <- plot + geom_vline(xintercept=Q_18M, linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_vline(xintercept=(Q_18M*2), linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_vline(xintercept=(Q_18M*3), linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_vline(xintercept=(Q_18M*4), linetype="dashed", 
                          color = "red", size=0.5); plot

plot <- plot + annotate(geom="text", x=0, y=(Q_24M/2), label="Q1",
                        color="red"); plot
plot <- plot + annotate(geom="text", x=nrow(df.18M)-10, y=nrow(df.24M), label="Q5",
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
plot <- plot + geom_hline(yintercept=(Q_18M*3), linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_hline(yintercept=(Q_18M*4), linetype="dashed", 
                          color = "red", size=0.5); plot

plot <- plot + geom_vline(xintercept=Q_12M, linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_vline(xintercept=(Q_12M*2), linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_vline(xintercept=(Q_12M*3), linetype="dashed", 
                          color = "red", size=0.5); plot
plot <- plot + geom_vline(xintercept=(Q_12M*4), linetype="dashed", 
                          color = "red", size=0.5); plot

plot <- plot + annotate(geom="text", x=0, y=(Q_18M/2), label="Q1",
                        color="red"); plot
plot <- plot + annotate(geom="text", x=nrow(df.12M)-10, y=nrow(df.18M), label="Q5",
                        color="red"); plot
plotList[[3]] <- plot


plotFileName <- "weight_loss_ranking_plots.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 5, height = 5)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]],
               ncol = 1, nrow = 1)
}
dev.off()
