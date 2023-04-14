#Author: Alicia Sorgen
#Date: 2022 June 09
#Description: Linear modeling for taxa by time

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
params <- c(params, "MetaPhlAn2")
params <- c(params, "Tertile")
params <- c(params, 0)
params <- c(params, 1)
params <- c(params, 6)
params <- c(params, 12)
params <- c(params, 18)
params <- c(params, 24)
end_month <- params[length(params)]

levels <- vector()
levels <- c(levels, "phylum")
# levels <- c(levels, "class")
# levels <- c(levels, "order")
# levels <- c(levels, "family")
levels <- c(levels, "genus")
# levels <- c(levels, "species")

moduleRoot <- "pValuePlots"

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
library(ggrepel); message("ggrepel: Version ", packageVersion("ggrepel"))

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
  rm(gitRoot)
  
  root <- paste0("~/BioLockJ_pipelines/")
  pipeline <- paste0(ANALYSIS,"_analysis_")
  
  if (any(dir(root) %like% pipeline) == TRUE) {
    root <- paste0(root,"/",str_subset(dir(root), pipeline), "/")
  } else {
    today <- as.character(format(Sys.Date(), "%Y%b%d"))
    root <- paste0(root,ANALYSIS,"_analysis_", today, "/")
    dir.create(root, showWarnings = FALSE)
  }
  rm(pipeline)
  
  if (any(dir(root) == "input") == FALSE) {
    file.copy(gitInput,
              root,
              recursive = TRUE)
    
  }
  rm(gitInput)
  module <- moduleRoot
  
  if (args[2] == "MetaPhlAn2" | args[2] == "Kraken2") {
    module <- paste0(module, "_", args[2])
  } 
  
  if (args[3] %in% c("Quintiles", "Quartiles", "Thirds", "Half")) {
    module <- paste0(module, "_", args[3])
  } 
  
  if (exists("end_month") == TRUE) {
    module <- paste0(module, "_BLto", end_month, "M")
  }
  
  if (any(dir(root) %like% module) == TRUE) {
    moduleDir <- paste0(root,str_subset(dir(root), module), "/")
  } else {
    moduleDir <- paste0(root, module, "/")
    dir.create(moduleDir, showWarnings = FALSE)
  }
  rm(module, root)
  
  scriptDir <- paste0(moduleDir, "script/")
  if (any(dir(moduleDir) == "script") == FALSE) {
    dir.create(scriptDir, showWarnings = FALSE)
  }
  rm(scriptDir)
  
  outputDir <- paste0(moduleDir, "output/")
  if (any(dir(moduleDir) == "output") == FALSE) {
    dir.create(outputDir, showWarnings = FALSE)
  }
  
  files <- list.files(outputDir, recursive = TRUE, full.names = TRUE)
  file.remove(files)
  rm(outputDir, files)
  
  resourcesDir <- paste0(moduleDir, "resources/")
  if (any(dir(moduleDir) == "resources") == FALSE) {
    dir.create(resourcesDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/functions.R")
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
    
    script = paste0(gitScripts,"/", moduleRoot, ".R"); script
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
    rm(script)
    
  }
  rm(resourcesDir, gitScripts)
  
  setwd(paste0(moduleDir, "script/"))
  rm(moduleDir)
  
}
rm(params, moduleRoot, ANALYSIS, end_month)

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
classifier <- args[2]
division <- args[3]
included <- args[4:length(args)]

prevModule <- paste0(classifier, "_Taxa_over_time_BLto", included[length(included)], "M")
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### p value box plot #####
for (level in levels) {
  
  message(paste0("\n***** Starting ", level, " ", division, " p value boxplots *****\n"))
  
  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  
  inputLevel <- paste0(inputDir, level, "/", division, "/")
  file.path <- paste0(inputLevel, level, "_by_Timepoint_", division, "_LM_pValues_BLto", included[length(included)], "M_", classifier, ".tsv")
  qFrame <- read.table(file.path, sep="\t",header = TRUE, check.names = FALSE)
  
  qFrame$log <- ifelse(qFrame$direction == "increase", -log10(qFrame$adjPValues),
                       log10(qFrame$adjPValues))
  qFrame <- na.omit(qFrame)
  
  # Label taxa as significant or not based on adjusted p-value
  qFrame$Significance <- ifelse(qFrame$adjPValues < 0.05, "Significant",
                                "Not Significant")
  
  f.test <- fisher.test(qFrame$groupNumbers, qFrame$Significance, workspace = 2e8)
  f.pVal <- roundP(f.test$p.value)
  
  title.lab <- paste0(classifier, " ", str_to_title(level), " level taxonomic changes over time (BL to ", included[length(included)], " months)")
  subtitle.lab <- paste0("Fisher's Exact Test (", f.pVal, ")")
  
  hLine1 <- -log10(0.05)
  hLine2 <- log10(0.05)
  
  labels1 <- getWeightGroupLabels1(division, qFrame, hLine1, hLine2)
  labels2 <- getWeightGroupLabels2(division, qFrame, hLine1, hLine2)
  quants1 <- getWeightGroupLocation1(division, qFrame, hLine1, hLine2)
  quants2 <- getWeightGroupLocation2(division, qFrame, hLine1, hLine2)
  
  
  ## Generate plot ##
  
  plot <- ggboxplot(
    qFrame, x = "groupNumbers", y = "log",
    add = "jitter",
    palette =c("#018571", "#4DAC26", "#5E3C99", "#0571B0", "#D01C8B"),
    fill = "groupNumbers",
    # color = "groupNumbers",
    # facet.by = "groupNumbers",
    scales = "free",
    # label="qbugName",
    outlier.shape = 1
  ); plot
  
  # Add title & subtitle
  plot <- plot +
    labs(title = title.lab,
         subtitle = subtitle.lab); plot
  
  # Edit x- and y-axis labels
  plot <- plot +
    labs(x = "Weight Loss Group", y = "p value (log10)"); plot
  
  # Add dashed line for significant p values increasing
  plot <- plot + geom_hline(yintercept=hLine1, linetype="dashed", color = "red"); plot
  
  # Add dashed line for significant p values decreasing
  plot <- plot + geom_hline(yintercept=hLine2, linetype="dashed", color = "red"); plot
  
  # Increase text size of labels
  plot <- plot +
    theme(text = element_text(size = 15)); plot
  
  # Remove legend
  plot <- plot + theme(legend.position = "none"); plot
  
  # Add text indicating the number of significant taxa increasing over time
  plot <- plot + geom_text(data = quants1, aes(x = groupNumbers, y = quant, label = labels1), size = 5); plot
  
  # Add text indicating the number of significant taxa decreasing over time
  plot <- plot + geom_text(data = quants2, aes(x = groupNumbers, y = quant, label = labels2), size = 5); plot
  
  # Add taxa names showing significant change 
  plot <- plot + geom_text_repel(min.segment.length = Inf, aes(label=ifelse(log>hLine1 | log < hLine2,as.character(qbugName),'')),
                                 # box.padding   = 0.35,
                                 # point.padding = 0.5,
                                 force = 0.5,
                                 nudge_x = 0.25,
                                 size = 2.5,
                                 direction = "y", 
                                 hjust = 0,
                                 segment.color = 'grey50'); plot
  
  
  # Add up arrow and annotation
  plot <- plot + annotate("segment", x=0.05, y=0.25, xend=0.05, yend=2,
                          col="black", arrow=arrow(length=unit(0.3, "cm"))); plot
  plot <- plot + annotate("text", x=0.1, y=1.1, hjust = 0, label = "Increases in abundance"); plot
  
  # Add down arrow and annotation
  plot <- plot + annotate("segment", x=0.05, y=-0.25, xend=0.05, yend=-2,
                          col="black", arrow=arrow(length=unit(0.3, "cm"))); plot
  plot <- plot + annotate("text", x=0.1, y=-1.1, hjust = 0, label = "Decreases in abundance"); plot
  
  # Print the plot
  file.path <- paste0(outputLevel, level, paste0("_by_Timepoint_", division, "_LM_pValue_boxplot_BLto", included[length(included)], "M_", classifier, ".pdf"))
  pdf(file.path, width = 14, height = 7)
  print(plot)
  dev.off()
  
} # for (level in levels)







##### p value scatter plot #####
for (level in levels) {
  
  message(paste0("\n***** Starting ", level, " ", division, " p value scatter plots *****\n"))
  
  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  
  inputLevel <- paste0(inputDir, level, "/", division, "/")
  file.path <- paste0(inputLevel, level, "_by_Timepoint_", division, "_LM_pValues_BLto", included[length(included)], "M_", classifier, ".tsv")
  qFrame <- read.table(file.path, sep="\t",header = TRUE, check.names = FALSE)
  
  qFrame$log <- ifelse(qFrame$direction == "increase", -log10(qFrame$adjPValues),
                       log10(qFrame$adjPValues))
  
  title.lab <- paste0(classifier, " ", str_to_title(level), " level taxonomic changes over time (BL to ", included[length(included)], " months)")
  subtitle.lab <- paste0("Fisher's Exact Test (", f.pVal, ")")
  
  hLine1 <- -log10(0.05)
  hLine2 <- log10(0.05)
  
  labels1 <- getWeightGroupLabels1(division, qFrame, hLine1, hLine2)
  labels2 <- getWeightGroupLabels2(division, qFrame, hLine1, hLine2)
  quants1 <- getWeightGroupLocation1(division, qFrame, hLine1, hLine2)
  quants2 <- getWeightGroupLocation2(division, qFrame, hLine1, hLine2)
  
  
  ## Generate plot ##
  
  pos <- position_jitter(width = 0.2, seed = 2)
  plot<- ggplot(qFrame, aes(x= groupNumbers, y= log,
                            add = "jitter",
                            # palette =c("#018571", "#4DAC26", "#5E3C99", "#0571B0", "#D01C8B"),
                            color = groupNumbers,
                            label=ifelse(log>hLine1 | log < hLine2,as.character(qbugName),''))) +
    geom_point(position = pos) +
    scale_color_gradient(low="blue", high="red"); plot
  
  if (division == "Quintiles") {
    # plot <- plot +     scale_x_continuous(breaks = 1:5, labels = c("Quintile 1", "Quintile 2", "Quintile 3", "Quintile 4", "Quintile 5"),
    #                                       expand = expansion(mult = 0.1))
    plot <- plot +     scale_x_continuous(limits = c(0, 6), breaks = 0:6, labels = c("", "Quintile 1", "Quintile 2", "Quintile 3", "Quintile 4", "Quintile 5", "")); plot
  }
  if (division == "Quartiles") {
    # plot <- plot +     scale_x_continuous(breaks = 1:4, labels = c("Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4"),
    #                                       expand = expansion(mult = 0.1))
    plot <- plot +     scale_x_continuous(limits = c(0, 5), breaks = 0:5, labels = c("", "Quartile 1", "Quartile 2", "Quartile 3", "Quartile 4", ""))
  }
  if (division == "Thirds") {
    # plot <- plot +     scale_x_continuous(breaks = 1:3, labels = c("Bottom 33%", "Middle 33%", "Top 33%"),
    #                                       expand = expansion(mult = 0.05))
    plot <- plot +     scale_x_continuous(limits = c(0, 4), breaks = 0:4, labels = c("", "Bottom 33%", "Middle 33%", "Top 33%", ""))
  }
  if (division == "Half") {
    # plot <- plot +     scale_x_continuous(breaks = 1:2, labels = c("Bottom 50%", "Top 50%"),
    #                                       expand = expansion(mult = 0.1))
    plot <- plot +     scale_x_continuous(limits = c(0, 3), breaks = 0:3, labels = c("", "Bottom 50%", "Top 50%", ""))
  }
  
  
  # Add title & subtitle
  plot <- plot +
    labs(title = title.lab,
         subtitle = subtitle.lab); plot
  
  # Edit x- and y-axis labels
  plot <- plot +
    labs(x = "Weight Loss Group", y = "p value (log10)"); plot
  
  # Add dashed line for significant p values increasing
  plot <- plot + geom_hline(yintercept=hLine1, linetype="dashed", color = "red"); plot
  
  # Add dashed line for significant p values decreasing
  plot <- plot + geom_hline(yintercept=hLine2, linetype="dashed", color = "red"); plot
  
  # Increase text size of labels
  plot <- plot +
    theme(text = element_text(size = 15)); plot
  
  # Remove legend
  plot <- plot + theme(legend.position = "none"); plot
  
  # Add text indicating the number of significant taxa increasing over time
  plot <- plot + geom_text(data = quants1, aes(x = groupNumbers, y = quant, label = labels1), size = 5); plot
  
  # Add text indicating the number of significant taxa decreasing over time
  plot <- plot + geom_text(data = quants2, aes(x = groupNumbers, y = quant, label = labels2), size = 5); plot
  
  # Add taxa names showing significant change 
  plot <- plot + geom_text_repel(min.segment.length = Inf, aes(label=ifelse(log>hLine1 | log < hLine2,as.character(qbugName),'')),
                                 # box.padding   = 0.35,
                                 # point.padding = 0.5,
                                 force = 0.5,
                                 nudge_x = 0.25,
                                 size = 2.5,
                                 direction = "y", 
                                 hjust = 0,
                                 segment.color = 'grey50'); plot
  
  
  # Add up arrow and annotation
  plot <- plot + annotate("segment", x=0.05, y=0.25, xend=0.05, yend=2,
                          col="black", arrow=arrow(length=unit(0.3, "cm"))); plot
  plot <- plot + annotate("text", x=0.1, y=1.1, hjust = 0, label = "Increases in abundance"); plot
  
  # Add down arrow and annotation
  plot <- plot + annotate("segment", x=0.05, y=-0.25, xend=0.05, yend=-2,
                          col="black", arrow=arrow(length=unit(0.3, "cm"))); plot
  plot <- plot + annotate("text", x=0.1, y=-1.1, hjust = 0, label = "Decreases in abundance"); plot
  
  # Print the plot
  file.path <- paste0(outputLevel, level, paste0("_by_Timepoint_", division, "_LM_pValue_scatterplot_BLto", included[length(included)], "M_", classifier, ".pdf"))
  pdf(file.path, width = 14, height = 7)
  print(plot)
  dev.off()
  
} # for (level in levels)



