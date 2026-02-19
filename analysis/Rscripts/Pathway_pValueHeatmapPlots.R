#Author: Alicia Sorgen
#Date: 2023 Jan 26
#Description: Create heatmap of taxa changes over time by tertile

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/git/gut-microbiome-bariatric-weight-outcomes")
params <- c(params, "HumanN2")
params <- c(params, "Tertile")
# params <- c(params, 0)
# params <- c(params, 1)
# params <- c(params, 6)
# params <- c(params, 12)
# params <- c(params, 18)
# params <- c(params, 24)
# end_month <- params[length(params)]


moduleRoot <- "Pathway_pValueHeatmapPlots"

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
    # module <- paste0(args[2], "_", module)
    module <- gsub("Taxa", args[2], module)
  } 
  
  if (args[3] %in% c("Quintiles", "Quartiles", "Thirds", "Half")) {
    module <- paste0(module, "_", args[3])
  } 
  
  if (exists("end_month") == TRUE) {
    module <- paste0(module, "_BLto", end_month, "M")
    rm(end_month)
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
rm(params, moduleRoot, ANALYSIS)

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
classifier <- args[2]; message("Classifier: ", classifier)
division <- args[3]; message("Division: ", division)
# included <- args[4:length(args)]

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Set up script variables #####
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

##### Read in files #####
months <- c(1,6,12,18,24)
message("\n>>>>> Starting pathway analysis <<<<<\n")
merged.df <- data.frame()

for (month in months) {
  
  message("\n>>>>> Starting ", month, " merge")
  prevModule <- paste0(classifier, "_Pathway_over_time_BLto", month, "M")
  inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
  
  inputLevel <- paste0(inputDir, division, "/")
  file.path <- paste0(inputLevel, "Pathway_by_Timepoint_", division, "_LM_pValues_BLto", month, "M_", classifier, ".tsv")
  qFrame <- read.table(file.path, sep="\t",header = TRUE, check.names = FALSE)
  qFrame$Time <- paste0("BL - ", month, "M")
  
  merged.df <- rbind(merged.df, qFrame)
} # for (month in months)

sigPathway <- unique(merged.df$qpathName[which(merged.df$adjPValues < 0.05)])

df <- merged.df[merged.df$qpathName %in% sigPathway,]

plotList <- list()
index <- 1
for (month in months) {
  
  message("\n>>>>> Starting ", month, " heatmap")
  
  df2 <- df[df$Time == paste0("BL - ", month, "M"),]
  df2$groupNumbers <- gsub("1", "Tertile 1", df2$groupNumbers)
  df2$groupNumbers <- gsub("2", "Tertile 2", df2$groupNumbers)
  df2$groupNumbers <- gsub("3", "Tertile 3", df2$groupNumbers)
  df2$groupNumbers <- as.factor(df2$groupNumbers)
  
  df2$logP <- ifelse(df2$direction == "increase", -log10(df2$adjPValues),
                     log(df2$adjPValues))

  breakMarks <- c(-100, log10(0.001), log10(0.01), log10(0.05),
                  -log10(0.05), -log10(0.01), -log10(0.001), 100)
  labelNames <- c("Decrease p < 0.001", " Decrease p < 0.01", "Decrease p < 0.05", "p > 0.05", "Increase p < 0.05", "Increase p < 0.01", "Increase p < 0.001")
  
  
  df3 <- df2 %>%
    mutate(qpathName=factor(qpathName, levels=rev(sort(unique(qpathName))))) %>%
    mutate(pValfactor = cut(logP, breaks = breakMarks,
                            labels = labelNames)) %>%
    mutate(pValfactor = factor(as.character(pValfactor), levels = rev(levels(pValfactor))))
  colPalette <- c("#66cc66", "#90ee90", "#ccffcc", "#fffafa", "#ffcccc", "#ff9999", "#ff3300")
  
  xCols <- which(rev(labelNames) %in% unique(df3$pValfactor))
  title.lab <- paste0("BL - ", month, "M")
  
  # Heatmap 
  plot <- ggplot(df3, aes(groupNumbers, qpathName, fill= pValfactor)) + 
    geom_tile(color = "black")+
    # geom_tile(data = df2[c(signif),], color = "black", size = 2)+
    scale_y_discrete(expand=c(0, 0)) +
    # scale_fill_distiller(palette = "RdPu"); plot
    # scale_fill_gradient2(low="red", mid = "white", high="green", limits = c(-lim, lim))+
    scale_fill_manual(values=colPalette[xCols], na.value = "grey90"); plot
  
  x.lab <- "Weight Loss Group"
  y.lab <- "Pathway"
  plot <- plot + labs(x=x.lab, y = y.lab, title = title.lab); plot
  
  plot <- plot + labs(fill = "Abundance Change"); plot
  
  plot <- plot + theme(plot.background = element_blank(), # remove plot background
                       panel.border = element_blank() # remove plot border
  ); plot
  
  plotList[[index]] <- plot
  index <- index + 1
} # for (month in months)


message("\n>>>>> Outputting plot")
# Output plot
file.path <- paste0(outputDir, "Pathway_Heatmap.pdf")
pdf(file.path, width = 10, height = 10)
for (i in 1:length(plotList)) {
  grid.arrange(
    # plotList[[1]], plotList[[2]], 
    plotList[[i]], 
    ncol = 1, nrow = 1
  )
  
}
dev.off()

