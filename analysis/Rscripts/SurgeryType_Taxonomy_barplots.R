#Author: Alicia Sorgen
#Date: 2022 Feb 09
#Description: Generate basic taxonomic relative abundance graphs

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/git/gut-microbiome-bariatric-weight-outcomes")
params <- c(params, "MetaPhlAn2")
# params <- c(params, 0)
# params <- c(params, 1)
# params <- c(params, 6)
# params <- c(params, 12)
# params <- c(params, 18)
# params <- c(params, 24)
# end <- params[length(params)]

levels <- vector()
levels <- c(levels, "phylum")
# levels <- c(levels, "class")
# levels <- c(levels, "order")
# levels <- c(levels, "family")
levels <- c(levels, "genus")
# levels <- c(levels, "species")

moduleRoot <- "SurgeryType_Taxonomy_barplots"

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
library(ggplot2); message("ggplot2:", packageVersion("ggplot2"))
library(gridExtra); message("gridExtra:", packageVersion("gridExtra"))
library(scales); message("scales:", packageVersion("scales"))
library(dplyr); message("dplyr:", packageVersion("dplyr"))

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
    module <- paste0(args[2], "_", module)
  } 
  
  if (exists("end_month") == TRUE) {
    module <- paste0(module, "_BLto", end_month)
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
classifier <- args[2]
# included <- args[3:length(args)]

prevModule <- paste0(classifier, "_TaxaMetaMerge")
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")

logCountFile <- paste0("_LogNormalizedCounts_", classifier, ".tsv")
relCountFile <- paste0("_RelativeAbundanceCounts_", classifier, ".tsv")
rawCountFile <- paste0("_rawCounts_", classifier, ".tsv")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Set up script variables #####
palette=c("blue", "red","cyan", "lightgreen","orange","yellow","grey","pink","purple", "green4", "hotpink", "goldenrod")



###### Working script #####
plotList <- list()
index <- 1

for (level in levels) {
  
  message(paste0("\n", level, "\n"))
  
  # outputLevel <- paste0(output, level, "/")
  # dir.create(outputLevel, showWarnings = FALSE)
  
  logCounts<-read.table(paste0(inputDir,level,logCountFile),sep="\t",header = TRUE, row.names = 1,check.names = FALSE)
  relCounts<-read.table(paste0(inputDir,level,relCountFile),sep="\t",header = TRUE, row.names = 1,check.names = FALSE)
  rawCounts<-read.table(paste0(inputDir,level,rawCountFile),sep="\t",header = TRUE, row.names = 1,check.names = FALSE)
  startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1
  endMetadataIndex <- which(colnames(logCounts)=="ResponderStatus")
  
  metaData <- relCounts[,1:endMetadataIndex]
  taxaData <- relCounts[,startAbundanceIndex:ncol(relCounts)]
  taxaData[is.na(taxaData)] <- 0

  allNames <- colnames(taxaData)
  otherNames <- vector()
  j <- 1
  for (i in 1:ncol(taxaData)) {


    if (max(taxaData[i], na.rm = TRUE) < 1) {
      otherNames[j] <- colnames(taxaData)[i]
      j<-j+1
    }

  }

  otherNames <- c(otherNames, colnames(relCounts)[which(colnames(relCounts) %like% "Virus")])
  otherNames <- c(otherNames, colnames(relCounts)[which(colnames(relCounts) %like% "noname")])
  otherNames <- c(otherNames, colnames(relCounts)[which(colnames(relCounts) %like% "unclassified")])
  otherNames <- c(otherNames, colnames(relCounts)[which(colnames(relCounts) %like% "virus")])
  
  
  ok_taxa <- allNames[!allNames %in% otherNames]
  rawCountsEdit <- taxaData[ , ok_taxa]
  ok_taxa_counts <- as.data.frame(colSums(rawCountsEdit))
  ok_taxa_counts$Taxa <- rownames(ok_taxa_counts)
  ok_taxa_counts <- ok_taxa_counts[order(ok_taxa_counts[,1], decreasing = TRUE),]
  
  top <- 10
  topNames <- ok_taxa_counts$Taxa[1:top]
  topNames <- na.omit(topNames)
  topTaxa <- taxaData[,topNames]
  bottomNames <- ok_taxa_counts$Taxa[top+1:nrow(ok_taxa_counts)]
  otherNames <- c(otherNames, bottomNames)
  otherNames <- na.omit(otherNames)
  rawCountsOther <- taxaData[,otherNames]
  Other <- rowSums(rawCountsOther, na.rm = TRUE)
  
  rawCounts2 <- cbind(metaData, topTaxa)
  rawCounts2 <- cbind(rawCounts2, Other)
  
  #### Gather table for plot prep ####
  df <- rawCounts2 %>%
    gather("Taxa", "raw", startAbundanceIndex:ncol(rawCounts2))
  df <- df[df$Surgery %in% c("RYGB", "SG"),]
  aggregate(df$raw, by=list(Surgery=df$Surgery), FUN=sum)
  
  df2 <- df %>%
    group_by(Surgery, Timepoint, Taxa) %>%
    summarise(counts = mean(raw))
  
  #### Edit timepoint names ####
  df2$Timepoint <- gsub(pattern = "ONE", replacement = "Month 01", df2$Timepoint)
  df2$Timepoint <- gsub(pattern = "SIX", replacement = "Month 06", df2$Timepoint)
  df2$Timepoint <- gsub(pattern = "TWELVE", replacement = "Month 12", df2$Timepoint)
  df2$Timepoint <- gsub(pattern = "EIGHTEEN", replacement = "Month 18", df2$Timepoint)
  df2$Timepoint <- gsub(pattern = "TWENTY_FOUR", replacement = "Month 24", df2$Timepoint)
  
  ### Change facet labels
  variable_names <- list(
    "BL" = "Baseline" ,
    "Month 01" = "1 month post-op",
    "Month 06" = "6 months post-op",
    "Month 12" = "12 months post-op", 
    "Month 18" = "18 months post-op",
    "Month 24" = "24 months post-op"
  )

  variable_labeller <- function(variable,value){
    return(variable_names[value])
  }
  
  #### Generate plot ####
  plot <- ggplot(df2, aes(x = Surgery, y = counts, fill = Taxa))+
    geom_col(position="fill")+
    scale_y_continuous(labels = scales::percent)+
    facet_wrap(~Timepoint,  ncol=3, labeller = variable_labeller)+
    scale_fill_manual(values = palette); plot
  
  Level <- str_to_title(level)
  title.lab <- paste0(Level, " - Top ", top, " taxa")
  
  #### Change plot labels
  plot <- plot + labs(title = title.lab,
              x = "Type of Surgery",
              y = "Relative Abundance"); plot
  
  #### Change legend title
  plot <- plot + guides(fill=guide_legend(title = Level)); plot
  
  #### Increase x axis labels
  plot <- plot + theme(axis.text=element_text(size=10)); plot

  plot <- plot + theme(strip.text.x = element_text(size = 10)); plot
  
  plotList[[index]] <- plot
  index=index+1
  

} # for (level in levels)

file.path <- paste0(outputDir, "SurgeryType_taxonomy_boxplots.pdf")
pdf(file.path, width = 10, height = 5)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
}
dev.off()
