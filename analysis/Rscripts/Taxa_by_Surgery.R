#Author: Alicia Sorgen
#Date: 2022 August 24
#Description: Modeling to compare taxonomic differences between surgery type at each timepoint.

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/git/gut-microbiome-bariatric-weight-outcomes")
params <- c(params, "MetaPhlAn2")
params <- c(params, 0)
params <- c(params, 1)
params <- c(params, 6)
params <- c(params, 12)
params <- c(params, 18)
params <- c(params, 24)
# endM <- params[length(params)]

levels <- vector()
levels <- c(levels, "phylum")
# levels <- c(levels, "class")
# levels <- c(levels, "order")
# levels <- c(levels, "family")
levels <- c(levels, "genus")
# levels <- c(levels, "species")

moduleRoot <- "Taxa_by_Surgery"

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
    module <- paste0(args[2], "_", module)
  } 
  
  if (exists("endM") == TRUE) {
    module <- paste0(module, "_BLto", end, "M")
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
included <- args[3:length(args)]

prevModule <- paste0(classifier, "_TaxaMetaMerge")
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")

logCountFile <- paste0("_LogNormalizedCounts_", classifier, ".tsv")

month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")


##### Taxa by Surgery - Wilcoxon (grouped by Timepoint) #####
for (level in levels) {
  
  message(paste0("\n***** Starting ", level, " analysis *****\n"))
  
  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  
  file.path <- paste0(inputDir,level, logCountFile)
  logCounts<-read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
  # logCounts <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
  
  startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1
  
  # Parse out metadata from counts
  myT<-logCounts[,startAbundanceIndex:ncol(logCounts)]
  
  PatientID <- logCounts$PatientID
  Timepoint <- logCounts$time
  Surgery <- logCounts$Surgery

  dFrame <- data.frame()

  
  for (i in 1:ncol(myT)){
    
    bug<-myT[,i]
    
    if ( mean(bug > 0, na.rm = TRUE) > 0.1 ){
      
      df<-data.frame(bug,PatientID, Timepoint, Surgery)
      df <- df[df$Timepoint %in% included,]
      df <- na.omit(df)
      df$Timepoint <- as.factor(df$Timepoint)

      stats <- df %>%
        group_by(Timepoint) %>%
        wilcox_test(bug ~ Surgery) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance()
      
      stats$.y.<-colnames(myT)[i]
      
      averages <- df %>%
        group_by(Timepoint, Surgery) %>%
        get_summary_stats(bug, type = "mean_sd")
      # averages$variable<-colnames(myT)[i]
      
      averages$Mean_SD <- paste(round(averages$mean, 2), "±", round(averages$sd, 2))
      
      stats$RYGB_mean <- averages$mean[which(averages$Surgery == "RYGB")]
      stats$SG_mean <- averages$mean[which(averages$Surgery == "SG")]
      
      stats$RYGB_sd <- averages$sd[which(averages$Surgery == "RYGB")]
      stats$SG_sd <- averages$sd[which(averages$Surgery == "SG")]
      
      stats$RYGB_Mean_SD <- averages$Mean_SD[which(averages$Surgery == "RYGB")]
      stats$SG_Mean_SD <- averages$Mean_SD[which(averages$Surgery == "SG")]
      
      stats$HigherAbundance <- ifelse(stats$RYGB_mean > stats$SG_mean, "RYGB",
                                      ifelse(stats$RYGB_mean < stats$SG_mean, "SG",
                                             "same"))
      
      dFrame <- rbind(dFrame, stats)
      
    } #if (mean(bug>0, na.rm = TRUE)>0.1)
    
  } #for (i in 1:ncol(myT))
  
  
  pVal.FileName <- paste0("_by_Surgery_wilcox_", classifier, ".tsv")
  file.path <- paste0(outputLevel, level, pVal.FileName)
  write.table(dFrame, file.path,sep="\t",quote = FALSE, row.names = FALSE)
  
  
  
  if (sum(dFrame$p.adj < 0.05) > 0) {
    
    dFrame <- dFrame[order(dFrame$p.adj),]
    
    file.path <- paste0(outputLevel, "Significant_", level, "_by_Surgery_wilcox_", classifier, ".pdf")
    pdf(file.path)
    
    par(mfrow=c(2,2))
    
    for (i in 1:nrow(dFrame)) {
      
      if (dFrame$p.adj[i] < 0.05) {
        
        taxa <- dFrame$.y.[i]
        timepoint <- dFrame$Timepoint[i]
        
        Timepoint <- logCounts$time
        Surgery <- logCounts$Surgery
        bugCol <- logCounts[,which(colnames(logCounts) == taxa)]
        
        df <- data.frame(Timepoint, Surgery, bugCol)
        df <- df[df$Timepoint == timepoint, ]
        
        title.lab <- paste0(taxa, "\n", 
                            month.labs[which(included == timepoint)], " (adj p = ", format(dFrame$p.adj[i], length = 3), ")")
        
        
        boxplot( df$bugCol ~ df$Surgery,
                 xlab = "Surgery Type",
                 ylab = taxa,
                 main = title.lab,
                 cex.main = 1
          )


      }
      
    }
    
    dev.off()
    
  }
  
  
} # for (level in levels)

