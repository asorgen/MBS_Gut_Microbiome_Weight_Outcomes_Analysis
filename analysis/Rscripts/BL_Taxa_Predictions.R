#Author: Alicia Sorgen
#Date: 2022 July 11
#Description: Modeling to determine if baseline taxa predicts weight outcomes.

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
# end <- params[length(params)]

levels <- vector()
levels <- c(levels, "phylum")
# levels <- c(levels, "class")
# levels <- c(levels, "order")
# levels <- c(levels, "family")
levels <- c(levels, "genus")
# levels <- c(levels, "species")

moduleRoot <- "BL_Taxa_Predictions"

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
  
  if (exists("end") == TRUE) {
    module <- paste0(module, "_BLto", end)
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

##### Analysis #####
for (level in levels) {
  
  message(paste0("\n***** Starting ", level, " *****\n"))
  
  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  
  file.path <- paste0(inputDir,level, logCountFile)
  logCounts<-read.table(file.path, sep="\t", header = TRUE, row.names = 1, check.names = FALSE)
  
  plotList <- list()
  dFrame <- data.frame()
  index <- 1
  
  for (j in 1:length(included)) {
    
    # Month X rows
    xM = which(logCounts$time == included[j])
    
    startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1
    
    # Parse out metadata from counts
    myT<-logCounts[,startAbundanceIndex:ncol(logCounts)]
    
    PatientID <- logCounts$PatientID
    Timepoint <- logCounts$time
    Percent_Loss_kg <- logCounts$Percent_Loss_kg
    
    for (i in 1:ncol(myT)) {
      
      bug <- myT[,i]
      
      if (mean(bug>0, na.rm = TRUE)>0.1){
        
        df <- data.frame(bug, PatientID, Timepoint, Percent_Loss_kg)
        
        # Month X abundance
        xM_bug <- df[xM, "bug"]
        names(xM_bug) <- df[xM, "PatientID"]
        
        # Add info to each row
        df$Bug_at_xM <- xM_bug[df$PatientID]
        
        df.x <- df[!(df$Timepoint == 0),]
        # df.x <- df
        df.x <- na.omit(df.x)
        stat.df <- df.x %>%
          group_by(Timepoint) %>%
          cor_test(Bug_at_xM, Percent_Loss_kg, method = "kendall") %>%
          adjust_pvalue(method = "BH")
        stat.df$var1 <- paste0(month.labs[j], " ", colnames(myT)[i], " abundance")
        stat.df$var2 <- paste0("Weight loss (%): BL to ", stat.df$Timepoint, " month(s)")
        stat.df$n <- table(df.x$Timepoint)
        dFrame <- rbind(dFrame, stat.df)
        
        
        stat.df <- df.x %>%
          group_by(Timepoint) %>%
          cor_test(Bug_at_xM, Percent_Loss_kg, method = "spearman") %>%
          adjust_pvalue(method = "BH")
        stat.df$var1 <- paste0(month.labs[j], " ", colnames(myT)[i], " abundance")
        stat.df$var2 <- paste0("Weight loss (%): BL to ", stat.df$Timepoint, " month(s)")
        stat.df$n <- table(df.x$Timepoint)
        dFrame <- rbind(dFrame, stat.df)
        
        
        for (month in included) {
          
          if (month != 0) {
            
            df2 <- df[df$Timepoint == month,]
            
            plot <- ggscatter(df2, x = "Bug_at_xM", y = "Percent_Loss_kg",
                              shape = 21, size = 2.5, # Points shape and size
                              add = "reg.line",  # Add regression line
                              conf.int = TRUE, # Add confidence interval
                              cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                              cor.coeff.args = list(method = "kendall"),
                              add.params = list(fill = "lightgray")
            ); plot
            
            x.lab <- paste0(month.labs[j], " ", colnames(myT)[i], " abundance")
            y.lab <- paste0("Weight loss (%): BL to ", month, " month(s)")
            plot <- plot + labs(x=x.lab, y = y.lab); plot
            
            plotList[[index]] <- plot
            index <- index + 1
            
            
            
          } # if (month != 0)
          
        } # for (month in included)
        
      } # if (mean(bug>0, na.rm = TRUE)>0.1)
      
    } # for (i in 1:ncol(myT))
    
  } # for (i in 1:length(included))
  
  FileName <- paste0(level, "_ken_spear_predictions_BLto", included[length(included)], "M_v2.tsv")
  file.path <- paste0(outputLevel, FileName)
  write.table(dFrame, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  
  plotFileName <- paste0(level, "_kendall_prediction_plots_BLto", included[length(included)], "M.pdf")
  file.path <- paste0(outputLevel, plotFileName)
  pdf(file.path, width = 5, height = 5)
  for (i in 1:length(plotList)) {
    grid.arrange(plotList[[i]],
                 ncol = 1, nrow = 1)
  }
  dev.off()
  
} # for (level in levels)
