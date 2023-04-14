#Author: Alicia Sorgen
#Date: 2022 July 11
#Description: Modeling to determine if baseline taxa predicts weight outcomes.

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
params <- c(params, "MetaPhlAn2")
# params <- c(params, 0)
# params <- c(params, 1)
# params <- c(params, 6)
# params <- c(params, 12)
# params <- c(params, 18)
# params <- c(params, 24)
# endM <- params[length(params)]

levels <- vector()
levels <- c(levels, "phylum")
# levels <- c(levels, "class")
# levels <- c(levels, "order")
# levels <- c(levels, "family")
levels <- c(levels, "genus")
# levels <- c(levels, "species")

moduleRoot <- "Taxa_Changes_1M"

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
  
  if (exists("endM") == TRUE) {
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

prevModule <- paste0(classifier, "_TaxaMetaMerge")
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
logCountFile <- paste0("_LogNormalizedCounts_", classifier, ".tsv")
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Analysis #####
included <- c(0, 1, 6, 12, 18, 24)
for (level in levels) {
  
  message(paste0("\n***** Starting ", level, " *****\n"))
  
  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  
  file.path <- paste0(inputDir,level, logCountFile)
  logCounts<-read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
  
  startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1
  
  
  # Parse out metadata from counts
  myT<-logCounts[,startAbundanceIndex:ncol(logCounts)]
  
  PatientID <- logCounts$PatientID
  Timepoint <- logCounts$time
  Percent_Loss_kg <- logCounts$Percent_Loss_kg
  
  plotList <- list()
  dFrame <- data.frame()
  index <- 1
  
  for (i in 1:ncol(myT)) {
    
    bug <- myT[,i]
    
    if (mean(bug>0, na.rm = TRUE)>0.1){
      
      df <- data.frame(bug, PatientID, Timepoint)

      df <- df[(df$Timepoint %in% c(0, 1)),]
      df <- spread(df, Timepoint, bug)
      df <- na.omit(df)
      df$bugChange <- df$`1` - df$`0`
      df$bugPercentChange <- ( df$bugChange / df$`0` ) * 100
      df[is.na(df)] <- 0
      
      for (j in 2:length(included)) {
        
        month <- included[j]
        
        # month rows
        m1 = which(logCounts$time == month)
        # 1 month abundance
        PWL <- logCounts[m1, "Percent_Loss_kg"]
        names(PWL) <- logCounts[m1, "PatientID"]
        
        # Add PWL info to each row
        df$PWL <- PWL[df$PatientID]
        
        plot <- ggscatter(df, x = "bugPercentChange", y = "PWL",
                          shape = 21, size = 2.5, # Points shape and size
                          add = "reg.line",  # Add regression line
                          conf.int = TRUE, # Add confidence interval
                          # cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                          # cor.coeff.args = list(method = "kendall"),
                          add.params = list(fill = "lightgray")
        ); plot
        
        x.lab <- paste0(colnames(myT)[i], " change (%): BL to 1 month")
        y.lab <- paste0("Weight loss (%): BL to ", month, " month")
        plot <- plot + labs(x=x.lab, y = y.lab); plot
        
        plotList[[index]] <- plot
        index <- index + 1
        
        stat.df <- df %>%
          cor_test(PWL, bugPercentChange, method = "kendall") 
        stat.df$var1 <- y.lab
        stat.df$var2 <- x.lab
        
        stat.df2 <- df %>%
          cor_test(PWL, bugPercentChange, method = "spearman") 
        
        stat.df$Spearman_r <- stat.df2$cor
        dFrame <- rbind(dFrame, stat.df)
        
      }
      
      
      
      
    } # if (mean(bug>0, na.rm = TRUE)>0.1)

  } # for (i in 1:ncol(myT))

  plotFileName <- paste0(level, "_kendall_Taxa_Weight_change_plots_BLto1M.pdf")
  file.path <- paste0(outputLevel, plotFileName)
  pdf(file.path, width = 5, height = 5)
  for (i in 1:length(plotList)) {
    grid.arrange(plotList[[i]],
                 ncol = 1, nrow = 1)
  }
  dev.off()
  
  final <- data.frame()
  for (x in unique(dFrame$var1)) {
    
    dFrame.Filt <- dFrame[dFrame$var1 == x,]
    dFrame.Filt$pAdj <- p.adjust(dFrame.Filt$p, method = "BH")
    final <- rbind(final, dFrame.Filt)
    
  }
  
  FileName <- paste0(level, "_kendall_Taxa_Weight_change_results_BLto1M.tsv")
  file.path <- paste0(outputLevel, FileName)
  write.table(final, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
} # for (level in levels)


##### Analysis (grouped by Surgery) #####
for (level in levels) {
  
  message(paste0("\n***** Starting ", level, " *****\n"))
  
  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  
  file.path <- paste0(inputDir,level, logCountFile)
  logCounts<-read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
  
  plotList <- list()
  dFrame <- data.frame()
  index <- 1
  
  for (surgType in c("RYGB", "SG")) {
    
    SurgFilt.df <- logCounts[logCounts$Surgery == surgType,]
    
    startAbundanceIndex <- which(colnames(SurgFilt.df)=="ResponderStatus")+1
    
    
    # Parse out metadata from counts
    myT<-SurgFilt.df[,startAbundanceIndex:ncol(SurgFilt.df)]
    
    PatientID <- SurgFilt.df$PatientID
    Timepoint <- SurgFilt.df$time
    Percent_Loss_kg <- SurgFilt.df$Percent_Loss_kg
    
    
    for (i in 1:ncol(myT)) {
      
      bug <- myT[,i]
      
      if (mean(bug>0, na.rm = TRUE)>0.1){
        
        df <- data.frame(bug, PatientID, Timepoint)
        
        df <- df[(df$Timepoint %in% c(0, 1)),]
        df <- spread(df, Timepoint, bug)
        df <- na.omit(df)
        df$bugChange <- df$`1` - df$`0`
        df$bugPercentChange <- ( df$bugChange / df$`0` ) * 100
        df[is.na(df)] <- 0
        
        for (j in 2:length(included)) {
          
          month <- included[j]
          
          # month rows
          m1 = which(SurgFilt.df$time == month)
          # 1 month abundance
          PWL <- SurgFilt.df[m1, "Percent_Loss_kg"]
          names(PWL) <- SurgFilt.df[m1, "PatientID"]
          
          # Add PWL info to each row
          df$PWL <- PWL[df$PatientID]
          
          plot <- ggscatter(df, x = "bugPercentChange", y = "PWL",
                            shape = 21, size = 2.5, # Points shape and size
                            add = "reg.line",  # Add regression line
                            conf.int = TRUE, # Add confidence interval
                            # cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                            # cor.coeff.args = list(method = "kendall"),
                            add.params = list(fill = "lightgray")
          ); plot
          
          x.lab <- paste0(colnames(myT)[i], " change (%): BL to 1 month")
          y.lab <- paste0("Weight loss (%): BL to ", month, " month")
          plot <- plot + labs(x=x.lab, y = y.lab,
                              title = surgType); plot
          
          plotList[[index]] <- plot
          index <- index + 1
          
          stat.df <- df %>%
            cor_test(PWL, bugPercentChange, method = "kendall") 
          stat.df$var1 <- y.lab
          stat.df$var2 <- x.lab
          
          stat.df2 <- df %>%
            cor_test(PWL, bugPercentChange, method = "spearman") 
          
          stat.df$Spearman_r <- stat.df2$cor
          stat.df$Surgery <- surgType
          dFrame <- rbind(dFrame, stat.df)
          
        }
        
        
        
        
      } # if (mean(bug>0, na.rm = TRUE)>0.1)
      
    } # for (i in 1:ncol(myT))
    
    plotFileName <- paste0(level, "_kendall_Taxa_Weight_change_plots_BLto1M_by_Surgery.pdf")
    file.path <- paste0(outputLevel, plotFileName)
    pdf(file.path, width = 5, height = 5)
    for (i in 1:length(plotList)) {
      grid.arrange(plotList[[i]],
                   ncol = 1, nrow = 1)
    }
    dev.off()
    
    final <- data.frame()
    for (x in unique(dFrame$var1)) {
      
      dFrame.Filt <- dFrame[dFrame$var1 == x,]
      
      for (y in unique(dFrame.Filt$Surgery)) {
        
        dFrame.Filt2 <- dFrame.Filt[dFrame.Filt$Surgery == y,]
        dFrame.Filt2$pAdj <- p.adjust(dFrame.Filt2$p, method = "BH")
        final <- rbind(final, dFrame.Filt2)
        
      }
      
    }
    
    FileName <- paste0(level, "_kendall_Taxa_Weight_change_results_BLto1M_by_Surgery.tsv")
    file.path <- paste0(outputLevel, FileName)
    write.table(final, file.path, sep="\t",quote = FALSE, row.names = FALSE)
    
  }
} # for (level in levels)
