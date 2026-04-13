#Author: Alicia Sorgen
#Date: 2023 Nov 28
#Description: Linear modeling for taxa by time

##### Edits for script #####
rm(list=ls())
set.seed(1989)

ANALYSIS <- "MetaPhlAn2_microbiome"
params <- vector()
params <- c(params, "~/UNCC/Projects/Bariatric_Surgery/Git_Repositories/gut-microbiome-bariatric-weight-outcomes")
# params <- c(params, 0)
# params <- c(params, 1)
# params <- c(params, 6)
# params <- c(params, 12)
# params <- c(params, 18)
params <- c(params, 24)


levels <- vector()
levels <- c(levels, "phylum")
# levels <- c(levels, "class")
# levels <- c(levels, "order")
# levels <- c(levels, "family")
levels <- c(levels, "genus")
# levels <- c(levels, "species")

moduleRoot <- "4.3_Taxa_over_time_by_Surgery_MLM"

##### Libraries #####
R <- sessionInfo()
message(R$R.version$version.string); rm(R)

library(stringr); message(">>> stringr: Version ", packageVersion("stringr"), "\n")
library(ggpubr); message(">>> ggpubr: Version ", packageVersion("ggpubr"), "\n")
library(tidyr); message(">>> tidyr: Version ", packageVersion("tidyr"), "\n")
library(rstatix); message(">>> rstatix: Version ", packageVersion("rstatix"), "\n")
library(nlme); message(">>> nlme: Version ", packageVersion("nlme"), "\n")
library(data.table); message(">>> data.table: Version ", packageVersion("data.table"), "\n")
library(gridExtra); message(">>> gridExtra: Version ", packageVersion("gridExtra"), "\n")

##### Set up working environment #####
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- params
}; rm(params)

if (args[1] == "BLJ") {
  message("\n************* Running in BioLockJ *************")
} else {
  message("\n************* Running locally *************")
  gitRoot    <- args[1]
  inputEnv <- Sys.getenv("INPUT_ROOT"); gitInput <- if (nchar(inputEnv) > 0) inputEnv else file.path(gitRoot, "..", "Data")
  gitScripts <- file.path(gitRoot, "analysis", "Rscripts")

  resultsEnv <- Sys.getenv("RESULTS_ROOT"); root <- if (nchar(resultsEnv) > 0) resultsEnv else file.path(gitRoot, "Results")
  dir.create(root, showWarnings = FALSE, recursive = TRUE)

  if (length(list.files(file.path(root, "input"), recursive = TRUE)) == 0) {
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
endM <- args[length(args)]

prevModule <- paste0("TaxaMetaMerge")
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
files <- dir(inputDir)
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Models #####
# level <- levels
for (level in levels) {
  
  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  months <- c(0, 1, 6, 12, 18, 24)
  
  logCountFile <- files[which(files %like% paste0(level, "_LogNormalizedCounts"))]
  file.path <-  paste0(inputDir,logCountFile)
  logCounts <- read.table(file.path, sep="\t",header = TRUE, check.names = FALSE)
  startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1
  
  for (surgType in c("RYGB", "SG")) {
    logCounts2 <- logCounts[logCounts$Surgery == surgType,]
    
    #----- Mixed linear model -----
    results.df.list <- list()
    df.index <- 1
    for (t in 1:length(months)) {
      if (t != length(months)) {
        otherMonths <- c((t+1):length(months))
        for (t1 in otherMonths) {
          timeFrame <- months[t:t1]
          log.df <- logCounts2[ logCounts2$time %in% timeFrame, ]
          startAbundanceIndex <- which(colnames(log.df)=="ResponderStatus")+1
          
          ## Parse out metadata from counts
          myT<-log.df[,startAbundanceIndex:ncol(log.df)]
          
          PatientID <- log.df$PatientID
          Timepoint <- log.df$time
          Weight_kg <- log.df$Weight_kg
          
          #----- Mixed Linear Modeling -----
          bugName <- vector()
          pVal_Timepoint <- vector()
          pVal_Weight_kg <- vector()
          
          index <- 1
          
          for (i in 1:ncol(myT)){
            
            bug<-myT[,i]
            
            if (mean(bug>0, na.rm = TRUE)>0.1){
              
              df <- data.frame(bug, PatientID, Timepoint, Weight_kg)
              df <- na.omit(df)
              df$Timepoint <- as.factor(df$Timepoint)
              
              # message(paste0(i, ". ", colnames(myT)[i]))
              tryCatch({
                mlm <- lme( bug ~ Timepoint + Weight_kg, method = "REML", random = ~1 | PatientID, data = df )
                fit<-anova(mlm)
                
                pVal_Timepoint[index] <- fit$"p-value"[2]
                pVal_Weight_kg[index] <- fit$"p-value"[3]
                bugName[index]<-colnames(myT)[i]
                
                index<-index+1
                
              }, error = function(e){})
            } # if (mean(bug>0, na.rm = TRUE)>0.1)
            
          } # for (i in 1:ncol(myT))
          # message("Loop done.")
          
          dFrame<-data.frame(pVal_Timepoint, pVal_Weight_kg)
          dFrame$Adj_pVal_Timepoint <- p.adjust(dFrame$pVal_Timepoint,method = "BH")
          dFrame$Adj_pVal_Weight_kg <- p.adjust(dFrame$pVal_Weight_kg,method = "BH")
          
          names(dFrame) <- paste0(names(dFrame), "_", timeFrame[1], "to", timeFrame[length(timeFrame)])
          dFrame <- data.frame(bugName, dFrame)
          #----- Mixed Linear Model END -----
          
          results.df.list[[df.index]] <- dFrame
          df.index <- df.index + 1
          
        } # for (t1 in otherMonths)
      } # if (t != length(months))
    } # for (t in 1:length(months))
    
    
    merged.results <- results.df.list[[1]]
    if (length(results.df.list) > 1) {
      for (l in 2:length(results.df.list)) {
        
        merged.results <- merge(merged.results, results.df.list[[l]], by = "bugName")
        
      } # for (l in 2:length(results.df.list))
    }
    #-----  end  --------- 
    
    #---- Print table ----
    write.table(merged.results, paste0(outputLevel, level,  "~Timepoint+Weight_kg_MixedLinearModelResults_", surgType, ".tsv"),sep="\t",row.names = FALSE,quote = FALSE)
    #-----  end  --------- 
    
    
    #----- Print results plots -----
    file.path <- paste0(outputLevel, level, "~Timepoint+Weight_kg_MixedLinearModel_PLOTS_", surgType, ".pdf")
    pdf(file.path)
    
    par(mfrow=c(2,2))
    
    for( i in 1:nrow(merged.results)) {
      BugName <- merged.results$bugName[i]
      
      for (t in 1:length(months)) {
        if (t != length(months)) {
          otherMonths <- c((t+1):length(months))
          for (t1 in otherMonths) {
            timeFrame <- months[t:t1]
            log.df <- logCounts2[ logCounts2$time %in% timeFrame, ]
            bugCol <- log.df[ , which(colnames(log.df) == BugName) ]
            Month <- log.df$time
            Weight <- log.df$Weight_kg
            PatientID <- log.df$PatientID
            
            df <- data.frame(PatientID, Month, Weight, bugCol)
            df <- na.omit(df)
            
            colEnding <- paste0(timeFrame[1], "to", timeFrame[length(timeFrame)])
            Timepoint_p <- merged.results[ i , which(colnames(merged.results) == paste0("Adj_pVal_Timepoint_", colEnding))]
            Weight_p <- merged.results[ i , which(colnames(merged.results) == paste0("Adj_pVal_Weight_kg_", colEnding))]
            
            patient_n <- length(unique(df$PatientID))
            n <- nrow(df)
            
            title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Weight_p, length = 3), "; n = ", n)
            plot(  df$Weight, df$bugCol,
                   ylab =BugName ,
                   xlab="Weight (kg)",
                   main=title.lab,
                   cex.main = 1 
            )
            
            title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Timepoint_p, length = 3), "; n = ", n)
            plot(  df$Month, df$bugCol, xaxt = 'n',
                   ylab =BugName ,
                   xlab="Timepoint (months)",
                   main=title.lab,
                   cex.main = 1 
            )
            axis(1, at = seq(0, max(df$Month), by = 1))
          } # for (t1 in otherMonths)
        } # if (t != length(months))
      } # for (t in 1:length(months))
      
    } # for( i in 1:nrow(merged.results))
    dev.off()
    #-----  end  --------- 
    
    #---- Print histograms ----
    file.path <- paste0(outputLevel, level,  "~Timepoint+Weight_kg_MixedLinearModel_pValueHistograms_", surgType, ".pdf")
    pdf(file.path)
    
    for (t in 1:length(months)) {
      if (t != length(months)) {
        otherMonths <- c((t+1):length(months))
        for (t1 in otherMonths) {
          timeFrame <- months[t:t1]
          colEnding <- paste0(timeFrame[1], "to", timeFrame[length(timeFrame)])
          
          pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Timepoint_", colEnding))]
          n=length(pValues)
          title.lab <- paste0(colEnding, "\nMLM Timepoint p values (n = ", n, ")")
          hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
          
          pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Weight_kg_", colEnding))]
          n=length(pValues)
          title.lab <- paste0(colEnding, "\nMLM Weight_kg p values (n = ", n, ")")
          hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
          
        } # for (t1 in otherMonths)
      } # if (t != length(months))
    } # for (t in 1:length(months))
    
    dev.off()
    #-----  end  --------- 
    
  } # for (surgType in c("RYGB", "SG"))
  
} # for (level in levels)

