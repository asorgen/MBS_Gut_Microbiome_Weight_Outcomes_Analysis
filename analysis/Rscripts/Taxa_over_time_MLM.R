#Author: Alicia Sorgen
#Date: 2023 Nov 28
#Description: Linear modeling for taxa by time

##### Edits for script #####
rm(list=ls())
set.seed(1989)

ANALYSIS <- "MetaPhlAn2_microbiome"
params <- vector()
params <- c(params, "~/git/gut-microbiome-bariatric-weight-outcomes")
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

moduleRoot <- "Taxa_over_time_MLM"

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
  gitRoot <- args[1]
  gitInput <- file.path(gitRoot, "analysis", "input")
  gitScripts <- file.path(gitRoot, "analysis", "Rscripts"); rm(gitRoot)
  
  root <- paste0("~/BioLockJ_pipelines/")
  pipeline <- paste0(ANALYSIS,"_analysis_")
  
  if (any(dir(root) %like% pipeline) == TRUE) {
    root <- paste0(root,"/",str_subset(dir(root), pipeline), "/")
  } else {
    today <- as.character(format(Sys.Date(), "%Y%b%d"))
    root <- paste0(root,ANALYSIS,"_analysis_", today, "/")
    dir.create(root, showWarnings = FALSE)
  }; rm(pipeline, ANALYSIS)
  
  if (any(dir(root) == "input") == FALSE) {
    file.copy(gitInput,
              root,
              recursive = TRUE)
  }; rm(gitInput)
  
  module <- moduleRoot
  
  if (exists("end_month") == TRUE) {
    mod1 <- paste0(sapply(strsplit(module, "_Weight"), "[", 1), "_"); mod1
    mod2 <- paste0("_Weight", sapply(strsplit(module, "_Weight"), "[", 2)); mod2
    module <- paste0(mod1, end_month, "M", mod2, "_BLto24M"); module
  }
  
  if (any(dir(root) %like% module) == TRUE) {
    moduleDir <- paste0(root,str_subset(dir(root), module), "/")
  } else {
    moduleDir <- paste0(root, module, "/")
    dir.create(moduleDir, showWarnings = FALSE)
  }; rm(root, module)
  
  scriptDir <- paste0(moduleDir, "script/")
  if (any(dir(moduleDir) == "script") == FALSE) {
    dir.create(scriptDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/", moduleRoot, ".R"); script
    file.copy(script,
              scriptDir,
              recursive = TRUE)
  }; rm(scriptDir, moduleRoot)
  
  outputDir <- paste0(moduleDir, "output/")
  if (any(dir(moduleDir) == "output") == FALSE) {
    dir.create(outputDir, showWarnings = FALSE)
  }
  
  files <- list.files(outputDir, recursive = TRUE, full.names = TRUE); rm(outputDir)
  # file.remove(files)
  rm(files)
  
  resourcesDir <- paste0(moduleDir, "resources/")
  if (any(dir(moduleDir) == "resources") == FALSE) {
    dir.create(resourcesDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/functions.R")
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
  }; rm(gitScripts, resourcesDir)
  setwd(paste0(moduleDir, "script/")); rm(moduleDir)
  
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
endM <- args[length(args)]

prevModule <- paste0("TaxaMetaMerge")
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
files <- dir(inputDir)
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Models #####
comparisons <- c("1_v_0", "6_v_0", "12_v_0", "18_v_0", "24_v_0",
                 "6_v_1", "12_v_1", "18_v_1", "24_v_1",
                 "12_v_6", "18_v_6", "24_v_6",
                 "18_v_12", "24_v_12",
                 "24_v_18")
comparisons <- paste0("Timepoint", comparisons)

for (level in levels) {
  
  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  months <- c(0, 1, 6, 12, 18, 24)
  
  logCountFile <- files[which(files %like% paste0(level, "_LogNormalizedCounts"))]
  file.path <-  paste0(inputDir,logCountFile)
  logCounts <- read.table(file.path, sep="\t",header = TRUE, check.names = FALSE)
  startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1

  #----- Mixed linear model -----
  results.df.list <- list()
  df.index <- 1
  for (t in 1:length(months)) {
    if (t != length(months)) {
      otherMonths <- c((t+1):length(months))
      for (t1 in otherMonths) {
        timeFrame <- months[t:t1]
        log.df <- logCounts[ logCounts$time %in% timeFrame, ]
        startAbundanceIndex <- which(colnames(log.df)=="ResponderStatus")+1
        
        ## Parse out metadata from counts
        myT<-log.df[,startAbundanceIndex:ncol(log.df)]
        
        if (level == "genus") {
          myT <- myT[,-(which(colnames(myT) %like% "ales_noname"))]
          myT <- myT[,-(which(colnames(myT) %like% "Candidatus_Saccharibacteria"))]
          termstoRemove <- c("Bacteroidetes_noname")
          myT <- myT[,-(which(colnames(myT) %in% termstoRemove))]
          # colnames(myT)
        }
        PatientID <- log.df$PatientID
        Timepoint <- log.df$time
        Surgery <- log.df$Surgery
        Weight_kg <- log.df$Weight_kg
        
        #----- Mixed Linear Modeling -----
        bugName <- vector()
        pVal_Timepoint <- vector()
        pVal_Surgery <- vector()
        pVal_Weight_kg <- vector()
        
        p_Timepoint1 <- vector()
        p_Timepoint6 <- vector()
        p_Timepoint12 <- vector()
        p_Timepoint18 <- vector()
        p_Timepoint24 <- vector()
        
        s_Timepoint1 <- vector()
        s_Timepoint6 <- vector()
        s_Timepoint12 <- vector()
        s_Timepoint18 <- vector()
        s_Timepoint24 <- vector()
        
        index <- 1
        
        for (i in 1:ncol(myT)){
          
          bug<-myT[,i]
          
          if (mean(bug>0, na.rm = TRUE)>0.1){
            
            df <- data.frame(bug, PatientID, Timepoint, Surgery, Weight_kg)
            df <- na.omit(df)
            df$Timepoint <- as.factor(df$Timepoint)
            
            # message(paste0(i, ". ", colnames(myT)[i]))
            tryCatch({
              mlm <- lme( bug ~ Timepoint + Surgery +  Weight_kg, method = "REML", random = ~1 | PatientID, data = df )
              fit<-anova(mlm)
              sm <- summary(mlm)
              
              pVal_Timepoint[index] <- fit$"p-value"[2]
              pVal_Surgery[index] <- fit$"p-value"[3]
              pVal_Weight_kg[index] <- fit$"p-value"[4]
              
              if ("Timepoint1" %in% names(sm$tTable[,1])) {p_Timepoint1[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint1"),5]}
              if ("Timepoint6" %in% names(sm$tTable[,1])) {p_Timepoint6[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint6"),5]}
              if ("Timepoint12" %in% names(sm$tTable[,1])) {p_Timepoint12[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint12"),5]}
              if ("Timepoint18" %in% names(sm$tTable[,1])) {p_Timepoint18[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint18"),5]}
              if ("Timepoint24" %in% names(sm$tTable[,1])) {p_Timepoint24[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint24"),5]}
              
              if ("Timepoint1" %in% names(sm$tTable[,1])) {s_Timepoint1[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint1"),1]}
              if ("Timepoint6" %in% names(sm$tTable[,1])) {s_Timepoint6[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint6"),1]}
              if ("Timepoint12" %in% names(sm$tTable[,1])) {s_Timepoint12[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint12"),1]}
              if ("Timepoint18" %in% names(sm$tTable[,1])) {s_Timepoint18[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint18"),1]}
              if ("Timepoint24" %in% names(sm$tTable[,1])) {s_Timepoint24[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint24"),1]}
              
              bugName[index]<-colnames(myT)[i]
              
              index<-index+1
              
            }, error = function(e){})
          } # if (mean(bug>0, na.rm = TRUE)>0.1)
          
        } # for (i in 1:ncol(myT))
        # message("Loop done.")
        
        dFrame<-data.frame(pVal_Timepoint, pVal_Surgery, pVal_Weight_kg)
        dFrame$Adj_pVal_Timepoint <- p.adjust(dFrame$pVal_Timepoint,method = "BH")
        dFrame$Adj_pVal_Surgery <- p.adjust(dFrame$pVal_Surgery,method = "BH")
        dFrame$Adj_pVal_Weight_kg <- p.adjust(dFrame$pVal_Weight_kg,method = "BH")
        
        if (length(p_Timepoint1) > 0) {
          dFrame <- data.frame(dFrame, p_Timepoint1, s_Timepoint1)
          dFrame$Adj_p_Timepoint1 <- p.adjust(dFrame$p_Timepoint1,method = "BH")
        }
        
        if (length(p_Timepoint6) > 0) {
          dFrame <- data.frame(dFrame, p_Timepoint6, s_Timepoint6)
          dFrame$Adj_p_Timepoint6 <- p.adjust(dFrame$p_Timepoint6,method = "BH")
        }
        
        if (length(p_Timepoint12) > 0) {
          dFrame <- data.frame(dFrame, p_Timepoint12, s_Timepoint12)
          dFrame$Adj_p_Timepoint12 <- p.adjust(dFrame$p_Timepoint12,method = "BH")
        }
        
        if (length(p_Timepoint18) > 0) {
          dFrame <- data.frame(dFrame, p_Timepoint18, s_Timepoint18)
          dFrame$Adj_p_Timepoint18 <- p.adjust(dFrame$p_Timepoint18,method = "BH")
        }
        
        if (length(p_Timepoint24) > 0) {
          dFrame <- data.frame(dFrame, p_Timepoint24, s_Timepoint24)
          dFrame$Adj_p_Timepoint24 <- p.adjust(dFrame$p_Timepoint24,method = "BH")
        }
        
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
  write.table(merged.results, paste0(outputLevel, level,  "~Timepoint+Surgery+Weight_kg_MixedLinearModelResults.tsv"),sep="\t",row.names = FALSE,quote = FALSE)
  #-----  end  --------- 
  
  
  #----- Print results plots -----
  file.path <- paste0(outputLevel, level, "~Timepoint+Surgery+Weight_kg_MixedLinearModel_PLOTS.pdf")
  pdf(file.path)
  
  par(mfrow=c(2,2))
  
  for( i in 1:nrow(merged.results)) {
    BugName <- merged.results$bugName[i]
    
    for (t in 1:length(months)) {
      if (t != length(months)) {
        otherMonths <- c((t+1):length(months))
        for (t1 in otherMonths) {
          timeFrame <- months[t:t1]
          log.df <- logCounts[ logCounts$time %in% timeFrame, ]
          bugCol <- log.df[ , which(colnames(log.df) == BugName) ]
          Month <- log.df$time
          Surgery <- log.df$Surgery
          Weight <- log.df$Weight_kg
          PatientID <- log.df$PatientID
          
          df <- data.frame(PatientID, Month, Surgery, Weight, bugCol)
          df <- na.omit(df)
          
          colEnding <- paste0(timeFrame[1], "to", timeFrame[length(timeFrame)])
          Timepoint_p <- merged.results[ i , which(colnames(merged.results) == paste0("Adj_pVal_Timepoint_", colEnding))]
          Surgery_p <- merged.results[ i , which(colnames(merged.results) == paste0("Adj_pVal_Surgery_", colEnding))]
          Weight_p <- merged.results[ i , which(colnames(merged.results) == paste0("Adj_pVal_Weight_kg_", colEnding))]
          
          patient_n <- length(unique(df$PatientID))
          n <- nrow(df)
          
          title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Weight_p, length = 3), "; n = ", n)
          plot(  df$Weight, df$bugCol,
                 ylab =BugName ,
                 xlab="Weight (kg)",
                 main=title.lab,
                 cex.main = 1 ,
                 col = ifelse( df$Surgery == "SG", "red", "blue" )
          )
          legend("topright", legend=c("SG", "RYGB"),
                 col=c("red", "blue"), pch = 21, cex=0.8)
          
          title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Timepoint_p, length = 3), "; n = ", n)
          plot(  df$Month, df$bugCol, xaxt = 'n',
                 ylab =BugName ,
                 xlab="Timepoint (months)",
                 main=title.lab,
                 cex.main = 1 ,
                 col = ifelse( df$Surgery == "SG", "red", "blue" )
          )
          axis(1, at = seq(0, max(df$Month), by = 1))
          legend("top", legend=c("SG", "RYGB"),
                 col=c("red", "blue"), pch = 21, cex=0.8)
          
          title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Surgery_p, length = 3), "; n = ", n)
          boxplot( df$bugCol ~  df$Surgery ,
                   xlab="Surgery Type",
                   ylab =BugName,
                   main = title.lab,
                   cex.main = 1
          )
          
        } # for (t1 in otherMonths)
      } # if (t != length(months))
    } # for (t in 1:length(months))
    
  } # for( i in 1:nrow(merged.results))
  dev.off()
  #-----  end  --------- 
  
  #---- Print histograms ----
  file.path <- paste0(outputLevel, level,  "~Timepoint+Surgery+Weight_kg_MixedLinearModel_pValueHistograms.pdf")
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
        
        pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Surgery_", colEnding))]
        n=length(pValues)
        title.lab <- paste0(colEnding, "\nMLM Surgery p values (n = ", n, ")")
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

  
  
  ##----- Heatmap set up --------------
  results.df.list <- list()
  df.index <- 1
  for (t in 1:length(months)) {
    if (t != length(months)) {
      
      otherMonths <- months[-t]
      factorOrder <- c(months[t], otherMonths)
      
      # otherMonths <- c((t+1):length(months))
      # for (t1 in otherMonths) {
      
      timeFrame <- months[1:length(months)]
      log.df <- logCounts[ logCounts$time %in% timeFrame, ]
      startAbundanceIndex <- which(colnames(log.df)=="ResponderStatus")+1
      
      ## Parse out metadata from counts
      myT<-log.df[,startAbundanceIndex:ncol(log.df)]
      
      if (level == "genus") {
        myT <- myT[,-(which(colnames(myT) %like% "ales_noname"))]
        myT <- myT[,-(which(colnames(myT) %like% "Candidatus_Saccharibacteria"))]
        termstoRemove <- c("Bacteroidetes_noname")
        myT <- myT[,-(which(colnames(myT) %in% termstoRemove))]
        # colnames(myT)
      }
      PatientID <- log.df$PatientID
      Timepoint <- log.df$time
      Surgery <- log.df$Surgery
      Weight_kg <- log.df$Weight_kg
      
      #----- Mixed Linear Modeling -----
      bugName <- vector()
      pVal_Timepoint <- vector()
      pVal_Surgery <- vector()
      pVal_Weight_kg <- vector()
      
      p_Timepoint0 <- vector()
      p_Timepoint1 <- vector()
      p_Timepoint6 <- vector()
      p_Timepoint12 <- vector()
      p_Timepoint18 <- vector()
      p_Timepoint24 <- vector()
      
      s_Timepoint0 <- vector()
      s_Timepoint1 <- vector()
      s_Timepoint6 <- vector()
      s_Timepoint12 <- vector()
      s_Timepoint18 <- vector()
      s_Timepoint24 <- vector()
      
      index <- 1
      
      for (i in 1:ncol(myT)){
        
        bug<-myT[,i]
        
        if (mean(bug>0, na.rm = TRUE)>0.1){
          
          df <- data.frame(bug, PatientID, Timepoint, Surgery, Weight_kg)
          df <- na.omit(df)
          
          smstat <- df %>%
            group_by(Timepoint) %>%
            get_summary_stats(bug, type = "mean_se")
          smstat$Timepoint <- as.numeric(smstat$Timepoint)
          
          df$Timepoint <- as.factor(df$Timepoint)
          df$Timepoint <- factor(df$Timepoint, levels = factorOrder)
          
          # message(paste0(i, ". ", colnames(myT)[i]))
          tryCatch({
            mlm <- lme( bug ~ Timepoint + Surgery +  Weight_kg, method = "REML", random = ~1 | PatientID, data = df )
            fit<-anova(mlm)
            sm <- summary(mlm)
            
            pVal_Timepoint[index] <- fit$"p-value"[2]
            pVal_Surgery[index] <- fit$"p-value"[3]
            pVal_Weight_kg[index] <- fit$"p-value"[4]
            
            if ("Timepoint0" %in% names(sm$tTable[,1])) {p_Timepoint0[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint0"),5]}
            if ("Timepoint1" %in% names(sm$tTable[,1])) {p_Timepoint1[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint1"),5]}
            if ("Timepoint6" %in% names(sm$tTable[,1])) {p_Timepoint6[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint6"),5]}
            if ("Timepoint12" %in% names(sm$tTable[,1])) {p_Timepoint12[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint12"),5]}
            if ("Timepoint18" %in% names(sm$tTable[,1])) {p_Timepoint18[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint18"),5]}
            if ("Timepoint24" %in% names(sm$tTable[,1])) {p_Timepoint24[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint24"),5]}
            
            if ("Timepoint0" %in% names(sm$tTable[,1])) {s_Timepoint0[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint0"),1]}
            if ("Timepoint1" %in% names(sm$tTable[,1])) {s_Timepoint1[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint1"),1]}
            if ("Timepoint6" %in% names(sm$tTable[,1])) {s_Timepoint6[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint6"),1]}
            if ("Timepoint12" %in% names(sm$tTable[,1])) {s_Timepoint12[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint12"),1]}
            if ("Timepoint18" %in% names(sm$tTable[,1])) {s_Timepoint18[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint18"),1]}
            if ("Timepoint24" %in% names(sm$tTable[,1])) {s_Timepoint24[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint24"),1]}
            
            bugName[index]<-colnames(myT)[i]
            
            index<-index+1
            
          }, error = function(e){})
        } # if (mean(bug>0, na.rm = TRUE)>0.1)
        
      } # for (i in 1:ncol(myT))
      # message("Loop done.")
      
      dFrame <- data.frame(bugName)
      
      # dFrame<-data.frame(pVal_Timepoint, pVal_Surgery, pVal_Weight_kg)
      # dFrame$Adj_pVal_Timepoint <- p.adjust(dFrame$pVal_Timepoint,method = "BH")
      # dFrame$Adj_pVal_Surgery <- p.adjust(dFrame$pVal_Surgery,method = "BH")
      # dFrame$Adj_pVal_Weight_kg <- p.adjust(dFrame$pVal_Weight_kg,method = "BH")
      
      if (length(p_Timepoint0) > 0) {
        dFrame <- data.frame(dFrame, p_Timepoint0, s_Timepoint0)
        dFrame$Adj_p_Timepoint0 <- p.adjust(dFrame$p_Timepoint0,method = "BH")
      }
      
      if (length(p_Timepoint1) > 0) {
        dFrame <- data.frame(dFrame, p_Timepoint1, s_Timepoint1)
        dFrame$Adj_p_Timepoint1 <- p.adjust(dFrame$p_Timepoint1,method = "BH")
      }
      
      if (length(p_Timepoint6) > 0) {
        dFrame <- data.frame(dFrame, p_Timepoint6, s_Timepoint6)
        dFrame$Adj_p_Timepoint6 <- p.adjust(dFrame$p_Timepoint6,method = "BH")
      }
      
      if (length(p_Timepoint12) > 0) {
        dFrame <- data.frame(dFrame, p_Timepoint12, s_Timepoint12)
        dFrame$Adj_p_Timepoint12 <- p.adjust(dFrame$p_Timepoint12,method = "BH")
      }
      
      if (length(p_Timepoint18) > 0) {
        dFrame <- data.frame(dFrame, p_Timepoint18, s_Timepoint18)
        dFrame$Adj_p_Timepoint18 <- p.adjust(dFrame$p_Timepoint18,method = "BH")
      }
      
      if (length(p_Timepoint24) > 0) {
        dFrame <- data.frame(dFrame, p_Timepoint24, s_Timepoint24)
        dFrame$Adj_p_Timepoint24 <- p.adjust(dFrame$p_Timepoint24,method = "BH")
      }
      
      names(dFrame)[2:ncol(dFrame)] <- paste0(names(dFrame)[2:ncol(dFrame)], "_v_", months[t])
      # dFrame <- data.frame(bugName, dFrame)
      #----- Mixed Linear Model END -----
      
      results.df.list[[df.index]] <- dFrame
      df.index <- df.index + 1
      
      # } # for (t1 in otherMonths)
    } # if (t != length(months))
  } # for (t in 1:length(months))
  
  merged.results <- results.df.list[[1]]
  if (length(results.df.list) > 1) {
    for (l in 2:length(results.df.list)) {
      
      merged.results <- merge(merged.results, results.df.list[[l]], by = "bugName")
      
    } # for (l in 2:length(results.df.list))
  }
  
  #---- Print table ----
  write.table(merged.results, paste0(outputLevel, level,  "~Timepoint+Surgery+Weight_kg_MixedLinearModelResults_Heatmap.tsv"),sep="\t",row.names = FALSE,quote = FALSE)
  #-----  end  --------- 
  
  merge2 <- merged.results %>% gather("Key", "Value", 2:ncol(merged.results))
  merge2$New <- ifelse(merge2$Key %like% "s_", "t-table Value",
                       ifelse(merge2$Key %like% "Adj_", "Adjusted p-value",
                              "p-value"))
  merge2$Key <- gsub("Adj_", "", merge2$Key)
  merge2$Key <- gsub("s_", "", merge2$Key)
  merge2$Key <- gsub("p_", "", merge2$Key)
  
  spread.df <- spread(merge2, New, Value)
  spread.df <- spread.df[spread.df$Key %in% comparisons,]
  spread.df$Key <- gsub("Timepoint", "", spread.df$Key)
  #---- Print table ----
  write.table(spread.df, paste0(outputLevel, level,  "~Timepoint+Surgery+Weight_kg_MixedLinearModelResults_Heatmap_Long.tsv"),sep="\t",row.names = FALSE,quote = FALSE)
  #-----  end  --------- 
  
  spread.df$start <- as.numeric(sapply(strsplit(spread.df$Key, "_v_"), "[", 2))
  spread.df$end <- as.numeric(sapply(strsplit(spread.df$Key, "_v_"), "[", 1))
  spread.df$signif <- sigStars(spread.df$`Adjusted p-value`)
  
  line.index <- 1
  line.List <- list()
  for (taxa in unique(spread.df$bugName)) {
    
    spread.df2 <- spread.df[spread.df$bugName == taxa,]
    spread.df3 <- spread.df2[spread.df2$`Adjusted p-value` < 0.05,]
    spread.df3 <- spread.df3[order(spread.df3$start, spread.df3$end), ]
    
    if (nrow(spread.df3) > 0) {
      timepoint <- logCounts$time
      bug<-logCounts[,which(colnames(logCounts) == taxa)]
      df <- data.frame(bug, timepoint)
      df <- na.omit(df)
      
      smstat <- df %>%
        group_by(timepoint) %>%
        get_summary_stats(bug, type = "mean_se")
      smstat$max <- smstat$mean + smstat$se
      max <- max(smstat$max)
      nudge <- max(smstat$max)*.05

      plot <- ggplot(smstat, aes(x = timepoint, y = mean, group = 1)) +
        geom_line(color = "blue", size = 2) +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, color = "black") +
        theme_bw() +
        labs(x = "Time point (month)", y = "Log10 abundance", title = taxa); plot
      
      for (p in 1:nrow(spread.df3)) {
        max <- max + nudge
        plot <- plot + geom_signif(y_position = c(max), xmin = spread.df3[p,"start"], 
                                   xmax = spread.df3[p,"end"], annotation = spread.df3[p,"signif"],
                                   tip_length = 0)
        
      }
      line.List[[line.index]] <- plot
      line.index <- line.index + 1
    } # if (nrow(spread.df3) > 0)
    
    
  } # for (taxa in unique(spread.df$bugName))
  
  pdf(file.path(outputLevel, paste0(level, "_linegraph.pdf")),height = 5, width = 10)
  for (i in 1:length(line.List)) {
    print(line.List[[i]])
  }
  dev.off()
  
  ##----- Heatmap set up classified taxa only--------------
  results.df.list <- list()
  df.index <- 1
  for (t in 1:length(months)) {
    if (t != length(months)) {
      
      otherMonths <- months[-t]
      factorOrder <- c(months[t], otherMonths)
      
      # otherMonths <- c((t+1):length(months))
      # for (t1 in otherMonths) {
      
      timeFrame <- months[1:length(months)]
      log.df <- logCounts[ logCounts$time %in% timeFrame, ]
      startAbundanceIndex <- which(colnames(log.df)=="ResponderStatus")+1
      
      ## Parse out metadata from counts
      myT<-log.df[,startAbundanceIndex:ncol(log.df)] #163
      
      if (level == "genus") {
        myT <- myT[,-(which(colnames(myT) %like% "_noname"))] #148
        myT <- myT[,-(which(colnames(myT) %like% "_unclassified"))] #139
        # colnames(myT)
      }
      PatientID <- log.df$PatientID
      Timepoint <- log.df$time
      Surgery <- log.df$Surgery
      Weight_kg <- log.df$Weight_kg
      
      #----- Mixed Linear Modeling -----
      bugName <- vector()
      pVal_Timepoint <- vector()
      pVal_Surgery <- vector()
      pVal_Weight_kg <- vector()
      
      p_Timepoint0 <- vector()
      p_Timepoint1 <- vector()
      p_Timepoint6 <- vector()
      p_Timepoint12 <- vector()
      p_Timepoint18 <- vector()
      p_Timepoint24 <- vector()
      
      s_Timepoint0 <- vector()
      s_Timepoint1 <- vector()
      s_Timepoint6 <- vector()
      s_Timepoint12 <- vector()
      s_Timepoint18 <- vector()
      s_Timepoint24 <- vector()
      
      index <- 1
      
      for (i in 1:ncol(myT)){
        
        bug<-myT[,i]
        
        if (mean(bug>0, na.rm = TRUE)>0.1){
          
          df <- data.frame(bug, PatientID, Timepoint, Surgery, Weight_kg)
          df <- na.omit(df)
          
          df$Timepoint <- as.factor(df$Timepoint)
          df$Timepoint <- factor(df$Timepoint, levels = factorOrder)
          
          
          tryCatch({
            mlm <- lme( bug ~ Timepoint + Surgery +  Weight_kg, method = "REML", random = ~1 | PatientID, data = df )
            fit<-anova(mlm)
            sm <- summary(mlm)
            
            pVal_Timepoint[index] <- fit$"p-value"[2]
            pVal_Surgery[index] <- fit$"p-value"[3]
            pVal_Weight_kg[index] <- fit$"p-value"[4]
            
            if ("Timepoint0" %in% names(sm$tTable[,1])) {p_Timepoint0[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint0"),5]}
            if ("Timepoint1" %in% names(sm$tTable[,1])) {p_Timepoint1[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint1"),5]}
            if ("Timepoint6" %in% names(sm$tTable[,1])) {p_Timepoint6[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint6"),5]}
            if ("Timepoint12" %in% names(sm$tTable[,1])) {p_Timepoint12[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint12"),5]}
            if ("Timepoint18" %in% names(sm$tTable[,1])) {p_Timepoint18[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint18"),5]}
            if ("Timepoint24" %in% names(sm$tTable[,1])) {p_Timepoint24[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint24"),5]}
            
            if ("Timepoint0" %in% names(sm$tTable[,1])) {s_Timepoint0[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint0"),1]}
            if ("Timepoint1" %in% names(sm$tTable[,1])) {s_Timepoint1[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint1"),1]}
            if ("Timepoint6" %in% names(sm$tTable[,1])) {s_Timepoint6[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint6"),1]}
            if ("Timepoint12" %in% names(sm$tTable[,1])) {s_Timepoint12[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint12"),1]}
            if ("Timepoint18" %in% names(sm$tTable[,1])) {s_Timepoint18[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint18"),1]}
            if ("Timepoint24" %in% names(sm$tTable[,1])) {s_Timepoint24[index] <- sm$tTable[which(names(sm$tTable[,1]) == "Timepoint24"),1]}
            
            bugName[index]<-colnames(myT)[i]
            
            index<-index+1
            
          }, error = function(e){})
        } # if (mean(bug>0, na.rm = TRUE)>0.1)
        
      } # for (i in 1:ncol(myT))

      dFrame <- data.frame(bugName)
      
      if (length(p_Timepoint0) > 0) {
        dFrame <- data.frame(dFrame, p_Timepoint0, s_Timepoint0)
        dFrame$Adj_p_Timepoint0 <- p.adjust(dFrame$p_Timepoint0,method = "BH")
      }
      
      if (length(p_Timepoint1) > 0) {
        dFrame <- data.frame(dFrame, p_Timepoint1, s_Timepoint1)
        dFrame$Adj_p_Timepoint1 <- p.adjust(dFrame$p_Timepoint1,method = "BH")
      }
      
      if (length(p_Timepoint6) > 0) {
        dFrame <- data.frame(dFrame, p_Timepoint6, s_Timepoint6)
        dFrame$Adj_p_Timepoint6 <- p.adjust(dFrame$p_Timepoint6,method = "BH")
      }
      
      if (length(p_Timepoint12) > 0) {
        dFrame <- data.frame(dFrame, p_Timepoint12, s_Timepoint12)
        dFrame$Adj_p_Timepoint12 <- p.adjust(dFrame$p_Timepoint12,method = "BH")
      }
      
      if (length(p_Timepoint18) > 0) {
        dFrame <- data.frame(dFrame, p_Timepoint18, s_Timepoint18)
        dFrame$Adj_p_Timepoint18 <- p.adjust(dFrame$p_Timepoint18,method = "BH")
      }
      
      if (length(p_Timepoint24) > 0) {
        dFrame <- data.frame(dFrame, p_Timepoint24, s_Timepoint24)
        dFrame$Adj_p_Timepoint24 <- p.adjust(dFrame$p_Timepoint24,method = "BH")
      }
      
      names(dFrame)[2:ncol(dFrame)] <- paste0(names(dFrame)[2:ncol(dFrame)], "_v_", months[t])
      
      
      # dFrame <- data.frame(bugName, dFrame)
      #----- Mixed Linear Model END -----
      
      results.df.list[[df.index]] <- dFrame
      df.index <- df.index + 1
      
      # } # for (t1 in otherMonths)
    } # if (t != length(months))
  } # for (t in 1:length(months))
  
  merged.results <- results.df.list[[1]]
  if (length(results.df.list) > 1) {
    for (l in 2:length(results.df.list)) {
      
      merged.results <- merge(merged.results, results.df.list[[l]], by = "bugName")
      
    } # for (l in 2:length(results.df.list))
  }
  
  
  if (level == "genus") {
    #---- Print table ----
    write.table(merged.results, paste0(outputLevel, level,  "~Timepoint+Surgery+Weight_kg_MixedLinearModelResults_Heatmap_Classified.tsv"),sep="\t",row.names = FALSE,quote = FALSE)
    #-----  end  --------- 
    
    merge2 <- merged.results %>% gather("Key", "Value", 2:ncol(merged.results))
    merge2$New <- ifelse(merge2$Key %like% "s_", "t-table Value",
                         ifelse(merge2$Key %like% "Adj_", "Adjusted p-value",
                                "p-value"))
    merge2$Key <- gsub("Adj_", "", merge2$Key)
    merge2$Key <- gsub("s_", "", merge2$Key)
    merge2$Key <- gsub("p_", "", merge2$Key)
    
    spread.df <- spread(merge2, New, Value)
    spread.df <- spread.df[spread.df$Key %in% comparisons,]
    #---- Print table ----
    write.table(spread.df, paste0(outputLevel, level,  "~Timepoint+Surgery+Weight_kg_MixedLinearModelResults_Heatmap__Classified_Long.tsv"),sep="\t",row.names = FALSE,quote = FALSE)
    #-----  end  --------- 
  }
  
  
} # for (level in levels)





