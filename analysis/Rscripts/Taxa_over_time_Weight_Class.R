#Author: Alicia Sorgen
#Date: 2022 June 09
#Description: Linear modeling for taxa by time

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/git/gut-microbiome-bariatric-weight-outcomes")
params <- c(params, "MetaPhlAn2")
# params <- c(params, 0)
params <- c(params, 1)
# params <- c(params, 6)
# params <- c(params, 12)
# params <- c(params, 18)
# params <- c(params, 24)
end_month <- params[length(params)]

levels <- vector()
levels <- c(levels, "phylum")
# levels <- c(levels, "class")
# levels <- c(levels, "order")
# levels <- c(levels, "family")
levels <- c(levels, "genus")
# levels <- c(levels, "species")

moduleRoot <- "Taxa_over_time_Weight_Class"

##### Libraries #####
R <- sessionInfo()
message(R$R.version$version.string)

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
    module <- paste0(args[2], "_", module); module
  } 
  
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
endM <- args[length(args)]

prevModule <- paste0(classifier, "_TaxaMetaMerge")
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
logCountFile <- paste0("_LogNormalizedCounts_", classifier, ".tsv")
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")
# divisions <- c("Quintile", "Quartile", "Tertile", "Half")
divisions <- c("Tertile")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Models (Weight loss groups) #####
# level <- levels
for (level in levels) {
  
  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  months <- c(0, 1, 6, 12, 18, 24)
  
  # division <- divisions
  for (division in divisions) {
    
    message(paste0("\n***** Starting ", level, " multivariate (", division, ") *****\n"))
    
    outputWLdivisions <- paste0(outputLevel, division, "/")
    dir.create(outputWLdivisions, showWarnings = FALSE)
    
    file.path <-  paste0(inputDir,level, logCountFile)
    logCounts <- read.table(file.path, sep="\t",header = TRUE, check.names = FALSE)
    startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1
    logCounts$time <- as.factor(logCounts$time)
    logCounts$time <- factor(logCounts$time, levels = months)
    
    # logCounts <- logCounts[ !is.na( logCounts[,startAbundanceIndex+1]), ] # remove rows with bug NAs
    # logCounts <- logCounts[ !is.na( logCounts$Weight_kg), ] # remove rows with Weight NAs
    # logCounts <- logCounts[ !is.na( logCounts$Percent_Loss_kg ), ] # remove rows with percent change NAs
    
    results.df.list <- list()
    df.index <- 1
    for ( t in 2:length(months) ) {
      
      timeFrame <- months[1:t]
      log.df <- logCounts[ logCounts$time %in% timeFrame, ]
      startAbundanceIndex <- which(colnames(log.df)=="ResponderStatus")+1
      divMonth <- paste0(division, "Assignment_", endM, "M")
      
      ## Parse out metadata from counts
      myT<-log.df[,startAbundanceIndex:ncol(log.df)]

      PatientID <- log.df$PatientID
      Timepoint <- log.df$time
      Surgery <- log.df$Surgery
      Weight_kg <- log.df$Weight_kg
      Division <- log.df[,which(colnames(log.df) == divMonth)]
      
      bugName <- vector()
      pVal_Timepoint <- vector()
      pVal_Surgery <- vector()
      pVal_Weight_kg <- vector()
      
      index <- 1
      
      for (i in 1:ncol(myT)){
        
        bug<-myT[,i]
        
        if (mean(bug>0, na.rm = TRUE)>0.1){
          
          df <- data.frame(bug, PatientID, Timepoint, Surgery, Weight_kg, Division)
          df <- na.omit(df)
          df$Timepoint <- as.factor(df$Timepoint)

          mlm <- lme( bug ~ Timepoint + Surgery +  Weight_kg, method = "REML", random = ~1 | PatientID, data = df )
          fit<-anova(mlm)

          pVal_Timepoint[index] <- fit$"p-value"[2]
          pVal_Surgery[index] <- fit$"p-value"[3]
          pVal_Weight_kg[index] <- fit$"p-value"[4]

          bugName[index]<-colnames(myT)[i]

          index<-index+1
          
        } # if (mean(bug>0, na.rm = TRUE)>0.1)
        
      } # for (i in 1:ncol(myT))
      
      dFrame<-data.frame(pVal_Timepoint, pVal_Surgery, pVal_Weight_kg)
      dFrame$Adj_pVal_Timepoint <- p.adjust(dFrame$pVal_Timepoint,method = "BH")
      dFrame$Adj_pVal_Surgery <- p.adjust(dFrame$pVal_Surgery,method = "BH")
      dFrame$Adj_pVal_Weight_kg <- p.adjust(dFrame$pVal_Weight_kg,method = "BH")
      
      names(dFrame) <- paste0(names(dFrame), "_BLto", timeFrame[length(timeFrame)],"M")
      dFrame <- data.frame(bugName, dFrame)
      
      results.df.list[[df.index]] <- dFrame
      df.index <- df.index + 1
      
    } # for ( t in 2:length(months) )
    
    merged.results <- results.df.list[[1]]
    if (length(results.df.list) > 1) {
      for (l in 2:length(results.df.list)) {
        
        merged.results <- merge(merged.results, results.df.list[[l]], by = "bugName")
        
      } # for (l in 2:length(results.df.list))
      
    }
    write.table(merged.results, paste0(outputLevel, level,  "_by_Timepoint_Surgery_Weight_kg__mixedLM_BLto", months[length(months)], "M_", classifier, "_RESULTS.tsv"),sep="\t",row.names = FALSE,quote = FALSE)
    
    
    file.path <- paste0(outputWLdivisions, level, "_by_Timepoint_Surgery_Weight_kg__mixedLM_", division, "_BLto", months[length(months)], "M_", classifier, "_PLOTS.pdf")
    pdf(file.path)

    par(mfrow=c(2,2))

    divPValues <- vector()
    groupNumbers <- vector()
    qbugName <- vector()
    direction <- vector()
    LengthTime <- vector()
    qIndex <- 1

    for( i in 1:nrow(merged.results)) {
      BugName <- merged.results$bugName[i]
      
      for ( t in 2:length(months) ) {
        
        timeFrame <- months[1:t]
        log.df <- logCounts[ logCounts$time %in% timeFrame, ]
        bugCol <- log.df[ , which(colnames(log.df) == BugName) ]
        Month <- log.df$time
        Surgery <- log.df$Surgery
        Weight <- log.df$Weight_kg
        PatientID <- log.df$PatientID
        divCol <- log.df[,which(colnames(log.df) == divMonth)]
        df <- data.frame(PatientID, Month, Surgery, Weight, bugCol)
        df <- na.omit(df)
        
        colEnding <- paste0("BLto", timeFrame[length(timeFrame)],"M")
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
        plot(  df$Month, df$bugCol,
                   ylab =BugName ,
                   xlab="Timepoint (months)",
                   main=title.lab,
                   cex.main = 1 ,
                   col = ifelse( df$Surgery == "SG", "red", "blue" )
            )
        legend("topright", legend=c("SG", "RYGB"),
               col=c("red", "blue"), pch = 21, cex=0.8)

        title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Surgery_p, length = 3), "; n = ", n)
        boxplot( df$bugCol ~  df$Surgery ,
                 xlab="Surgery Type",
                 ylab =BugName,
                 main = title.lab,
                 cex.main = 1
        )

        if (division == "Quintile") {groups <- 1:5}
        if (division == "Quartile") {groups <- 1:4}
        if (division == "Tertile") {groups <- 1:3}
        if (division == "Half") {groups <- 1:2}
          
        qFrameFinal <- data.frame()
          
          for( j in groups ) {
            
            lm.df <- data.frame(PatientID, Month, Surgery, Weight, bugCol, divCol)
            lm.df <- na.omit(lm.df)
            lm.df <- lm.df[lm.df$divCol == j,]
  
            myLm <- lm(  lm.df$bugCol ~ lm.df$Month)
            p = anova(myLm)$"Pr(>F)"[1]
            divPValues[qIndex] <- p
            groupNumbers[qIndex] <- j
            qbugName[qIndex] <- BugName
  
            title.lab <- paste0( division, " ", j, " (", length(unique(lm.df$PatientID)), " patients)\n", BugName , "\nTimepoint (p = ", format(p, length = 3), ")")
            plot(  lm.df$Month, lm.df$bugCol ,
                   ylab =BugName ,
                   xlab="Timepoint (months)",
                   main=title.lab,
                   cex.main = 1,
                   col = ifelse( lm.df$Surgery == "SG", "red", "blue" )
            )
            legend("topright", legend=c("SG", "RYGB"),
                   col=c("red", "blue"), pch = 21, cex=0.8)
  
            averages <- lm.df %>%
              group_by(Month) %>%
              get_summary_stats(bugCol, type = "mean_sd")
  
            BL_avg <- averages$mean[1]
            end_avg <- averages$mean[nrow(averages)]
  
            direction[qIndex] <- ifelse(BL_avg > end_avg, "decrease",
                                        ifelse(BL_avg < end_avg, "increase",
                                               "same"))
            LengthTime[qIndex] <- colEnding 
  
            qIndex <- qIndex + 1
  
          } # for( j in groups )
            
        qFrame <- data.frame( qbugName, LengthTime, groupNumbers,divPValues, direction)
        qFrame$adjPValues <- p.adjust(qFrame$divPValues, method="BH")
        
        qFrameFinal <- rbind(qFrameFinal, qFrame)

        } # for ( t in 2:length(months) )

    } # for( i in 1:nrow(merged.results))
    dev.off()


    
    qFrameFinal <- qFrameFinal[ order( qFrameFinal$divPValues),]

    file.path <- paste0(outputWLdivisions, level, "_by_Timepoint_", division, "_LM_pValues_BLto", months[length(months)], "M_", classifier, ".tsv")
    write.table(qFrameFinal, file=file.path, sep="\t", row.names=FALSE)


    file.path <- paste0(outputLevel, level,  "_by_Timepoint_Surgery_Weight_kg__MLM_pValHistogram_BLto", months[length(months)], "M_", classifier, ".pdf")
    pdf(file.path)
    
    for ( t in 2:length(months) ) {
      
      timeFrame <- months[1:t]
      colEnding <- paste0("BLto", timeFrame[length(timeFrame)],"M")
      
      pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Timepoint_", colEnding))]
      n=length(pValues)
      title.lab <- paste0(classifier, " ", colEnding, "\nMLM Timepoint p values (n = ", n, ")")
      hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
      
      pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Surgery_", colEnding))]
      n=length(pValues)
      title.lab <- paste0(classifier, " ", colEnding, "\nMLM Surgery p values (n = ", n, ")")
      hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
      
      pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Weight_kg_", colEnding))]
      n=length(pValues)
      title.lab <- paste0(classifier, " ", colEnding, "\nMLM Weight_kg p values (n = ", n, ")")
      hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
      
      
    } # for ( t in 2:length(months) )

    dev.off()
    
    dFrame.Final <- data.frame()
    for ( t in 1:length(months) ) {
      
      timeFrame <- months[t]
      log.df <- logCounts[ logCounts$time %in% timeFrame, ]
      Surgery <- log.df$Surgery
      Weight_kg <- log.df$Weight_kg
      PatientID <- log.df$PatientID
      Division <- log.df[,which(colnames(log.df) == divMonth)]
      
      myT<-log.df[,startAbundanceIndex:ncol(log.df)]
      
      pVal_Division <- vector()
      bugName <- vector()
      index <- 1
      
      for (i in 1:ncol(myT)){
        
        bug<-myT[,i]
        
        if (mean(bug>0, na.rm = TRUE)>0.1){
          
          df <- data.frame(bug, PatientID, Surgery, Weight_kg, Division)
          df <- na.omit(df)
          
          lm <- lm( bug ~ Division, data = df )
          fit<-anova(lm)
          
          pVal_Division[index] <- fit$`Pr(>F)`[1]
          
          bugName[index]<-colnames(myT)[i]
          Timepoint[index] <- t
          index<-index+1
          
        } # if (mean(bug>0, na.rm = TRUE)>0.1)
        
      } # for (i in 1:ncol(myT))
      
      dFrame<-data.frame(bugName, pVal_Division)
      dFrame$Adj_pVal_Division <- p.adjust(dFrame$pVal_Division,method = "BH")
      dFrame$Timepoint <- months[t]
      
      dFrame.Final <- rbind(dFrame.Final, dFrame)
    } # for ( t in 1:length(months) )
    
    file.path <- paste0(outputWLdivisions, level, "_", classifier, "_by_", division,"_Timepoint_LM_Results.tsv")
    write.table(dFrame.Final, file.path, sep="\t",quote = FALSE, row.names = FALSE)
    
    
    tukey.df <- data.frame()
    plotList <- list()
    index <- 1
    for (i in unique(dFrame.Final$bugName)) {
      
      dFrame2 <- dFrame.Final[dFrame.Final$bugName == i,]
      # logCounts2 <- logCounts[ !is.na( logCounts[which(colnames(logCounts) == divMonth)]), ]
      # logCounts2 <- logCounts2[ !is.na( logCounts2[which(colnames(logCounts2) == i)]), ]
      logCounts2 <- logCounts[logCounts$time %in% months,]
      
      Division <- logCounts2[ , which(colnames(logCounts2) == divMonth)]
      Time <- logCounts2$time
      bug <- logCounts2[, which(colnames(logCounts2) == i)]
      df2 <- data.frame(Division, Time, bug)
      df2 <- na.omit(df2)
      x.lab <- division
      y.lab <- i
      
      plot <- ggboxplot(
        df2, x = "Division", y = "bug", color = "black",
        fill = "Division", palette = c("blue", "pink", "red"),
        # facet.by = "month",
        scales = "free", add = "jitter"
      ); plot
      
      month.labsP <- paste0(month.labs, "\n", roundP(dFrame2$Adj_pVal_Division[1:nrow(dFrame2)]))
      plot <- facet(plot, facet.by = "Time"
                    , panel.labs = list(Time = month.labsP[1:length(months)])
      ); plot
      
      plot <- plot + labs(x=x.lab, y = y.lab)+ theme(legend.position = "none"); plot
      
      df2$Division <- as.factor(df2$Division)
      stat.test <- df2 %>%
        group_by(Time) %>%
        tukey_hsd(bug ~ Division)
      
      tukey.df <- rbind(tukey.df, stat.test)
      
      stat.test <- stat.test %>%
        add_xy_position(x = "Division", fun = "max")
      
      plot <- plot +
        stat_pvalue_manual(
          stat.test,
          bracket.nudge.y = 0.07,
          step.increase = .1,
          size = 7,
          hide.ns = TRUE,
          label = "p.adj.signif"
        ); plot
      
      plotList[[index]] <- plot
      index <- index + 1
    }
    
    file.path <- paste0(outputWLdivisions, level, "_", classifier, "_by_", division, "_Timepoint_BoxPlot.pdf")
    pdf(file.path, width = 10, height = 7)
    for (x in 1:length(plotList)) {
      print(plotList[[x]])
    }
    dev.off()
    
    
    file.path <- paste0(outputWLdivisions, level, "_", classifier, "_by_", division,"_Timepoint_Tukey_Results.tsv")
    write.table(tukey.df, file.path, sep="\t",quote = FALSE, row.names = FALSE)
    
    
  } # for (division in divisions)
  
} # for (level in levels)

