#Author: Alicia Sorgen
#Date: 2022 June 09
#Description: Linear modeling for taxa by time

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
params <- c(params, "MetaPhlAn2")
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

moduleRoot <- "Taxa_over_time"

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
    module <- paste0(args[2], "_", module)
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
# divisions <- c("Quintile", "Quartile", "Tertile", "Half")
divisions <- c("Tertile")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Univariate linear model #####
for (level in levels) {

  message(paste0("\n***** Starting ", level, " univariate analysis *****\n"))

  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)

  file.path <- paste0(inputDir,level, logCountFile)
  logCounts <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

  startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1

  # Parse out metadata from counts
  myT<-logCounts[,startAbundanceIndex:ncol(logCounts)]

  PatientID <- logCounts$PatientID
  Timepoint <- logCounts$time
  bugName <- vector()
  pValue.MLM <- vector()
  pValue.LM <- vector()
  p_BLv1M <- vector()
  p_BLv6M <- vector()
  p_BLv12M <- vector()
  p_BLv18M <- vector()
  p_BLv24M <- vector()

  index <- 1

  for (i in 1:ncol(myT)){

    bug<-myT[,i]

    if (mean(bug>0, na.rm = TRUE)>0.1){

      df<-data.frame(bug,PatientID, Timepoint)
      df <- df[df$Timepoint %in% included,]
      df <- na.omit(df)
      df$Timepoint <- as.factor(df$Timepoint)

      mlm <- lme(bug ~ Timepoint, method="REML", random=~1 | PatientID, data=df)
      smry <- summary(mlm)
      fit<-anova(mlm)
      pValue.MLM[index]<-fit$`p-value`[2]

      lm <- lm(df$bug ~ df$Timepoint)
      smry <- summary(lm)
      fit<-anova(lm)
      pValue.LM[index]<-fit$`Pr(>F)`[1]
      p_BLv1M[index] <- smry$coefficients[2,4]
      if (length(included) > 2) {  p_BLv6M[index] <- smry$coefficients[3,4]}
      if (length(included) > 3) {  p_BLv12M[index] <- smry$coefficients[4,4]}
      if (length(included) > 4) {  p_BLv18M[index] <- smry$coefficients[5,4]}
      if (length(included) > 5) {  p_BLv24M[index] <- smry$coefficients[6,4]}
      bugName[index]<-colnames(myT)[i]

      index<-index+1

    } #if (mean(bug>0, na.rm = TRUE)>0.1)

  } #for (i in 1:ncol(myT))


  taxa.df<-data.frame(bugName, pValue.MLM, pValue.LM, p_BLv1M)
  taxa.df$Adj_pValue.MLM <- p.adjust(taxa.df$pValue.MLM,method = "BH")
  taxa.df$MLM_sig <- sigStars(taxa.df$Adj_pValue.MLM)
  taxa.df$Adj_pValue.LM <- p.adjust(taxa.df$pValue.LM,method = "BH")
  taxa.df$LM_sig <- sigStars(taxa.df$Adj_pValue.LM)
  taxa.df$Adj_p_BLv1M<- p.adjust(taxa.df$p_BLv1M,method = "BH")
  taxa.df$BLv1M_sig <- sigStars(taxa.df$Adj_p_BLv1M)

  if (length(included) > 2) {
    taxa.df$p_BLv6M <- p_BLv6M
    taxa.df$Adj_p_BLv6M <- p.adjust(taxa.df$p_BLv6M,method = "BH")
    taxa.df$BLv6M_sig <- sigStars(taxa.df$Adj_p_BLv6M)
    }

  if (length(included) > 3) {
    taxa.df$p_BLv12M <- p_BLv12M
    taxa.df$Adj_p_BLv12M <- p.adjust(taxa.df$p_BLv12M,method = "BH")
    taxa.df$BLv12M_sig <- sigStars(taxa.df$Adj_p_BLv12M)
  }

  if (length(included) > 4) {
    taxa.df$p_BLv18M <- p_BLv18M
    taxa.df$Adj_p_BLv18M <- p.adjust(taxa.df$p_BLv18M,method = "BH")
    taxa.df$BLv18M_sig <- sigStars(taxa.df$Adj_p_BLv18M)
  }

  if (length(included) > 5) {
    taxa.df$p_BLv24M <- p_BLv24M
    taxa.df$Adj_p_BLv24M <- p.adjust(taxa.df$p_BLv24M,method = "BH")
    taxa.df$BLv24M_sig <- sigStars(taxa.df$Adj_p_BLv24M)
  }

  pVal.FileName <- paste0("_by_Timepoint_mixedLM_and_LM_BLto", included[length(included)], "M_", classifier, ".tsv")
  file.path <- paste0(outputLevel, level, pVal.FileName)
  write.table(taxa.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)




  file.path <- paste0(outputLevel, level, "_by_Timepoint_mixedLM_and_LM_pValHistogram_BLto", included[length(included)], "M_", classifier, ".pdf")

  pdf(file.path)

  pValues <- taxa.df[,which(colnames(taxa.df) == "pValue.MLM")]
  n=length(pValues)
  title.lab <- paste0(classifier, ": ", level, " level\nMixed LM (bug ~ Timepoint) p values (n = ", n, ")")
  hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))

  pValues <- taxa.df[,which(colnames(taxa.df) == "pValue.LM")]
  n=length(pValues)
  title.lab <- paste0(classifier, ": ", level, " level\nLM (bug ~ Timepoint) p values (n = ", n, ")")
  hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))

  pValues <- taxa.df[,which(colnames(taxa.df) == "p_BLv1M")]
  n=length(pValues)
  title.lab <- paste0(classifier, ": ", level, " level\nLinear model BL v. 1 month p values (n = ", n, ")")
  hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))

  if (length(included) > 2) {
    pValues <- taxa.df[,which(colnames(taxa.df) == "p_BLv6M")]
    n=length(pValues)
    title.lab <- paste0(classifier, ": ", level, " level\nLinear model BL v. 6 month p values (n = ", n, ")")
    hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
  }

  if (length(included) > 3) {
    pValues <- taxa.df[,which(colnames(taxa.df) == "p_BLv12M")]
    n=length(pValues)
    title.lab <- paste0(classifier, ": ", level, " level\nLinear model BL v. 12 month p values (n = ", n, ")")
    hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
  }

  if (length(included) > 4) {
    pValues <- taxa.df[,which(colnames(taxa.df) == "p_BLv18M")]
    n=length(pValues)
    title.lab <- paste0(classifier, ": ", level, " level\nLinear model BL v. 18 month p values (n = ", n, ")")
    hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
  }

  if (length(included) > 5) {
    pValues <- taxa.df[,which(colnames(taxa.df) == "p_BLv24M")]
    n=length(pValues)
    title.lab <- paste0(classifier, ": ", level, " level\nLinear model BL v. 24 month p values (n = ", n, ")")
    hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
  }

  dev.off()



} # for (level in levels)

##### Univariate linear model plots #####
for (level in levels) {
  
  message(paste0("\n***** Starting ", level, " univariate plots *****\n"))
  # level <- "phylum"
  outputLevel <- paste0(outputDir, level, "/")
  
  file.path <-  paste0(inputDir,level, logCountFile)
  logCounts <- read.table(file.path, sep="\t",header = TRUE, check.names = FALSE)
  logCounts2 <- logCounts[logCounts$time %in% included,]
  
  
  file.path <- paste0(outputLevel, level, pVal.FileName)
  pVal.df <- read.table(file.path, sep="\t",header = TRUE, check.names = FALSE)
  pVal.df$pFinal <- roundP(pVal.df$Adj_pValue.MLM)
  
  ordered <- pVal.df[order(pVal.df$Adj_pValue.MLM),]
  plotList <- list()
  
  for (row in 1:nrow(ordered)) {
    
    bugName <- ordered$bugName[row]
    
    bugCol <- logCounts2[,which(colnames(logCounts2) == bugName)]
    month <- logCounts2$time
    df <- data.frame(bugCol, month)
    df <- na.omit(df)
    
    stat.test <- df %>%
      t_test(bugCol ~ month, ref.group = "0") %>%
      adjust_pvalue(method = "BH") %>%
      add_significance()
    
    stat.test <- stat.test %>%
      add_xy_position(x = "month", fun = "max")
    
    # stat.test$y.position[2] <- stat.test$y.position[1] + 0.3
    # stat.test$y.position[3] <- stat.test$y.position[2] + 0.3
    
    
    plot <- ggboxplot(
      df, x = "month", y = "bugCol",
      # fill = "#2EDF52",
      # facet.by = c("Surgery"),
      scales = "free", add = "jitter"
    )
    
    plot <- plot +
      labs(x = "Post-op Time (months)", y = "Abundance (log)")
    
    plot <- plot +
      labs(title = paste0(bugName, " (", ordered$pFinal[row], ")"))

    plot <- plot +
      labs(caption = "Pairwise comparisons calculated using t-tests")
    
    plot <- plot +
      stat_pvalue_manual(
        stat.test,
        bracket.nudge.y = 0.5,
        size = 6,
        hide.ns = TRUE,
        label = "p.adj.signif",
        step.increase = 0.07
      )
    plot
    
    plotList[[row]] <- plot
    
  } # for (row in 1:nrow(ordered))
  
  
  
  # Output plot
  file.path <- paste0(outputLevel, level, "_by_Timepoint_MLM_BLto", included[length(included)], "M_Boxplots_", classifier, ".pdf")
  pdf(file.path, width = 7, height = 4)
  for (i in 1:length(plotList)) {
    grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
  }
  dev.off()
  
} # for (level in levels)

##### RYGB v SG at each timepoint (t-test) #####
for (level in levels) {
  
  outputLevel <- paste0(outputDir, level, "/")
  
  file.path <-  paste0(inputDir,level, logCountFile)
  logCounts <- read.table(file.path, sep="\t",header = TRUE, check.names = FALSE)
  
  df <- logCounts[logCounts$time %in% included,]
  startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1
  myT<-df[,startAbundanceIndex:ncol(df)]
  
  
  month <- df$time
  surgery <- df$Surgery
  
  # index <- 1
  t.testResults <- data.frame()
  plotList <- list()
  index <- 1
  
  for (i in 1:ncol(myT)){
    
    bug <- myT[,i]
    
    if (mean(bug>0, na.rm = TRUE)>0.1){
      
      df2<-data.frame(bug, surgery, month)
      df2$month <- as.factor(df2$month)
      df2 <- na.omit(df2)
      
      stat.test <- df2 %>%
        group_by(month) %>%
        t_test(bug ~ surgery) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance()
      
      stat.test$.y.<-colnames(myT)[i]
      
      t.testResults <- rbind(t.testResults, stat.test)
      
      stat.test <- stat.test %>%
        add_xy_position(x = "surgery", fun = "max")
      
      plot <- ggboxplot(
        df2, x = "surgery", y = "bug", color = "month",
        # fill = "#2EDF52",
        # facet.by = "month",
        scales = "free", add = "jitter"
      )
      
      plot <- facet(plot, facet.by = "month", panel.labs = list(month = month.labs[1:length(included)]))
      
      plot <- plot +
        labs(x = "Surgery Type", y = "Log Abundance")
      
      plot <- plot +
        labs(title = colnames(myT)[i])
      
      caption.lab <- "* significance bars determined by t-test"
      plot <- plot +
        labs(caption = caption.lab)
      
      plot <- plot +
        stat_pvalue_manual(
          stat.test,
          bracket.nudge.y = 0.05,
          size = 5,
          hide.ns = TRUE,
          label = "p.adj.signif"
        )
      
      plot <- plot + theme(legend.position = "none")
      
      plotList[[index]] <- plot
      index <- index + 1
      
    } #if (mean(bug>0, na.rm = TRUE)>0.1)
    
  } # for (i in 1:ncol(myT))
  
  write.table(t.testResults, paste0(outputLevel, level, "_by_Surgery_ttest_BLto", included[length(included)], "M_Results_", classifier, ".tsv"),sep="\t",row.names = FALSE,quote = FALSE)
  
  # Output plot
  file.path <- paste0(outputLevel, level, "_by_Surgery_ttest_BLto", included[length(included)], "M_Boxplots_", classifier, ".pdf")
  pdf(file.path, width = 7, height = 5)
  for (i in 1:length(plotList)) {
    grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
  }
  dev.off()
  
} # for (level in levels)

##### Multivariate models (Weight loss groups) #####
for (level in levels) {
  
  outputLevel <- paste0(outputDir, level, "/")
  
  for (division in divisions) {
    
    message(paste0("\n***** Starting ", level, " multivariate (", division, ") *****\n"))
    
    outputWLdivisions <- paste0(outputLevel, division, "/")
    dir.create(outputWLdivisions, showWarnings = FALSE)
    
    file.path <-  paste0(inputDir,level, logCountFile)
    logCounts <- read.table(file.path, sep="\t",header = TRUE, check.names = FALSE)
    startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1
    
    
    logCounts <- logCounts[ !is.na( logCounts[,startAbundanceIndex+1]), ] # remove rows with bug NAs
    logCounts <- logCounts[ !is.na( logCounts$Weight_kg), ] # remove rows with Weight NAs
    logCounts <- logCounts[ !is.na( logCounts$Percent_Loss_kg ), ] # remove rows with percent change NAs
    logCounts <- logCounts[ logCounts$time %in% included, ]
    
    startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1
    divMonth <- which(colnames(logCounts) == paste0(division, "Assignment_", included[length(included)], "M"))
    
    ## Parse out metadata from counts
    myT<-logCounts[,startAbundanceIndex:ncol(logCounts)]
    
    PatientID <- logCounts$PatientID
    Timepoint <- logCounts$time
    Surgery <- logCounts$Surgery
    Weight_kg <- logCounts$Weight_kg
    
    bugName <- vector()
    pVal_Timepoint <- vector()
    pVal_Surgery <- vector()
    pVal_Weight_kg <- vector()
    
    index <- 1
    
    for (i in 1:ncol(myT)){
      
      bug<-myT[,i]
      
      if (mean(bug>0, na.rm = TRUE)>0.1){
        
        df<-data.frame(bug,PatientID, Timepoint, Surgery, Weight_kg)
        df <- na.omit(df)
        df$Timepoint <- as.factor(df$Timepoint)
        
        mlm <- lme( bug ~ Timepoint + Surgery +  Weight_kg, method = "REML", random = ~1 | PatientID, data = df )
        fit<-anova(mlm)
        
        pVal_Timepoint[index] <- fit$"p-value"[2]
        pVal_Surgery[index] <- fit$"p-value"[3]
        pVal_Weight_kg[index] <- fit$"p-value"[4]
        
        bugName[index]<-colnames(myT)[i]
        
        index<-index+1
        
      } #if (mean(bug>0, na.rm = TRUE)>0.1)
      
    } #for (i in 1:ncol(myT))
    
    
    dFrame<-data.frame(bugName, pVal_Timepoint, pVal_Surgery, pVal_Weight_kg)
    
    dFrame$Adj_pVal_Timepoint <- p.adjust(dFrame$pVal_Timepoint,method = "BH")
    
    dFrame$Adj_pVal_Surgery <- p.adjust(dFrame$pVal_Surgery,method = "BH")
    
    dFrame$Adj_pVal_Weight_kg <- p.adjust(dFrame$pVal_Weight_kg,method = "BH")
    
    
    write.table(dFrame, paste0(outputLevel, level,  "_by_Timepoint_Surgery_Weight_kg__mixedLM_BLto", included[length(included)], "M_", classifier, "_RESULTS.tsv"),sep="\t",row.names = FALSE,quote = FALSE)
    
    
    file.path <- paste0(outputWLdivisions, level, "_by_Timepoint_Surgery_Weight_kg__mixedLM_", division, "_BLto", included[length(included)], "M_", classifier, "_PLOTS.pdf")
    pdf(file.path)
    
    par(mfrow=c(2,2))
    
    divPValues <- vector()
    groupNumbers <- vector()
    qbugName <- vector()
    direction <- vector()
    qIndex <- 1
    
    for( i in 1:nrow(dFrame)) {
      dataCol <- which( dFrame$bugName[i] == names(logCounts) )
      
      title.lab <- paste0( dFrame$bugName[i] , "\nWeight (adj p = ", format(dFrame$Adj_pVal_Weight_kg[i], length = 3), ")")
      plot(  logCounts$Weight_kg, logCounts[,dataCol],
             ylab =dFrame$bugName[i] ,
             xlab="Weight (kg)",
             main=title.lab,
             cex.main = 1 ,
             col = ifelse( logCounts$Surgery == "SG", "red", "blue" )
      )
      legend("topright", legend=c("SG", "RYGB"),
             col=c("red", "blue"), pch = 21, cex=0.8)
      
      title.lab <- paste0( dFrame$bugName[i] , "\nTimepoint (adj p = ", format(dFrame$Adj_pVal_Timepoint[i], length = 3), ")")
      plot(  logCounts$time, logCounts[,dataCol],
             ylab =dFrame$bugName[i] ,
             xlab="Timepoint (months)",
             main=title.lab,
             cex.main = 1 ,
             col = ifelse( logCounts$Surgery == "SG", "red", "blue" )
      )
      legend("topright", legend=c("SG", "RYGB"),
             col=c("red", "blue"), pch = 21, cex=0.8)
      
      title.lab <- paste0( dFrame$bugName[i] , "\nSurgery (adj p = ", format(dFrame$Adj_pVal_Surgery[i], length = 3), ")")
      boxplot( logCounts[,dataCol] ~  logCounts$Surgery ,
               xlab="Surgery Type",
               ylab =dFrame$bugName[i],
               main = title.lab,
               cex.main = 1
      )
      
      if (division == "Quintile") {groups <- 1:5}
      if (division == "Quartile") {groups <- 1:4}
      if (division == "Tertile") {groups <- 1:3}
      if (division == "Half") {groups <- 1:2}
      
      for( j in groups ) {
        bugCol <- logCounts[,dataCol]
        timeCol <- logCounts$time
        divCol <- logCounts[,divMonth]
        surgCol <- logCounts$Surgery
        id <- logCounts$PatientID
        lm.df <- data.frame(id, bugCol, timeCol, divCol, surgCol)
        lm.df <- na.omit(lm.df)
        
        lm.df <- lm.df[lm.df$divCol == j,]
        
        myLm <- lm(  lm.df$bugCol ~ lm.df$timeCol)
        p = anova(myLm)$"Pr(>F)"[1]
        divPValues[qIndex] <- p
        groupNumbers[qIndex] <- j
        qbugName[qIndex] <- dFrame$bugName[i]
        
        title.lab <- paste0( division, " ", j, " (", length(unique(lm.df$id)), " patients)\n", dFrame$bugName[i] , "\nTimepoint (adj p = ", format(p, length = 3), ")")
        plot(  lm.df$timeCol, lm.df$bugCol ,
               ylab =dFrame$bugName[i] ,
               xlab="Timepoint (months)",
               main=title.lab,
               cex.main = 1,
               col = ifelse( lm.df$surgCol == "SG", "red", "blue" )
        )
        legend("topright", legend=c("SG", "RYGB"),
               col=c("red", "blue"), pch = 21, cex=0.8)
        
        averages <- lm.df %>%
          group_by(timeCol) %>%
          get_summary_stats(bugCol, type = "mean_sd")
        
        BL_avg <- averages$mean[1]
        end_avg <- averages$mean[nrow(averages)]
        
        direction[qIndex] <- ifelse(BL_avg > end_avg, "decrease",
                                    ifelse(BL_avg < end_avg, "increase",
                                           "same"))
        
        qIndex <- qIndex + 1
        
      }
    }
    dev.off()
    
    
    qFrame <- data.frame( qbugName, groupNumbers,divPValues, direction)
    qFrame <- qFrame[ order( qFrame$divPValues),]
    qFrame$adjPValues <- p.adjust(qFrame$divPValues, method="BH")
    
    file.path <- paste0(outputWLdivisions, level, "_by_Timepoint_", division, "_LM_pValues_BLto", included[length(included)], "M_", classifier, ".tsv")
    write.table(qFrame, file=file.path, sep="\t", row.names=FALSE)
    
    
    file.path <- paste0(outputLevel, level,  "_by_Timepoint_Surgery_Weight_kg__MLM_pValHistogram_BLto", included[length(included)], "M_", classifier, ".pdf")
    pdf(file.path)
    
    pValues <- dFrame[,which(colnames(dFrame) == "pVal_Timepoint")]
    n=length(pValues)
    title.lab <- paste0(classifier, "\nMLM Timepoint p values (n = ", n, ")")
    hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
    
    pValues <- dFrame[,which(colnames(dFrame) == "pVal_Surgery")]
    n=length(pValues)
    title.lab <- paste0(classifier, "\nMLM Surgery p values (n = ", n, ")")
    hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
    
    pValues <- dFrame[,which(colnames(dFrame) == "pVal_Weight_kg")]
    n=length(pValues)
    title.lab <- paste0(classifier, "\nMLM Weight_kg p values (n = ", n, ")")
    hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
    
    dev.off()
    
    
  } # for (division in divisions)
  
} # for (level in levels)

##### Univariate taxa change linear model #####
for (level in levels) {
  
  message(paste0("\n***** Starting ", level, " change univariate analysis *****\n"))
  
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
  bugName <- vector()
  pValue.MLM <- vector()
  pValue.LM <- vector()
  p_BLv1M <- vector()
  p_BLv6M <- vector()
  p_BLv12M <- vector()
  p_BLv18M <- vector()
  p_BLv24M <- vector()
  
  index <- 1
  
  for (i in 1:ncol(myT)){
    
    bug<-myT[,i]
    
    if (mean(bug>0, na.rm = TRUE)>0.1){
      
      df<-data.frame(bug,PatientID, Timepoint)
      # df <- df[df$Timepoint %in% included,]
      # df <- na.omit(df)

      # baseline rows
      bl = which(df$Timepoint==0)
      
      # baseline abundance
      BL_bug <- df[bl, "bug"]
      names(BL_bug) <- df[bl, "PatientID"]
      
      # Add bug info to each row
      df$BL_bug <- BL_bug[df$PatientID]
      
      # Remove baseline rows (won't show any change)
      df <- df[df$Timepoint != 0,]
      df$bugChange <- df$bug - df$BL_bug
      df$bugPercentChange <- ( df$bugChange / df$BL_bug ) * 100
      
      df$Timepoint <- as.factor(df$Timepoint)
      
      df <- df[!(is.na(df$bugChange)),]
      
      mlm <- lme(bugChange ~ Timepoint, method="REML", random=~1 | PatientID, data=df)
      smry <- summary(mlm)
      fit<-anova(mlm)
      pValue.MLM[index]<-fit$`p-value`[2]
      
      lm <- lm(df$bugChange ~ df$Timepoint)
      smry <- summary(lm)
      fit<-anova(lm)
      pValue.LM[index]<-fit$`Pr(>F)`[1]
      # p_BLv1M[index] <- smry$coefficients[2,4]
      # if (length(included) > 2) {  p_BLv6M[index] <- smry$coefficients[3,4]}
      # if (length(included) > 3) {  p_BLv12M[index] <- smry$coefficients[4,4]}
      # if (length(included) > 4) {  p_BLv18M[index] <- smry$coefficients[5,4]}
      # if (length(included) > 5) {  p_BLv24M[index] <- smry$coefficients[6,4]}
      bugName[index]<-colnames(myT)[i]
      
      index<-index+1
      
    } #if (mean(bug>0, na.rm = TRUE)>0.1)
    
  } #for (i in 1:ncol(myT))
  
  
  # taxa.df<-data.frame(bugName, pValue.MLM, pValue.LM, p_BLv1M)
  taxa.df<-data.frame(bugName, pValue.MLM, pValue.LM)
  taxa.df$Adj_pValue.MLM <- p.adjust(taxa.df$pValue.MLM,method = "BH")
  taxa.df$MLM_sig <- sigStars(taxa.df$Adj_pValue.MLM)
  taxa.df$Adj_pValue.LM <- p.adjust(taxa.df$pValue.LM,method = "BH")
  taxa.df$LM_sig <- sigStars(taxa.df$Adj_pValue.MLM)
  # taxa.df$Adj_p_BLv1M<- p.adjust(taxa.df$p_BLv1M,method = "BH")
  # taxa.df$BLv1M_sig <- sigStars(taxa.df$Adj_pValue.MLM)
  
  # if (length(included) > 2) {
  #   taxa.df$p_BLv6M <- p_BLv6M
  #   taxa.df$Adj_p_BLv6M <- p.adjust(taxa.df$p_BLv6M,method = "BH")
  #   taxa.df$BLv6M_sig <- sigStars(taxa.df$Adj_pValue.MLM)
  # }
  # 
  # if (length(included) > 3) {
  #   taxa.df$p_BLv12M <- p_BLv12M
  #   taxa.df$Adj_p_BLv12M <- p.adjust(taxa.df$p_BLv12M,method = "BH")
  #   taxa.df$BLv12M_sig <- sigStars(taxa.df$Adj_pValue.MLM)
  # }
  # 
  # if (length(included) > 4) {
  #   taxa.df$p_BLv18M <- p_BLv18M
  #   taxa.df$Adj_p_BLv18M <- p.adjust(taxa.df$p_BLv18M,method = "BH")
  #   taxa.df$BLv18M_sig <- sigStars(taxa.df$Adj_pValue.MLM)
  # }
  # 
  # if (length(included) > 5) {
  #   taxa.df$p_BLv24M <- p_BLv24M
  #   taxa.df$Adj_p_BLv24M <- p.adjust(taxa.df$p_BLv24M,method = "BH")
  #   taxa.df$BLv24M_sig <- sigStars(taxa.df$Adj_pValue.MLM)
  # }
  
  file.path <- paste0(outputLevel, level, "Change_by_Timepoint_mixedLM_and_LM_Results_BLto", included[length(included)], "M_", classifier, ".tsv")
  write.table(taxa.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)
  
  
  
  
  file.path <- paste0(outputLevel, level, "Change_by_Timepoint_MLM_and_LM_pValHistogram_BLto", included[length(included)], "M_", classifier, ".pdf")
  
  pdf(file.path)
  
  pValues <- taxa.df[,which(colnames(taxa.df) == "pValue.MLM")]
  n=length(pValues)
  title.lab <- paste0(classifier, ": ", level, " level\nTimepoint MLM p values (n = ", n, ")")
  hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
  
  pValues <- taxa.df[,which(colnames(taxa.df) == "pValue.LM")]
  n=length(pValues)
  title.lab <- paste0(classifier, ": ", level, " level\nTimepoint LM p values (n = ", n, ")")
  hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
  
  # pValues <- taxa.df[,which(colnames(taxa.df) == "p_BLv1M")]
  # n=length(pValues)
  # title.lab <- paste0(classifier, ": ", level, " level\nLinear model BL v. 1 month p values (n = ", n, ")")
  # hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
  # 
  # if (length(included) > 2) {
  #   pValues <- taxa.df[,which(colnames(taxa.df) == "p_BLv6M")]
  #   n=length(pValues)
  #   title.lab <- paste0(classifier, ": ", level, " level\nLinear model BL v. 6 month p values (n = ", n, ")")
  #   hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
  # }
  # 
  # if (length(included) > 3) {
  #   pValues <- taxa.df[,which(colnames(taxa.df) == "p_BLv12M")]
  #   n=length(pValues)
  #   title.lab <- paste0(classifier, ": ", level, " level\nLinear model BL v. 12 month p values (n = ", n, ")")
  #   hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
  # }
  # 
  # if (length(included) > 4) {
  #   pValues <- taxa.df[,which(colnames(taxa.df) == "p_BLv18M")]
  #   n=length(pValues)
  #   title.lab <- paste0(classifier, ": ", level, " level\nLinear model BL v. 18 month p values (n = ", n, ")")
  #   hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
  # }
  # 
  # if (length(included) > 5) {
  #   pValues <- taxa.df[,which(colnames(taxa.df) == "p_BLv24M")]
  #   n=length(pValues)
  #   title.lab <- paste0(classifier, ": ", level, " level\nLinear model BL v. 24 month p values (n = ", n, ")")
  #   hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
  # }
  
  dev.off()
  
  
  
} # for (level in levels)


##### Top v Bottom WL Groups by timepoint #####
for (level in levels) {
  
  outputLevel <- paste0(outputDir, level, "/")
  
  for (division in divisions) {
    
    message(paste0("\n***** Starting ", level, " level ", division, " top v bottom *****\n"))
    
    outputWLdivisions <- paste0(outputLevel, division, "/")
    dir.create(outputWLdivisions, showWarnings = FALSE)
    
    file.path <-  paste0(inputDir,level, logCountFile)
    # file.path <-  paste0(outputDir,level, inputFileName)
    logCounts <- read.table(file.path, sep="\t",header = TRUE, check.names = FALSE)
    
    df <- logCounts[logCounts$time %in% included,]
    startAbundanceIndex <- which(colnames(df)=="ResponderStatus")+1
    myT<-df[,startAbundanceIndex:ncol(df)]
    
    divMonth <- which(colnames(df) == paste0(division, "Assignment_", included[length(included)], "M"))
    Group1 <- df$PatientID[which(df[,divMonth] == 1)]
    Group2 <- df$PatientID[which(df[,divMonth] == 2)]
    
    if (division == "Half") {
      
      WLGroup <- ifelse(df$PatientID %in% Group1, "Bottom WL Group",
                        ifelse(df$PatientID %in% Group2, "Top WL Group",
                               NA))
      
    }
    
    if (division == "Tertile") {
      
      Group3 <- df$PatientID[which(df[,divMonth] == 3)]
      
      WLGroup <- ifelse(df$PatientID %in% Group1, "Bottom WL Group",
                        ifelse(df$PatientID %in% Group2, "Group2",
                               ifelse(df$PatientID %in% Group3, "Top WL Group",
                                      NA)))
      
    }
    if (division == "Quartile") {
      
      Group3 <- df$PatientID[which(df[,divMonth] == 3)]
      Group4 <- df$PatientID[which(df[,divMonth] == 4)]
      
      WLGroup <- ifelse(df$PatientID %in% Group1, "Bottom WL Group",
                        ifelse(df$PatientID %in% Group2, "Group2",
                               ifelse(df$PatientID %in% Group3, "Group3",
                                      ifelse(df$PatientID %in% Group4, "Top WL Group",
                                             NA))))
    }
    
    if (division == "Quintile") {
      
      Group3 <- df$PatientID[which(df[,divMonth] == 3)]
      Group4 <- df$PatientID[which(df[,divMonth] == 4)]
      Group5 <- df$PatientID[which(df[,divMonth] == 5)]
      
      WLGroup <- ifelse(df$PatientID %in% Group1, "Bottom WL Group",
                        ifelse(df$PatientID %in% Group2, "Group2",
                               ifelse(df$PatientID %in% Group3, "Group3",
                                      ifelse(df$PatientID %in% Group4, "Group4",
                                             ifelse(df$PatientID %in% Group5, "Top WL Group",
                                                    NA)))))
    }
    
    month <- df$time
    
    # index <- 1
    t.testResults <- data.frame()
    plotList <- list()
    index <- 1
    
    for (i in 1:ncol(myT)){
      
      bug <- myT[,i]
      
      if (mean(bug>0, na.rm = TRUE)>0.1){
        
        df2<-data.frame(bug, WLGroup, month)
        df2 <- df2[df2$WLGroup %in% c("Bottom WL Group", "Top WL Group"),]
        df2$month <- as.factor(df2$month)
        df2 <- na.omit(df2)
        
        # message(">> t-test <<")
        stat.test <- df2 %>%
          group_by(month) %>%
          t_test(bug ~ WLGroup) %>%
          adjust_pvalue(method = "BH") %>%
          add_significance()
        
        stat.test$.y.<-colnames(myT)[i]
        
        t.testResults <- rbind(t.testResults, stat.test)
        
        stat.test <- stat.test %>%
          add_xy_position(x = "WLGroup", fun = "max")
        # message(">> t-test COMPLETE <<")
        
        
        plot <- ggboxplot(
          df2, x = "WLGroup", y = "bug", color = "month",
          # fill = "#2EDF52",
          # facet.by = "month",
          scales = "free", add = "jitter"
        )
        
        plot <- facet(plot, facet.by = "month", panel.labs = list(month = month.labs[1:length(included)]))
        
        plot <- plot +
          labs(x = paste0("Weight Loss ", division, "s"), y = "Log Abundance")
        
        plot <- plot +
          labs(title = colnames(myT)[i])
        
        plot <- plot +
          stat_pvalue_manual(
            stat.test,
            bracket.nudge.y = 0.05,
            size = 5,
            hide.ns = TRUE,
            label = "p.adj.signif"
          )
        
        plot <- plot + theme(legend.position = "none")
        
        plotList[[index]] <- plot
        index <- index + 1
      } # if (mean(bug>0, na.rm = TRUE)>0.1)
      
    } # for (i in 1:ncol(myT))
    
    write.table(t.testResults, paste0(outputWLdivisions, level, "_by_Top_Bottom_", division, "_ttest_BLto", included[length(included)], "M_Results_", classifier, ".tsv"),sep="\t",row.names = FALSE,quote = FALSE)
    
    # Output plot
    file.path <- paste0(outputWLdivisions, level, "_by_Top_Bottom_", division, "_ttest_BLto", included[length(included)], "M_Boxplots_", classifier, ".pdf")
    pdf(file.path, width = 10, height = 5)
    for (i in 1:length(plotList)) {
      grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
    }
    dev.off()
    
    
    
  } # for (division in divisions)
} # for (level in levels)
