#Author: Alicia Sorgen
#Date: 2022 July 11
#Description: Modeling to determine if baseline taxa predicts weight outcomes.

# Set up ------------------------------------------------------------------
rm(list=ls())
set.seed(1989)
H1 = function(comment) {
   delim = "#"
   side = "#"
   fill = " "
   maxLen <- 90
   lineLen <- 60
   message("")
   message(paste0(strrep(delim, maxLen)))
   if (nchar(comment) > lineLen) {
      cut <- strsplit(comment, " ")
      line <- ""
      for (word in cut[[1]]) {
         if ((nchar(line) + 1 + nchar(word)) > lineLen) {
            edge1 <- round((maxLen - nchar(line))/2, 0) - 2
            edge2 <- maxLen - edge1 - nchar(line) - 4
            line_print <- paste0(side, strrep(fill, edge1), " ", line, " ", strrep(fill, edge2), side)
            message(line_print)
            line <- word
         } else {
            line <- paste(line, word)
         }
      }
      edge1 <- round((maxLen - nchar(line))/2, 0) - 2
      edge2 <- maxLen - edge1 - nchar(line) - 4
      line_print <- paste0(side, strrep(fill, edge1), " ", line, " ", strrep(fill, edge2), side)
      message(line_print)
   } else {
      edge1 <- round((maxLen - nchar(comment))/2, 0) - 2
      edge2 <- maxLen - edge1 - nchar(comment) - 4
      line_print <- paste0(side, strrep(fill, edge1), " ", comment, " ", strrep(fill, edge2), side)
      message(line_print)
   }
   message(paste0(strrep(delim, maxLen)))
   message("")
}

# Script Edits --------------------------------------------------------------------------------


params <- vector()
params <- c(params, getOption("mbs.pipe_root", default = stop("Set options(mbs.pipe_root = ...) in your local .Rprofile (gitignored)")))
params <- c(params, "MetaPhlAn2")
# params <- c(params, 0)
# params <- c(params, 1)
params <- c(params, 6)
# params <- c(params, 12)
# params <- c(params, 18)
# params <- c(params, 24)
endM <- params[length(params)]

levels <- vector()
levels <- c(levels, "phylum")
# levels <- c(levels, "class")
# levels <- c(levels, "order")
# levels <- c(levels, "family")
levels <- c(levels, "genus")
# levels <- c(levels, "species")

moduleRoot <- "Taxa_Changes"


# Packages ------------------------------------------------------------------------------------


library(stringr)
library(ggpubr)
library(tidyr)
library(rstatix)
library(nlme)
library(data.table)
library(gridExtra)

message("\n\n---------- [ R Version ] ----------\n")
R <- sessionInfo()
message(R$R.version$version.string)

message("\n\n---------- [ Packages ] ----------\n")
for (i in 1:length(R$otherPkgs)) {
  pkg <- R$otherPkgs[[i]]["Package"]
  v <- R$otherPkgs[[i]]["Version"]
  message(pkg, ":", v)
}
rm(pkg, R, v)



# Set working environment ---------------------------------------------------------
H1("Working Environment")
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
   args <- params
}

message("\n************* Running locally *************")
proj_root    <- args[1]
message("Project root directory: ", proj_root, "\n")
inputEnv <- Sys.getenv("INPUT_ROOT"); input_root <- if (nchar(inputEnv) > 0) inputEnv else file.path(dirname(proj_root), "Data")
script_root <- file.path(proj_root, basename(proj_root), "analysis", "Rscripts")

resultsEnv <- Sys.getenv("RESULTS_ROOT"); pipeRoot <- if (nchar(resultsEnv) > 0) resultsEnv else file.path(proj_root, gsub('Analysis', 'Results', basename(proj_root)))

module <- paste0("3.0_Taxa_Changes_", args[3], "M")
moduleDir <- file.path(pipeRoot, module)

outputDir <- file.path(moduleDir, "output/")

rm(proj_root, inputEnv, resultsEnv, module, params)



# Set functions ------------------------------------------------
H1("Functions")
message("Loading functions from: ", script_root)
funcScript <- file.path(script_root, "functions.R")
source(funcScript); rm(funcScript)


# Input ---------------------------------------------------------------------------------------


classifier <- args[2]
predictionMonth <- args[3]

prevModule <- paste0("TaxaMetaMerge")
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
logCountFile <- paste0("_LogNormalizedCounts_", classifier, ".tsv")
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")


# Output --------------------------------------------------------------------------------------


outputDir = file.path(moduleDir,"output/")


# Analysis ------------------------------------------------------------------------------------


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

      df <- df[(df$Timepoint %in% c(0, predictionMonth)),]
      df <- spread(df, Timepoint, bug)
      df <- na.omit(df)
      df$bugChange <- df[,which(colnames(df) == predictionMonth)] - df$`0`
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
                          # conf.int = TRUE, # Add confidence interval
                          # cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                          # cor.coeff.args = list(method = "kendall"),
                          add.params = list(fill = "lightgray")
        ); plot
        
        x.lab <- paste0(colnames(myT)[i], " change (%): BL to ", predictionMonth, " month(s)")
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

  plotFileName <- paste0(level, "_kendall_Taxa_Weight_change_plots_BLto", predictionMonth,"M.pdf")
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
  
  FileName <- paste0(level, "_kendall_Taxa_Weight_change_results_BLto", predictionMonth,"M.tsv")
  file.path <- paste0(outputLevel, FileName)
  write.table(final, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
} # for (level in levels)



# Analysis (grouped by surgery) ---------------------------------------------------------------


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
        
        df <- df[(df$Timepoint %in% c(0, predictionMonth)),]
        df <- spread(df, Timepoint, bug)
        df <- na.omit(df)
        df$bugChange <- df[,which(colnames(df) == predictionMonth)] - df$`0`
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
          
          x.lab <- paste0(colnames(myT)[i], " change (%): BL to ", predictionMonth," month(s)")
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
    
    plotFileName <- paste0(level, "_kendall_Taxa_Weight_change_plots_BLto", predictionMonth,"M_by_Surgery.pdf")
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
    
    FileName <- paste0(level, "_kendall_Taxa_Weight_change_results_BLto", predictionMonth,"M_by_Surgery.tsv")
    file.path <- paste0(outputLevel, FileName)
    write.table(final, file.path, sep="\t",quote = FALSE, row.names = FALSE)
    
  }
} # for (level in levels)
