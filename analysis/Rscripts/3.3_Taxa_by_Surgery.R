#Author: Alicia Sorgen
#Date: 2022 August 24
#Description: Modeling to compare taxonomic differences between surgery type at each timepoint.

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


# Libraries ----------------------------------------------------------------
H1("Libraries")
R <- sessionInfo()
message(R$R.version$version.string); rm(R)

library(stringr); message("stringr: Version ", packageVersion("stringr"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))


# Script-specific edits -------------------------------------------------
params <- vector()
params <- c(params, getOption("mbs.pipe_root", default = stop("Set options(mbs.pipe_root = ...) in your local .Rprofile (gitignored)")))
params <- c(params, "MetaPhlAn2")
params <- c(params, 0)
params <- c(params, 1)
params <- c(params, 6)
params <- c(params, 12)
params <- c(params, 18)
params <- c(params, 24)

module <- "3.3_Taxa_by_Surgery"

levels <- vector()
levels <- c(levels, "phylum")
# levels <- c(levels, "class")
# levels <- c(levels, "order")
# levels <- c(levels, "family")
levels <- c(levels, "genus")
# levels <- c(levels, "species")


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
# dir.create(pipeRoot, showWarnings = FALSE, recursive = TRUE)

# if (length(list.files(file.path(pipeRoot, "input"), recursive = TRUE)) == 0) {
#    dir.create(file.path(pipeRoot, "input"), showWarnings = FALSE, recursive = TRUE)
#    invisible(file.copy(list.files(input_root, full.names = TRUE, include.dirs = TRUE),
#                        file.path(pipeRoot, "input"), recursive = TRUE))
# }

moduleDir <- file.path(pipeRoot, module)
# dir.create(moduleDir, showWarnings = FALSE)

outputDir <- file.path(moduleDir, "output/")
# dir.create(outputDir, showWarnings = FALSE)
unlink(list.files(outputDir, full.names = TRUE, recursive = TRUE))

rm(proj_root, inputEnv, resultsEnv, module, params)


# Set functions ------------------------------------------------
H1("Functions")
message("Loading functions from: ", script_root)
funcScript <- file.path(script_root, "functions.R")
source(funcScript); rm(funcScript)


##### Set up input #####
classifier <- args[2]
included <- args[3:length(args)]

prevModule <- paste0("TaxaMetaMerge")
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

      stats <- tryCatch({
        df %>%
          group_by(Timepoint) %>%
          wilcox_test(bug ~ Surgery) %>%
          adjust_pvalue(method = "BH") %>%
          add_significance()
      }, error = function(e) NULL)

      if (is.null(stats)) next

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

