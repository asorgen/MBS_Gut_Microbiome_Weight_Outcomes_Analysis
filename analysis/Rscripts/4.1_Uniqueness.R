#Author: Alicia Sorgen
#Date: 2023 June 19
#Description: 
# getwd()

# Script set up -----------------------------------------------------------------------------------------------------------------------
rm(list=ls())
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
set.seed(1989)


# ANALYSIS <- "MetaPhlAn2_microbiome"
# params <- vector()
# params <- c(params, "~/UNCC/Projects/Bariatric_Surgery/Git_Repositories/gut-microbiome-bariatric-weight-outcomes")
# params <- c(params, "MetaPhlAn2")

moduleRoot <- "4.1_Uniqueness"
module <- moduleRoot

taxaLevels <- vector()
# taxaLevels <- c(taxaLevels, "phylum")
# taxaLevels <- c(taxaLevels, "class")
# taxaLevels <- c(taxaLevels, "order")
# taxaLevels <- c(taxaLevels, "family")
taxaLevels <- c(taxaLevels, "genus")
# taxaLevels <- c(taxaLevels, "species")

##### Libraries #####
H1("Libraries")
R <- sessionInfo()
message(R$R.version$version.string)

library(stringr); message("stringr: Version ", packageVersion("stringr"))
library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
library(scales); message("scales: Version ", packageVersion("scales"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))
library(vegan); message("vegan: Version ", packageVersion("vegan"))



##### Set up working environment #####
params <- vector()
params <- c(params, "~/UNCC/Projects/Bariatric_Surgery/Git_Repositories/MBS_Gut_Microbiome_Weight_Outcomes_Analysis/MBS_Gut_Microbiome_Weight_Outcomes_Analysis")
params <- c(params, "MetaPhlAn2")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- params
}

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
  outputRoot <- root
}

if (args[1] == "BLJ") {
  pipeRoot  <- dirname(dirname(getwd()))
  moduleDir <- dirname(getwd())
  outputRoot <- pipeRoot
}




##### Set up functions file #####
H1("Functions")
funcScript <- if (args[1] == "BLJ") file.path(moduleDir, "resources", "functions.R") else file.path(gitScripts, "functions.R")
source(funcScript)

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]
rm(str, funcScript)

##### Set up input #####
H1("Input")
classifier <- args[2]

prevModule <- paste0("Diversity_Metrics")
inputDir <- paste0(outputRoot,"/",str_subset(dir(outputRoot), prevModule),"/output")
message(inputDir)

file <- paste0("_LogNormalizedCounts_", classifier, "_alpha.tsv")
files <- list()
index <- 0
for (taxaLevel in taxaLevels) {
  
  index <- index + 1
  fileName <- paste0(taxaLevel, file)
  message(fileName)
  
  filePath <- file.path(inputDir, taxaLevel, fileName)
  
  logcounts_and_meta = read.table(filePath, header = T, sep = "\t")
  logcounts_and_meta$time  <- ifelse(logcounts_and_meta$Timepoint == "BL", 0,
                                     ifelse(logcounts_and_meta$Timepoint == "ONE", 1,
                                            ifelse(logcounts_and_meta$Timepoint == "SIX", 6,
                                                   ifelse(logcounts_and_meta$Timepoint == "TWELVE", 12,
                                                          ifelse(logcounts_and_meta$Timepoint == "EIGHTEEN", 18,
                                                                 24)))))
  files[[index]] <- logcounts_and_meta
  names(files)[index] <- taxaLevel

}

##### Compute uniqueness metrics #####
H1("Compute uniqueness metrics")
for (taxaLevel in taxaLevels) {
  df <- files[[taxaLevel]]

  taxaStart <- which(colnames(df) == "ResponderStatus") + 1
  taxaEnd   <- which(colnames(df) == "ShannonIndex") - 1
  taxaMat   <- as.matrix(df[, taxaStart:taxaEnd])
  taxaMat[is.na(taxaMat)] <- 0

  # Bray-Curtis uniqueness: mean BC distance to all other samples
  # Samples with zero total abundance get NaN distance — convert to NA
  bray_dist <- as.matrix(vegdist(taxaMat, method = "bray"))
  bray_dist[is.nan(bray_dist)] <- NA
  diag(bray_dist) <- NA
  df$Bray_Uniqueness <- rowMeans(bray_dist, na.rm = TRUE)

  # Kendall uniqueness: 1 - mean Kendall correlation to all other samples
  # Use Spearman (fast BLAS) as a close approximation to Kendall rank similarity
  kend_cor <- cor(t(taxaMat), method = "spearman")
  diag(kend_cor) <- NA
  df$Kendall_Uniqueness <- 1 - rowMeans(kend_cor, na.rm = TRUE)

  files[[taxaLevel]] <- df
}

##### Set up output #####
H1("Set up output")
outputDir = file.path(moduleDir,"output/")
message(outputDir)
for (taxaLevel in taxaLevels) {
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  dir.create(outputDirLevel, showWarnings = FALSE)
  dir.create(outputDirLevelPlots, showWarnings = FALSE)
  dir.create(outputDirLevelResults, showWarnings = FALSE)

}

##### Write updated _alpha.tsv with uniqueness columns #####
H1("Write updated alpha file with uniqueness metrics")
for (taxaLevel in taxaLevels) {
  df <- files[[taxaLevel]]
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outFile <- paste0(outputDirLevel, taxaLevel, "_LogNormalizedCounts_", classifier, "_alpha.tsv")
  write.table(df, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote: ", outFile)
}
##### Script variables #####
H1("Script variables")
surgeryPalette <- c("steelblue", "gold")
outcomePalette <- c("steelblue", "tomato")
colorPalette <- c("pink", "purple", "orange", "green", "grey", "steelblue", "tomato")

month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")
timepoints <- c(0, 1, 6, 12, 18, 24)

#####- Bray-Curtis Uniqueness ~ Timepoint -- Tukey -------------------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  x.lab <- "Timepoint"
  y.lab <- "Bray-Curtis Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT <- myT[is.finite(myT$Bray_Uniqueness), ]
  myT$time <- factor(myT$time, levels = timepoints)

  stat.test <- myT %>%
    # group_by(Surgery) %>%
    tukey_hsd(Bray_Uniqueness ~ time) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_BrayCurtis_Uniqueness_by_Timepoint_Tukey_Results.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  stat.test <- stat.test %>%
    add_xy_position(x = "time", fun = "max")
  
  plot <- ggboxplot(
    myT, x = "time", y = "Bray_Uniqueness", color = "black",
    fill = "time", palette = colorPalette,
    # facet.by = "month",
    add = "jitter",
    scales = "free"
  )
  
  # plot <- facet(plot, facet.by = "Surgery")
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.01,
      size = 5,
      hide.ns = TRUE,
      step.increase = 0.025,
      label = "p.adj.signif"
    )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_BrayCurtis_Uniqueness_by_Timepoint_Tukey_BoxPlot.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)

#####- Bray-Curtis Uniqueness ~ Timepoint (+ Surgery) -- Tukey -------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  x.lab <- "Timepoint"
  y.lab <- "Bray-Curtis Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT <- myT[!is.na(myT$Surgery),]
  myT$time <- factor(myT$time, levels = timepoints)
  
  stat.test <- myT %>%
    group_by(Surgery) %>%
    tukey_hsd(Bray_Uniqueness ~ time) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_BrayCurtis_Uniqueness_by_Timepoint_by_Surgery_Tukey_Results.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  stat.test <- stat.test %>%
    add_xy_position(x = "time", fun = "max")
  
  plot <- ggboxplot(
    myT, x = "time", y = "Bray_Uniqueness", color = "black",
    fill = "time", palette = colorPalette,
    # facet.by = "month",
    add = "jitter",
    scales = "free"
  )
  
  plot <- facet(plot, facet.by = "Surgery")
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.01,
      size = 5,
      hide.ns = TRUE,
      step.increase = 0.01,
      label = "p.adj.signif"
    )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_BrayCurtis_Uniqueness_by_Timepoint_by_Surgery_Tukey_BoxPlot.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)

#####- Kendall Uniqueness ~ Timepoint -- Tukey -------------------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  x.lab <- "Timepoint"
  y.lab <- "Kendall Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT$time <- factor(myT$time, levels = timepoints)
  
  stat.test <- myT %>%
    # group_by(Surgery) %>%
    tukey_hsd(Kendall_Uniqueness ~ time) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Kendall_Uniqueness_by_Timepoint_Tukey_Results.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  stat.test <- stat.test %>%
    add_xy_position(x = "time", fun = "max")
  
  plot <- ggboxplot(
    myT, x = "time", y = "Kendall_Uniqueness", color = "black",
    fill = "time", palette = colorPalette,
    # facet.by = "month",
    add = "jitter",
    scales = "free"
  )
  
  # plot <- facet(plot, facet.by = "Surgery")
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.01,
      size = 5,
      hide.ns = TRUE,
      step.increase = 0.025,
      label = "p.adj.signif"
    )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_Kendall_Uniqueness_by_Timepoint_Tukey_BoxPlot.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)

#####- Kendall Uniqueness ~ Timepoint (+ Surgery) -- Tukey -------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  x.lab <- "Timepoint"
  y.lab <- "Kendall Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT <- myT[!is.na(myT$Surgery),]
  myT$time <- factor(myT$time, levels = timepoints)
  
  stat.test <- myT %>%
    group_by(Surgery) %>%
    tukey_hsd(Kendall_Uniqueness ~ time) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Kendall_Uniqueness_by_Timepoint_by_Surgery_Tukey_Results.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  stat.test <- stat.test %>%
    add_xy_position(x = "time", fun = "max")
  
  plot <- ggboxplot(
    myT, x = "time", y = "Kendall_Uniqueness", color = "black",
    fill = "time", palette = colorPalette,
    # facet.by = "month",
    add = "jitter",
    scales = "free"
  )
  
  plot <- facet(plot, facet.by = "Surgery")
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.01,
      size = 5,
      hide.ns = TRUE,
      step.increase = 0.01,
      label = "p.adj.signif"
    )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_Kendall_Uniqueness_by_Timepoint_by_Surgery_Tukey_BoxPlot.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)


#####- Uniqueness Percentiles ----------------------------------------------------------------######
files2 <- list()
for (taxaLevel in taxaLevels) {
  
  
  myT <- files[[which(names(files) == taxaLevel)]]
  rownames(myT) <- myT$SampleID
  
  taxaStart <- which(colnames(myT) == "ResponderStatus") + 1
  
  metaTable <- myT[,1:(taxaStart - 1)]
  taxaTable <- myT[,taxaStart:ncol(myT)]
  
  dissimilarityMatrix_bray = generateDissimilarityMatrix(countsTable = taxaTable,dissimilarityMetric = "bray")
  dissimilarityMatrix_kendallcorr = generateCorrelationDissimilarity(countsTable = taxaTable, corrMethod = "kendall")
  
  # myT <- files[[which(names(files) == taxaLevel)]]
  myT$time <- factor(myT$time, levels = timepoints)
  
  my_list_of_percentiles = seq(0,1,by=0.1)

  bray_uniqueness_percentiles = list()
  kendall_uniqueness_percentiles = list()
  
  for(i in 1:length(my_list_of_percentiles)){
    
    bray_uniqueness_percentiles[[i]] = getDissimilarityValue_at_a_percentile(dissimilarityMatrix_bray,
                                                                             percentile = my_list_of_percentiles[i])
    
    kendall_uniqueness_percentiles[[i]] = getDissimilarityValue_at_a_percentile(dissimilarityMatrix_kendallcorr,
                                                                                percentile = my_list_of_percentiles[i])
    names(bray_uniqueness_percentiles)[i] <- paste0("Bray_Percentile_", my_list_of_percentiles[i])
    names(kendall_uniqueness_percentiles)[i] <- paste0("Kendall_Percentile_", my_list_of_percentiles[i])
  }
  
  myT1 <- myT
  for(i in length(bray_uniqueness_percentiles):1){
    myT1 <- data.frame(bray_uniqueness_percentiles[[i]], myT1)
    colnames(myT1)[1] <- names(bray_uniqueness_percentiles)[i]
  }
  
  for(i in length(kendall_uniqueness_percentiles):1){
    myT1 <- data.frame(kendall_uniqueness_percentiles[[i]], myT1)
    colnames(myT1)[1] <- names(kendall_uniqueness_percentiles)[i]
  }
  
  files[[which(taxaLevels == taxaLevel)]] <- myT1
  
  
  my_list_of_percentiles = seq(0.9,1,by=0.01)
  
  bray_uniqueness_percentiles = list()
  kendall_uniqueness_percentiles = list()
  
  for(i in 1:length(my_list_of_percentiles)){
    
    bray_uniqueness_percentiles[[i]] = getDissimilarityValue_at_a_percentile(dissimilarityMatrix_bray,
                                                                             percentile = my_list_of_percentiles[i])
    
    kendall_uniqueness_percentiles[[i]] = getDissimilarityValue_at_a_percentile(dissimilarityMatrix_kendallcorr,
                                                                                percentile = my_list_of_percentiles[i])
    names(bray_uniqueness_percentiles)[i] <- paste0("Bray_Percentile_", my_list_of_percentiles[i])
    names(kendall_uniqueness_percentiles)[i] <- paste0("Kendall_Percentile_", my_list_of_percentiles[i])
  }
  
  myT1 <- myT
  for(i in length(bray_uniqueness_percentiles):1){
    myT1 <- data.frame(bray_uniqueness_percentiles[[i]], myT1)
    colnames(myT1)[1] <- names(bray_uniqueness_percentiles)[i]
  }
  
  for(i in length(kendall_uniqueness_percentiles):1){
    myT1 <- data.frame(kendall_uniqueness_percentiles[[i]], myT1)
    colnames(myT1)[1] <- names(kendall_uniqueness_percentiles)[i]
  }
  
  files2[[which(taxaLevels == taxaLevel)]] <- myT1
  names(files2)[index] <- taxaLevel
  
} # for (taxaLevel in taxaLevels)


#####- Uniqueness Percentiles ~ Timepoint -- Tukey -------------------------------------------######
for (taxaLevel in taxaLevels) {
  
  
  myT <- files[[which(names(files) == taxaLevel)]]
  rownames(myT) <- myT$SampleID
  myT$time <- factor(myT$time, levels = timepoints)
  
  my_list_of_percentiles = seq(0,1,by=0.1)
  
  for (beta in c("Bray", "Kendall")) {
    
    betaName <- ifelse(beta == "Bray", "Bray-Curtis", "Kendall")
    plotList <- list()
    statResults <- data.frame()
    index <- 0
    
    for(i in my_list_of_percentiles){
      
      index <- index + 1
      myT2 <- myT
      colnames(myT2)[which(colnames(myT) == paste0(beta,"_Percentile_", i)) ] <- "DataCol"
      
      stat.test <- myT2 %>%
        # group_by(Surgery) %>%
        tukey_hsd(DataCol ~ time) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance()
      stat.test$Percentile <- i
      
      statResults <- rbind(statResults, stat.test)
      
      stat.test <- stat.test %>%
        add_xy_position(x = "time", fun = "max")
      
      
      x.lab <- "Timepoint"
      y.lab <- paste0(betaName, " Uniqueness (", i*100, "th percentile)")
      title.lab <- str_to_title(taxaLevel)
      subtitle.lab <- ""
      
      plot <- ggboxplot(
        myT2, x = "time", y = "DataCol", color = "black",
        fill = "time", palette = colorPalette,
        # facet.by = "month",
        # add = "jitter",
        scales = "free"
      )
      
      # plot <- facet(plot, facet.by = "Surgery")
      
      plot <- plot + labs(x=x.lab, y = y.lab)
      if (index == 1) {
        plot <- plot + labs(title = title.lab)
      }
      # plot <- plot + labs(subtitle = subtitle.lab)
      # plot <- plot + labs(caption = caption.lab)
      plot <- plot + theme(legend.position = "none")
      
      plot <- plot +
        stat_pvalue_manual(
          stat.test,
          bracket.nudge.y = 0.01,
          size = 5,
          hide.ns = TRUE,
          step.increase = 0.025,
          label = "p.adj.signif"
        )
      
      plotList[[index]] <- plot
      
      
    } # for(i in my_list_of_percentiles)
    
    file.path <- paste0(taxaLevel, "_", beta,"_Uniqueness_Percentiles_by_Timepoint_Tukey_Results.tsv")
    write.table(statResults, file.path, sep="\t",quote = FALSE, row.names = FALSE)
    
    fileName <-   paste0(taxaLevel, "_", beta,"_Uniqueness_Percentiles_by_Timepoint_Tukey_BoxPlot.pdf")
    file.path <- paste0(outputDirLevelPlots, fileName)
    pdf(file.path, width = 10, height = 7)
    
    grid.arrange(plotList[[1]], plotList[[2]], plotList[[3]],
                 plotList[[4]], plotList[[5]], plotList[[6]],
                 ncol = 3, nrow = 2)
    grid.arrange(plotList[[7]], plotList[[8]], plotList[[9]],
                 plotList[[10]], plotList[[11]],
                 ncol = 3, nrow = 2)
    
    dev.off()
    
  } # for (beta in c("Bray", "Kendall"))
  
  
  
} # for (taxaLevel in taxaLevels)


#####- Uniqueness Percentiles ~ Timepoint (+ Surgery) -- Tukey -------------------------------------------######
for (taxaLevel in taxaLevels) {
  
  
  myT <- files[[which(names(files) == taxaLevel)]]
  rownames(myT) <- myT$SampleID
  myT$time <- factor(myT$time, levels = timepoints)
  
  my_list_of_percentiles = seq(0,1,by=0.1)
  
  for (beta in c("Bray", "Kendall")) {
    
    betaName <- ifelse(beta == "Bray", "Bray-Curtis", "Kendall")
    plotList <- list()
    statResults <- data.frame()
    index <- 0
    
    for(i in my_list_of_percentiles){
      
      index <- index + 1
      myT2 <- myT
      myT2 <- myT2[!is.na(myT2$Surgery),]
      colnames(myT2)[which(colnames(myT) == paste0(beta,"_Percentile_", i)) ] <- "DataCol"
      
      stat.test <- myT2 %>%
        group_by(Surgery) %>%
        tukey_hsd(DataCol ~ time) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance()
      stat.test$Percentile <- i
      
      statResults <- rbind(statResults, stat.test)
      
      stat.test <- stat.test %>%
        add_xy_position(x = "time", fun = "max")
      
      
      x.lab <- "Timepoint"
      y.lab <- paste0(betaName, " Uniqueness (", i*100, "th percentile)")
      title.lab <- str_to_title(taxaLevel)
      subtitle.lab <- ""
      
      plot <- ggboxplot(
        myT2, x = "time", y = "DataCol", color = "black",
        fill = "time", palette = colorPalette,
        # facet.by = "month",
        # add = "jitter",
        scales = "free"
      )
      
      plot <- facet(plot, facet.by = "Surgery")
      
      plot <- plot + labs(x=x.lab, y = y.lab)
      if (index == 1) {
        plot <- plot + labs(title = title.lab)
      }
      # plot <- plot + labs(subtitle = subtitle.lab)
      # plot <- plot + labs(caption = caption.lab)
      plot <- plot + theme(legend.position = "none")
      
      plot <- plot +
        stat_pvalue_manual(
          stat.test,
          bracket.nudge.y = 0.01,
          size = 5,
          hide.ns = TRUE,
          step.increase = 0.025,
          label = "p.adj.signif"
        )
      
      plotList[[index]] <- plot
      
      
    } # for(i in my_list_of_percentiles)
    
    file.path <- paste0(taxaLevel, "_", beta,"_Uniqueness_Percentiles_by_Timepoint_by_Surgery_Tukey_Results.tsv")
    write.table(statResults, file.path, sep="\t",quote = FALSE, row.names = FALSE)
    
    fileName <-   paste0(taxaLevel, "_", beta,"_Uniqueness_Percentiles_by_Timepoint_by_Surgery_Tukey_BoxPlot.pdf")
    file.path <- paste0(outputDirLevelPlots, fileName)
    pdf(file.path, width = 15, height = 7)
    
    grid.arrange(plotList[[1]], plotList[[2]], plotList[[3]],
                 plotList[[4]], plotList[[5]], plotList[[6]],
                 ncol = 3, nrow = 2)
    grid.arrange(plotList[[7]], plotList[[8]], plotList[[9]],
                 plotList[[10]], plotList[[11]],
                 ncol = 3, nrow = 2)
    
    dev.off()
    
  } # for (beta in c("Bray", "Kendall"))
  
  
  
} # for (taxaLevel in taxaLevels)


#####- Uniqueness Percentiles ~ Timepoint -- Mixed Linear Model ------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  myT <- files[[which(names(files) == taxaLevel)]]
  rownames(myT) <- myT$SampleID
  myT$time <- factor(myT$time, levels = timepoints)
  
  my_list_of_percentiles = seq(0,1,by=0.1)
  
  pVal_Timepoint <- vector()
  pVal_Surgery <- vector()
  pVal_Weight_kg <- vector()
  Uniqueness <- vector()
  Percentile <- vector()
  index <- 0
  
  for (beta in c("Bray", "Kendall")) {
    
    betaName <- ifelse(beta == "Bray", "Bray-Curtis", "Kendall")
    
    for(i in my_list_of_percentiles){
      
      index <- index + 1
      # myT2 <- myT
      DataCol <- myT[,which(colnames(myT) == paste0(beta,"_Percentile_", i))]
      Timepoint <- myT$time
      Surgery <- myT$Surgery
      Weight_kg <- myT$Weight_kg
      PatientID <- myT$PatientID
      
      df <- data.frame(PatientID, DataCol, Timepoint, Surgery, Weight_kg)
      df$Timepoint <- factor(df$Timepoint, levels = timepoints)
      df <- na.omit(df)
      
      # mlm <- lme( DataCol ~ Timepoint + Surgery +  Weight_kg, method = "REML", random = ~1 | PatientID, data = df )
      mlm <- lme( DataCol ~ Timepoint, method = "REML", random = ~1 | PatientID, data = df )
      fit<-anova(mlm)
      
      pVal_Timepoint[index] <- fit$"p-value"[2]
      # pVal_Surgery[index] <- fit$"p-value"[3]
      # pVal_Weight_kg[index] <- fit$"p-value"[4]
      Uniqueness[index] <- betaName
      Percentile[index] <- i
      
    } # for(i in my_list_of_percentiles)
    

  } # for (beta in c("Bray", "Kendall"))
  
  # dFrame <- data.frame(Uniqueness, Percentile, pVal_Timepoint, pVal_Surgery, pVal_Weight_kg)
  dFrame <- data.frame(Uniqueness, Percentile, pVal_Timepoint)
  dFrame$Adj_pVal_Timepoint <- NA
  # dFrame$Adj_pVal_Surgery <- NA
  # dFrame$Adj_pVal_Weight_kg <- NA
  
  for (beta in unique(dFrame$Uniqueness)) {
    dFrame$Adj_pVal_Timepoint[which(dFrame$Uniqueness == beta)] <- p.adjust(dFrame$pVal_Timepoint[which(dFrame$Uniqueness == beta)],method = "BH")
    # dFrame$Adj_pVal_Surgery[which(dFrame$Uniqueness == beta)] <- p.adjust(dFrame$pVal_Surgery[which(dFrame$Uniqueness == beta)],method = "BH")
    # dFrame$Adj_pVal_Weight_kg[which(dFrame$Uniqueness == beta)] <- p.adjust(dFrame$pVal_Weight_kg[which(dFrame$Uniqueness == beta)],method = "BH")
  } # for (beta in unique(dFrame$Uniqueness))
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Uniqueness_Percentiles_by_Timepoint_MixedLinearModel_Results.tsv")
  write.table(dFrame, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
} # for (taxaLevel in taxaLevels)



#####- Uniqueness Percentiles ~ Timepoint -- Mixed Linear Model ------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  myT <- files[[which(names(files) == taxaLevel)]]
  rownames(myT) <- myT$SampleID
  myT$time <- factor(myT$time, levels = timepoints)
  
  my_list_of_percentiles = seq(0,1,by=0.1)
  
  pVal_Timepoint <- vector()
  pVal_Surgery <- vector()
  pVal_Weight_kg <- vector()
  Uniqueness <- vector()
  Percentile <- vector()
  index <- 0
  
  for (beta in c("Bray", "Kendall")) {
    
    betaName <- ifelse(beta == "Bray", "Bray-Curtis", "Kendall")
    
    for(i in my_list_of_percentiles){
      
      index <- index + 1
      # myT2 <- myT
      DataCol <- myT[,which(colnames(myT) == paste0(beta,"_Percentile_", i))]
      Timepoint <- myT$time
      Surgery <- myT$Surgery
      Weight_kg <- myT$Weight_kg
      PatientID <- myT$PatientID
      
      df <- data.frame(PatientID, DataCol, Timepoint, Surgery, Weight_kg)
      df$Timepoint <- factor(df$Timepoint, levels = timepoints)
      df <- na.omit(df)
      
      # mlm <- lme( DataCol ~ Timepoint + Surgery +  Weight_kg, method = "REML", random = ~1 | PatientID, data = df )
      mlm <- lme( DataCol ~ Timepoint, method = "REML", random = ~1 | PatientID, data = df )
      fit<-anova(mlm)
      
      pVal_Timepoint[index] <- fit$"p-value"[2]
      # pVal_Surgery[index] <- fit$"p-value"[3]
      # pVal_Weight_kg[index] <- fit$"p-value"[4]
      Uniqueness[index] <- betaName
      Percentile[index] <- i
      
    } # for(i in my_list_of_percentiles)
    
    
  } # for (beta in c("Bray", "Kendall"))
  
  # dFrame <- data.frame(Uniqueness, Percentile, pVal_Timepoint, pVal_Surgery, pVal_Weight_kg)
  dFrame <- data.frame(Uniqueness, Percentile, pVal_Timepoint)
  dFrame$Adj_pVal_Timepoint <- NA
  # dFrame$Adj_pVal_Surgery <- NA
  # dFrame$Adj_pVal_Weight_kg <- NA
  
  for (beta in unique(dFrame$Uniqueness)) {
    dFrame$Adj_pVal_Timepoint[which(dFrame$Uniqueness == beta)] <- p.adjust(dFrame$pVal_Timepoint[which(dFrame$Uniqueness == beta)],method = "BH")
    # dFrame$Adj_pVal_Surgery[which(dFrame$Uniqueness == beta)] <- p.adjust(dFrame$pVal_Surgery[which(dFrame$Uniqueness == beta)],method = "BH")
    # dFrame$Adj_pVal_Weight_kg[which(dFrame$Uniqueness == beta)] <- p.adjust(dFrame$pVal_Weight_kg[which(dFrame$Uniqueness == beta)],method = "BH")
  } # for (beta in unique(dFrame$Uniqueness))
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Uniqueness_Percentiles_by_Timepoint_MixedLinearModel_Results.tsv")
  write.table(dFrame, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
} # for (taxaLevel in taxaLevels)



#####- Uniqueness Percentiles ~ Timepoint (+ Surgery) -- Mixed Linear Model ------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  myT <- files[[which(names(files) == taxaLevel)]]
  rownames(myT) <- myT$SampleID
  myT$time <- factor(myT$time, levels = timepoints)
  
  my_list_of_percentiles = seq(0,1,by=0.1)
  
  pVal_Timepoint <- vector()
  SurgeryType <- vector()
  Uniqueness <- vector()
  Percentile <- vector()
  index <- 0
  
  for (beta in c("Bray", "Kendall")) {
    
    betaName <- ifelse(beta == "Bray", "Bray-Curtis", "Kendall")
    
    for(i in my_list_of_percentiles){
      
      # myT2 <- myT
      DataCol <- myT[,which(colnames(myT) == paste0(beta,"_Percentile_", i))]
      Timepoint <- myT$time
      Surgery <- myT$Surgery
      Weight_kg <- myT$Weight_kg
      PatientID <- myT$PatientID
      
      df <- data.frame(PatientID, DataCol, Timepoint, Surgery, Weight_kg)
      df$Timepoint <- factor(df$Timepoint, levels = timepoints)
      df <- na.omit(df)
      
      for (surg in unique(df$Surgery)) {
        
        index <- index + 1
        df2 <- df[df$Surgery == surg,]
        mlm <- lme( DataCol ~ Timepoint, method = "REML", random = ~1 | PatientID, data = df2 )
        fit<-anova(mlm)
        
        pVal_Timepoint[index] <- fit$"p-value"[2]
        Uniqueness[index] <- betaName
        Percentile[index] <- i
        SurgeryType[index] <- surg
        
      } # for (surg in unique(df$Surgery))
      
    } # for(i in my_list_of_percentiles)
    
    
  } # for (beta in c("Bray", "Kendall"))
  
  # dFrame <- data.frame(Uniqueness, Percentile, pVal_Timepoint, pVal_Surgery, pVal_Weight_kg)
  dFrame <- data.frame(Uniqueness, Percentile, SurgeryType, pVal_Timepoint)
  dFrame$Adj_pVal_Timepoint <- NA

  for (beta in unique(dFrame$Uniqueness)) {
    for (surg in unique(dFrame$SurgeryType)) {
      dFrame$Adj_pVal_Timepoint[which(dFrame$Uniqueness == beta & dFrame$SurgeryType == surg)] <- p.adjust(dFrame$pVal_Timepoint[which(dFrame$Uniqueness == beta & dFrame$SurgeryType == surg)],method = "BH")
      
    } # for (surg in unique(dFrame$Surgery))
  } # for (beta in unique(dFrame$Uniqueness))
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Uniqueness_Percentiles_by_Timepoint_by_Surgery_MixedLinearModel_Results.tsv")
  write.table(dFrame, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
} # for (taxaLevel in taxaLevels)



#####- Uniqueness Percentiles ~ Timepoint -- Mixed Linear Model p-value plots ----------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  myT <- files[[which(names(files) == taxaLevel)]]
  rownames(myT) <- myT$SampleID
  myT$time <- factor(myT$time, levels = timepoints)
  
  my_list_of_percentiles = seq(0,1,by=0.1)
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Uniqueness_Percentiles_by_Timepoint_MixedLinearModel_Results.tsv")
  mlmResults = read.table(file.path, header = T, sep = "\t")
  
  for (beta in c("Bray-Curtis", "Kendall")) {
    
    df <- mlmResults[mlmResults$Uniqueness == beta,]
    # df$Adj_pVal_Timepoint[which(df$Adj_pVal_Timepoint == 0)] <- 2e-16
    df$log10p <- log10(df$Adj_pVal_Timepoint) * -1
    # df$log10p[is.infinite(df$log10p)] <- 0    
    
    x.lab <- "Percentiles"
    y.lab <- "-log10 p-values"
    title.lab <- paste(str_to_title(taxaLevel), beta, "Uniqueness")
    subtitle.lab <- "Mixed Linear Model: Uniqueness Percentile ~ Timepoint (all)"
    
    plot <- ggplot(df, aes(Percentile, log10p))+
      geom_point(shape = 21, size = 2.5)
    
    plot <- plot + labs(x=x.lab, y = y.lab)
    plot <- plot + labs(title = title.lab)
    plot <- plot + labs(subtitle = subtitle.lab)
    # plot <- plot + labs(caption = caption.lab)
    # plot <- plot + theme(legend.position = "none")
    
    plot
    fileName <- paste0(taxaLevel, "_", beta, "_Uniqueness_Percentiles_by_Timepoint_MixedLinearModel_Plot.pdf")
    file.path <- paste0(outputDirLevelPlots, fileName)
    pdf(file.path, width = 7, height = 5)
    print(plot)
    dev.off()
    
  } # for (beta in c("Bray", "Kendall"))
  

} # for (taxaLevel in taxaLevels)



#####- Uniqueness Percentiles ~ Timepoint (+ Surgery) -- Mixed Linear Model p-value plots ----######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  myT <- files[[which(names(files) == taxaLevel)]]
  rownames(myT) <- myT$SampleID
  myT$time <- factor(myT$time, levels = timepoints)
  
  my_list_of_percentiles = seq(0,1,by=0.1)
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Uniqueness_Percentiles_by_Timepoint_by_Surgery_MixedLinearModel_Results.tsv")
  mlmResults = read.table(file.path, header = T, sep = "\t")
  
  for (beta in c("Bray-Curtis", "Kendall")) {
    
    df <- mlmResults[mlmResults$Uniqueness == beta,]
    # df$Adj_pVal_Timepoint[which(df$Adj_pVal_Timepoint == 0)] <- 2e-16
    df$log10p <- log10(df$Adj_pVal_Timepoint) * -1
    # df$log10p[is.infinite(df$log10p)] <- 0    
    
    x.lab <- "Percentiles"
    y.lab <- "-log10 p-values"
    title.lab <- paste(str_to_title(taxaLevel), beta, "Uniqueness")
    subtitle.lab <- "Mixed Linear Model: Uniqueness Percentile ~ Timepoint (all)"
    
    plot <- ggplot(df, aes(Percentile, log10p))+
      geom_point(shape = 21, size = 2.5)
    
    plot <- facet(plot, facet.by = "SurgeryType")
    plot <- plot + geom_hline(yintercept=-log10(0.05), linetype = "dotted", col = "red")
    
    plot <- plot + labs(x=x.lab, y = y.lab)
    plot <- plot + labs(title = title.lab)
    plot <- plot + labs(subtitle = subtitle.lab)
    # plot <- plot + labs(caption = caption.lab)
    # plot <- plot + theme(legend.position = "none")
    
    plot
    fileName <- paste0(taxaLevel, "_", beta, "_Uniqueness_Percentiles_by_Timepoint_by_Surgery_MixedLinearModel_Plot.pdf")
    file.path <- paste0(outputDirLevelPlots, fileName)
    pdf(file.path, width = 8, height = 5)
    print(plot)
    dev.off()
    
  } # for (beta in c("Bray", "Kendall"))
  
  
} # for (taxaLevel in taxaLevels)



#####- Uniqueness Percentiles (90-100) ~ Timepoint -- Mixed Linear Model ---------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  myT <- files2[[which(names(files2) == taxaLevel)]]
  rownames(myT) <- myT$SampleID
  myT$time <- factor(myT$time, levels = timepoints)
  
  my_list_of_percentiles = seq(0.9,1,by=0.01)
  
  pVal_Timepoint <- vector()
  pVal_Surgery <- vector()
  pVal_Weight_kg <- vector()
  Uniqueness <- vector()
  Percentile <- vector()
  index <- 0
  
  for (beta in c("Bray", "Kendall")) {
    
    betaName <- ifelse(beta == "Bray", "Bray-Curtis", "Kendall")
    
    for(i in my_list_of_percentiles){
      
      index <- index + 1
      # myT2 <- myT
      DataCol <- myT[,which(colnames(myT) == paste0(beta,"_Percentile_", i))]
      Timepoint <- myT$time
      Surgery <- myT$Surgery
      Weight_kg <- myT$Weight_kg
      PatientID <- myT$PatientID
      
      df <- data.frame(PatientID, DataCol, Timepoint, Surgery, Weight_kg)
      df$Timepoint <- factor(df$Timepoint, levels = timepoints)
      df <- na.omit(df)
      
      # mlm <- lme( DataCol ~ Timepoint + Surgery +  Weight_kg, method = "REML", random = ~1 | PatientID, data = df )
      mlm <- lme( DataCol ~ Timepoint, method = "REML", random = ~1 | PatientID, data = df )
      fit<-anova(mlm)
      
      pVal_Timepoint[index] <- fit$"p-value"[2]
      # pVal_Surgery[index] <- fit$"p-value"[3]
      # pVal_Weight_kg[index] <- fit$"p-value"[4]
      Uniqueness[index] <- betaName
      Percentile[index] <- i
      
    } # for(i in my_list_of_percentiles)
    
    
  } # for (beta in c("Bray", "Kendall"))
  
  # dFrame <- data.frame(Uniqueness, Percentile, pVal_Timepoint, pVal_Surgery, pVal_Weight_kg)
  dFrame <- data.frame(Uniqueness, Percentile, pVal_Timepoint)
  dFrame$Adj_pVal_Timepoint <- NA
  # dFrame$Adj_pVal_Surgery <- NA
  # dFrame$Adj_pVal_Weight_kg <- NA
  
  for (beta in unique(dFrame$Uniqueness)) {
    dFrame$Adj_pVal_Timepoint[which(dFrame$Uniqueness == beta)] <- p.adjust(dFrame$pVal_Timepoint[which(dFrame$Uniqueness == beta)],method = "BH")
    # dFrame$Adj_pVal_Surgery[which(dFrame$Uniqueness == beta)] <- p.adjust(dFrame$pVal_Surgery[which(dFrame$Uniqueness == beta)],method = "BH")
    # dFrame$Adj_pVal_Weight_kg[which(dFrame$Uniqueness == beta)] <- p.adjust(dFrame$pVal_Weight_kg[which(dFrame$Uniqueness == beta)],method = "BH")
  } # for (beta in unique(dFrame$Uniqueness))
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Uniqueness_Percentiles_90_100_by_Timepoint_MixedLinearModel_Results.tsv")
  write.table(dFrame, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
} # for (taxaLevel in taxaLevels)



#####- Uniqueness Percentiles (90-100) ~ Timepoint -- Mixed Linear Model p-value plots -------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  myT <- files2[[which(names(files2) == taxaLevel)]]
  rownames(myT) <- myT$SampleID
  myT$time <- factor(myT$time, levels = timepoints)
  
  my_list_of_percentiles = seq(0.9,1,by=0.01)
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Uniqueness_Percentiles_90_100_by_Timepoint_MixedLinearModel_Results.tsv")
  mlmResults = read.table(file.path, header = T, sep = "\t")
  
  for (beta in c("Bray-Curtis", "Kendall")) {
    
    df <- mlmResults[mlmResults$Uniqueness == beta,]
    # df$Adj_pVal_Timepoint[which(df$Adj_pVal_Timepoint == 0)] <- 2e-16
    df$log10p <- log10(df$Adj_pVal_Timepoint) * -1
    # df$log10p[is.infinite(df$log10p)] <- 0    
    
    x.lab <- "Percentiles"
    y.lab <- "-log10 p-values"
    title.lab <- paste(str_to_title(taxaLevel), beta, "Uniqueness")
    subtitle.lab <- "Mixed Linear Model: Uniqueness Percentile ~ Timepoint (all)"
    
    plot <- ggplot(df, aes(Percentile, log10p))+
      geom_point(shape = 21, size = 2.5)
    
    plot <- plot + labs(x=x.lab, y = y.lab)
    plot <- plot + labs(title = title.lab)
    plot <- plot + labs(subtitle = subtitle.lab)
    # plot <- plot + labs(caption = caption.lab)
    # plot <- plot + theme(legend.position = "none")
    
    plot
    fileName <- paste0(taxaLevel, "_", beta, "_Uniqueness_Percentiles_90_100_by_Timepoint_MixedLinearModel_Plot.pdf")
    file.path <- paste0(outputDirLevelPlots, fileName)
    pdf(file.path, width = 7, height = 5)
    print(plot)
    dev.off()
    
  } # for (beta in c("Bray", "Kendall"))
  
  
} # for (taxaLevel in taxaLevels)




#####- Bray-Curtis Uniqueness ~ Weight_kg -- Kendall -----------------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  x.lab <- "Weight (kg)"
  y.lab <- "Bray-Curtis Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT$time <- factor(myT$time, levels = timepoints)
  
  stat.test <- myT %>%
    # group_by(Surgery) %>%
    cor_test(Bray_Uniqueness, Weight_kg, method = "kendall") %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_BrayCurtis_Uniqueness_by_Weight_kg_KendallResults.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  # stat.test <- stat.test %>%
  #   add_xy_position(x = "Weight_kg", fun = "max")
  
  # plot <- ggboxplot(
  #   myT, x = "Weight_kg", y = "Bray_Uniqueness", color = "black",
  #   fill = "Weight_kg", palette = colorPalette,
  #   # facet.by = "month",
  #   add = "jitter",
  #   scales = "free"
  # )
  
  plot <- ggscatter(myT, x = "Weight_kg", y = "Bray_Uniqueness",
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    add.params = list(color = "red"), # Customize reg. line
                    # conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    # cor.coef.coord = c(min, max),
                    cor.coeff.args = list(method = "kendall")
  )
  
  # plot <- facet(plot, facet.by = "Surgery")
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  # plot <- plot +
  #   stat_pvalue_manual(
  #     stat.test,
  #     bracket.nudge.y = 0.01,
  #     size = 5,
  #     hide.ns = TRUE,
  #     step.increase = 0.025,
  #     label = "p.adj.signif"
  #   )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_BrayCurtis_Uniqueness_by_Weight_kg_Kendall_Scatter.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)


#####- Bray-Curtis Uniqueness ~ Weight_kg (+ Timepoint) -- Kendall -----------------------------######
for (taxaLevel in taxaLevels) {
  
  metricType <- "raw"
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  i <- 1
  
  if (metricType == "loss") {
    i <- 2
  }
  
  x.lab <- "Weight (kg)"
  y.lab <- "Bray-Curtis Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT <- myT[myT$time %in% timepoints[i:length(timepoints)],]
  myT$time <- factor(myT$time, levels = timepoints[i:length(timepoints)])
  
  stat.test <- myT %>%
    group_by(time) %>%
    cor_test(Bray_Uniqueness, Weight_kg, method = "kendall") %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_BrayCurtis_Uniqueness_by_Weight_kg__by_Timepoint_KendallResults.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  # stat.test <- stat.test %>%
  #   add_xy_position(x = "Weight_kg", fun = "max")
  
  # plot <- ggboxplot(
  #   myT, x = "Weight_kg", y = "Bray_Uniqueness", color = "black",
  #   fill = "Weight_kg", palette = colorPalette,
  #   # facet.by = "month",
  #   add = "jitter",
  #   scales = "free"
  # )
  
  plot <- ggscatter(myT, x = "Weight_kg", y = "Bray_Uniqueness",
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    add.params = list(color = "red"), # Customize reg. line
                    # conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    # cor.coef.coord = c(min, max),
                    cor.coeff.args = list(method = "kendall")
  )
  
  facet.Labels <- month.labs[i:length(timepoints)]
  # times <- timepoints[i:length(timepoints)]
  plot <- facet(plot, facet.by = "time"
                , panel.labs = list(time = facet.Labels)
  ); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  # plot <- plot +
  #   stat_pvalue_manual(
  #     stat.test,
  #     bracket.nudge.y = 0.01,
  #     size = 5,
  #     hide.ns = TRUE,
  #     step.increase = 0.025,
  #     label = "p.adj.signif"
  #   )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_BrayCurtis_Uniqueness_by_Weight_kg_by_Timepoint_Kendall_Scatter.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)

#####- Bray-Curtis Uniqueness ~ Percent_Loss_kg -- Kendall -----------------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  x.lab <- "Weight loss from BL (%)"
  y.lab <- "Bray-Curtis Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT$time <- factor(myT$time, levels = timepoints)
  
  stat.test <- myT %>%
    # group_by(Surgery) %>%
    cor_test(Bray_Uniqueness, Percent_Loss_kg, method = "kendall") %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_BrayCurtis_Uniqueness_by_Percent_Loss_kg_KendallResults.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  # stat.test <- stat.test %>%
  #   add_xy_position(x = "Percent_Loss_kg", fun = "max")
  
  # plot <- ggboxplot(
  #   myT, x = "Percent_Loss_kg", y = "Bray_Uniqueness", color = "black",
  #   fill = "Percent_Loss_kg", palette = colorPalette,
  #   # facet.by = "month",
  #   add = "jitter",
  #   scales = "free"
  # )
  
  plot <- ggscatter(myT, x = "Percent_Loss_kg", y = "Bray_Uniqueness",
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    add.params = list(color = "red"), # Customize reg. line
                    # conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    # cor.coef.coord = c(min, max),
                    cor.coeff.args = list(method = "kendall")
  )
  
  # plot <- facet(plot, facet.by = "Surgery")
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  # plot <- plot +
  #   stat_pvalue_manual(
  #     stat.test,
  #     bracket.nudge.y = 0.01,
  #     size = 5,
  #     hide.ns = TRUE,
  #     step.increase = 0.025,
  #     label = "p.adj.signif"
  #   )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_BrayCurtis_Uniqueness_by_Percent_Loss_kg_Kendall_Scatter.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)


#####- Bray-Curtis Uniqueness ~ Percent_Loss_kg (+ Timepoint) -- Kendall -----------------------------######
for (taxaLevel in taxaLevels) {
  
  metricType <- "loss"
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  i <- 1
  
  if (metricType == "loss") {
    i <- 2
  }
  
  x.lab <- "Weight loss from BL (%)"
  y.lab <- "Bray-Curtis Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT <- myT[myT$time %in% timepoints[i:length(timepoints)],]
  myT$time <- factor(myT$time, levels = timepoints[i:length(timepoints)])
  
  stat.test <- myT %>%
    group_by(time) %>%
    cor_test(Bray_Uniqueness, Percent_Loss_kg, method = "kendall") %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_BrayCurtis_Uniqueness_by_Percent_Loss_kg__by_Timepoint_KendallResults.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  # stat.test <- stat.test %>%
  #   add_xy_position(x = "Percent_Loss_kg", fun = "max")
  
  # plot <- ggboxplot(
  #   myT, x = "Percent_Loss_kg", y = "Bray_Uniqueness", color = "black",
  #   fill = "Percent_Loss_kg", palette = colorPalette,
  #   # facet.by = "month",
  #   add = "jitter",
  #   scales = "free"
  # )
  
  plot <- ggscatter(myT, x = "Percent_Loss_kg", y = "Bray_Uniqueness",
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    add.params = list(color = "red"), # Customize reg. line
                    # conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    # cor.coef.coord = c(min, max),
                    cor.coeff.args = list(method = "kendall")
  )
  
  facet.Labels <- month.labs[i:length(timepoints)]
  # times <- timepoints[i:length(timepoints)]
  plot <- facet(plot, facet.by = "time"
                , panel.labs = list(time = facet.Labels)
  ); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  # plot <- plot +
  #   stat_pvalue_manual(
  #     stat.test,
  #     bracket.nudge.y = 0.01,
  #     size = 5,
  #     hide.ns = TRUE,
  #     step.increase = 0.025,
  #     label = "p.adj.signif"
  #   )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_BrayCurtis_Uniqueness_by_Percent_Loss_kg_by_Timepoint_Kendall_Scatter.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)

#####- Bray-Curtis Uniqueness ~ PEWL_kg -- Kendall -----------------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  x.lab <- "Excess weight loss from BL (%)"
  y.lab <- "Bray-Curtis Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT$time <- factor(myT$time, levels = timepoints)
  
  stat.test <- myT %>%
    # group_by(Surgery) %>%
    cor_test(Bray_Uniqueness, PEWL_kg, method = "kendall") %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_BrayCurtis_Uniqueness_by_PEWL_kg_KendallResults.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  # stat.test <- stat.test %>%
  #   add_xy_position(x = "PEWL_kg", fun = "max")
  
  # plot <- ggboxplot(
  #   myT, x = "PEWL_kg", y = "Bray_Uniqueness", color = "black",
  #   fill = "PEWL_kg", palette = colorPalette,
  #   # facet.by = "month",
  #   add = "jitter",
  #   scales = "free"
  # )
  
  plot <- ggscatter(myT, x = "PEWL_kg", y = "Bray_Uniqueness",
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    add.params = list(color = "red"), # Customize reg. line
                    # conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    # cor.coef.coord = c(min, max),
                    cor.coeff.args = list(method = "kendall")
  )
  
  # plot <- facet(plot, facet.by = "Surgery")
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  # plot <- plot +
  #   stat_pvalue_manual(
  #     stat.test,
  #     bracket.nudge.y = 0.01,
  #     size = 5,
  #     hide.ns = TRUE,
  #     step.increase = 0.025,
  #     label = "p.adj.signif"
  #   )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_BrayCurtis_Uniqueness_by_PEWL_kg_Kendall_Scatter.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)


#####- Bray-Curtis Uniqueness ~ PEWL_kg (+ Timepoint) -- Kendall -----------------------------######
for (taxaLevel in taxaLevels) {
  
  metricType <- "loss"
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  i <- 1
  
  if (metricType == "loss") {
    i <- 2
  }
  
  x.lab <- "Excess weight loss from BL (%)"
  y.lab <- "Bray-Curtis Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT <- myT[myT$time %in% timepoints[i:length(timepoints)],]
  myT$time <- factor(myT$time, levels = timepoints[i:length(timepoints)])
  
  stat.test <- myT %>%
    group_by(time) %>%
    cor_test(Bray_Uniqueness, PEWL_kg, method = "kendall") %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_BrayCurtis_Uniqueness_by_PEWL_kg__by_Timepoint_KendallResults.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  # stat.test <- stat.test %>%
  #   add_xy_position(x = "PEWL_kg", fun = "max")
  
  # plot <- ggboxplot(
  #   myT, x = "PEWL_kg", y = "Bray_Uniqueness", color = "black",
  #   fill = "PEWL_kg", palette = colorPalette,
  #   # facet.by = "month",
  #   add = "jitter",
  #   scales = "free"
  # )
  
  plot <- ggscatter(myT, x = "PEWL_kg", y = "Bray_Uniqueness",
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    add.params = list(color = "red"), # Customize reg. line
                    # conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    # cor.coef.coord = c(min, max),
                    cor.coeff.args = list(method = "kendall")
  )
  
  facet.Labels <- month.labs[i:length(timepoints)]
  # times <- timepoints[i:length(timepoints)]
  plot <- facet(plot, facet.by = "time"
                , panel.labs = list(time = facet.Labels)
  ); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  # plot <- plot +
  #   stat_pvalue_manual(
  #     stat.test,
  #     bracket.nudge.y = 0.01,
  #     size = 5,
  #     hide.ns = TRUE,
  #     step.increase = 0.025,
  #     label = "p.adj.signif"
  #   )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_BrayCurtis_Uniqueness_by_PEWL_kg_by_Timepoint_Kendall_Scatter.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)


#####- Kendall Uniqueness ~ Weight_kg -- Kendall -----------------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  x.lab <- "Weight (kg)"
  y.lab <- "Kendall Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT$time <- factor(myT$time, levels = timepoints)
  
  stat.test <- myT %>%
    # group_by(Surgery) %>%
    cor_test(Kendall_Uniqueness, Weight_kg, method = "kendall") %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Kendall_Uniqueness_by_Weight_kg_KendallResults.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  # stat.test <- stat.test %>%
  #   add_xy_position(x = "Weight_kg", fun = "max")
  
  # plot <- ggboxplot(
  #   myT, x = "Weight_kg", y = "Kendall_Uniqueness", color = "black",
  #   fill = "Weight_kg", palette = colorPalette,
  #   # facet.by = "month",
  #   add = "jitter",
  #   scales = "free"
  # )
  
  plot <- ggscatter(myT, x = "Weight_kg", y = "Kendall_Uniqueness",
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    add.params = list(color = "red"), # Customize reg. line
                    # conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    # cor.coef.coord = c(min, max),
                    cor.coeff.args = list(method = "kendall")
  )
  
  # plot <- facet(plot, facet.by = "Surgery")
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  # plot <- plot +
  #   stat_pvalue_manual(
  #     stat.test,
  #     bracket.nudge.y = 0.01,
  #     size = 5,
  #     hide.ns = TRUE,
  #     step.increase = 0.025,
  #     label = "p.adj.signif"
  #   )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_Kendall_Uniqueness_by_Weight_kg_Kendall_Scatter.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)


#####- Kendall Uniqueness ~ Weight_kg (+ Timepoint) -- Kendall -----------------------------######
for (taxaLevel in taxaLevels) {
  
  metricType <- "raw"
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  i <- 1
  
  if (metricType == "loss") {
    i <- 2
  }
  
  x.lab <- "Weight (kg)"
  y.lab <- "Kendall Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT <- myT[myT$time %in% timepoints[i:length(timepoints)],]
  myT$time <- factor(myT$time, levels = timepoints[i:length(timepoints)])
  
  stat.test <- myT %>%
    group_by(time) %>%
    cor_test(Kendall_Uniqueness, Weight_kg, method = "kendall") %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Kendall_Uniqueness_by_Weight_kg__by_Timepoint_KendallResults.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  # stat.test <- stat.test %>%
  #   add_xy_position(x = "Weight_kg", fun = "max")
  
  # plot <- ggboxplot(
  #   myT, x = "Weight_kg", y = "Kendall_Uniqueness", color = "black",
  #   fill = "Weight_kg", palette = colorPalette,
  #   # facet.by = "month",
  #   add = "jitter",
  #   scales = "free"
  # )
  
  plot <- ggscatter(myT, x = "Weight_kg", y = "Kendall_Uniqueness",
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    add.params = list(color = "red"), # Customize reg. line
                    # conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    # cor.coef.coord = c(min, max),
                    cor.coeff.args = list(method = "kendall")
  )
  
  facet.Labels <- month.labs[i:length(timepoints)]
  # times <- timepoints[i:length(timepoints)]
  plot <- facet(plot, facet.by = "time"
                , panel.labs = list(time = facet.Labels)
  ); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  # plot <- plot +
  #   stat_pvalue_manual(
  #     stat.test,
  #     bracket.nudge.y = 0.01,
  #     size = 5,
  #     hide.ns = TRUE,
  #     step.increase = 0.025,
  #     label = "p.adj.signif"
  #   )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_Kendall_Uniqueness_by_Weight_kg_by_Timepoint_Kendall_Scatter.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)

#####- Kendall Uniqueness ~ Percent_Loss_kg -- Kendall -----------------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  x.lab <- "Weight loss from BL (%)"
  y.lab <- "Kendall Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT$time <- factor(myT$time, levels = timepoints)
  
  stat.test <- myT %>%
    # group_by(Surgery) %>%
    cor_test(Kendall_Uniqueness, Percent_Loss_kg, method = "kendall") %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Kendall_Uniqueness_by_Percent_Loss_kg_KendallResults.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  # stat.test <- stat.test %>%
  #   add_xy_position(x = "Percent_Loss_kg", fun = "max")
  
  # plot <- ggboxplot(
  #   myT, x = "Percent_Loss_kg", y = "Kendall_Uniqueness", color = "black",
  #   fill = "Percent_Loss_kg", palette = colorPalette,
  #   # facet.by = "month",
  #   add = "jitter",
  #   scales = "free"
  # )
  
  plot <- ggscatter(myT, x = "Percent_Loss_kg", y = "Kendall_Uniqueness",
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    add.params = list(color = "red"), # Customize reg. line
                    # conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    # cor.coef.coord = c(min, max),
                    cor.coeff.args = list(method = "kendall")
  )
  
  # plot <- facet(plot, facet.by = "Surgery")
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  # plot <- plot +
  #   stat_pvalue_manual(
  #     stat.test,
  #     bracket.nudge.y = 0.01,
  #     size = 5,
  #     hide.ns = TRUE,
  #     step.increase = 0.025,
  #     label = "p.adj.signif"
  #   )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_Kendall_Uniqueness_by_Percent_Loss_kg_Kendall_Scatter.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)


#####- Kendall Uniqueness ~ Percent_Loss_kg (+ Timepoint) -- Kendall -----------------------------######
for (taxaLevel in taxaLevels) {
  
  metricType <- "loss"
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  i <- 1
  
  if (metricType == "loss") {
    i <- 2
  }
  
  x.lab <- "Weight loss from BL (%)"
  y.lab <- "Kendall Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT <- myT[myT$time %in% timepoints[i:length(timepoints)],]
  myT$time <- factor(myT$time, levels = timepoints[i:length(timepoints)])
  
  stat.test <- myT %>%
    group_by(time) %>%
    cor_test(Kendall_Uniqueness, Percent_Loss_kg, method = "kendall") %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Kendall_Uniqueness_by_Percent_Loss_kg__by_Timepoint_KendallResults.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  # stat.test <- stat.test %>%
  #   add_xy_position(x = "Percent_Loss_kg", fun = "max")
  
  # plot <- ggboxplot(
  #   myT, x = "Percent_Loss_kg", y = "Kendall_Uniqueness", color = "black",
  #   fill = "Percent_Loss_kg", palette = colorPalette,
  #   # facet.by = "month",
  #   add = "jitter",
  #   scales = "free"
  # )
  
  plot <- ggscatter(myT, x = "Percent_Loss_kg", y = "Kendall_Uniqueness",
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    add.params = list(color = "red"), # Customize reg. line
                    # conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    # cor.coef.coord = c(min, max),
                    cor.coeff.args = list(method = "kendall")
  )
  
  facet.Labels <- month.labs[i:length(timepoints)]
  # times <- timepoints[i:length(timepoints)]
  plot <- facet(plot, facet.by = "time"
                , panel.labs = list(time = facet.Labels)
  ); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  # plot <- plot +
  #   stat_pvalue_manual(
  #     stat.test,
  #     bracket.nudge.y = 0.01,
  #     size = 5,
  #     hide.ns = TRUE,
  #     step.increase = 0.025,
  #     label = "p.adj.signif"
  #   )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_Kendall_Uniqueness_by_Percent_Loss_kg_by_Timepoint_Kendall_Scatter.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)

#####- Kendall Uniqueness ~ PEWL_kg -- Kendall -----------------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  x.lab <- "Excess weight loss from BL (%)"
  y.lab <- "Kendall Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT$time <- factor(myT$time, levels = timepoints)
  
  stat.test <- myT %>%
    # group_by(Surgery) %>%
    cor_test(Kendall_Uniqueness, PEWL_kg, method = "kendall") %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Kendall_Uniqueness_by_PEWL_kg_KendallResults.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  # stat.test <- stat.test %>%
  #   add_xy_position(x = "PEWL_kg", fun = "max")
  
  # plot <- ggboxplot(
  #   myT, x = "PEWL_kg", y = "Kendall_Uniqueness", color = "black",
  #   fill = "PEWL_kg", palette = colorPalette,
  #   # facet.by = "month",
  #   add = "jitter",
  #   scales = "free"
  # )
  
  plot <- ggscatter(myT, x = "PEWL_kg", y = "Kendall_Uniqueness",
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    add.params = list(color = "red"), # Customize reg. line
                    # conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    # cor.coef.coord = c(min, max),
                    cor.coeff.args = list(method = "kendall")
  )
  
  # plot <- facet(plot, facet.by = "Surgery")
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  # plot <- plot +
  #   stat_pvalue_manual(
  #     stat.test,
  #     bracket.nudge.y = 0.01,
  #     size = 5,
  #     hide.ns = TRUE,
  #     step.increase = 0.025,
  #     label = "p.adj.signif"
  #   )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_Kendall_Uniqueness_by_PEWL_kg_Kendall_Scatter.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)


#####- Kendall Uniqueness ~ PEWL_kg (+ Timepoint) -- Kendall -----------------------------######
for (taxaLevel in taxaLevels) {
  
  metricType <- "loss"
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  i <- 1
  
  if (metricType == "loss") {
    i <- 2
  }
  
  x.lab <- "Excess weight loss from BL (%)"
  y.lab <- "Kendall Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  subtitle.lab <- ""
  
  myT <- files[[which(names(files) == taxaLevel)]]
  myT <- myT[myT$time %in% timepoints[i:length(timepoints)],]
  myT$time <- factor(myT$time, levels = timepoints[i:length(timepoints)])
  
  stat.test <- myT %>%
    group_by(time) %>%
    cor_test(Kendall_Uniqueness, PEWL_kg, method = "kendall") %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevelResults, taxaLevel, "_Kendall_Uniqueness_by_PEWL_kg__by_Timepoint_KendallResults.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  # stat.test <- stat.test %>%
  #   add_xy_position(x = "PEWL_kg", fun = "max")
  
  # plot <- ggboxplot(
  #   myT, x = "PEWL_kg", y = "Kendall_Uniqueness", color = "black",
  #   fill = "PEWL_kg", palette = colorPalette,
  #   # facet.by = "month",
  #   add = "jitter",
  #   scales = "free"
  # )
  
  plot <- ggscatter(myT, x = "PEWL_kg", y = "Kendall_Uniqueness",
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    add.params = list(color = "red"), # Customize reg. line
                    # conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    # cor.coef.coord = c(min, max),
                    cor.coeff.args = list(method = "kendall")
  )
  
  facet.Labels <- month.labs[i:length(timepoints)]
  # times <- timepoints[i:length(timepoints)]
  plot <- facet(plot, facet.by = "time"
                , panel.labs = list(time = facet.Labels)
  ); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab)
  plot <- plot + labs(title = title.lab)
  # plot <- plot + labs(subtitle = subtitle.lab)
  # plot <- plot + labs(caption = caption.lab)
  plot <- plot + theme(legend.position = "none")
  
  # plot <- plot +
  #   stat_pvalue_manual(
  #     stat.test,
  #     bracket.nudge.y = 0.01,
  #     size = 5,
  #     hide.ns = TRUE,
  #     step.increase = 0.025,
  #     label = "p.adj.signif"
  #   )
  
  plot
  
  fileName <-   paste0(taxaLevel, "_Kendall_Uniqueness_by_PEWL_kg_by_Timepoint_Kendall_Scatter.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)


#####- Patient Kendall Uniqueness over time --------------------------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  x.lab <- "Timepoint"
  y.lab <- "Kendall Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  # subtitle.lab <- ""
  plotList <- list()
  index <- 1
  
  myT <- files[[which(names(files) == taxaLevel)]]
  time <- factor(myT$time, levels = timepoints)
  PatientID <- myT$PatientID
  Kendall_Uniqueness <- myT$Kendall_Uniqueness
  Surgery <- myT$Surgery
  
  df <- data.frame(PatientID, Kendall_Uniqueness, time, Surgery)
  
  ## Filter out patients with uniqueness at all timepoints
  patients <- vector()
  for (id in unique(df$PatientID)) {
    count <- length(which(df$PatientID == id))
    if (count == 6) {
      patients <- c(patients, id)
    }
  }
  
  df.1 <- df[df$PatientID %in% patients,]
  
  ## Plot patient uniqueness over time
  plot <- ggplot(data = df.1, aes(x = time, y = Kendall_Uniqueness, group = PatientID, color = PatientID)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  # plot <- facet(plot, facet.by = "Surgery")
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  # plot <- plot + labs(subtitle = subtitle.lab); plot
  # plot <- plot + labs(caption = caption.lab); plot
  # plot <- plot + theme(legend.position = "none"); plot
  plot <- plot + guides(color = "none"); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  ## Find top / bottom 10 most unique patients at BL
  df.2 <- df.1[df.1$time == 0,]
  df.2 <- df.2[order(df.2$Kendall_Uniqueness),]
  
  bottomIDs <- df.2$PatientID[1:10]
  topIDs <- df.2$PatientID[(nrow(df.2)-9):nrow(df.2)]
  
  
  ## Plot patient uniqueness over time from bottom 10 unique patients
  df.3 <- df[df$PatientID %in% bottomIDs,]
  df.3$Group <- "Bottom 10"
  title.lab <- paste0(str_to_title(taxaLevel), " - Bottom 10 unique patients at BL")
  
  plot <- ggplot(data = df.3, aes(x = time, y = Kendall_Uniqueness, group = PatientID, color = PatientID)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  # plot <- plot + theme(legend.position = "none"); plot
  plot <- plot + guides(color = "none"); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  ## Plot patient uniqueness over time from top 10 unique patients
  df.4 <- df[df$PatientID %in% topIDs,]
  df.4$Group <- "Top 10"
  title.lab <- paste0(str_to_title(taxaLevel), " - Top 10 unique patients at BL")
  
  plot <- ggplot(data = df.4, aes(x = time, y = Kendall_Uniqueness, group = PatientID, color = PatientID)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  # plot <- plot + theme(legend.position = "none"); plot
  plot <- plot + guides(color = "none"); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  ## Merge top & bottom uniqueness group tables
  df.5 <- rbind(df.3, df.4)
  title.lab <- paste0(str_to_title(taxaLevel), " - Top & Bottom 10 unique patients at BL")
  
  plot <- ggplot(data = df.5, aes(x = time, y = Kendall_Uniqueness, group = PatientID, color = Group)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  
  ## Find top / bottom 15 most unique patients at BL
  df.2 <- df.1[df.1$time == 0,]
  df.2 <- df.2[order(df.2$Kendall_Uniqueness),]
  
  bottomIDs <- df.2$PatientID[1:15]
  topIDs <- df.2$PatientID[(nrow(df.2)-14):nrow(df.2)]
  
  
  ## Plot patient uniqueness over time from bottom 15 unique patients
  df.3 <- df[df$PatientID %in% bottomIDs,]
  df.3$Group <- "Bottom 15"
  title.lab <- paste0(str_to_title(taxaLevel), " - Bottom 15 unique patients at BL")
  
  plot <- ggplot(data = df.3, aes(x = time, y = Kendall_Uniqueness, group = PatientID, color = PatientID)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  # plot <- plot + theme(legend.position = "none"); plot
  plot <- plot + guides(color = "none"); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  ## Plot patient uniqueness over time from top 15 unique patients
  df.4 <- df[df$PatientID %in% topIDs,]
  df.4$Group <- "Top 15"
  title.lab <- paste0(str_to_title(taxaLevel), " - Top 15 unique patients at BL")
  
  plot <- ggplot(data = df.4, aes(x = time, y = Kendall_Uniqueness, group = PatientID, color = PatientID)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  # plot <- plot + theme(legend.position = "none"); plot
  plot <- plot + guides(color = "none"); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  ## Merge top & bottom uniqueness group tables
  df.5 <- rbind(df.3, df.4)
  title.lab <- paste0(str_to_title(taxaLevel), " - Top & Bottom 15 unique patients at BL")
  
  plot <- ggplot(data = df.5, aes(x = time, y = Kendall_Uniqueness, group = PatientID, color = Group)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  fileName <-   paste0(taxaLevel, "_PatientID_Kendall_Uniqueness_by_Timepoint.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  for (i in 1:length(plotList)) {
    grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
  }  
  dev.off()
  
} # for (taxaLevel in taxaLevels)


#####- Patient Bray-Curtis Uniqueness over time --------------------------------------------------######
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  outputDirLevelPlots <- paste0(outputDirLevel, "Plots/")
  outputDirLevelResults <- paste0(outputDirLevel, "Results/")
  
  x.lab <- "Timepoint"
  y.lab <- "Bray-Curtis Uniqueness"
  title.lab <- str_to_title(taxaLevel)
  # subtitle.lab <- ""
  plotList <- list()
  index <- 1
  
  myT <- files[[which(names(files) == taxaLevel)]]
  time <- factor(myT$time, levels = timepoints)
  PatientID <- myT$PatientID
  Bray_Uniqueness <- myT$Bray_Uniqueness
  Surgery <- myT$Surgery
  
  df <- data.frame(PatientID, Bray_Uniqueness, time, Surgery)
  
  ## Filter out patients with uniqueness at all timepoints
  patients <- vector()
  for (id in unique(df$PatientID)) {
    count <- length(which(df$PatientID == id))
    if (count == 6) {
      patients <- c(patients, id)
    }
  }
  
  df.1 <- df[df$PatientID %in% patients,]
  
  ## Plot patient uniqueness over time
  plot <- ggplot(data = df.1, aes(x = time, y = Bray_Uniqueness, group = PatientID, color = PatientID)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  # plot <- facet(plot, facet.by = "Surgery")
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  # plot <- plot + labs(subtitle = subtitle.lab); plot
  # plot <- plot + labs(caption = caption.lab); plot
  # plot <- plot + theme(legend.position = "none"); plot
  plot <- plot + guides(color = "none"); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  ## Find top / bottom 10 most unique patients at BL
  df.2 <- df.1[df.1$time == 0,]
  df.2 <- df.2[order(df.2$Bray_Uniqueness),]
  
  bottomIDs <- df.2$PatientID[1:10]
  topIDs <- df.2$PatientID[(nrow(df.2)-9):nrow(df.2)]
  
  
  ## Plot patient uniqueness over time from bottom 10 unique patients
  df.3 <- df[df$PatientID %in% bottomIDs,]
  df.3$Group <- "Bottom 10"
  title.lab <- paste0(str_to_title(taxaLevel), " - Bottom 10 unique patients at BL")
  
  plot <- ggplot(data = df.3, aes(x = time, y = Bray_Uniqueness, group = PatientID, color = PatientID)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  # plot <- plot + theme(legend.position = "none"); plot
  plot <- plot + guides(color = "none"); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  ## Plot patient uniqueness over time from top 10 unique patients
  df.4 <- df[df$PatientID %in% topIDs,]
  df.4$Group <- "Top 10"
  title.lab <- paste0(str_to_title(taxaLevel), " - Top 10 unique patients at BL")
  
  plot <- ggplot(data = df.4, aes(x = time, y = Bray_Uniqueness, group = PatientID, color = PatientID)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  # plot <- plot + theme(legend.position = "none"); plot
  plot <- plot + guides(color = "none"); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  ## Merge top & bottom uniqueness group tables
  df.5 <- rbind(df.3, df.4)
  title.lab <- paste0(str_to_title(taxaLevel), " - Top & Bottom 10 unique patients at BL")
  
  plot <- ggplot(data = df.5, aes(x = time, y = Bray_Uniqueness, group = PatientID, color = Group)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  
  ## Find top / bottom 15 most unique patients at BL
  df.2 <- df.1[df.1$time == 0,]
  df.2 <- df.2[order(df.2$Bray_Uniqueness),]
  
  bottomIDs <- df.2$PatientID[1:15]
  topIDs <- df.2$PatientID[(nrow(df.2)-14):nrow(df.2)]
  
  
  ## Plot patient uniqueness over time from bottom 15 unique patients
  df.3 <- df[df$PatientID %in% bottomIDs,]
  df.3$Group <- "Bottom 15"
  title.lab <- paste0(str_to_title(taxaLevel), " - Bottom 15 unique patients at BL")
  
  plot <- ggplot(data = df.3, aes(x = time, y = Bray_Uniqueness, group = PatientID, color = PatientID)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  # plot <- plot + theme(legend.position = "none"); plot
  plot <- plot + guides(color = "none"); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  ## Plot patient uniqueness over time from top 15 unique patients
  df.4 <- df[df$PatientID %in% topIDs,]
  df.4$Group <- "Top 15"
  title.lab <- paste0(str_to_title(taxaLevel), " - Top 15 unique patients at BL")
  
  plot <- ggplot(data = df.4, aes(x = time, y = Bray_Uniqueness, group = PatientID, color = PatientID)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  # plot <- plot + theme(legend.position = "none"); plot
  plot <- plot + guides(color = "none"); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  ## Merge top & bottom uniqueness group tables
  df.5 <- rbind(df.3, df.4)
  title.lab <- paste0(str_to_title(taxaLevel), " - Top & Bottom 15 unique patients at BL")
  
  plot <- ggplot(data = df.5, aes(x = time, y = Bray_Uniqueness, group = PatientID, color = Group)) +
    geom_line() +
    geom_point(aes(shape = Surgery)); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  plot <- plot + labs(title = title.lab); plot
  
  plotList[[index]] <- plot
  index <- index + 1
  
  fileName <-   paste0(taxaLevel, "_PatientID_Bray_Uniqueness_by_Timepoint.pdf")
  file.path <- paste0(outputDirLevelPlots, fileName)
  pdf(file.path, width = 10, height = 7)
  for (i in 1:length(plotList)) {
    grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
  }  
  dev.off()
  
} # for (taxaLevel in taxaLevels)

