#Author: Alicia Sorgen
#Date: 2022 September 30
#Description: 

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/UNCC/Projects/Bariatric_Surgery/Git_Repositories/gut-microbiome-bariatric-weight-outcomes")
params <- c(params, "MetaPhlAn2")

moduleRoot <- "4.0_Diversity_Metrics"

taxaLevels <- vector()
taxaLevels <- c(taxaLevels, "phylum")
# taxaLevels <- c(taxaLevels, "class")
# taxaLevels <- c(taxaLevels, "order")
# taxaLevels <- c(taxaLevels, "family")
taxaLevels <- c(taxaLevels, "genus")
# taxaLevels <- c(taxaLevels, "species")

##### Libraries #####
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
}

if (args[1] == "BLJ") {
  pipeRoot  <- dirname(dirname(getwd()))
  moduleDir <- dirname(getwd())
}



##### Set up functions file #####
funcScript <- if (args[1] == "BLJ") file.path(moduleDir, "resources", "functions.R") else file.path(gitScripts, "functions.R")
source(funcScript)

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]
rm(str, funcScript)

##### Set up input #####
classifier <- args[2]

prevModule <- paste0("TaxaMetaMerge")
inputDir <- paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Script variables #####
surgeryPalette <- c("steelblue", "gold")
outcomePalette <- c("steelblue", "tomato")
colorPalette <- c("pink", "purple", "orange", "green", "grey", "steelblue", "tomato")

month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")
##### Calculate Alpha Diversity #####
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  dir.create(outputDirLevel, showWarnings = FALSE)
  
  file.path <- paste0(inputDir, taxaLevel,"_LogNormalizedCounts_", classifier, ".tsv")
  myT <- read.table(file.path, sep='\t', header = TRUE)
  rownames(myT) <- myT$SampleID
  
  taxaStart <- which(colnames(myT) == "ResponderStatus") + 1
  
  taxaTable <- myT[,taxaStart:ncol(myT)]
  
  myT$ShannonIndex <- diversity(taxaTable, index = "shannon")
  myT$SimpsonIndex <- diversity(taxaTable, index = "simpson")
  myT$Richness     <- specnumber(taxaTable)
  myT$Evenness     <- ifelse(myT$Richness > 1, myT$ShannonIndex / log(myT$Richness), 0)
  
  file.path <- paste0(outputDirLevel, taxaLevel,"_LogNormalizedCounts_", classifier, "_alpha.tsv")
  write.table(myT, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
} # for (taxaLevel in taxaLevels) {

##### Shannon Diversity ~ Surgery Box Plots -- Wilcoxon (grouped by time) #####
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")

  file.path <- paste0(outputDirLevel, taxaLevel,"_LogNormalizedCounts_", classifier, "_alpha.tsv")
  myT <- read.table(file.path, sep='\t', header = TRUE)
  rownames(myT) <- myT$SampleID
  
  myT2 <- myT[!(is.na(myT$time)),]
  
  x.lab <- "Surgery Type"
  y.lab <- "Shannon Index"
  title.lab <- ""
  subtitle.lab <- ""
  caption.lab <- "Wilcoxon rank sum test"
  
  
  stat.test <- myT2 %>%
    group_by(time) %>%
    wilcox_test(ShannonIndex ~ Surgery) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_ShannonIndex_by_Surgery_Timepoint_Wilcoxon_Results.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  
  stat.test <- stat.test %>%
    add_xy_position(x = "time", fun = "mean_se")
  
  plot <- ggboxplot(
    myT2, x = "Surgery", y = "ShannonIndex", color = "black",
    fill = "Surgery", palette = surgeryPalette,
    # facet.by = "month",
    scales = "free", add = "jitter"
  ); plot
  
  plot <- facet(plot, facet.by = "time"
                , panel.labs = list(time = month.labs)
  ); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  # plot <- plot + labs(title = title.lab); plot
  # plot <- plot + labs(subtitle = subtitle.lab); plot
  plot <- plot + labs(caption = caption.lab); plot

  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.05,
      size = 10,
      hide.ns = TRUE,
      label = "p.adj.signif"
    ); plot
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_ShannonIndex_by_Surgery_Timepoint_Wilcoxon_BoxPlot.pdf")
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)


##### Shannon Diversity ~ Surgery Box Plots -- t-test (grouped by time) #####
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  
  file.path <- paste0(outputDirLevel, taxaLevel,"_LogNormalizedCounts_", classifier, "_alpha.tsv")
  myT <- read.table(file.path, sep='\t', header = TRUE)
  rownames(myT) <- myT$SampleID
  
  myT2 <- myT[!(is.na(myT$time)),]
  
  x.lab <- "Surgery Type"
  y.lab <- "Shannon Index"
  title.lab <- ""
  subtitle.lab <- ""
  caption.lab <- "t-test"
  
  
  stat.test <- myT2 %>%
    group_by(time) %>%
    t_test(ShannonIndex ~ Surgery) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_ShannonIndex_by_Surgery_Timepoint_ttest_Results.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  
  stat.test <- stat.test %>%
    add_xy_position(x = "time", fun = "mean_se")
  
  plot <- ggboxplot(
    myT2, x = "Surgery", y = "ShannonIndex", color = "black",
    fill = "Surgery", palette = surgeryPalette,
    # facet.by = "month",
    scales = "free", add = "jitter"
  ); plot
  
  plot <- facet(plot, facet.by = "time"
                , panel.labs = list(time = month.labs)
  ); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  # plot <- plot + labs(title = title.lab); plot
  # plot <- plot + labs(subtitle = subtitle.lab); plot
  plot <- plot + labs(caption = caption.lab); plot
  
  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.05,
      size = 10,
      hide.ns = TRUE,
      label = "p.adj.signif"
    ); plot
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_ShannonIndex_by_Surgery_Timepoint_ttest_BoxPlot.pdf")
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)



##### Shannon Diversity ~ Site Box Plots -- Wilcoxon (grouped by time) #####
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  
  file.path <- paste0(outputDirLevel, taxaLevel,"_LogNormalizedCounts_", classifier, "_alpha.tsv")
  myT <- read.table(file.path, sep='\t', header = TRUE)
  rownames(myT) <- myT$SampleID
  
  myT2 <- myT[!(is.na(myT$time)),]
  
  x.lab <- "Site Type"
  y.lab <- "Shannon Index"
  title.lab <- ""
  subtitle.lab <- ""
  caption.lab <- "Wilcoxon rank sum test"
  
  
  stat.test <- myT2 %>%
    group_by(time) %>%
    wilcox_test(ShannonIndex ~ Site) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_ShannonIndex_by_Site_Timepoint_Wilcoxon_Results.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  
  stat.test <- stat.test %>%
    add_xy_position(x = "time", fun = "mean_se")
  
  plot <- ggboxplot(
    myT2, x = "Site", y = "ShannonIndex", color = "black",
    fill = "Site", palette = surgeryPalette,
    # facet.by = "month",
    scales = "free", add = "jitter"
  ); plot
  
  plot <- facet(plot, facet.by = "time"
                , panel.labs = list(time = month.labs)
  ); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  # plot <- plot + labs(title = title.lab); plot
  # plot <- plot + labs(subtitle = subtitle.lab); plot
  plot <- plot + labs(caption = caption.lab); plot
  
  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.05,
      size = 10,
      hide.ns = TRUE,
      label = "p.adj.signif"
    ); plot
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_ShannonIndex_by_Site_Timepoint_Wilcoxon_BoxPlot.pdf")
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)



##### Bray-Curtis ~ Surgery PCoA #####
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")

  file.path <- paste0(inputDir, taxaLevel,"_LogNormalizedCounts_", classifier, ".tsv")
  myT <- read.table(file.path, sep='\t', header = TRUE)
  # rownames(myT) <- myT$SampleID
  
  taxaStart <- which(colnames(myT) == "ResponderStatus") + 1
  
  taxaTable <- myT[,taxaStart:ncol(myT)]
  taxaTable$Surgery <- myT$Surgery
  taxaTable$Timepoint <- myT$time
  taxaTable <- na.omit(taxaTable)
  end <- which(colnames(taxaTable) == "Surgery") -1

  myMDS <- capscale(taxaTable[,1:end] ~ 1, distance="bray")
  percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
  pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
  
  col2=adjustcolor(surgeryPalette[factor(myT$Surgery)], alpha.f = 1)
  
  adon.results<-adonis2(taxaTable[,1:end] ~ taxaTable$Surgery, method="bray",perm=999)
  # print(adon.results)
  # capture.output(adon.results, file = paste0(outputDir, level, "_",studyName,"_adonis_Surgery.txt"))
  Title <- paste0(str_to_title(taxaLevel), " PERMANOVA p = ", adon.results$`Pr(>F)`[1])
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_BrayCurtis_by_Surgery_PCoA.pdf")
  pdf(file.path)
  par(mar=c(5, 5, 5, 5), xpd=FALSE, mfrow = c(1,1))
  
  all.combindations <- combn(1:5, 2)
  for (combination in 1:ncol(all.combindations)) {
    PCoA_a <- all.combindations[1,combination]
    PCoA_b <- all.combindations[2,combination]
    
    pcoa12 <- ordiplot(myMDS, choices = c(PCoA_a,PCoA_b), display = "sites", type = "none", cex.lab = 2, 
                       xlab = pcoa_p[PCoA_a], ylab = pcoa_p[PCoA_b], main = Title)
    points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
    for (n in 1:length(unique(na.omit(taxaTable$Surgery)))) {
      ordiellipse(pcoa12, taxaTable$Surgery, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                  col = surgeryPalette[n], show.groups = levels(factor(taxaTable$Surgery))[n], label = T, 
                  font = 2, cex = 1)
    }
    legend("topright", unique(taxaTable$Surgery),
           col = surgeryPalette, cex = 1.5, pch = 16, bty = "n")
    
  }
  
  
  dev.off()
  
  
} # for (taxaLevel in taxaLevels)

##### Bray-Curtis ~ Site PCoA #####
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  
  file.path <- paste0(inputDir, taxaLevel,"_LogNormalizedCounts_", classifier, ".tsv")
  myT <- read.table(file.path, sep='\t', header = TRUE)
  # rownames(myT) <- myT$SampleID
  
  taxaStart <- which(colnames(myT) == "ResponderStatus") + 1
  
  taxaTable <- myT[,taxaStart:ncol(myT)]
  taxaTable$Site <- myT$Site
  taxaTable$Timepoint <- myT$time
  taxaTable <- na.omit(taxaTable)
  end <- which(colnames(taxaTable) == "Site") -1
  
  myMDS <- capscale(taxaTable[,1:end] ~ 1, distance="bray")
  percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
  pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
  
  col2=adjustcolor(colorPalette[factor(myT$Site)], alpha.f = 1)
  
  adon.results<-adonis2(taxaTable[,1:end] ~ taxaTable$Site, method="bray",perm=999)
  # print(adon.results)
  # capture.output(adon.results, file = paste0(outputDir, level, "_",studyName,"_adonis_Site.txt"))
  Title <- paste0(str_to_title(taxaLevel), " PERMANOVA p = ", adon.results$`Pr(>F)`[1])
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_BrayCurtis_by_Site_PCoA.pdf")
  pdf(file.path)
  par(mar=c(5, 5, 5, 5), xpd=FALSE, mfrow = c(1,1))
  
  all.combindations <- combn(1:5, 2)
  for (combination in 1:ncol(all.combindations)) {
    PCoA_a <- all.combindations[1,combination]
    PCoA_b <- all.combindations[2,combination]
    
    pcoa12 <- ordiplot(myMDS, choices = c(PCoA_a,PCoA_b), display = "sites", type = "none", cex.lab = 2, 
                       xlab = pcoa_p[PCoA_a], ylab = pcoa_p[PCoA_b], main = Title)
    points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
    for (n in 1:length(unique(na.omit(taxaTable$Site)))) {
      ordiellipse(pcoa12, taxaTable$Site, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                  col = colorPalette[n], show.groups = levels(factor(taxaTable$Site))[n], label = T, 
                  font = 2, cex = 1)
    }
    legend("topright", unique(taxaTable$Site),
           col = colorPalette, cex = 1.5, pch = 16, bty = "n")
    
  }
  
  
  dev.off()
  
  
} # for (taxaLevel in taxaLevels)

##### Bray-Curtis ~ TertileAssignment_12M PCoA #####
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  
  file.path <- paste0(inputDir, taxaLevel,"_LogNormalizedCounts_", classifier, ".tsv")
  myT <- read.table(file.path, sep='\t', header = TRUE)
  # rownames(myT) <- myT$SampleID
  
  taxaStart <- which(colnames(myT) == "ResponderStatus") + 1
  
  taxaTable <- myT[,taxaStart:ncol(myT)]
  taxaTable$TertileAssignment_12M <- myT$TertileAssignment_12M
  taxaTable$TertileAssignment_12M <- as.factor(taxaTable$TertileAssignment_12M)
  
  # taxaTable$Timepoint <- myT$time
  taxaTable <- na.omit(taxaTable)
  end <- which(colnames(taxaTable) == "TertileAssignment_12M") -1
  
  myMDS <- capscale(taxaTable[,1:end] ~ 1, distance="bray")
  percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
  pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
  
  col2=adjustcolor(colorPalette[factor(taxaTable$TertileAssignment_12M)], alpha.f = 1)
  
  adon.results<-adonis2(taxaTable[,1:end] ~ taxaTable$TertileAssignment_12M, method="bray",perm=999)
  # print(adon.results)
  capture.output(adon.results, file = paste0(outputDirLevel, taxaLevel, "_", classifier, "_PERMANOVA_TertileAssignment_12M.txt"))
  Title <- paste0("Patient Tertiles 12 months post-op\n", str_to_title(taxaLevel), " PERMANOVA p = ", adon.results$`Pr(>F)`[1])
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_BrayCurtis_by_TertileAssignment_12M_PCoA.pdf")
  pdf(file.path)
  par(mar=c(5, 5, 5, 5), xpd=FALSE, mfrow = c(1,1))
  
  all.combindations <- combn(1:5, 2)
  for (combination in 1:ncol(all.combindations)) {
    PCoA_a <- all.combindations[1,combination]
    PCoA_b <- all.combindations[2,combination]
    
    pcoa12 <- ordiplot(myMDS, choices = c(PCoA_a,PCoA_b), display = "sites", type = "none", cex.lab = 2, 
                       xlab = pcoa_p[PCoA_a], ylab = pcoa_p[PCoA_b], main = Title)
    points(pcoa12, "sites", col = taxaTable$TertileAssignment_12M, pch = 16, cex = 1.5)
    # points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
    # for (n in 1:length(unique(na.omit(taxaTable$TertileAssignment_12M)))) {
    #   ordiellipse(pcoa12, taxaTable$TertileAssignment_12M, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
    #               col = colorPalette[n], show.groups = levels(factor(taxaTable$TertileAssignment_12M))[n], label = T, 
    #               font = 2, cex = 1)
    # }
    legend("topright", legend = levels(taxaTable$TertileAssignment_12M),
           col = 1:length(levels(taxaTable$TertileAssignment_12M)), 
           cex = 1.5, pch = 16, bty = "n")
    
  }
  
  
  dev.off()
  
  
} # for (taxaLevel in taxaLevels)

##### Bray-Curtis ~ TertileAssignment_24M PCoA #####
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  
  file.path <- paste0(inputDir, taxaLevel,"_LogNormalizedCounts_", classifier, ".tsv")
  myT <- read.table(file.path, sep='\t', header = TRUE)
  # rownames(myT) <- myT$SampleID
  
  taxaStart <- which(colnames(myT) == "ResponderStatus") + 1
  
  taxaTable <- myT[,taxaStart:ncol(myT)]
  taxaTable$TertileAssignment_24M <- myT$TertileAssignment_24M
  taxaTable$TertileAssignment_24M <- as.factor(taxaTable$TertileAssignment_24M)
  
  # taxaTable$Timepoint <- myT$time
  taxaTable <- na.omit(taxaTable)
  end <- which(colnames(taxaTable) == "TertileAssignment_24M") -1
  
  myMDS <- capscale(taxaTable[,1:end] ~ 1, distance="bray")
  percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
  pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
  
  col2=adjustcolor(colorPalette[factor(taxaTable$TertileAssignment_24M)], alpha.f = 1)
  
  adon.results<-adonis2(taxaTable[,1:end] ~ taxaTable$TertileAssignment_24M, method="bray",perm=999)
  # print(adon.results)
  capture.output(adon.results, file = paste0(outputDirLevel, taxaLevel, "_", classifier, "_PERMANOVA_TertileAssignment_24M.txt"))
  Title <- paste0("Patient Tertiles 24 months post-op\n", str_to_title(taxaLevel), " PERMANOVA p = ", adon.results$`Pr(>F)`[1])
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_BrayCurtis_by_TertileAssignment_24M_PCoA.pdf")
  pdf(file.path)
  par(mar=c(5, 5, 5, 5), xpd=FALSE, mfrow = c(1,1))
  
  all.combindations <- combn(1:5, 2)
  for (combination in 1:ncol(all.combindations)) {
    PCoA_a <- all.combindations[1,combination]
    PCoA_b <- all.combindations[2,combination]
    
    pcoaPlot <- ordiplot(myMDS, choices = c(PCoA_a,PCoA_b), display = "sites", type = "none", cex.lab = 2, 
                       xlab = pcoa_p[PCoA_a], ylab = pcoa_p[PCoA_b], main = Title)
    points(pcoaPlot, "sites", col = taxaTable$TertileAssignment_24M, pch = 16, cex = 1.5)
    # points(pcoaPlot, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
    # for (n in 1:length(unique(na.omit(taxaTable$TertileAssignment_24M)))) {
    #   ordiellipse(pcoaPlot, taxaTable$TertileAssignment_24M, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
    #               col = colorPalette[n], show.groups = levels(factor(taxaTable$TertileAssignment_24M))[n], label = T, 
    #               font = 2, cex = 1)
    # }
    legend("topright", legend = levels(taxaTable$TertileAssignment_24M),
           col = 1:length(levels(taxaTable$TertileAssignment_24M)), 
           cex = 1.5, pch = 16, bty = "n")
    
  }
  
  
  dev.off()
  
  
} # for (taxaLevel in taxaLevels)

##### Bray-Curtis ~ Surgery PCoA #####
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  
  file.path <- paste0(inputDir, taxaLevel,"_LogNormalizedCounts_", classifier, ".tsv")
  myT <- read.table(file.path, sep='\t', header = TRUE)
  # rownames(myT) <- myT$SampleID
  
  taxaStart <- which(colnames(myT) == "ResponderStatus") + 1
  
  taxaTable <- myT[,taxaStart:ncol(myT)]
  taxaTable$Surgery <- myT$Surgery
  taxaTable$Timepoint <- myT$time
  taxaTable$Timepoint <- as.factor(taxaTable$Timepoint)
  taxaTable <- na.omit(taxaTable)
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_BrayCurtis_by_Surgery_by_Timepoint_PCoA.pdf")
  pdf(file.path)
  
  Month <- vector()
  p <- vector()
  R2 <- vector()
  Fval <- vector()
  
  for (TP in levels(taxaTable$Timepoint)) {
    
    taxaTable2 <- taxaTable[taxaTable$Timepoint == TP,]
    end <- which(colnames(taxaTable2) == "Surgery") -1
    myMDS <- capscale(taxaTable2[,1:end] ~ 1, distance="bray")
    percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
    pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
    
    col2=adjustcolor(surgeryPalette[factor(taxaTable2$Surgery)], alpha.f = 1)
    
    adon.results<-adonis2(taxaTable2[,1:end] ~ taxaTable2$Surgery, method="bray",perm=999)
    # print(adon.results)
    # capture.output(adon.results, file = paste0(outputDir, level, "_",studyName,"_adonis_Surgery.txt"))
    Title <- paste0(str_to_title(taxaLevel), " Bray-Curtis dissimilarity\nMonth ", TP, " - PERMANOVA p = ", adon.results$`Pr(>F)`[1])
    
    par(mar=c(5, 5, 5, 5), xpd=FALSE, mfrow = c(1,1))
    
    all.combindations <- combn(1:5, 2)
    for (combination in 1:ncol(all.combindations)) {
      PCoA_a <- all.combindations[1,combination]
      PCoA_b <- all.combindations[2,combination]
      
      pcoa12 <- ordiplot(myMDS, choices = c(PCoA_a,PCoA_b), display = "sites", type = "none", cex.lab = 2, 
                         xlab = pcoa_p[PCoA_a], ylab = pcoa_p[PCoA_b], main = Title)
      points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
      for (n in 1:length(unique(na.omit(taxaTable2$Surgery)))) {
        ordiellipse(pcoa12, taxaTable2$Surgery, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                    col = surgeryPalette[n], show.groups = levels(factor(taxaTable2$Surgery))[n], label = T, 
                    font = 2, cex = 1)
      }
      legend("topright", unique(taxaTable2$Surgery),
             col = surgeryPalette, cex = 1.5, pch = 16, bty = "n")
      
    }
    
    Month <- c(Month, TP)
    p <- c(p, adon.results$`Pr(>F)`[1])
    R2 <- c(R2, adon.results$R2[1])
    Fval <- c(Fval, adon.results$F[1])
    
  } # for (TP in levels(taxaTable$Timepoint))
  
  dev.off()
  
  dFrame <- data.frame(Month, p, R2, Fval)
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_BrayCurtis_by_Surgery_by_Timepoint_PERMANOVAResults.tsv")
  write.table(dFrame, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  
} # for (taxaLevel in taxaLevels)


##### Bray-Curtis ~ TertileAssignment_24M & Timepoint PCoA #####
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  
  file.path <- paste0(inputDir, taxaLevel,"_LogNormalizedCounts_", classifier, ".tsv")
  myT <- read.table(file.path, sep='\t', header = TRUE)
  # rownames(myT) <- myT$SampleID
  
  taxaStart <- which(colnames(myT) == "ResponderStatus") + 1
  
  taxaTable <- myT[,taxaStart:ncol(myT)]
  taxaTable$TertileAssignment_24M <- myT$TertileAssignment_24M
  taxaTable$TertileAssignment_24M <- gsub("1", "Tertile 1", taxaTable$TertileAssignment_24M)
  taxaTable$TertileAssignment_24M <- gsub("2", "Tertile 2", taxaTable$TertileAssignment_24M)
  taxaTable$TertileAssignment_24M <- gsub("3", "Tertile 3", taxaTable$TertileAssignment_24M)
  taxaTable$TertileAssignment_24M <- as.factor(taxaTable$TertileAssignment_24M)
  taxaTable$Timepoint <- myT$time
  taxaTable$Timepoint <- as.factor(taxaTable$Timepoint)
  
  taxaTable <- na.omit(taxaTable)
  end <- which(colnames(taxaTable) == "TertileAssignment_24M") -1
  
  Month <- vector()
  p <- vector()
  R2 <- vector()
  Fval <- vector()
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_BrayCurtis_by_TertileAssignment_24M_by_Timepoint_PCoA.pdf")
  pdf(file.path,width = 10, height = 5)
  par(mar=c(5, 5, 5, 5), xpd=FALSE, mfrow = c(2,3))
  
  for (month in unique(taxaTable$Timepoint)) {
    
    taxaTable2 <- taxaTable[taxaTable$Timepoint == month,]
    myMDS <- capscale(taxaTable2[,1:end] ~ 1, distance="bray")
    percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
    pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
    
    
    adon.results<-adonis2(taxaTable2[,1:end] ~ taxaTable2$TertileAssignment_24M, method="bray",perm=999)
    print(adon.results)
    # capture.output(adon.results, file = paste0(outputDirLevel, taxaLevel, "_", classifier, "_PERMANOVA_TertileAssignment_24M.txt"))
    
    Month <- c(Month, month)
    p <- c(p, adon.results$`Pr(>F)`[1])
    R2 <- c(R2, adon.results$R2[1])
    Fval <- c(Fval, adon.results$F[1])
    
    Title <- paste0("Weight Loss Tertiles (24-mo)\n", str_to_title(taxaLevel), " PERMANOVA p = ", adon.results$`Pr(>F)`[1])
    Subtitle <- month.labs[which(unique(taxaTable$Timepoint) == month)]
    
    
    all.combindations <- combn(1:2, 2)
    for (combination in 1:ncol(all.combindations)) {
      PCoA_a <- all.combindations[1,combination]
      PCoA_b <- all.combindations[2,combination]
      
      pcoaPlot <- ordiplot(myMDS, choices = c(PCoA_a,PCoA_b), display = "sites", type = "none", cex.lab = 2,
                           xlab = pcoa_p[PCoA_a], ylab = pcoa_p[PCoA_b], main = paste0(Title, "\n", Subtitle))
      colorPalette <- c("purple", "orange", "yellowgreen", "grey", "steelblue", "tomato", "pink")
      col2=adjustcolor(colorPalette[factor(taxaTable2$TertileAssignment_24M)], alpha.f = 1)
      # points(pcoaPlot, "sites", col = taxaTable2$TertileAssignment_24M, pch = 16, cex = 1.5)
      points(pcoaPlot, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
      for (n in 1:length(levels(taxaTable2$TertileAssignment_24M))) {
        ordiellipse(pcoaPlot, taxaTable2$TertileAssignment_24M, kind = "se", conf = 0.95, lwd = 2, draw = "lines",
                    col = col2[n], show.groups = levels(factor(taxaTable2$TertileAssignment_24M))[n], label = F,
                    font = 2, cex = 1)
      }
      legend("topright", 
             legend = levels(taxaTable2$TertileAssignment_24M),
             col = col2[1:length(levels(taxaTable2$TertileAssignment_24M))],
             cex = 0.75, pch = 16,
             horiz = FALSE)
      
    }
    
    
    
  }
  
  dev.off()
  
  dFrame <- data.frame(Month, p, R2, Fval)
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_BrayCurtis_by_TertileAssignment_24M_by_Timepoint_PERMANOVAResults.tsv")
  write.table(dFrame, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  
} # for (taxaLevel in taxaLevels)



##### Shannon Diversity ~ Tertile Box Plots -- ANOVA (grouped by time) #####
for (taxaLevel in taxaLevels) {
  
  outputDirLevel <- paste0(outputDir, taxaLevel, "/")
  
  file.path <- paste0(outputDirLevel, taxaLevel,"_LogNormalizedCounts_", classifier, "_alpha.tsv")
  myT <- read.table(file.path, sep='\t', header = TRUE)
  rownames(myT) <- myT$SampleID
  
  PatientID <- myT$PatientID
  Tertile <- myT$TertileAssignment_24M
  Timepoint <- myT$time
  ShannonIndex <- myT$ShannonIndex
  myT2 <- data.frame(PatientID, Tertile, Timepoint, ShannonIndex)
  myT2$Tertile <- as.factor(myT2$Tertile)
  
  myT2 <- na.omit(myT2)
  
  x.lab <- "Tertile"
  y.lab <- "Shannon Index"
  title.lab <- ""
  subtitle.lab <- "ANOVA"
  caption.lab <- "pairwise comparisons visa t-test"
  
  stat.test <- myT2 %>%
    group_by(Timepoint) %>%
    anova_test(ShannonIndex ~ Tertile) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_ShannonIndex_by_Tertile_Timepoint_LinearModel_Results.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  stat.test <- myT2 %>%
    group_by(Timepoint) %>%
    tukey_hsd(ShannonIndex ~ Tertile)
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_ShannonIndex_by_Tertile_Timepoint_Tukey_Results.tsv")
  write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  
  stat.test <- stat.test %>%
    add_xy_position(x = "Tertile", fun = "max")
  
  
  plot <- ggboxplot(
    myT2, x = "Tertile", y = "ShannonIndex", color = "black",
    fill = "Tertile", palette = c("blue", "pink", "red"),
    # facet.by = "month",
    scales = "free", add = "jitter"
  ); plot
  
  plot <- facet(plot, facet.by = "Timepoint"
                , panel.labs = list(Timepoint = month.labs)
  ); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  # plot <- plot + labs(title = title.lab); plot
  # plot <- plot + labs(subtitle = subtitle.lab); plot
  # plot <- plot + labs(caption = caption.lab); plot
  
  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.05,
      size = 7,
      hide.ns = TRUE,
      label = "p.adj.signif"
    ); plot
  
  file.path <- paste0(outputDirLevel, taxaLevel, "_", classifier, "_ShannonIndex_by_Tertile_Timepoint_BoxPlot.pdf")
  pdf(file.path, width = 10, height = 7)
  print(plot)
  dev.off()
  
} # for (taxaLevel in taxaLevels)



