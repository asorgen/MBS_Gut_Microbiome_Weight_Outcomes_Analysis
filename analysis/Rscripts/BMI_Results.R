#Author: Alicia Sorgen
#Date: 2022 June 13
#Description: BMI analysis

##### Edits for script #####
rm(list=ls())
params <- vector()

ANALYSIS <- "microbiome_n124"
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
params <- c(params, "weight_update_BLonly_excluded.txt")

moduleRoot <- paste0("BMI_Results")


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
  
  if (length(args) > 1) {
    if (args[2] == "MetaPhlAn2" | args[2] == "Kraken2") {
      module <- paste0(args[2], module)
    } 
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
inputDir = paste0(pipeRoot, "/input/metadataTables/")
inputFile <- args[2]

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Read in tables #####
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
myTable$time = as.numeric(gsub("BIO-.-...-","",myTable$SampleID))
myTable$Surgery <- myTable$TypeofSurgery
myTable$Surgery <- ifelse(myTable$Surgery == "Gastric Bypass", "RYGB",
                          ifelse(myTable$Surgery == "Sleeve Gastrectomy", "SG",
                                 NA))
myTable$Timepoint <- as.factor(myTable$Timepoint)
myTable$Timepoint <- factor(myTable$Timepoint, levels = c("BL", "ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR"))
row.names(myTable) = myTable$SampleID
# BL.df <- myTable[myTable$Timepoint == "BL",]
myTable <- myTable[!is.na(myTable$Surgery),]


##### RYGB patient BMI over time - Mixed Linear Model #####
Filt <- myTable[myTable$Surgery == "RYGB",]
Filt <- na.omit(Filt)

mlm <- lme(BMI_kgm2 ~ Timepoint, method = 'REML', random = ~1 | PatientID, data = Filt)
smry <- summary(mlm)
lm.df <- data.frame(smry$tTable)
lm.df$Comparison <- row.names(lm.df)
lm.df$Comparison <- gsub(pattern = "Timepoint", replacement = "BL v ", lm.df$Comparison)
lm.df$Metric <- "BMI"
lm.df$Adj_pvalue <- p.adjust(lm.df$p.value, method = "BH")
lm.df$p_value <- roundP(lm.df$Adj_pvalue)

file.path <- paste0(outputDir,"RYGB_BMI_LM_changes_over_time.tsv")
write.table(lm.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)

##### SG patient BMI over time - Mixed Linear Model #####
Filt <- myTable[myTable$Surgery == "SG",]
Filt <- na.omit(Filt)

mlm <- lme(BMI_kgm2 ~ Timepoint, method = 'REML', random = ~1 | PatientID, data = Filt)
smry <- summary(mlm)
lm.df <- data.frame(smry$tTable)
lm.df$Comparison <- row.names(lm.df)
lm.df$Comparison <- gsub(pattern = "Timepoint", replacement = "BL v ", lm.df$Comparison)
lm.df$Metric <- "BMI"
lm.df$Adj_pvalue <- p.adjust(lm.df$p.value, method = "BH")
lm.df$p_value <- roundP(lm.df$Adj_pvalue)


file.path <- paste0(outputDir, "SG_BMI_LM_changes_over_time.tsv")
write.table(lm.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)

myTable$Timepoint <- factor(myTable$Timepoint)

# model <- aov(BMI_kgm2 ~ Surgery * Timepoint + Error(PatientID / Timepoint), data = myTable)
##### BMI differences in surgery type at each timepoint - wilcox #####
months <- c(0, 1, 6, 12, 18, 24)

Timepoint <- vector()
Avg_BMI_RYGB <- vector()
SD_BMI_RYGB <- vector()
Avg_BMI_SG <- vector()
SD_BMI_SG <- vector()
p_value <- vector()
index <- 1

for (month in months) {
  
  Timepoint[index] <- month
  Filt <- myTable[myTable$time == month,]
  
  Avg_BMI_RYGB[index] <- mean(Filt$BMI_kgm2[Filt$Surgery=="RYGB"], na.rm = TRUE)
  SD_BMI_RYGB[index] <- sd(Filt$BMI_kgm2[Filt$Surgery=="RYGB"], na.rm = TRUE)
  
  Avg_BMI_SG[index] <- mean(Filt$BMI_kgm2[Filt$Surgery=="SG"], na.rm = TRUE)
  SD_BMI_SG[index] <- sd(Filt$BMI_kgm2[Filt$Surgery=="SG"], na.rm = TRUE)
  
  wilcox.test <- wilcox.test(Filt$BMI_kgm2 ~ Filt$Surgery)
  p_value[index] <- wilcox.test$p.value
  
  index <- index + 1
}

result.df <- data.frame(Timepoint, Avg_BMI_RYGB, SD_BMI_RYGB, Avg_BMI_SG, SD_BMI_SG, p_value)

result.df$BH_p <- p.adjust(result.df$p_value, method = "BH")
result.df$p_value_rounded <- roundP(result.df$BH_p)
result.df$significance <- sigStars(result.df$BH_p)

m <- c(0, 1, 6, 12, 18, 24)
bio.conc <- vector()
index <- 1
for (i in 1:nrow(result.df)) {
  
  RYGB.mean <- paste0("(", round(result.df$Avg_BMI_RYGB[i], 2), " +/- ", round(result.df$SD_BMI_RYGB[i], 2), " kg/m2)")
  SG.mean <- paste0("(", round(result.df$Avg_BMI_SG[i], 2), " +/- ", round(result.df$SD_BMI_SG[i], 2), " kg/m2)")
  
  if (result.df$p_value[i] < 0.05) {
    bio.conc[index] <- paste0("There was a significant difference in BMI in RYGB ", RYGB.mean, " and SG patients ", SG.mean, " at ", m[i], " month(s) (", result.df$p_value_rounded[i], ", Wilcoxon).")
    print(bio.conc[index])
  } else {
    bio.conc[index] <- paste0("There was NO significant difference in BMI in RYGB ", RYGB.mean, " and SG patients ", SG.mean, " at ", m[i], " month(s) (", result.df$p_value_rounded[i], ", Wilcoxon).")
    print(bio.conc[index])
  }
  
  
  index <- index + 1
  
}
result.df$Conclusion <- bio.conc
file.path <- paste0(outputDir,"SurgeryType_BMI_wilcox_at_each_timepoint.tsv")
write.table(result.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)

##### BMI loss over time between surgery types - wilcox #####
# Find baseline rows
bl <- which(myTable$Timepoint == "BL")

# Extract baseline BMI
blb <- myTable[bl, "BMI_kgm2"]
names(blb) = myTable[bl, "PatientID"]

# Add info to each row
myTable$Baseline_BMI = blb[myTable$PatientID]

# Calculate average loss from baseline
myTable$Loss_from_BL_BMI <- myTable$Baseline_BMI - myTable$BMI_kgm2

months <- c(1, 6, 12, 18, 24)
Timepoint <- vector()
Avg_Loss_RYGB <- vector()
SD_Loss_RYGB <- vector()
Avg_Loss_SG <- vector()
SD_Loss_SG <- vector()
p_value <- vector()
index <- 1

for (month in months) {
  
  Timepoint[index] <- month
  Filt <- myTable[myTable$time == month,]
  
  Avg_Loss_RYGB[index] <- mean(Filt$Loss_from_BL_BMI[Filt$Surgery=="RYGB"], na.rm = TRUE)
  SD_Loss_RYGB[index] <- sd(Filt$Loss_from_BL_BMI[Filt$Surgery=="RYGB"], na.rm = TRUE)
  
  Avg_Loss_SG[index] <- mean(Filt$Loss_from_BL_BMI[Filt$Surgery=="SG"], na.rm = TRUE)
  SD_Loss_SG[index] <- sd(Filt$Loss_from_BL_BMI[Filt$Surgery=="SG"], na.rm = TRUE)
  
  wilcox.test <- wilcox.test(Filt$Loss_from_BL_BMI ~ Filt$Surgery)
  p_value[index] <- wilcox.test$p.value
  
  index <- index + 1
}

result.df <- data.frame(Timepoint, Avg_Loss_RYGB, SD_Loss_RYGB, Avg_Loss_SG, SD_Loss_SG, p_value)
result.df$BH_p <- p.adjust(result.df$p_value, method = "BH")
result.df$p_value_rounded <- roundP(result.df$BH_p)
result.df$significance <- sigStars(result.df$BH_p)

bio.conc <- vector()
index <- 1

for (i in 1:nrow(result.df)) {
  
  RYGB.mean <- paste0("(", round(result.df$Avg_Loss_RYGB[i], 2), " +/- ", round(result.df$SD_Loss_RYGB[i], 2), " kg/m2)")
  SG.mean <- paste0("(", round(result.df$Avg_Loss_SG[i], 2), " +/- ", round(result.df$SD_Loss_SG[i], 2), " kg/m2)")
  
  if (result.df$p_value[i] < 0.05) {
    bio.conc[index] <- paste0("There was a significant difference in BMI loss between RYGB ", RYGB.mean, " and SG patients ", SG.mean, " at ", m[i], " month(s) (", result.df$p_value_rounded[i], ", Wilcoxon).")
    print(bio.conc[index])
  } else {
    bio.conc[index] <- paste0("There was NO significant difference in BMI loss between RYGB ", RYGB.mean, " and SG patients ", SG.mean, " at ", m[i], " month(s) (", result.df$p_value_rounded[i], ", Wilcoxon).")
    print(bio.conc[index])
  }
  
  
  index <- index + 1
  
}

result.df$Conclusion <- bio.conc
file.path <- paste0(outputDir,"SurgeryType_BMI_Loss_from_BL_wilcox_at_each_timepoint.tsv")
write.table(result.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)

##### BMI differences in location at each timepoint - wilccox #####
months <- c(0, 1, 6, 12, 18, 24)

Timepoint <- vector()
Avg_Fargo <- vector()
SD_Fargo <- vector()
Avg_Cleveland <- vector()
SD_Cleveland <- vector()
p_value <- vector()
index <- 1

for (month in months) {
  
  Timepoint[index] <- month
  Filt <- myTable[myTable$time == month,]
  
  Avg_Fargo[index] <- mean(Filt$BMI_kgm2[Filt$Site=="Fargo"], na.rm = TRUE)
  SD_Fargo[index] <- sd(Filt$BMI_kgm2[Filt$Site=="Fargo"], na.rm = TRUE)
  
  Avg_Cleveland[index] <- mean(Filt$BMI_kgm2[Filt$Site=="Cleveland"], na.rm = TRUE)
  SD_Cleveland[index] <- sd(Filt$BMI_kgm2[Filt$Site=="Cleveland"], na.rm = TRUE)
  
  wilcox.test <- wilcox.test(Filt$BMI_kgm2 ~ Filt$Site)
  p_value[index] <- wilcox.test$p.value
  
  index <- index + 1
}

result.df <- data.frame(Timepoint, Avg_Fargo, SD_Fargo, Avg_Cleveland, SD_Cleveland, p_value)
result.df$BH_p <- p.adjust(result.df$p_value, method = "BH")
result.df$p_value_rounded <- roundP(result.df$BH_p)
result.df$significance <- sigStars(result.df$BH_p)

m <- c(0, 1, 6, 12, 18, 24)
bio.conc <- vector()
index <- 1
for (i in 1:nrow(result.df)) {
  
  Fargo.mean <- paste0("(", round(result.df$Avg_Fargo[i], 2), " +/- ", round(result.df$SD_Fargo[i], 2), " kg/m2)")
  Cleveland.mean <- paste0("(", round(result.df$Avg_Cleveland[i], 2), " +/- ", round(result.df$SD_Cleveland[i], 2), " kg/m2)")
  
  if (result.df$p_value[i] < 0.05) {
    bio.conc[index] <- paste0("There was a significant difference in BMI between Fargo ", Fargo.mean, " and Cleveland patients ", Cleveland.mean, " at ", m[i], " month(s) (", result.df$p_value_rounded[i], ", Wilcoxon).")
    print(bio.conc[index])
  } else {
    bio.conc[index] <- paste0("There was NO significant difference in BMI between Fargo ", Fargo.mean, " and Cleveland patients ", Cleveland.mean, " at ", m[i], " month(s) (", result.df$p_value_rounded[i], ", Wilcoxon).")
    print(bio.conc[index])
  }
  
  
  index <- index + 1
  
}

result.df$Conclusion <- bio.conc
file.path <- paste0(outputDir,"Site_BMI_wilcox_at_each_timepoint.tsv")
write.table(result.df, file.path,sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, "updated_weight.tsv")
write.table(myTable, file.path,sep="\t",quote = FALSE, row.names = FALSE)


##### Surgery BMI_kgm2 over time boxplots #####
metricName <- "BMI_kgm2"
metric <- myTable[,which(colnames(myTable) == metricName)]
y.lab <- "BMI (kg/m^2)"

Surgery <- myTable$Surgery
Timepoint <- factor(myTable$time)
x.lab <- "Time (months) post-surgery"

df <- data.frame(Surgery, metric, Timepoint)

stats.test <- df %>%
  group_by(Surgery) %>%
  tukey_hsd(metric ~ Timepoint)

stats.test$pValue <- roundP(stats.test$p.adj)

stats.test <- stats.test %>%
  add_xy_position(x = "Timepoint")

stats.test2 <- stats.test[stats.test$group1 == "0",]

plot <- ggboxplot(df, x = "Timepoint", y = "metric", 
                  add = "jitter"); plot

plot <- facet(plot, facet.by = "Surgery"); plot

plot <- plot + labs(y = y.lab,
                    x = x.lab)

plotA <- plot +
  stat_pvalue_manual(
    stats.test2,
    bracket.nudge.y = 0.05,
    size = 8,
    hide.ns = TRUE,
    label = "p.adj.signif"
  ); plot

plotB <- plot +
  stat_pvalue_manual(
    stats.test,
    bracket.nudge.y = 0.05,
    size = 8,
    hide.ns = TRUE,
    label = "p.adj.signif"
  ); plot

plotC <- plot +
  stat_pvalue_manual(
    stats.test,
    bracket.nudge.y = 0.09,
    # step.increase = 0.05,
    size = 5,
    hide.ns = TRUE,
    label = "pValue"
  ); plot

plotD <- plot +
  stat_pvalue_manual(
    stats.test2,
    bracket.nudge.y = 0.05,
    # step.increase = 0.05,
    size = 5,
    hide.ns = TRUE,
    label = "pValue"
  ); plot

file.path <- paste0(outputDir, metricName, "_over_Timepoint_Surgery_facet_Tukey_BoxPlot.pdf")
pdf(file.path, width = 10, height = 7)
print(plotA)
print(plotB)
print(plotC)
print(plotD)
dev.off()


##### Surgery Weight_kg over time boxplots #####
metricName <- "Weight_kg"
metric <- myTable[,which(colnames(myTable) == metricName)]
y.lab <- "Weight (kg)"

Surgery <- myTable$Surgery
Timepoint <- factor(myTable$time)
x.lab <- "Time (months) post-surgery"

df <- data.frame(Surgery, metric, Timepoint)

stats.test <- df %>%
  group_by(Surgery) %>%
  tukey_hsd(metric ~ Timepoint)

stats.test$pValue <- roundP(stats.test$p.adj)

stats.test <- stats.test %>%
  add_xy_position(x = "Timepoint")

stats.test2 <- stats.test[stats.test$group1 == "0",]

plot <- ggboxplot(df, x = "Timepoint", y = "metric", 
                  add = "jitter"); plot

plot <- facet(plot, facet.by = "Surgery"); plot

plot <- plot + labs(y = y.lab,
                    x = x.lab)

plotA <- plot +
  stat_pvalue_manual(
    stats.test2,
    bracket.nudge.y = 0.05,
    size = 8,
    hide.ns = TRUE,
    label = "p.adj.signif"
  ); plot

plotB <- plot +
  stat_pvalue_manual(
    stats.test,
    bracket.nudge.y = 0.05,
    size = 8,
    hide.ns = TRUE,
    label = "p.adj.signif"
  ); plot

plotC <- plot +
  stat_pvalue_manual(
    stats.test,
    bracket.nudge.y = 0.09,
    step.increase = 0.05,
    size = 5,
    hide.ns = TRUE,
    label = "pValue"
  ); plot

file.path <- paste0(outputDir, metricName, "_over_Timepoint_Surgery_facet_Tukey_BoxPlot.pdf")
pdf(file.path, width = 10, height = 7)
print(plotA)
print(plotB)
print(plotC)
dev.off()


##### BMI_kgm2 over time boxplots #####
metricName <- "BMI_kgm2"
metric <- myTable[,which(colnames(myTable) == metricName)]
y.lab <- "BMI (kg/m2)"

Surgery <- myTable$Surgery
Timepoint <- factor(myTable$time)
x.lab <- "Time (months) post-surgery"

df <- data.frame(Surgery, metric, Timepoint)

stats.test <- df %>%
  # group_by(Surgery) %>%
  tukey_hsd(metric ~ Timepoint)

stats.test$pValue <- roundP(stats.test$p.adj)

stats.test <- stats.test %>%
  add_xy_position(x = "Timepoint")

stats.test2 <- stats.test[stats.test$group1 == "0",]

plot <- ggboxplot(df, x = "Timepoint", y = "metric", 
                  add = "jitter"); plot

# plot <- facet(plot, facet.by = "Surgery"); plot

plot <- plot + labs(y = y.lab,
                    x = x.lab)

plotA <- plot +
  stat_pvalue_manual(
    stats.test2,
    bracket.nudge.y = 0.05,
    size = 8,
    hide.ns = TRUE,
    label = "p.adj.signif"
  ); plot

plotB <- plot +
  stat_pvalue_manual(
    stats.test,
    bracket.nudge.y = 0.05,
    size = 8,
    hide.ns = TRUE,
    label = "p.adj.signif"
  ); plot

plotC <- plot +
  stat_pvalue_manual(
    stats.test,
    bracket.nudge.y = 0.09,
    step.increase = 0.05,
    size = 5,
    hide.ns = TRUE,
    label = "pValue"
  ); plot

file.path <- paste0(outputDir, metricName, "_over_Timepoint_Tukey_BoxPlot.pdf")
pdf(file.path, width = 10, height = 7)
print(plotA)
print(plotB)
print(plotC)
dev.off()



##### Weight_kg over time boxplots #####
metricName <- "Weight_kg"
metric <- myTable[,which(colnames(myTable) == metricName)]
y.lab <- "Weight (kg)"

Surgery <- myTable$Surgery
Timepoint <- factor(myTable$time)
x.lab <- "Time (months) post-surgery"

df <- data.frame(Surgery, metric, Timepoint)

stats.test <- df %>%
  # group_by(Surgery) %>%
  tukey_hsd(metric ~ Timepoint)

stats.test$pValue <- roundP(stats.test$p.adj)

stats.test <- stats.test %>%
  add_xy_position(x = "Timepoint")

stats.test2 <- stats.test[stats.test$group1 == "0",]

plot <- ggboxplot(df, x = "Timepoint", y = "metric", 
                  add = "jitter"); plot

# plot <- facet(plot, facet.by = "Surgery"); plot

plot <- plot + labs(y = y.lab,
                    x = x.lab)

plotA <- plot +
  stat_pvalue_manual(
    stats.test2,
    bracket.nudge.y = 0.05,
    size = 8,
    hide.ns = TRUE,
    label = "p.adj.signif"
  ); plot

plotB <- plot +
  stat_pvalue_manual(
    stats.test,
    bracket.nudge.y = 0.05,
    size = 8,
    hide.ns = TRUE,
    label = "p.adj.signif"
  ); plot

plotC <- plot +
  stat_pvalue_manual(
    stats.test,
    bracket.nudge.y = 0.09,
    step.increase = 0.05,
    size = 5,
    hide.ns = TRUE,
    label = "pValue"
  ); plot

file.path <- paste0(outputDir, metricName, "_over_Timepoint_Tukey_BoxPlot.pdf")
pdf(file.path, width = 10, height = 7)
print(plotA)
print(plotB)
print(plotC)
dev.off()


