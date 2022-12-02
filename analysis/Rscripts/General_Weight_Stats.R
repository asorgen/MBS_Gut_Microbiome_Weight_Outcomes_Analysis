#Author: Alicia Sorgen
#Date: 2022 June 13
#Description: General weight comparison stats

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
params <- c(params, 0)
params <- c(params, 1)
params <- c(params, 6)
params <- c(params, 12)
params <- c(params, 18)
params <- c(params, 24)
# end_month <- params[length(params)]

moduleRoot <- "General_Weight_Stats"


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
  # file.remove(files)
  
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
prevModule <- "ExcessWeightLoss"
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
inputFile <- "updated_weight2.tsv"
included <- args[2:length(args)]

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Read in tables #####
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

# Order Timepoint as factor
myTable$Timepoint <- as.factor(myTable$Timepoint)
myTable$Timepoint <- factor(myTable$Timepoint, levels = c("BL", "ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR"))
# length(unique(myTable$PatientID))

weightMetrics <- c("Weight_kg", "Loss_from_BL_kg", "Loss_from_BL_BMI", "BMI_kgm2", "Excess_kg", "PEWL_kg", "Percent_Loss_kg")
weightLabels <- c("Weight (kg)", "Weight loss (kg)", "BMI loss", "BMI", "Excess weight (kg)", "Excess weight loss (%)", "Weight loss (%)")

##### Sex by Weight Metric - Wilcoxon (grouped by Timepoint) #####
stats_RAW <- data.frame()
averages_RAW <- data.frame()
finalTable <- data.frame()

Timepoint <- myTable$Timepoint
Sex <- myTable$Sex
Sex <- gsub(pattern = "F", replacement = "Female", Sex)
Sex <- gsub(pattern = "M", replacement = "Male", Sex)

for (i in weightMetrics) {
  
  Weight <- myTable[,which(colnames(myTable) == i)]
  
  df <- data.frame(Timepoint, Sex, Weight)
  
  if (i %in% c("Loss_from_BL_kg", "Loss_from_BL_BMI", "Excess_Loss_from_BL_kg", "Excess_Loss_from_BL_BMI", "PEWL_kg", "Percent_Loss_kg")) {
    df <- df[!(df$Timepoint == "BL"),]
  }
  
  df <- na.omit(df)
  
  stats <- df %>%
    group_by(Timepoint) %>%
    wilcox_test(Weight ~ Sex) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  stats$.y. <- i
  
  stats_RAW <- rbind(stats_RAW, stats)
  
  averages <- df %>%
    group_by(Timepoint, Sex) %>%
    get_summary_stats(Weight, type = "mean_sd")
  
  averages$variable <- i
  
  averages_RAW <- rbind(averages_RAW, averages)
  
  averages$Mean_SD <- paste(round(averages$mean, 2), "±", round(averages$sd, 2))
  averages <- averages[,-(3:6)]
  average2 <- spread(averages, Sex, Mean_SD)
  average3 <- average2[,-1]
  summary <- cbind(stats, average3)

  finalTable <- rbind(finalTable, summary)
  
}

finalTable$p_value <- roundP(finalTable$p.adj)
finalTable$Timepoint <- gsub(pattern = "BL", replacement = "Baseline", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "ONE", replacement = "Postop 1 month", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "SIX", replacement = "Postop 6 months", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "TWELVE", replacement = "Postop 12 months", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "EIGHTEEN", replacement = "Postop 18 months", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "TWENTY_FOUR", replacement = "Postop 24 months", finalTable$Timepoint)
finalTable <- finalTable[,-(3:10)]

file.path <- paste0(outputDir, "Sex_raw_wilcoxon_results.tsv")
write.table(stats_RAW, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, "Sex_raw_averages.tsv")
write.table(averages_RAW, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, "Sex_final_wilcoxon_summary.tsv")
write.table(finalTable, file.path, sep="\t",quote = FALSE, row.names = FALSE)


##### Ethnicity by Weight Metric - Kruskal-Wallis (grouped by Timepoint) #####
stats_RAW <- data.frame()
averages_RAW <- data.frame()
finalTable <- data.frame()

Timepoint <- myTable$Timepoint
Ethnicity <- myTable$Ethnicity

for (i in weightMetrics) {
  
  Weight <- myTable[,which(colnames(myTable) == i)]
  
  df <- data.frame(Timepoint, Ethnicity, Weight)
  
  if (i %in% c("Loss_from_BL_kg", "Loss_from_BL_BMI", "Excess_Loss_from_BL_kg", "Excess_Loss_from_BL_BMI", "PEWL_kg", "Percent_Loss_kg")) {
    df <- df[!(df$Timepoint == "BL"),]
  }
  
  df <- na.omit(df)
  
  stats <- df %>%
    group_by(Timepoint) %>%
    kruskal_test(Weight ~ Ethnicity) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  stats$.y. <- i
  
  stats_RAW <- rbind(stats_RAW, stats)
  
  averages <- df %>%
    group_by(Timepoint, Ethnicity) %>%
    get_summary_stats(Weight, type = "mean_sd")
  
  averages$variable <- i
  
  averages_RAW <- rbind(averages_RAW, averages)
  
  averages$Mean_SD <- paste(round(averages$mean, 2), "±", round(averages$sd, 2))
  averages <- averages[,-(3:6)]
  average2 <- spread(averages, Ethnicity, Mean_SD)
  average3 <- average2[,-1]
  summary <- cbind(stats, average3)
  
  finalTable <- rbind(finalTable, summary)
  
}

finalTable$p_value <- roundP(finalTable$p.adj)
finalTable$Timepoint <- gsub(pattern = "BL", replacement = "Baseline", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "ONE", replacement = "Postop 1 month", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "SIX", replacement = "Postop 6 months", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "TWELVE", replacement = "Postop 12 months", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "EIGHTEEN", replacement = "Postop 18 months", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "TWENTY_FOUR", replacement = "Postop 24 months", finalTable$Timepoint)
finalTable <- finalTable[,-(3:9)]

file.path <- paste0(outputDir, "Ethnicity_raw_kruskal_results.tsv")
write.table(stats_RAW, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, "Ethnicity_raw_averages.tsv")
write.table(averages_RAW, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, "Ethnicity_final_kruskal_summary.tsv")
write.table(finalTable, file.path, sep="\t",quote = FALSE, row.names = FALSE)



##### Surgery by Weight Metric - Wilcoxon (grouped by Timepoint) #####
stats_RAW <- data.frame()
averages_RAW <- data.frame()
finalTable <- data.frame()

Timepoint <- myTable$Timepoint
Surgery <- myTable$Surgery

for (i in weightMetrics) {
  
  Weight <- myTable[,which(colnames(myTable) == i)]
  
  df <- data.frame(Timepoint, Surgery, Weight)
  
  if (i %in% c("Loss_from_BL_kg", "Loss_from_BL_BMI", "Excess_Loss_from_BL_kg", "Excess_Loss_from_BL_BMI", "PEWL_kg", "Percent_Loss_kg")) {
    df <- df[!(df$Timepoint == "BL"),]
  }
  
  df <- na.omit(df)
  
  stats <- df %>%
    group_by(Timepoint) %>%
    wilcox_test(Weight ~ Surgery) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  stats$.y. <- i
  
  stats_RAW <- rbind(stats_RAW, stats)
  
  averages <- df %>%
    group_by(Timepoint, Surgery) %>%
    get_summary_stats(Weight, type = "mean_sd")
  
  averages$variable <- i
  
  averages_RAW <- rbind(averages_RAW, averages)
  
  averages$Mean_SD <- paste(round(averages$mean, 2), "±", round(averages$sd, 2))
  averages <- averages[,-(3:6)]
  average2 <- spread(averages, Surgery, Mean_SD)
  average3 <- average2[,-1]
  summary <- cbind(stats, average3)
  
  finalTable <- rbind(finalTable, summary)
  
}

finalTable$p_value <- roundP(finalTable$p.adj)
finalTable$Timepoint <- gsub(pattern = "BL", replacement = "Baseline", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "ONE", replacement = "Postop 1 month", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "SIX", replacement = "Postop 6 months", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "TWELVE", replacement = "Postop 12 months", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "EIGHTEEN", replacement = "Postop 18 months", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "TWENTY_FOUR", replacement = "Postop 24 months", finalTable$Timepoint)
finalTable <- finalTable[,-(3:10)]

file.path <- paste0(outputDir, "Surgery_raw_wilcoxon_results.tsv")
write.table(stats_RAW, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, "Surgery_raw_averages.tsv")
write.table(averages_RAW, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, "Surgery_final_wilcoxon_summary.tsv")
write.table(finalTable, file.path, sep="\t",quote = FALSE, row.names = FALSE)



##### Site by Weight Metric - Wilcoxon (grouped by Timepoint) #####
stats_RAW <- data.frame()
averages_RAW <- data.frame()
finalTable <- data.frame()

Timepoint <- myTable$Timepoint
Site <- myTable$Site

for (i in weightMetrics) {
  
  Weight <- myTable[,which(colnames(myTable) == i)]
  
  df <- data.frame(Timepoint, Site, Weight)
  
  if (i %in% c("Loss_from_BL_kg", "Loss_from_BL_BMI", "Excess_Loss_from_BL_kg", "Excess_Loss_from_BL_BMI", "PEWL_kg", "Percent_Loss_kg")) {
    df <- df[!(df$Timepoint == "BL"),]
  }
  
  df <- na.omit(df)
  
  stats <- df %>%
    group_by(Timepoint) %>%
    wilcox_test(Weight ~ Site) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  stats$.y. <- i
  
  stats_RAW <- rbind(stats_RAW, stats)
  
  averages <- df %>%
    group_by(Timepoint, Site) %>%
    get_summary_stats(Weight, type = "mean_sd")
  
  averages$variable <- i
  
  averages_RAW <- rbind(averages_RAW, averages)
  
  averages$Mean_SD <- paste(round(averages$mean, 2), "±", round(averages$sd, 2))
  averages <- averages[,-(3:6)]
  average2 <- spread(averages, Site, Mean_SD)
  average3 <- average2[,-1]
  summary <- cbind(stats, average3)
  
  finalTable <- rbind(finalTable, summary)
  
}

finalTable$p_value <- roundP(finalTable$p.adj)
finalTable$Timepoint <- gsub(pattern = "BL", replacement = "Baseline", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "ONE", replacement = "Postop 1 month", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "SIX", replacement = "Postop 6 months", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "TWELVE", replacement = "Postop 12 months", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "EIGHTEEN", replacement = "Postop 18 months", finalTable$Timepoint)
finalTable$Timepoint <- gsub(pattern = "TWENTY_FOUR", replacement = "Postop 24 months", finalTable$Timepoint)
finalTable <- finalTable[,-(3:10)]

file.path <- paste0(outputDir, "Site_raw_wilcoxon_results.tsv")
write.table(stats_RAW, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, "Site_raw_averages.tsv")
write.table(averages_RAW, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, "Site_final_wilcoxon_summary.tsv")
write.table(finalTable, file.path, sep="\t",quote = FALSE, row.names = FALSE)


##### Weight Metric by Timepoint - Mixed linear model #####
PatientID <- myTable$PatientID
MLM.df <- data.frame()
Tukey.df <- data.frame()

for (i in weightMetrics) {
  
  Weight <- myTable[,which(colnames(myTable) == i)]
  Timepoint <- myTable$Timepoint
  
  df <- data.frame(PatientID, Timepoint, Weight)
  
  if (!(i %in% c("Loss_from_BL_kg", "Loss_from_BL_BMI", "Excess_Loss_from_BL_kg", "Excess_Loss_from_BL_BMI", "PEWL_kg", "Percent_Loss_kg"))) {
    
    df <- na.omit(df)
    
    mlm <- lme(Weight ~ Timepoint, method = 'REML', random = ~1 | PatientID, data = df)
    smry <- summary(mlm)
    lm.df <- data.frame(smry$tTable)
    lm.df$Comparison <- row.names(lm.df)
    lm.df$Comparison <- gsub(pattern = "Timepoint", replacement = "BL v ", lm.df$Comparison)
    lm.df$Metric <- i
    lm.df$Adj_pvalue <- p.adjust(lm.df$p.value, method = "BH")
    lm.df$p_value <- roundP(lm.df$Adj_pvalue)
    
    anova <- aov(mlm)
    tukey <- TukeyHSD(anova, conf.level = .95)
    tukey.df <- as.data.frame(tukey$Timepoint)
    tukey.df$Comparison <- row.names(tukey.df)
    tukey.df$Metric <- i
    tukey.df$p_value <- roundP(tukey.df$`p adj`)
    
    MLM.df <- rbind(MLM.df, lm.df)
    Tukey.df <- rbind(Tukey.df, tukey.df)
    
  }
  
  

}


file.path <- paste0(outputDir, "Weight_Metric_by_Timepoint_MLM.tsv")
write.table(MLM.df, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, "Weight_Metric_by_Timepoint_Tukey.tsv")
write.table(Tukey.df, file.path, sep="\t",quote = FALSE, row.names = FALSE)


##### Box plots for weight metrics over time #####
plotList <- list()
index <- 1

for (i in weightMetrics) {
  
  label <- weightLabels[which(weightMetrics == i)]
  
  Weight <- myTable[,which(colnames(myTable) == i)]
  
  df <- data.frame(PatientID, Timepoint, Weight)
  df <- na.omit(df)
  
  if ((i %in% c("Loss_from_BL_kg", "Loss_from_BL_BMI", "Excess_Loss_from_BL_kg", "Excess_Loss_from_BL_BMI", "PEWL_kg", "Percent_Loss_kg"))) {
    
    df <- df[!(df$Timepoint == "BL"),]
    
    df$Timepoint <- gsub(pattern = "ONE", replacement = "BL - 1M", df$Timepoint)
    df$Timepoint <- gsub(pattern = "SIX", replacement = "BL - 6M", df$Timepoint)
    df$Timepoint <- gsub(pattern = "TWELVE", replacement = "BL - 12M", df$Timepoint)
    df$Timepoint <- gsub(pattern = "EIGHTEEN", replacement = "BL - 18M", df$Timepoint)
    df$Timepoint <- gsub(pattern = "TWENTY_FOUR", replacement = "BL - 24M", df$Timepoint)
    
    df$Timepoint <- as.factor(df$Timepoint)
    df$Timepoint <- factor( df$Timepoint, levels = c("BL - 1M", "BL - 6M", "BL - 12M", "BL - 18M", "BL - 24M"))
    plot <- ggboxplot(
      df, x = "Timepoint", y = "Weight",
      scales = "free", add = "jitter"
    ); plot
    
    plot <- plot +
      labs(y = label); plot
    
  } else {
    
    df$Timepoint <- gsub(pattern = "BL", replacement = "Baseline", df$Timepoint)
    df$Timepoint <- gsub(pattern = "ONE", replacement = "1 month", df$Timepoint)
    df$Timepoint <- gsub(pattern = "SIX", replacement = "6 months", df$Timepoint)
    df$Timepoint <- gsub(pattern = "TWELVE", replacement = "12 months", df$Timepoint)
    df$Timepoint <- gsub(pattern = "EIGHTEEN", replacement = "18 months", df$Timepoint)
    df$Timepoint <- gsub(pattern = "TWENTY_FOUR", replacement = "24 months", df$Timepoint)
    
    df$Timepoint <- as.factor(df$Timepoint)
    df$Timepoint <- factor( df$Timepoint, levels = c("Baseline", "1 month", "6 months", "12 months", "18 months", "24 months"))
    plot <- ggboxplot(
      df, x = "Timepoint", y = "Weight",
      scales = "free", add = "jitter"
    ); plot
    
    plot <- plot +
      labs(y = label); plot
    
  }
  
  plotList[[index]] <- plot
  index <- index + 1
  
}

# Output plot
file.path <- paste0(outputDir, "Weight_Metric_Boxplots.pdf")
pdf(file.path, width = 7, height = 5)
for (i in 1:length(plotList)) {
  grid.arrange(plotList[[i]], ncol = 1, nrow = 1)
}
dev.off()




