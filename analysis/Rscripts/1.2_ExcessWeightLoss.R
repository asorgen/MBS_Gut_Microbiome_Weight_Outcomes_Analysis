#Author: Alicia Sorgen
#Date: 2022 June 13
#Description: Linear Modeling of ASA24 Metadata

##### Edits for script #####
rm(list=ls())
params <- vector()

ANALYSIS <- "microbiome_n124"
params <- c(params, "~/UNCC/Projects/Bariatric_Surgery/Git_Repositories/gut-microbiome-bariatric-weight-outcomes")

moduleRoot <- "1.2_ExcessWeightLoss"


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
  gitRoot    <- args[1]
  gitInput   <- file.path(gitRoot, "..", "Data")
  gitScripts <- file.path(gitRoot, "analysis", "Rscripts")

  root <- file.path(gitRoot, "Results")
  dir.create(root, showWarnings = FALSE, recursive = TRUE)

  if (!dir.exists(file.path(root, "input"))) {
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
prevModule <- "BMI_Results"
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
inputFile <- "updated_weight.tsv"

##### Set up output #####
outputDir = file.path(moduleDir,"output/")
statsOutput <- "EWL_Summary_table.tsv"
tableupdate <- "updated_weight2.tsv"
summaryFile0 <- "Responder_Summary.tsv"
summaryFile1 <- "Responder_Summary_by_SurgeryType.tsv"
summaryFile <- "Avg_PEWL_for_Responders_at_each_timepoint.tsv"
summaryFile2 <- "Avg_PEWL_at_each_timepoint.tsv"
ttestFile <- "PEWL_for_Responders_at_each_timepoint_ttest.tsv"
wilcoxFile <- "PEWL_for_Responders_at_each_timepoint_wilcox.tsv"
ttestFile_surgery <- "PEWL_for_SurgeryType_at_each_timepoint_ttest.tsv"
wilcoxFile_surgery <- "PEWL_for_SurgeryType_at_each_timepoint_wilcox.tsv"
summaryFile_surgery <- "Avg_PEWL_for_SurgeryType_at_each_timepoint.tsv"
summaryFile_all <- "Avg_PEWL_for_Combined_Responders_at_each_timepoint.tsv"
weightmetricFile <- "SurgeryType_t_wilcox_weight_metric_results.tsv"
surgTypeAvgFile <- "Avg_weight_metric_by_SurgeryType.tsv"

##### Read in tables #####
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

# Order Timepoint as factor
myTable$Timepoint <- as.factor(myTable$Timepoint)
myTable$Timepoint <- factor(myTable$Timepoint, levels = c("BL", "ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR"))
# Assign SampleID to as dataframe row names
row.names(myTable) = myTable$SampleID

# Convert the timepoint data to numeric values.
myTable$time = as.numeric(myTable$time)

# Baseline rows
bl = which(myTable$Timepoint=="BL")

# Baseline weight (kg)
BL_kg = myTable[bl, "Weight_kg"]
names(BL_kg) = myTable[bl, "PatientID"]

# Add info to each row
myTable$Baseline_kg = BL_kg[myTable$PatientID]

##### Calculate average loss relative to baseline.   #####
myTable$Loss_per_Month_kg = (myTable$Baseline_kg - myTable$Weight_kg) / myTable$time
myTable$Loss_per_Month_BMI = (myTable$Baseline_BMI - myTable$BMI_kgm2) / myTable$time


myTable$Loss_from_BL_kg = (myTable$Baseline_kg - myTable$Weight_kg)
summary.df <-myTable %>%
  group_by(Timepoint) %>%
  get_summary_stats(Loss_from_BL_kg, type = "mean_sd")
summary.df
final <- summary.df

myTable$Percent_Loss_kg <- (myTable$Loss_from_BL_kg / myTable$Baseline_kg) * 100
summary.df <-myTable %>%
  group_by(Timepoint) %>%
  get_summary_stats(Percent_Loss_kg, type = "mean_sd")
final <- rbind(final, summary.df)

summary.df <-myTable %>%
  group_by(Timepoint) %>%
  get_summary_stats(Loss_from_BL_BMI, type = "mean_sd")
final <- rbind(final, summary.df)

summary.df <-myTable %>%
  group_by(Timepoint) %>%
  get_summary_stats(Weight_kg, type = "mean_sd")
final <- rbind(final, summary.df)

summary.df <-myTable %>%
  group_by(Timepoint) %>%
  get_summary_stats(BMI_kgm2, type = "mean_sd")
final <- rbind(final, summary.df)

##### Calculate excess weight (kg) for each patient #####
myTable$Ideal_kg <- myTable$Baseline_m^2 * 25
myTable$Excess_kg <- myTable$Weight_kg - myTable$Ideal_kg
myTable$Excess_BMI <- myTable$BMI_kgm2 - 25
summary.df <-myTable %>%
  group_by(Timepoint) %>%
  get_summary_stats(Excess_kg, type = "mean_sd")
summary.df
final <- rbind(final, summary.df)

summary.df <-myTable %>%
  group_by(Timepoint) %>%
  get_summary_stats(Excess_BMI, type = "mean_sd")
final <- rbind(final, summary.df)

##### Calculate % excess weight (kg) for each patient #####
myTable$Percent_Excess_kg <- (myTable$Excess_kg / myTable$Weight_kg) * 100
summary.df <-myTable %>%
  group_by(Timepoint) %>%
  get_summary_stats(Percent_Excess_kg, type = "mean_sd")
summary.df
final <- rbind(final, summary.df)

##### Calculate excess weight loss from baseline #####
# Baseline excess weight (kg)
BL_Excess_kg = myTable[bl, "Excess_kg"]
names(BL_Excess_kg) = myTable[bl, "PatientID"]

# Add info to each row
myTable$Baseline_Excess_kg = BL_Excess_kg[myTable$PatientID]

# Baseline excess weight (BMI)
BL_Excess_BMI = myTable[bl, "Excess_BMI"]
names(BL_Excess_BMI) = myTable[bl, "PatientID"]

# Add info to each row
myTable$Baseline_Excess_BMI = BL_Excess_BMI[myTable$PatientID]

myTable$BL_Ideal_Diff_kg <- myTable$Baseline_kg - myTable$Ideal_kg
myTable$BL_Ideal_Diff_BMI <- myTable$Baseline_BMI - 25

myTable$Excess_Loss_from_BL_kg <- myTable$Baseline_Excess_kg - myTable$Excess_kg
summary.df <-myTable %>%
  group_by(Timepoint) %>%
  get_summary_stats(Excess_Loss_from_BL_kg, type = "mean_sd")
summary.df
final <- rbind(final, summary.df)

myTable$Excess_Loss_from_BL_BMI <- myTable$Baseline_Excess_BMI - myTable$Excess_BMI
summary.df <-myTable %>%
  group_by(Timepoint) %>%
  get_summary_stats(Excess_Loss_from_BL_BMI, type = "mean_sd")
summary.df
final <- rbind(final, summary.df)

##### Calculate percent excess weight loss (PEWL) #####
myTable$PEWL_kg <- (myTable$Loss_from_BL_kg / myTable$BL_Ideal_Diff_kg) * 100
summary.df <-myTable %>%
  group_by(Timepoint) %>%
  get_summary_stats(PEWL_kg, type = "mean_sd")
summary.df
final <- rbind(final, summary.df)

myTable$ResponderStatus <- ifelse(myTable$PEWL_kg >= 50, paste0(myTable$time, "-month responder"),
                                  ifelse(myTable$PEWL_kg < 50, paste0(myTable$time, "-month non-responder"), 
                                         NA))

smmry  <- as.data.frame(table(myTable$ResponderStatus))
response.summary <- as.data.frame(table(myTable$ResponderStatus, myTable$Surgery))


file.path <- paste0(outputDir, summaryFile0)
write.table(smmry, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, summaryFile1)
write.table(response.summary, file.path, sep="\t",quote = FALSE, row.names = FALSE)

final$final <- paste(round(final$mean, 2), "±", round(final$sd, 2))
file.path <- paste0(outputDir, statsOutput)
write.table(final, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, tableupdate)
write.table(myTable, file.path, sep="\t",quote = FALSE, row.names = FALSE)

##### Calculate average EWL between responders and non-responders at each timepoint #####
avg2 <- myTable %>%
  group_by(Timepoint) %>%
  get_summary_stats(PEWL_kg, type = "mean_sd") 

ttest.df <- data.frame()
wilcox.df <- data.frame()
summary <- data.frame()
summary_all <- data.frame()

months <- c("TWELVE", "EIGHTEEN", "TWENTY_FOUR")
included <- c("ONE", "SIX")

for (month in months) {
  
  PatientID <- myTable$PatientID
  Timepoint <- myTable$Timepoint
  PEWL_kg <- myTable$PEWL_kg
  ResponderStatus <- myTable$ResponderStatus
  
  df <- data.frame(PatientID, Timepoint, PEWL_kg)
  df <- na.omit(df)
  df <- spread(df, Timepoint, PEWL_kg)
  m <- which(colnames(df) == month)
  
  df$ResponderStatus <- ifelse(df[,m] >= 50, "Responder",
                               ifelse(df[,m] < 50, "Non-responder", 
                                      NA))
  df2 <- df[!is.na(df$ResponderStatus),]
  df2 <- df2 %>%
    gather("Timepoint", "PEWL_kg", 2:7)
  
  included <- c(included, month)
  
  df2 <- df2[df2$Timepoint %in% included,]
  df2$Timepoint <- as.factor(df2$Timepoint)
  df2$Timepoint <- factor(df2$Timepoint, levels = included)
  
  
  
  # t-test
  stats.summ <- df2 %>%
    group_by(Timepoint) %>%
    t_test(PEWL_kg ~ ResponderStatus) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  stats.summ$p_value <- roundP(stats.summ$p.adj)
  
  stats.summ$ResponderMonth <- month
  
  ttest.df <- rbind(ttest.df, stats.summ)
  
  
  
  # Wilcoxon
  stats.summ <- df2 %>%
    group_by(Timepoint) %>%
    wilcox_test(PEWL_kg ~ ResponderStatus) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  stats.summ$p_value <- roundP(stats.summ$p.adj)
  
  stats.summ$ResponderMonth <- month
  
  wilcox.df <- rbind(wilcox.df, stats.summ)
  
  
  averages <- df2 %>%
    group_by(Timepoint, ResponderStatus) %>%
    get_summary_stats(PEWL_kg, type = "mean_sd")
  
  averages$ResponderMonth <- month
  
  summary <- rbind(summary, averages)
  
  averages_all <- df2 %>%
    group_by(Timepoint) %>%
    get_summary_stats(PEWL_kg, type = "mean_sd")
  
  averages_all$ResponderMonth <- month
  
  summary_all <- rbind(summary_all, averages_all)
  
  
}

summary$final <- paste(round(summary$mean, 2), "±", round(summary$sd, 2))
summary_all$final <- paste(round(summary_all$mean, 2), "±", round(summary_all$sd, 2))

file.path <- paste0(outputDir, summaryFile2)
write.table(avg2, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, summaryFile_all)
write.table(summary_all, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, summaryFile)
write.table(summary, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, ttestFile)
write.table(ttest.df, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, wilcoxFile)
write.table(wilcox.df, file.path, sep="\t",quote = FALSE, row.names = FALSE)

variables <- unique(final$variable)
results_combined <- data.frame()
surgtype_combined <- data.frame()

for (v1 in variables) {
  
  Timepoint <- myTable$Timepoint
  Surgery <- myTable$Surgery
  Variable <- myTable[,v1]
  
  df <- data.frame(Timepoint, Surgery, Variable)
  
  changeMetrics <- c("Loss_from_BL_BMI", "Loss_from_BL_kg", "Percent_Loss_kg", "Excess_Loss_from_BL_kg", "Excess_Loss_from_BL_BMI", "PEWL_kg")
  
  if (v1 %in% changeMetrics) {
    df2 <- df[!(df$Timepoint == "BL"),]
  } else {
    df2 <- df
  }
  
  df2 <- na.omit(df2)
  
  # t-test
  results <- df2 %>%
    group_by(Timepoint) %>%
    t_test(Variable ~ Surgery) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  results$.y. <- v1
  results$Model <- "t-test"
  results <- results[,-(which(colnames(results) == "df"))]
  
  results_combined <- rbind(results_combined, results)
  
  
  # Wilcoxon
  results <- df2 %>%
    group_by(Timepoint) %>%
    wilcox_test(Variable ~ Surgery) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  results$.y. <- v1
  results$Model <- "Wilcoxon"
  results_combined <- rbind(results_combined, results)
  
  SurgeryType_averages <- df2 %>%
    group_by(Timepoint, Surgery) %>%
    get_summary_stats(Variable, type = "mean_sd")
  
  SurgeryType_averages$variable <- v1
  
  surgtype_combined <- rbind(surgtype_combined, SurgeryType_averages)
  
}

results_combined$p_final <- roundP(results_combined$p.adj)

file.path <- paste0(outputDir, weightmetricFile)
write.table(results_combined, file.path, sep="\t",quote = FALSE, row.names = FALSE)


surgtype_combined$final <- paste0(round(surgtype_combined$mean, 2), " ± ", round(surgtype_combined$sd, 2))

file.path <- paste0(outputDir, surgTypeAvgFile)
write.table(surgtype_combined, file.path, sep="\t",quote = FALSE, row.names = FALSE)

ttest.df <- data.frame()
wilcox.df <- data.frame()
summary <- data.frame()

months <- c("TWELVE", "EIGHTEEN", "TWENTY_FOUR")
included <- c("ONE", "SIX")

for (month in months) {
  
  PatientID <- myTable$PatientID
  Timepoint <- myTable$Timepoint
  PEWL_kg <- myTable$PEWL_kg
  Surgery <- myTable$Surgery
  
  df <- data.frame(PatientID, Timepoint, PEWL_kg, Surgery)
  # df <- na.omit(df)
  df <- spread(df, Timepoint, PEWL_kg)
  m <- which(colnames(df) == month)
  
  df2 <- df[!is.na(df$Surgery),]
  df2 <- df2 %>%
    gather("Timepoint", "PEWL_kg", 3:8)
  
  included <- c(included, month)
  
  df2 <- df2[df2$Timepoint %in% included,]
  df2$Timepoint <- as.factor(df2$Timepoint)
  df2$Timepoint <- factor(df2$Timepoint, levels = included)
  
  
  
  # t-test
  stats.summ <- df2 %>%
    group_by(Timepoint) %>%
    t_test(PEWL_kg ~ Surgery) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  stats.summ$p_value <- roundP(stats.summ$p.adj)
  
  stats.summ$ResponderMonth <- month
  
  ttest.df <- rbind(ttest.df, stats.summ)
  
  
  
  # Wilcoxon
  stats.summ <- df2 %>%
    group_by(Timepoint) %>%
    wilcox_test(PEWL_kg ~ Surgery) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  stats.summ$p_value <- roundP(stats.summ$p.adj)
  
  stats.summ$ResponderMonth <- month
  
  wilcox.df <- rbind(wilcox.df, stats.summ)
  
  
  averages <- df2 %>%
    group_by(Timepoint, Surgery) %>%
    get_summary_stats(PEWL_kg, type = "mean_sd")
  
  averages$ResponderMonth <- month
  
  summary <- rbind(summary, averages)
  
}

summary$final <- paste0(round(summary$mean, 2), " ± ", round(summary$sd, 2))



file.path <- paste0(outputDir, summaryFile_surgery)
write.table(summary, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, ttestFile_surgery)
write.table(ttest.df, file.path, sep="\t",quote = FALSE, row.names = FALSE)

file.path <- paste0(outputDir, wilcoxFile_surgery)
write.table(wilcox.df, file.path, sep="\t",quote = FALSE, row.names = FALSE)


##### Create summary table separated by surgery type #####
summary.df <-myTable %>%
  group_by(Timepoint, Surgery) %>%
  get_summary_stats(Loss_from_BL_kg, type = "mean_sd")
summary.df$Formatted <- paste0(round(summary.df$mean, 2), " ± ", round(summary.df$sd, 2))
final <- summary.df

summary.df <-myTable %>%
  group_by(Timepoint, Surgery) %>%
  get_summary_stats(Percent_Loss_kg, type = "mean_sd")
summary.df$Formatted <- paste0(round(summary.df$mean, 2), " ± ", round(summary.df$sd, 2))
final <- rbind(final, summary.df)

summary.df <-myTable %>%
  group_by(Timepoint, Surgery) %>%
  get_summary_stats(Loss_from_BL_BMI, type = "mean_sd")
summary.df$Formatted <- paste0(round(summary.df$mean, 2), " ± ", round(summary.df$sd, 2))
final <- rbind(final, summary.df)

summary.df <-myTable %>%
  group_by(Timepoint, Surgery) %>%
  get_summary_stats(Weight_kg, type = "mean_sd")
summary.df$Formatted <- paste0(round(summary.df$mean, 2), " ± ", round(summary.df$sd, 2))
final <- rbind(final, summary.df)

summary.df <-myTable %>%
  group_by(Timepoint, Surgery) %>%
  get_summary_stats(BMI_kgm2, type = "mean_sd")
summary.df$Formatted <- paste0(round(summary.df$mean, 2), " ± ", round(summary.df$sd, 2))
final <- rbind(final, summary.df)

summary.df <-myTable %>%
  group_by(Timepoint, Surgery) %>%
  get_summary_stats(Excess_kg, type = "mean_sd")
summary.df$Formatted <- paste0(round(summary.df$mean, 2), " ± ", round(summary.df$sd, 2))
final <- rbind(final, summary.df)

summary.df <-myTable %>%
  group_by(Timepoint, Surgery) %>%
  get_summary_stats(BL_Ideal_Diff_kg, type = "mean_sd")
summary.df$Formatted <- paste0(round(summary.df$mean, 2), " ± ", round(summary.df$sd, 2))
final <- rbind(final, summary.df)

summary.df <-myTable %>%
  group_by(Timepoint, Surgery) %>%
  get_summary_stats(BL_Ideal_Diff_BMI, type = "mean_sd")
summary.df$Formatted <- paste0(round(summary.df$mean, 2), " ± ", round(summary.df$sd, 2))
final <- rbind(final, summary.df)

summary.df <-myTable %>%
  group_by(Timepoint, Surgery) %>%
  get_summary_stats(PEWL_kg, type = "mean_sd")
summary.df$Formatted <- paste0(round(summary.df$mean, 2), " ± ", round(summary.df$sd, 2))
final <- rbind(final, summary.df)

Timepoint <- vector()
Surgery_type  <- vector()
n <- vector()
Weight_kg <- vector()
BMI_kgm2 <- vector()
Excess_kg <- vector()
Loss_from_BL_kg <- vector()
Loss_from_BL_BMI <- vector()
Percent_Loss_kg <- vector()
PEWL_kg <- vector()
index <- 1


for (var in unique(final$Timepoint)) {
  
  df2 <- final[final$Timepoint == var,]
  
  for (i in 1:2) {
    Timepoint[index] <- var
    
    df3 <- df2[df2$variable == "Weight_kg",]
    Weight_kg[index] <- df3$Formatted[i]
    
    df3 <- df2[df2$variable == "Weight_kg",]
    Weight_kg[index] <- df3$Formatted[i]
    
    df3 <- df2[df2$variable == "BMI_kgm2",]
    BMI_kgm2[index] <- df3$Formatted[i]
    
    df3 <- df2[df2$variable == "Excess_kg",]
    Excess_kg[index] <- df3$Formatted[i]
    
    df3 <- df2[df2$variable == "Loss_from_BL_kg",]
    Loss_from_BL_kg[index] <- df3$Formatted[i]
    
    df3 <- df2[df2$variable == "Loss_from_BL_BMI",]
    Loss_from_BL_BMI[index] <- df3$Formatted[i]
    
    df3 <- df2[df2$variable == "Percent_Loss_kg",]
    Percent_Loss_kg[index] <- df3$Formatted[i]
    
    df3 <- df2[df2$variable == "PEWL_kg",]
    PEWL_kg[index] <- df3$Formatted[i]
    
    Surgery_type[index] <- df3$Surgery[i]
    n[index] <- df3$n[i]
    
    index <- index + 1
    
  }
}

dFrame <- data.frame(Timepoint, Surgery_type, n, Weight_kg, BMI_kgm2, Excess_kg, Loss_from_BL_kg, Loss_from_BL_BMI, Percent_Loss_kg, PEWL_kg)

dFrame$Timepoint <- gsub(pattern = "BL", replacement = "Baseline", dFrame$Timepoint)
dFrame$Timepoint <- gsub(pattern = "ONE", replacement = "Postop 1 month", dFrame$Timepoint)
dFrame$Timepoint <- gsub(pattern = "SIX", replacement = "Postop 6 months", dFrame$Timepoint)
dFrame$Timepoint <- gsub(pattern = "TWELVE", replacement = "Postop 12 months", dFrame$Timepoint)
dFrame$Timepoint <- gsub(pattern = "EIGHTEEN", replacement = "Postop 18 months", dFrame$Timepoint)
dFrame$Timepoint <- gsub(pattern = "TWENTY_FOUR", replacement = "Postop 24 months", dFrame$Timepoint)

names(dFrame)[names(dFrame) == "Surgery_type"] <- "Surgery type"
names(dFrame)[names(dFrame) == "Weight_kg"] <- "Weight (kg)"
names(dFrame)[names(dFrame) == "BMI_kgm2"] <- "BMI (kg/m2)"
names(dFrame)[names(dFrame) == "Excess_kg"] <- "Excess weight (kg)"
names(dFrame)[names(dFrame) == "Loss_from_BL_kg"] <- "Weight loss (kg)"
names(dFrame)[names(dFrame) == "Loss_from_BL_BMI"] <- "BMI loss (kg/m2)"
names(dFrame)[names(dFrame) == "Percent_Loss_kg"] <- "Percent loss"
names(dFrame)[names(dFrame) == "PEWL_kg"] <- "%EWL"

file.path <- paste0(outputDir, "EWL_Summary_Surgery.tsv")
write.table(dFrame, file.path,sep="\t",quote = FALSE, row.names = FALSE)
