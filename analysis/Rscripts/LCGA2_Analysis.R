#' ---
#' title: "2-class GMM-1 Analysis"
#' author: "Alicia Sorgen"
#' date: Dec 11, 2023
#' ---
#' 
## ----include=FALSE--------------------------------------------
rm(list=ls())
set.seed(1989)


## ----Library, include=FALSE-----------------------------------
R <- sessionInfo()
message(R$R.version$version.string); rm(R)

library(MASS); message("MASS: Version ", packageVersion("MASS"))
library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
library(lcmm); message("lcmm: Version ", packageVersion("lcmm"))
library(stringr); message("stringr: Version ", packageVersion("stringr"))
# library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
# library(scales); message("scales: Version ", packageVersion("scales"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))
library(vegan); message("vegan: Version ", packageVersion("vegan"))


## ----Script-Edits---------------------------------------------
rm(list=ls())

ANALYSIS <- "MetaPhlAn2_microbiome"
date = "2024Jun11"
params <- vector()
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
params <- c(params, "LCGA") # args[2]
params <- c(params, 2) # args[3]
params <- c(params, "RYGB") # args[4]
params <- c(params, "Univariate") # args[5]
params <- c(params, "Model2") # args[6]


level <- "genus"


## ----Working-Environment, include=FALSE-----------------------
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
  rm(gitRoot)
  
  root <- paste0("~/BioLockJ_pipelines/")
  pipeline <- paste0(ANALYSIS,"_analysis_", date)
  
  if (any(dir(root) %like% pipeline) == TRUE) {
    root <- paste0(root,"/",str_subset(dir(root), pipeline), "/")
  } else {
    today <- as.character(format(Sys.Date(), "%Y%b%d"))
    root <- paste0(root,ANALYSIS,"_analysis_", today, "/")
    dir.create(root, showWarnings = FALSE)
  }
  rm(pipeline)
  
  if (any(dir(root) == "input") == FALSE) {
    file.copy(gitInput,
              root,
              recursive = TRUE)
    
  }
  rm(gitInput)

    module <- paste0(args[5], "_",
                     args[4], "_",
                     gsub("-", "", args[2]), "_", 
                     args[6], "_",
                     args[3], "Class_Analysis")
  
  if (any(dir(root) %like% module) == TRUE) {
    moduleDir <- paste0(root,str_subset(dir(root), module), "/")
  } else {
    moduleDir <- paste0(root, module, "/")
    dir.create(moduleDir, showWarnings = FALSE)
  }
  rm(module, root)
  
  scriptDir <- paste0(moduleDir, "script/")
  if (any(dir(moduleDir) == "script") == FALSE) {
    dir.create(scriptDir, showWarnings = FALSE)
  }
  rm(scriptDir)
  
  outputDir <- paste0(moduleDir, "output")
  
  # Use unlink to remove all files and subdirectories within outputDir
  unlink(outputDir, recursive = TRUE)
  
  if (any(dir(moduleDir) == "output") == FALSE) {
    dir.create(outputDir, showWarnings = FALSE)
  }
  
  # files <- list.files(outputDir, recursive = TRUE, full.names = TRUE)
  # file.remove(files)
  # dirs <- list.dirs(outputDir, recursive = TRUE, full.names = TRUE)
  # file.remove(dirs)
  
  rm(outputDir, files)
  
  resourcesDir <- paste0(moduleDir, "resources/")
  if (any(dir(moduleDir) == "resources") == FALSE) {
    dir.create(resourcesDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/functions.R")
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
    
    script = paste0(gitScripts,"/", moduleRoot, ".R"); script
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
    rm(script)
    
  }
  rm(resourcesDir, gitScripts)
  
  setwd(paste0(moduleDir, "script/"))
  rm(moduleDir)
  
}

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())


## ----Functions------------------------------------------------
message("Functions\n")
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]
rm(str)

#' 
#' ### Input
## ----Input-------------------------------------------------------------
message("Input\n")
prevModule <- str_subset(dir(pipeRoot), "WeightMetaMerge")
inputDir = paste0(pipeRoot,"/",prevModule,"/output/")
inputFile <- "metadata.tsv"

prevModule2 <- paste0("Diversity_Metrics")
inputDir2 = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule2),"/output/", level, "/")
logCountFile <- paste0(level, "_LogNormalizedCounts_Diversity.tsv")



## ----Output-------------------------------------------------------------
message("Output\n")
outputDir = file.path(moduleDir,"output/")

## ----Script-Variables-----------------------------------------
message("Variables\n")
months <- c(0,1,6,12,18,24)
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

# outcomes <- c("Weight_kg", "BMI_kgm2", "Loss_from_BL_kg", "Percent_Loss_kg", "PEWL_kg")
# outcomeLabels <- c("Weight (kg)", "BMI (kg/m^2)", "Weight Loss from BL (kg)", "Weight Loss from BL (%)", "Excess Weight Loss from BL (%)")
outcomes <- c("Percent_Loss_kg", "Weight_kg")
outcomeLabel <- c("Weight Loss from BL (%)", "Weight_kg")
repeatedMeasure <- c("time")
covariates <- c("Surgery")
subjects <- c("SampleID", "PatientID")

modelType <- args[2]
n_class <- as.numeric(args[3])
surg <- args[4]
lmeMethod <- args[5]
ModelNum <- args[6]

Palette <- c("orange", "blue", "green3", "black", "red")

## ----Read-Table, echo=FALSE-----------------------------------
message("Read table\n")
# Read in table
filePath <- paste0(inputDir, inputFile)
myTable <- read.table(filePath, sep="\t", header = TRUE, check.names = FALSE)

filePath <- paste0(inputDir2, logCountFile)
taxa.df <- read.table(filePath, sep="\t", header = TRUE, check.names = FALSE)

# Trim variables
df.1 <- myTable[ , which(colnames(myTable) %in% c(outcomes, repeatedMeasure, covariates, subjects))]


## ----Data-Summary---------------------------------------------
message("Summary\n") 
outputDirNew <- paste0(outputDir, "Data-Summary/")
dir.create(outputDirNew, showWarnings = FALSE)

n <- nrow(df.1)
n.subjects <- length(unique(df.1$PatientID))

# Percent weight
sum.stats <- df.1 %>%
  group_by(time) %>%
  get_summary_stats(Percent_Loss_kg, type = "mean_sd")

# knitr::kable(sum.stats[, c(1,3)], "pipe", caption = "RYGB Patients at each timepoint", booktabs = TRUE)

outputName <- paste0("Percent_Loss_kg_Summary.tsv")
outputPath <- file.path(outputDirNew, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)


# Weight
sum.stats <- df.1 %>%
  group_by(time, Surgery) %>%
  get_summary_stats(Weight_kg, type = "mean_sd")

outputName <- paste0("Weight_kg_Summary.tsv")
outputPath <- file.path(outputDirNew, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
weight.table <- sum.stats[,c("time", "Surgery", "n")]

## ----BMI-Figure-----------------------------------------------
message("BMI figures\n") 
outputDirNew <- paste0(outputDir, "BMI-Figure/")
dir.create(outputDirNew, showWarnings = FALSE)

# Trim variables
df.2 <- myTable[ , which(colnames(myTable) %in% c( "BMI_kgm2", repeatedMeasure, covariates, subjects))]


df.2$Baseline_kg <- assignAllbyTerm(df.2, "time", "0", "BMI_kgm2", "PatientID")

df.2$`Observed - BL` <- df.2$BMI_kgm2 - df.2$Baseline_kg

df.2$`(Observed - BL)/BL` <- df.2$`Observed - BL` / df.2$Baseline_kg

df.2$observed_weight_change <- df.2$`(Observed - BL)/BL` * 100



df.RYGB <- df.2[df.2$Surgery == "RYGB",]

df.RYGB.O <- df.RYGB %>%
  group_by(time) %>%
  get_summary_stats(BMI_kgm2, type = "full")

df.SG <- df.2[df.2$Surgery == "SG",]

df.SG.O <- df.SG %>%
  group_by(time) %>%
  get_summary_stats(BMI_kgm2, type = "full")

stat.test <- df.2 %>%
  group_by(Surgery) %>%
  wilcox_test(BMI_kgm2 ~ time)
outputName <- paste0("observed_BMI_time_Wilcoxon.tsv")
outputPath <- file.path(outputDirNew, outputName)
write.table(stat.test, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

stat.test <- df.2 %>%
  group_by(time) %>%
  wilcox_test(BMI_kgm2 ~ Surgery)
outputName <- paste0("observed_BMI_Surgery_Wilcoxon.tsv")
outputPath <- file.path(outputDirNew, outputName)
write.table(stat.test, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

sum.stats <- df.2 %>%
  group_by(time, Surgery) %>%
  get_summary_stats(BMI_kgm2, type = "full")

outputName <- paste0("observed_BMI_Surgery_desc_statistics.tsv")
outputPath <- file.path(outputDirNew, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

# Generate the line graph using ggplot2
plot <- ggplot() +
  geom_point(data = df.RYGB.O, aes(x = time, y = median), color = "blue", shape = 17, size = 2) +
  geom_errorbar(data=df.RYGB.O, mapping=aes(x=time, ymin=q1, ymax=q3), width=0.5, linewidth=0.5, color="blue") + 
  geom_point(data = df.SG.O, aes(x = time, y = median), color = "forestgreen", shape = 18, size = 3) +
  geom_errorbar(data=df.SG.O, mapping=aes(x=time, ymin=q1, ymax=q3), width=0.5, linewidth=0.5, color="forestgreen") + 
  geom_line(data = df.RYGB.O, aes(x = time, y = mean), color = "blue" , linewidth=1) +
  geom_line(data = df.SG.O, aes(x = time, y = mean), color = "forestgreen", linewidth=1) +
  labs(x = "Follow-up time (months)", y = "Body mass index") +
  theme_bw(); plot



testPlot <- ggplot() +
  geom_point(data = df.RYGB.O, aes(x = time, y = median, color = "Median and IQR"), shape = 15, size = 3) +
  geom_errorbar(data=df.RYGB.O, mapping=aes(x=time, ymin=q1, ymax=q3), width=0.5, linewidth=0.5, color="black") + 
  geom_point(data = df.SG.O, aes(x = time, y = median, color = "Median and IQR"), shape = 15, size = 3) +
  geom_errorbar(data=df.SG.O, mapping=aes(x=time, ymin=q1, ymax=q3), width=0.5, linewidth=0.5, color="black") + 
  geom_line(data = df.RYGB.O, aes(x = time, y = mean, color = "RYGB (mean)") , size=1) +
  geom_line(data = df.SG.O, aes(x = time, y = mean, color = "SG (mean)"), size=1) +
  labs(x = "Follow-up time (months)", y = "Weight Change (%)") +
  scale_colour_manual(name="", values=c("black","blue", "forestgreen"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("blank", "solid", "solid"),
                        shape = c(15, NA, NA))))  +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")); testPlot

# now extract the legend
legend <- get_legend(testPlot)

filename.1 <- "Observed_BMI_by_Surgery_Median_IQR.pdf"
filePath <- paste0(outputDirNew, filename.1)

library(cowplot); message("cowplot: Version ", packageVersion("cowplot"))

pdf(filePath, width = 10, height = 5)
ggdraw(plot_grid(plot_grid(plot, ncol=1, align='v'),
                 plot_grid(legend, ncol=1),
                 rel_widths=c(1, 0.5)))

dev.off()





## ----ExcessWeight, echo=FALSE, fig.height=5, fig.width=10, warning=FALSE----
message("ExcessWeight\n")
outputDirNew <- paste0(outputDir, "ExcessWeight/")
dir.create(outputDirNew, showWarnings = FALSE)

# Trim variables
df.2 <- myTable[ , which(colnames(myTable) %in% c( "BMI_kgm2", repeatedMeasure, covariates, subjects))]

df.2$Baseline_BMI <- assignAllbyTerm(df.2, "time", "0", "BMI_kgm2", "PatientID")
df.2$Excess_BMI <- df.2$BMI_kgm2 - 25

df.2$Baseline_Excess_BMI <- assignAllbyTerm(df.2, "time", "0", "Excess_BMI", "PatientID")

df.2$BL_Ideal_Diff_BMI <- df.2$Baseline_BMI - 25

df.2$Excess_Loss_from_BL <- df.2$Baseline_Excess_BMI - df.2$Excess_BMI

df.2$Loss_from_BL <- df.2$Baseline_BMI - df.2$BMI_kgm2

df.2$PEWL <- (df.2$Loss_from_BL / df.2$BL_Ideal_Diff_BMI) * 100

sum.stats <- df.2 %>%
  group_by(time, Surgery) %>%
  get_summary_stats(PEWL, type = "full")

outputName <- paste0("observed_excess_weight_change_Surgery_desc_statistics.tsv")
outputPath <- file.path(outputDirNew, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

sum.stats <- df.2 %>%
  group_by(time) %>%
  get_summary_stats(PEWL, type = "full")

outputName <- paste0("observed_excess_weight_change_desc_statistics.tsv")
outputPath <- file.path(outputDirNew, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

## ----Figure1, echo=FALSE, fig.height=5, fig.width=10, warning=FALSE----
message("Figure 1\n")
outputDirNew <- paste0(outputDir, "Figure1/")
dir.create(outputDirNew, showWarnings = FALSE)

# Trim variables
df.1 <- myTable[ , which(colnames(myTable) %in% c( "Weight_kg", repeatedMeasure, covariates, subjects))]

if (lmeMethod == "Univariate") {
  df.2 <- data.frame()
  for (x in unique(df.1$Surgery)) {
    
    df.surg <- df.1[df.1$Surgery == x,]
    model <- lme( Weight_kg ~ time, method = "ML", random = ~1 | PatientID, data = df.surg, na.action = na.omit )
    df.surg$predicted_values <- predict(model, newdata = data.frame(time = df.surg$time, PatientID = df.surg$PatientID))
    
    df.2 <- rbind(df.2, df.surg)
    
  }
} else {
  df.2 <- df.1
  model <- lme( Weight_kg ~ time + Surgery, method = "ML", random = ~1 | PatientID, data = df.2, na.action = na.omit )
  df.2$predicted_values <- predict(model, newdata = data.frame(time = df.2$time, PatientID = df.2$PatientID, Surgery = df.2$Surgery))
  
}


df.2$modeled <- ifelse(is.na(df.2$Weight_kg) == TRUE, df.2$predicted_values, 
                       df.2$Weight_kg)

df.2$Baseline_kg <- assignAllbyTerm(df.2, "time", "0", "Weight_kg", "PatientID")

df.2$`Observed - BL` <- df.2$Weight_kg - df.2$Baseline_kg
df.2$`Modeled - BL` <- df.2$modeled - df.2$Baseline_kg

df.2$`(Observed - BL)/BL` <- df.2$`Observed - BL` / df.2$Baseline_kg
df.2$`(Modeled - BL)/BL` <- df.2$`Modeled - BL` / df.2$Baseline_kg

df.2$observed_weight_change <- df.2$`(Observed - BL)/BL` * 100
df.2$modeled_weight_change <- df.2$`(Modeled - BL)/BL` * 100

main.DF <- df.2

df.RYGB <- df.2[df.2$Surgery == "RYGB",]

df.RYGB.O <- df.RYGB %>%
  group_by(time) %>%
  get_summary_stats(observed_weight_change, type = "full")

df.RYGB.M <- df.RYGB %>%
  group_by(time) %>%
  get_summary_stats(modeled_weight_change, type = "full")


df.SG <- df.2[df.2$Surgery == "SG",]

df.SG.O <- df.SG %>%
  group_by(time) %>%
  get_summary_stats(observed_weight_change, type = "full")

df.SG.M <- df.SG %>%
  group_by(time) %>%
  get_summary_stats(modeled_weight_change, type = "full")


# Generate the line graph using ggplot2
plot <- ggplot() +
  geom_point(data = df.RYGB.O, aes(x = time, y = median, color = "Median and IQR (observed)"), shape = 15, size = 3) +
  geom_errorbar(data=df.RYGB.O, mapping=aes(x=time, ymin=q1, ymax=q3), width=0.5, linewidth=0.5, color="black") + 
  geom_point(data = df.SG.O, aes(x = time, y = median, color = "Median and IQR (observed)"), shape = 15, size = 3) +
  geom_errorbar(data=df.SG.O, mapping=aes(x=time, ymin=q1, ymax=q3), width=0.5, linewidth=0.5, color="black") + 
  geom_line(data = df.RYGB.M, aes(x = time, y = mean, color = "RYGB (modeled)") , size=1) +
  geom_line(data = df.SG.M, aes(x = time, y = mean, color = "SG (modeled)"), size=1) +
  labs(x = "Follow-up time (months)", y = "Weight Change (%)") +
  scale_colour_manual(name="", values=c("black","blue", "forestgreen"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("blank", "solid", "solid"),
                        shape = c(15, NA, NA))))  +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")); plot

filename.1 <- "MLM_Modeled_and_Observed_Weight_Change_by_Surgery_Median_IQR.pdf"
filePath <- paste0(outputDirNew, filename.1)

pdf(filePath, width = 10, height = 5)
print(plot)
dev.off()


# Generate the line graph using ggplot2
plot <- ggplot() +
  geom_point(data = df.RYGB.O, aes(x = time, y = mean, color = "Mean and SE (observed)"), shape = 15, size = 3) +
  geom_errorbar(data=df.RYGB.O, mapping=aes(x=time, ymin=(mean-se), ymax=(mean+se)), width=0.5, linewidth=0.5, color="black") + 
  geom_point(data = df.SG.O, aes(x = time, y = mean, color = "Mean and SE (observed)"), shape = 15, size = 3) +
  geom_errorbar(data=df.SG.O, mapping=aes(x=time, ymin=(mean-se), ymax=(mean+se)), width=0.5, linewidth=0.5, color="black") + 
  geom_line(data = df.RYGB.M, aes(x = time, y = mean, color = "RYGB (modeled)") , size=1) +
  geom_line(data = df.SG.M, aes(x = time, y = mean, color = "SG (modeled)"), size=1) +
  labs(x = "Follow-up time (months)", y = "Weight Change (%)") +
  scale_colour_manual(name="", values=c("black","blue", "forestgreen"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("blank", "solid", "solid"),
                        shape = c(15, NA, NA))))  +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")); plot

filename.1 <- "MLM_Modeled_and_Observed_Weight_Change_by_Surgery_Mean_SE.pdf"
filePath <- paste0(outputDirNew, filename.1)

pdf(filePath, width = 10, height = 5)
print(plot)
dev.off()

df.4 <- df.2[df.2$time != 0,]

stat.test <- df.4 %>%
  group_by(time) %>%
  wilcox_test(observed_weight_change ~ Surgery)

outputName <- paste0("observed_weight_change_Surgery_Wilcoxon.tsv")
outputPath <- file.path(outputDirNew, outputName)
write.table(stat.test, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

sum.stats <- df.4 %>%
  group_by(time, Surgery) %>%
  get_summary_stats(observed_weight_change, type = "full")

outputName <- paste0("observed_weight_change_Surgery_desc_statistics.tsv")
outputPath <- file.path(outputDirNew, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

sum.stats <- df.4 %>%
  group_by(time, Surgery) %>%
  get_summary_stats(`Observed - BL`, type = "full")

outputName <- paste0("net_observed_weight_change_Surgery_desc_statistics.tsv")
outputPath <- file.path(outputDirNew, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

stat.test <- df.4 %>%
  group_by(Surgery) %>%
  wilcox_test(observed_weight_change ~ time)
outputName <- paste0("observed_weight_change_time_Wilcoxon.tsv")
outputPath <- file.path(outputDirNew, outputName)
write.table(stat.test, outputPath, sep="\t", quote = FALSE, row.names = FALSE)



stat.test <- df.4 %>%
  group_by(time) %>%
  wilcox_test(modeled_weight_change ~ Surgery)

outputName <- paste0("modeled_weight_change_Surgery_Wilcoxon.tsv")
outputPath <- file.path(outputDirNew, outputName)
write.table(stat.test, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

sum.stats <- df.4 %>%
  group_by(time, Surgery) %>%
  get_summary_stats(modeled_weight_change, type = "full")

outputName <- paste0("modeled_weight_change_Surgery_desc_statistics.tsv")
outputPath <- file.path(outputDirNew, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
rm(df.RYGB, df.RYGB.M, df.RYGB.O, df.4, df.SG, df.SG.M, df.SG.O, df.surg)
rm(legend, model, plot, stat.test, sum.stats, testPlot, weight.table)
## ----Model----------------------------------------------------
message("Model\n")
outputDirNew <- paste0(outputDir, "Model/")
dir.create(outputDirNew, showWarnings = FALSE)

plotList <- list()
gmmSummary <- data.frame()
index <- 1


if (surg == "all") {
  surgTypes <- "all"
} else { 
  surgTypes <- c("RYGB", "SG")
}

df.merged <- data.frame()
for (z in 1:length(surgTypes)) {
  surg <- surgTypes[z]
  if (surg == "all") {
    df.1a <- df.2
  } else { 
    df.1a <- df.2[df.2$Surgery == surg,]
  }
  
  df.1a$PatientNum <- as.numeric(factor(df.1a$PatientID))
  
  if (ModelNum %in% c("Model1")) {
    if (modelType == "LCGA") {
      model1 <- hlme( data = df.1a,
                      modeled_weight_change ~ time + Surgery,
                      subject = "PatientNum",
                      ng = 1
      )
      
      modelx <- gridsearch( rep = 100,
                            maxiter = 10,
                            minit = model1,
                            hlme( data = df.1a,
                                  modeled_weight_change ~ time + Surgery,
                                  var.time = "time",
                                  subject = "PatientNum",
                                  ng = n_class,
                                  mixture = ~ time
                            )
      )
      
    }
    
    if (modelType == "GMM-1") {
      model1 <- hlme( data = df.1a,
                      modeled_weight_change ~ time + Surgery,
                      subject = "PatientNum",
                      random = ~1,  # GMM-1
                      ng = 1
      )
      
      modelx <- gridsearch( rep = 100,
                            maxiter = 10,
                            minit = model1,
                            hlme( data = df.1a,
                                  modeled_weight_change ~ time + Surgery,
                                  var.time = "time",
                                  subject = "PatientNum",
                                  random = ~1,  # GMM-1
                                  ng = n_class,
                                  mixture = ~ time,
                                  nwg = TRUE # GMM-1 & GMM-2
                            )
      )
    }
    
    if (modelType == "GMM-2") {
      model1 <- hlme( data = df.1a,
                      modeled_weight_change ~ time + Surgery,
                      subject = "PatientNum",
                      random = ~1 + time, # GMM-2
                      ng = 1
      )
      
      modelx <- gridsearch( rep = 100,
                            maxiter = 10,
                            minit = model1,
                            hlme( data = df.1a,
                                  modeled_weight_change ~ time + Surgery,
                                  var.time = "time",
                                  subject = "PatientNum",
                                  random = ~1 + time,  # GMM-2
                                  ng = n_class,
                                  mixture = ~ time,
                                  nwg = TRUE # GMM-1 & GMM-2
                            )
      )
      
    }
  }
  if (ModelNum %in% c("Model2")) {
    if (modelType == "LCGA") {
      model1 <- hlme( data = df.1a,
                      modeled_weight_change ~ time,
                      subject = "PatientNum",
                      ng = 1
      )
      
      modelx <- gridsearch( rep = 100,
                            maxiter = 10,
                            minit = model1,
                            hlme( data = df.1a,
                                  modeled_weight_change ~ time,
                                  var.time = "time",
                                  subject = "PatientNum",
                                  ng = n_class,
                                  mixture = ~ time
                            )
      )
      
    }
    
    if (modelType == "GMM-1") {
      model1 <- hlme( data = df.1a,
                      modeled_weight_change ~ time,
                      subject = "PatientNum",
                      random = ~1,  # GMM-1
                      ng = 1
      )
      
      modelx <- gridsearch( rep = 100,
                            maxiter = 10,
                            minit = model1,
                            hlme( data = df.1a,
                                  modeled_weight_change ~ time,
                                  var.time = "time",
                                  subject = "PatientNum",
                                  random = ~1,  # GMM-1
                                  ng = n_class,
                                  mixture = ~ time,
                                  nwg = TRUE # GMM-1 & GMM-2
                            )
      )
      
    }
    
    if (modelType == "GMM-2") {
      model1 <- hlme( data = df.1a,
                      modeled_weight_change ~ time,
                      subject = "PatientNum",
                      random = ~1 + time, # GMM-2
                      ng = 1
      )
      
      modelx <- gridsearch( rep = 100,
                            maxiter = 10,
                            minit = model1,
                            hlme( data = df.1a,
                                  modeled_weight_change ~ time,
                                  var.time = "time",
                                  subject = "PatientNum",
                                  random = ~1 + time,  # GMM-2
                                  ng = n_class,
                                  mixture = ~ time,
                                  nwg = TRUE # GMM-1 & GMM-2
                            )
      )
      
    }
  }
  
  smry <- as.data.frame( summarytable( modelx))
  
  class.df <- modelx$pprob
  x <- which(colnames(smry) == "BIC")
  BIC_value <- round(smry[1,x], 2)
  
  ng <- unique(class.df$class)
  
  df.3 <- merge(df.1a, class.df, by = "PatientNum", all.x = TRUE)
  df.3 <- df.3[,-((which(colnames(df.3) == "class")+1):ncol(df.3))]
  if (surg == "SG") {
    df.3$class <- ifelse(df.3$class == 1, 2, 1)
  }
  df.merged <- rbind(df.merged, df.3)
  n <- length(unique(df.3$PatientID))
  
  plot <- ggplot() +
    # xlim(0, 36)+
    labs(x = "Follow-up time (months)", y = "Weight Change (%)",
         title = paste0(modelType, " ", n_class, "-class model - ", surg, " patients")) +
    theme_bw(); plot
  
  df.C.O.List <- list()
  df.C.M.List <- list()
  ng <- ng[order(ng)]
  
  for (i in ng) {
    
    df.C <- df.3[df.3$class == i,]
    
    percent <- (length(unique(df.C$PatientID))/n)*100
    label <- paste0("Group ", i, " (n = ", length(unique(df.C$PatientID)), ", ", round(percent,1), "%)")
    
    
    df.C.O <- df.C %>%
      group_by(time) %>%
      get_summary_stats(observed_weight_change, type = "full")
    df.C.O.List[[i]] <- df.C.O
    
    df.C.M <- df.C %>%
      group_by(time) %>%
      get_summary_stats(modeled_weight_change, type = "full")
    df.C.M.List[[i]] <- df.C.M
    
    # Generate the line graph using ggplot2
    plot <- plot + 
      geom_point(data = df.C.O, aes(x = time, y = median), color = Palette[i], shape = (14+i), size = 3) +
      geom_errorbar(data=df.C.O, mapping=aes(x=time, ymin=q1, ymax=q3), width=0.5, linewidth=0.5, color=Palette[i]) + 
      geom_line(data = df.C.M, aes(x = time, y = mean), color = Palette[i], size=1) +
      annotate("text", x = 25, y = df.C.O$median[nrow(df.C.M)], hjust = 0, label = label, color = Palette[i])+
      coord_cartesian(xlim = c(0, 25), # This focuses the x-axis on the range of interest
                      clip = 'off') +   # This keeps the labels from disappearing
      theme(legend.title = element_blank(),
            legend.spacing.y = unit(0, "mm"), 
            panel.border = element_blank(),
            axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
            axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
            # aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            plot.margin = unit(c(1,10,1,1), "lines") # This widens the right margin; plot
      )
    
  } # for (i in 1:ng) 
  
  testPlot <- ggplot() +
    geom_point(data = df.C.O.List[[1]], aes(x = time, y = median, color = "Group median and IQR (observed)"), shape = 15, size = 3) +
    geom_errorbar(data=df.C.O.List[[1]], mapping=aes(x=time, ymin=q1, ymax=q3), width=0.5, linewidth=0.5, color="black") + 
    geom_point(data = df.C.O.List[[2]], aes(x = time, y = median, color = "Group median and IQR (observed)"), shape = 15, size = 3) +
    geom_errorbar(data=df.C.O.List[[2]], mapping=aes(x=time, ymin=q1, ymax=q3), width=0.5, linewidth=0.5, color="black") + 
    geom_line(data = df.C.M.List[[1]], aes(x = time, y = mean, color = "Group trajectory (modeled)") , size=1) +
    geom_line(data = df.C.M.List[[2]], aes(x = time, y = mean, color = "Group trajectory (modeled)"), size=1) +
    labs(x = "Follow-up time (months)", y = "Weight Change (%)") +
    scale_colour_manual(name="", values=c("grey","grey"),
                        guide = guide_legend(override.aes = list(
                          linetype = c("blank", "solid"),
                          shape = c(15, NA))))  +
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.spacing.y = unit(0, "mm"), 
          panel.border = element_rect(colour = "black", fill=NA),
          aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black")); testPlot
  
  # now extract the legend
  legend <- get_legend(testPlot)
  plotList[[index]] <- plot
  index <- index + 1
}




filename.2 <- paste0(modelType, "_Modeled_and_Observed_Weight_Change.pdf")
filePath <- paste0(outputDirNew, filename.2)

pdf(filePath, width = 10, height = 10)
grid.arrange(plotList[[1]], plotList[[2]],
             ncol = 1, nrow = 2)

dev.off()

plotList <- list()
gmmSummary <- data.frame()
index <- 1


df.merged <- data.frame()

df.2$PatientNum <- as.numeric(factor(df.2$PatientID))

if (ModelNum %in% c("Model1")) {
  if (modelType == "LCGA") {
    model1 <- hlme( data = df.2,
                    modeled_weight_change ~ time + Surgery,
                    subject = "PatientNum",
                    ng = 1
    )
    
    modelx <- gridsearch( rep = 100,
                          maxiter = 10,
                          minit = model1,
                          hlme( data = df.2,
                                modeled_weight_change ~ time + Surgery,
                                var.time = "time",
                                subject = "PatientNum",
                                ng = n_class,
                                mixture = ~ time
                          )
    )
    
  }
  
  if (modelType == "GMM-1") {
    model1 <- hlme( data = df.2,
                    modeled_weight_change ~ time + Surgery,
                    subject = "PatientNum",
                    random = ~1,  # GMM-1
                    ng = 1
    )
    
    modelx <- gridsearch( rep = 100,
                          maxiter = 10,
                          minit = model1,
                          hlme( data = df.2,
                                modeled_weight_change ~ time + Surgery,
                                var.time = "time",
                                subject = "PatientNum",
                                random = ~1,  # GMM-1
                                ng = n_class,
                                mixture = ~ time,
                                nwg = TRUE # GMM-1 & GMM-2
                          )
    )
  }
  
  if (modelType == "GMM-2") {
    model1 <- hlme( data = df.2,
                    modeled_weight_change ~ time + Surgery,
                    subject = "PatientNum",
                    random = ~1 + time, # GMM-2
                    ng = 1
    )
    
    modelx <- gridsearch( rep = 100,
                          maxiter = 10,
                          minit = model1,
                          hlme( data = df.2,
                                modeled_weight_change ~ time + Surgery,
                                var.time = "time",
                                subject = "PatientNum",
                                random = ~1 + time,  # GMM-2
                                ng = n_class,
                                mixture = ~ time,
                                nwg = TRUE # GMM-1 & GMM-2
                          )
    )
    
  }
}
if (ModelNum %in% c("Model2")) {
  if (modelType == "LCGA") {
    model1 <- hlme( data = df.2,
                    modeled_weight_change ~ time,
                    subject = "PatientNum",
                    ng = 1
    )
    
    modelx <- gridsearch( rep = 100,
                          maxiter = 10,
                          minit = model1,
                          hlme( data = df.2,
                                modeled_weight_change ~ time,
                                var.time = "time",
                                subject = "PatientNum",
                                ng = n_class,
                                mixture = ~ time
                          )
    )
    
  }
  
  if (modelType == "GMM-1") {
    model1 <- hlme( data = df.2,
                    modeled_weight_change ~ time,
                    subject = "PatientNum",
                    random = ~1,  # GMM-1
                    ng = 1
    )
    
    modelx <- gridsearch( rep = 100,
                          maxiter = 10,
                          minit = model1,
                          hlme( data = df.2,
                                modeled_weight_change ~ time,
                                var.time = "time",
                                subject = "PatientNum",
                                random = ~1,  # GMM-1
                                ng = n_class,
                                mixture = ~ time,
                                nwg = TRUE # GMM-1 & GMM-2
                          )
    )
    
  }
  
  if (modelType == "GMM-2") {
    model1 <- hlme( data = df.2,
                    modeled_weight_change ~ time,
                    subject = "PatientNum",
                    random = ~1 + time, # GMM-2
                    ng = 1
    )
    
    modelx <- gridsearch( rep = 100,
                          maxiter = 10,
                          minit = model1,
                          hlme( data = df.2,
                                modeled_weight_change ~ time,
                                var.time = "time",
                                subject = "PatientNum",
                                random = ~1 + time,  # GMM-2
                                ng = n_class,
                                mixture = ~ time,
                                nwg = TRUE # GMM-1 & GMM-2
                          )
    )
    
  }
}

smry <- as.data.frame( summarytable( modelx))

class.df <- modelx$pprob
x <- which(colnames(smry) == "BIC")
BIC_value <- round(smry[1,x], 2)

ng <- unique(class.df$class)

df.3 <- merge(df.2, class.df, by = "PatientNum", all.x = TRUE)
df.3 <- df.3[,-((which(colnames(df.3) == "class")+1):ncol(df.3))]
df.3$class <- ifelse(df.3$class == 1, 2, 1)
df.merged <- rbind(df.merged, df.3)
n <- length(unique(df.3$PatientID))

plot <- ggplot() +
  # xlim(0, 36)+
  labs(x = "Follow-up time (months)", y = "Weight Change (%)",
       title = paste0(modelType, " ", n_class, "-class model - Combined patients")) +
  theme_bw(); plot

df.C.O.List <- list()
df.C.M.List <- list()
ng <- ng[order(ng)]

for (i in ng) {
  
  df.C <- df.3[df.3$class == i,]
  
  percent <- (length(unique(df.C$PatientID))/n)*100
  label <- paste0("Group ", i, " (n = ", length(unique(df.C$PatientID)), ", ", round(percent,1), "%)")
  
  
  df.C.O <- df.C %>%
    group_by(time) %>%
    get_summary_stats(observed_weight_change, type = "full")
  df.C.O.List[[i]] <- df.C.O
  
  df.C.M <- df.C %>%
    group_by(time) %>%
    get_summary_stats(modeled_weight_change, type = "full")
  df.C.M.List[[i]] <- df.C.M
  
  # Generate the line graph using ggplot2
  plot <- plot + 
    geom_point(data = df.C.O, aes(x = time, y = median), color = Palette[i], shape = (14+i), size = 3) +
    geom_errorbar(data=df.C.O, mapping=aes(x=time, ymin=q1, ymax=q3), width=0.5, linewidth=0.5, color=Palette[i]) + 
    geom_line(data = df.C.M, aes(x = time, y = mean), color = Palette[i], linewidth=1) +
    annotate("text", x = 25, y = df.C.O$median[nrow(df.C.M)], hjust = 0, label = label, color = Palette[i])+
    coord_cartesian(xlim = c(0, 25), # This focuses the x-axis on the range of interest
                    clip = 'off') +   # This keeps the labels from disappearing
    theme(legend.title = element_blank(),
          legend.spacing.y = unit(0, "mm"), 
          panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
          # aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          plot.margin = unit(c(1,10,1,1), "lines") # This widens the right margin; plot
    )
  
} # for (i in 1:ng) 

testPlot <- ggplot() +
  geom_point(data = df.C.O.List[[1]], aes(x = time, y = median, color = "Group median and IQR (observed)"), shape = 15, size = 3) +
  geom_errorbar(data=df.C.O.List[[1]], mapping=aes(x=time, ymin=q1, ymax=q3), width=0.5, linewidth=0.5, color="black") + 
  geom_point(data = df.C.O.List[[2]], aes(x = time, y = median, color = "Group median and IQR (observed)"), shape = 15, size = 3) +
  geom_errorbar(data=df.C.O.List[[2]], mapping=aes(x=time, ymin=q1, ymax=q3), width=0.5, linewidth=0.5, color="black") + 
  geom_line(data = df.C.M.List[[1]], aes(x = time, y = mean, color = "Group trajectory (modeled)") , linewidth=1) +
  geom_line(data = df.C.M.List[[2]], aes(x = time, y = mean, color = "Group trajectory (modeled)"), linewidth=1) +
  labs(x = "Follow-up time (months)", y = "Weight Change (%)") +
  scale_colour_manual(name="", values=c("grey","grey"),
                      guide = guide_legend(override.aes = list(
                        linetype = c("blank", "solid"),
                        shape = c(15, NA))))  +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")); testPlot

# now extract the legend
legend <- get_legend(testPlot)
plotList[[index]] <- plot
index <- index + 1




filename.2 <- paste0(modelType, "_Modeled_and_Observed_Weight_Change_All.pdf")
filePath <- paste0(outputDirNew, filename.2)

pdf(filePath, width = 10, height = 5)
for (p in 1:length(plotList)) {
  print(plotList[[p]])
}

dev.off()






## ----Figure2-code, fig.height=5, fig.width=10, include=FALSE----
message("Figure 2\n")
outputDirNew <- paste0(outputDir, "Figure2/")
dir.create(outputDirNew, showWarnings = FALSE)

df.2 <- df.merged[df.merged$time == 0,]

f.test <- fisher.test(df.2$Surgery, df.2$class, workspace = 2e8)

freq <- as.data.frame(table(df.2$Surgery, df.2$class))
names(freq)[names(freq) == "Var1"] <- "Surgery"
names(freq)[names(freq) == "Var2"] <- "Group"

n_RYGB <- sum(freq$Freq[which(freq$Surgery == "RYGB")])
n_SG <- sum(freq$Freq[which(freq$Surgery == "SG")])

freq$Percent <- ifelse(freq$Surgery == "RYGB", (freq$Freq / n_RYGB) * 100,
                       (freq$Freq / n_SG) * 100)

plot <- ggplot(freq, aes(fill=Group, y=Freq, x=Surgery)) + 
  geom_bar(position="stack", stat="identity", color = "black")+
  labs(y = "Number of patients")+
  scale_fill_manual(values=Palette[1:n_class]); plot

subtitle.lab <- paste0("Fisher test ", round(f.test$p.value, 3))
plot <- plot + labs(subtitle = subtitle.lab); plot

plot <- plot + geom_text(aes(label = paste0(round(Percent, 1), "%")), size = 3, hjust = 0.5, position = "stack", vjust = 1.5, colour = "black"); plot
# print(plot)
filename.2 <- paste0(modelType, "_", n_class, "_Groups_Surgery_Fisher.pdf")
filePath <- paste0(outputDirNew, filename.2)

pdf(filePath, width = 10, height = 5)
print(plot)
dev.off()


## ----Overall-MLM, fig.height=5, fig.width=10, warning=FALSE, include=FALSE----
message("Overall MLM\n")
outputDirNew <- paste0(outputDir, "Overall_MLM/")
dir.create(outputDirNew, showWarnings = FALSE)

months <- c(0, 1, 6, 12, 18, 24)

outputResults <- paste0(outputDirNew)
# dir.create(outputResults, showWarnings = FALSE)

startAbundanceIndex <- which(colnames(taxa.df)=="ResponderStatus")+1
SampleID <- taxa.df$SampleID
taxa.df2 <- taxa.df[,startAbundanceIndex:ncol(taxa.df)]
taxa.df2$SampleID <- SampleID

DF <- merge(df.merged, taxa.df2, by = "SampleID")


main.DF <- DF
outputName <- paste0("Modeled_Weight_Taxa.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(main.DF, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

results.df.list <- list()
df.index <- 1
for ( t in 2:length(months) ) {
  
  timeFrame <- months[1:t]
  log.df <- DF[ DF$time %in% timeFrame, ]
  startAbundanceIndex <- which(colnames(log.df)=="modeled_weight_change")+1

  ## Parse out metadata from counts
  myT<-log.df[,startAbundanceIndex:ncol(log.df)]
  
  toRemove <- which(colnames(myT) %like% "noname")
  myT <- myT[,-toRemove]
  
  toRemove <- which(colnames(myT) %like% "unclassified")
  myT <- myT[,-toRemove]

  toRemove <- which(colnames(myT) %like% "virus")
  if (length(toRemove) > 0) {
    myT <- myT[,-toRemove]
  }

  PatientID <- log.df$PatientID
  Timepoint <- log.df$time
  Surgery <- log.df$Surgery
  Observed_Weight <- log.df$Weight_kg
  Modeled_Weight <- log.df$modeled
  Group <- log.df$class
  
  #----- Mixed Linear Modeling -----
  bugName <- vector()
  pVal_Timepoint <- vector()
  pVal_Surgery <- vector()
  pVal_Observed_Weight <- vector()
  pVal_Modeled_Weight <- vector()
  
  index <- 1
  
  for (i in 1:ncol(myT)){
    
    bug<-myT[,i]
    
    if (mean(bug>0, na.rm = TRUE)>0.1){
      
      df <- data.frame(bug, PatientID, Timepoint, Surgery, Observed_Weight, Modeled_Weight, Group)
      df <- na.omit(df)
      df$Timepoint <- as.factor(df$Timepoint)
      
      tryCatch({
        mlm <- lme( bug ~ Timepoint + Surgery +  Observed_Weight, method = "ML", random = ~1 | PatientID, data = df )
        fit<-anova(mlm)
        
        pVal_Timepoint[index] <- fit$"p-value"[2]
        pVal_Surgery[index] <- fit$"p-value"[3]
        pVal_Observed_Weight[index] <- fit$"p-value"[4]
        
        mlm <- lme( bug ~ Timepoint + Surgery +  Modeled_Weight, method = "ML", random = ~1 | PatientID, data = df )
        fit<-anova(mlm)
        
        pVal_Modeled_Weight[index] <- fit$"p-value"[4]
        
        bugName[index]<-colnames(myT)[i]
        
        index<-index+1
        
      }, error = function(e){})
    } # if (mean(bug>0, na.rm = TRUE)>0.1)
    
  } # for (i in 1:ncol(myT))

  dFrame<-data.frame(pVal_Timepoint, pVal_Surgery, pVal_Observed_Weight, pVal_Modeled_Weight)
  dFrame$Adj_pVal_Timepoint <- p.adjust(dFrame$pVal_Timepoint,method = "BH")
  dFrame$Adj_pVal_Surgery <- p.adjust(dFrame$pVal_Surgery,method = "BH")
  dFrame$Adj_pVal_Observed_Weight <- p.adjust(dFrame$pVal_Observed_Weight,method = "BH")
  dFrame$Adj_pVal_Modeled_Weight <- p.adjust(dFrame$pVal_Modeled_Weight,method = "BH")
  
  names(dFrame) <- paste0(names(dFrame), "_BLto", timeFrame[length(timeFrame)],"M")
  dFrame <- data.frame(bugName, dFrame)
  #----- Mixed Linear Model END -----
  
  results.df.list[[df.index]] <- dFrame
  df.index <- df.index + 1
  
} # for ( t in 2:length(months) )

merged.results <- results.df.list[[1]]
if (length(results.df.list) > 1) {
  for (l in 2:length(results.df.list)) {

    merged.results <- merge(merged.results, results.df.list[[l]], by = "bugName")

  } # for (l in 2:length(results.df.list))
}

filename.3 <- paste0(level,  "_by_Timepoint+Surgery+Observed_Modeled_Weight_MixedLinearModelResults", modelType, "_", n_class,".tsv")
write.table(merged.results, paste0(outputResults, filename.3),sep="\t",row.names = FALSE,quote = FALSE)


## ----LinearModel-Analysis-p-values, fig.height=5, fig.width=10, warning=FALSE, include=FALSE----
message("Linear Model\n")
outputDirNew <- paste0(outputDir, "LinearModel_Analysis_pvalues/")
dir.create(outputDirNew, showWarnings = FALSE)


months <- c(0, 1, 6, 12, 18, 24)

outputResults <- paste0(outputDirNew, "ResultsTables/")
dir.create(outputResults, showWarnings = FALSE)

startAbundanceIndex <- which(colnames(taxa.df)=="ResponderStatus")+1
SampleID <- taxa.df$SampleID
taxa.df2 <- taxa.df[,startAbundanceIndex:ncol(taxa.df)]
taxa.df2$SampleID <- SampleID

DF <- merge(df.merged, taxa.df2, by = "SampleID")

filename.4 <- paste0(level,  "_by_Timepoint+Surgery+Observed_Modeled_Weight_MixedLinearModelPlots", modelType, "_", n_class,".pdf")
filePath <- paste0(outputDirNew, filename.4)
pdf(filePath)

par(mfrow=c(2,2))

groupPValues1 <- vector()
groupPValues2 <- vector()
groupPValues3 <- vector()
groupPValues4 <- vector()
groupNumbers <- vector()
qbugName <- vector()
direction <- vector()
LengthTime <- vector()
Num_Patients <- vector()
qIndex <- 1

for( i in 1:nrow(merged.results)) {
  BugName <- merged.results$bugName[i]
  
  for ( t in 2:length(months) ) {
    
    timeFrame <- months[1:t]
    log.df <- DF[ DF$time %in% timeFrame, ]
    bugCol <- log.df[ , which(colnames(log.df) == BugName) ]
    Month <- log.df$time
    Surgery <- log.df$Surgery
    Observed_Weight <- log.df$Weight_kg
    Modeled_Weight <- log.df$modeled
    PatientID <- log.df$PatientID
    Group <- log.df$class
    
    colEnding <- paste0("BLto", timeFrame[length(timeFrame)],"M")
    Timepoint_p <- merged.results[ i , which(colnames(merged.results) == paste0("Adj_pVal_Timepoint_", colEnding))]
    Surgery_p <- merged.results[ i , which(colnames(merged.results) == paste0("Adj_pVal_Surgery_", colEnding))]
    Observed_Weight_p <- merged.results[ i , which(colnames(merged.results) == paste0("Adj_pVal_Observed_Weight_", colEnding))]
    Modeled_Weight_p <- merged.results[ i , which(colnames(merged.results) == paste0("Adj_pVal_Modeled_Weight_", colEnding))]
    
    df <- data.frame(PatientID, Surgery, Observed_Weight, bugCol)
    df <- na.omit(df)
    patient_n <- length(unique(df$PatientID))
    n <- nrow(df)
    
    title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Observed_Weight_p, digits = 3), "; n = ", n)
    plot(  df$Observed_Weight, df$bugCol,
           ylab =paste0(BugName, " log10 Abundance") ,
           xlab="Observed Weight (kg)",
           main=title.lab,
           cex.main = 1 ,
           col = ifelse( df$Surgery == "SG", "red", "blue" )
    )
    legend("topright", legend=c("SG", "RYGB"),
           col=c("red", "blue"), pch = 21, cex=0.8)
    
    df <- data.frame(PatientID, Surgery, Modeled_Weight, bugCol)
    df <- na.omit(df)
    patient_n <- length(unique(df$PatientID))
    n <- nrow(df)
    title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Modeled_Weight_p, digits = 3), "; n = ", n)
    plot(  df$Modeled_Weight, df$bugCol,
           ylab =paste0(BugName, " log10 Abundance") ,
           xlab="Modeled Weight (kg)",
           main=title.lab,
           cex.main = 1 ,
           col = ifelse( df$Surgery == "SG", "red", "blue" )
    )
    legend("topright", legend=c("SG", "RYGB"),
           col=c("red", "blue"), pch = 21, cex=0.8)
    
    df <- data.frame(PatientID, Surgery, Month, bugCol)
    df <- na.omit(df)
    patient_n <- length(unique(df$PatientID))
    n <- nrow(df)
    
    title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Timepoint_p, digits = 3), "; n = ", n)
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
    
    df <- data.frame(PatientID, Surgery, bugCol)
    df <- na.omit(df)
    patient_n <- length(unique(df$PatientID))
    n <- nrow(df)
    title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Surgery_p, digits = 3), "; n = ", n)
    boxplot( df$bugCol ~  df$Surgery ,
             xlab="Surgery Type",
             ylab =BugName,
             main = title.lab,
             cex.main = 1
    )
    
    qFrameFinal <- data.frame()
    for( j in 1:n_class ) {
      
      lm.df <- data.frame(PatientID, Month, Surgery, bugCol, Group)
      lm.df <- na.omit(lm.df)
      lm.df <- lm.df[lm.df$Group == j,]
      # Num_Patients[qIndex] <- length(unique(lm.df$PatientID))
      
      if (mean(lm.df$bugCol>0, na.rm = TRUE)>0.1){
        tryCatch({
          myLm <- lm(  lm.df$bugCol ~ lm.df$Month) # add surgery type, surgery*time
          p = anova(myLm)$"Pr(>F)"[1]
          groupPValues1[qIndex] <- p
          
          title.lab <- paste0( "Group ", j, " (", length(unique(lm.df$PatientID)), " patients)\n",
                               BugName ,
                               "\nTimepoint (p = ", format(p, digits = 3), ")")
          plot(  lm.df$Month, lm.df$bugCol ,
                 ylab =BugName ,
                 xlab="Timepoint (months)",
                 main=title.lab,
                 cex.main = 1,
                 col = ifelse( lm.df$Surgery == "SG", "red", "blue" )
          )
          legend("topright", legend=c("SG", "RYGB"),
                 col=c("red", "blue"), pch = 21, cex=0.8)
          
        }, error = function(e){})
        tryCatch({
          # myLm <- lm(  lm.df$bugCol ~ lm.df$Month + lm.df$Surgery + lm.df$Month:lm.df$Surgery) # add surgery type, surgery*time
          myLm <- lm(  lm.df$bugCol ~ lm.df$Surgery + lm.df$Month + lm.df$Surgery:lm.df$Month) # add surgery type, surgery*time
          p = anova(myLm)$"Pr(>F)"[2]
          groupPValues2[qIndex] <- p
          
        }, error = function(e){})
        tryCatch({
          # mlm <- lme( bugCol ~ Month + Surgery + Month:Surgery, method = "REML", random = ~1 | PatientID, data = lm.df )
          mlm <- lme( bugCol ~ Surgery + Month + Surgery:Month, method = "REML", random = ~1 | PatientID, data = lm.df )
          fit<-anova(mlm)
          p = fit$"p-value"[3]
          groupPValues3[qIndex] <- p
          
        }, error = function(e){})
        
        tryCatch({
          mlm <- lme( bugCol ~ Month, method = "REML", random = ~1 | PatientID, data = lm.df )
          fit<-anova(mlm)
          p = fit$"p-value"[2]
          groupPValues4[qIndex] <- p
        }, error = function(e){})
        
        groupNumbers[qIndex] <- j
        qbugName[qIndex] <- BugName
        
        
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
        
      }
      
    } # for( j in 1:ng )
    
    # length(
    #   # qbugName
    #   # LengthTime
    #   Num_Patients
    #   # groupNumbers
    #   # groupPValues1
    #   # groupPValues2
    #   # groupPValues3
    #   # groupPValues4
    #   # direction
    # )
    qFrame <- data.frame( qbugName, LengthTime, groupNumbers,groupPValues1, groupPValues2, groupPValues3, groupPValues4, direction)
    qFrame$adjPValues1 <- p.adjust(qFrame$groupPValues1, method="BH")
    qFrame$adjPValues2 <- p.adjust(qFrame$groupPValues2, method="BH")
    qFrame$adjPValues3 <- p.adjust(qFrame$groupPValues3, method="BH")
    qFrame$adjPValues4 <- p.adjust(qFrame$groupPValues4, method="BH")
    
    qFrameFinal <- rbind(qFrameFinal, qFrame)
    
  } # for ( t in 2:length(months) )
  
} # for( i in 1:nrow(merged.results))
dev.off()



months <- c(0, 1, 6, 12, 18, 24)


startAbundanceIndex <- which(colnames(taxa.df)=="ResponderStatus")+1
SampleID <- taxa.df$SampleID
taxa.df2 <- taxa.df[,startAbundanceIndex:ncol(taxa.df)]
taxa.df2$SampleID <- SampleID

DF <- merge(df.merged, taxa.df2, by = "SampleID")

qFrameFinal <- qFrameFinal[ order( qFrameFinal$groupPValues1),]

filename.5 <- paste0(level, "_by_Timepoint_groupedby_", modelType, "_", n_class,"_LinearModelResults.tsv")
filePath.5 <- paste0(outputResults, filename.5)
write.table(qFrameFinal, file=filePath.5, sep="\t", row.names=FALSE)


## ----Figure3-code, fig.height=5, fig.width=7, warning=FALSE, include=FALSE----
message("Figure 3\n")
outputDirNew <- paste0(outputDir, "Figure3/")
dir.create(outputDirNew, showWarnings = FALSE)


# viruses <- c("Viruses_noname", "C2likevirus")
# unnamed <- c("noname", "unclassified")

qFrame <- read.table(filePath.5, sep="\t",header = TRUE, check.names = FALSE)

qFrame$adjPValues1 <- ifelse(qFrame$adjPValues1 == 0, 2e-16,
                             qFrame$adjPValues1)

qFrame$log <- ifelse(qFrame$direction == "increase", -log10(qFrame$adjPValues1),
                     log10(qFrame$adjPValues1))
qFrame <- na.omit(qFrame)

# Label taxa as significant or not based on adjusted p-value
qFrame$Significance <- ifelse(qFrame$adjPValues1 < 0.05, "Significant",
                              "Not Significant")
# qFrame <- qFrame[!( qFrame$qbugName %in% viruses ),]

timeFrames <- c("BLto1M", "BLto6M", "BLto12M", "BLto18M", "BLto24M")
filename.9 <- paste0(level, "_analysis1_pvalue_boxplots.pdf")
filePath <- paste0(outputDirNew, filename.9)
pdf(filePath, width = 7, height = 5)

for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] ) {
  
  df <- qFrame[ qFrame$LengthTime == timeFrame, ]
  endPoint <- gsub("BLto", "", timeFrame)
  
  f.test <- fisher.test(df$groupNumbers, df$Significance, workspace = 2e8)
  f.pVal <- paste0("p = ", format(f.test$p.value, digits = 3))
  # f.pVal <- roundP(f.test$p.value)
  
  num_patients <- unique(df$Num_Patients)
  title.lab <- paste0( str_to_title(level),  " level taxonomic changes over time (BL to ", endPoint, ")"); title.lab
  subtitle.lab <- paste0("Fisher's Exact Test (", f.pVal, ")")
  caption.lab <- "lm(  Taxa ~ Timepoint )"
  
  hLine1 <- -log10(0.05)
  hLine2 <- log10(0.05)
  
  Q2u <- length(which(df$log > hLine1 & df$groupNumbers == 2))
  Q2d <- length(which(df$log < hLine2 & df$groupNumbers == 2))
  Q1u <- length(which(df$log > hLine1 & df$groupNumbers == 1))
  Q1d <- length(which(df$log < hLine2 & df$groupNumbers == 1))
  
  labels1 <- c(Q1u, Q2u)
  labels2 <- c(Q1d, Q2d)
  
  getWeightGroupLocation1 <- function(df, hLine1) {
    
    df <- na.omit(df)
    quants1 <- data.table(df)[, list(quant = as.numeric(max(log))), by = groupNumbers]
    quants1 <- quants1[order(quants1$groupNumbers),]
    for (k in 1:nrow(quants1)) {
      if (quants1$quant[k] < hLine1) {
        quants1$quant[k] <- hLine1
      }
    }
    quants1$quant <- quants1$quant + 0.5
    
    return(quants1)
  }
  getWeightGroupLocation2 <- function(df, hLine2) {
    
    df <- na.omit(df)
    quants2 <- data.table(df)[, list(quant = as.numeric(min(log))), by = groupNumbers]
    quants2 <- quants2[order(quants2$groupNumbers),]
    for (k in 1:nrow(quants2)) {
      if (quants2$quant[k] > hLine2) {
        quants2$quant[k] <- hLine2
      }
    }
    
    quants2$quant <- quants2$quant - 0.5
    
    return(quants2)
  }
  
  quants1 <- getWeightGroupLocation1(df, hLine1)
  quants2 <- getWeightGroupLocation2(df, hLine2)
  
  ## Generate plot ##
  
  plot <- ggboxplot(
    df, x = "groupNumbers", y = "log",
    add = "jitter", shape = 1, 
    palette = c("orange", "blue", "green3", "black", "red"),
    fill = "groupNumbers",
    # color = "groupNumbers",
    # facet.by = "groupNumbers",
    scales = "free",
    # label="qbugName",
    outlier.shape = 1
  )
  
  # plot <- plot +
  #   xlim(0,4)
  
  # Add title & subtitle
  plot <- plot +
    labs(title = title.lab,
         subtitle = subtitle.lab,
         caption = caption.lab)
  
  plot <- plot +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5)
    )
  
  # Edit x- and y-axis labels
  plot <- plot +
    labs(x = paste0(modelType, " Groups"), y = "p value (log10)")
  
  # Add dashed line for significant p values increasing
  plot <- plot + geom_hline(yintercept=hLine1, linetype="dashed", color = "red")
  
  # Add dashed line for significant p values decreasing
  plot <- plot + geom_hline(yintercept=hLine2, linetype="dashed", color = "red")
  
  # Increase text size of labels
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Remove legend
  plot <- plot + theme(legend.position = "none")
  
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Add text indicating the number of significant taxa increasing over time
  plot <- plot + geom_text(data = quants1, aes(x = groupNumbers, y = quant, label = labels1), size = 5)
  
  # Add text indicating the number of significant taxa decreasing over time
  plot <- plot + geom_text(data = quants2, aes(x = groupNumbers, y = quant, label = labels2), size = 5)
  
  # Add up arrow and annotation
  plot2 <- plot + annotate("segment", 
                           x=0.45, xend=0.45,
                           y=0.25, yend=2,
                           col="black", arrow=arrow(length=unit(0.25, "cm")))
  
  plot2 <- plot2 + annotate("text", x=0.5, y=1.5, hjust = 0, label = "Increases in abundance")
  
  # Add down arrow and annotation
  plot2 <- plot2 + annotate("segment", 
                            x=0.45, xend=0.45,
                            y=-0.25, yend=-2,
                            col="black", arrow=arrow(length=unit(0.25, "cm")))
  plot2 <- plot2 + annotate("text", x=0.5, y=-1.5, hjust = 0, label = "Decreases in abundance")
  
  print(plot2)
  
} # for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] )

dev.off()


# knitr::knit_exit()

## ----Figure4-code, fig.height=5, fig.width=7, warning=FALSE, include=FALSE----
message("Figure 4\n")

outputDirNew <- paste0(outputDir, "Figure4/")
dir.create(outputDirNew, showWarnings = FALSE)

viruses <- c("Viruses_noname", "C2likevirus")

qFrame <- read.table(filePath.5, sep="\t",header = TRUE, check.names = FALSE)

qFrame$adjPValues2 <- ifelse(qFrame$adjPValues2 == 0, 2e-16,
                             qFrame$adjPValues2)

qFrame$log <- ifelse(qFrame$direction == "increase", -log10(qFrame$adjPValues2),
                     log10(qFrame$adjPValues2))
qFrame <- na.omit(qFrame)

# Label taxa as significant or not based on adjusted p-value
qFrame$Significance <- ifelse(qFrame$adjPValues2 < 0.05, "Significant",
                              "Not Significant")
qFrame <- qFrame[!( qFrame$qbugName %in% viruses ),]

timeFrames <- c("BLto1M", "BLto6M", "BLto12M", "BLto18M", "BLto24M")
filename.10 <- paste0(level, "_analysis2_pvalue_boxplots.pdf")
filePath <- paste0(outputDirNew, filename.10)
pdf(filePath, width = 7, height = 5)

for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] ) {
  
  # message(timeFrame)
  df <- qFrame[ qFrame$LengthTime == timeFrame, ]
  endPoint <- gsub("BLto", "", timeFrame)
  
  f.test <- fisher.test(df$groupNumbers, df$Significance, workspace = 2e8)
  f.pVal <- paste0("p = ", format(f.test$p.value, digits = 3))
  # f.pVal <- roundP(f.test$p.value)
  
  num_patients <- unique(df$Num_Patients)
  title.lab <- paste0( str_to_title(level),  " level taxonomic changes over time (BL to ", endPoint, ")"); title.lab
  subtitle.lab <- paste0("Fisher's Exact Test (", f.pVal, ")")
  
  hLine1 <- -log10(0.05)
  hLine2 <- log10(0.05)
  
  Q2u <- length(which(df$log > hLine1 & df$groupNumbers == 2))
  Q2d <- length(which(df$log < hLine2 & df$groupNumbers == 2))
  Q1u <- length(which(df$log > hLine1 & df$groupNumbers == 1))
  Q1d <- length(which(df$log < hLine2 & df$groupNumbers == 1))
  
  labels1 <- c(Q1u, Q2u)
  labels2 <- c(Q1d, Q2d)
  
  getWeightGroupLocation1 <- function(df, hLine1) {
    
    df <- na.omit(df)
    quants1 <- data.table(df)[, list(quant = as.numeric(max(log))), by = groupNumbers]
    quants1 <- quants1[order(quants1$groupNumbers),]
    for (k in 1:nrow(quants1)) {
      if (quants1$quant[k] < hLine1) {
        quants1$quant[k] <- hLine1
      }
    }
    quants1$quant <- quants1$quant + 0.5
    
    return(quants1)
  }
  getWeightGroupLocation2 <- function(df, hLine2) {
    
    df <- na.omit(df)
    quants2 <- data.table(df)[, list(quant = as.numeric(min(log))), by = groupNumbers]
    quants2 <- quants2[order(quants2$groupNumbers),]
    for (k in 1:nrow(quants2)) {
      if (quants2$quant[k] > hLine2) {
        quants2$quant[k] <- hLine2
      }
    }
    
    quants2$quant <- quants2$quant - 0.5
    
    return(quants2)
  }
  
  quants1 <- getWeightGroupLocation1(df, hLine1)
  quants2 <- getWeightGroupLocation2(df, hLine2)
  
  ## Generate plot ##
  
  plot <- ggboxplot(
    df, x = "groupNumbers", y = "log",
    add = "jitter", shape = 1, 
    palette = c("orange", "blue", "green3", "black", "red"),
    fill = "groupNumbers",
    # color = "groupNumbers",
    # facet.by = "groupNumbers",
    scales = "free",
    # label="qbugName",
    outlier.shape = 1
  )
  
  # plot <- plot +
  #   xlim(0,4)
  caption.lab <- "lm(  Taxa ~ Surgery + Timepoint + Surgery:Timepoint )"

  # Add title & subtitle
  plot <- plot +
    labs(title = title.lab,
         subtitle = subtitle.lab,
         caption = caption.lab)
  
  plot <- plot +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5)
    )
  
  # Edit x- and y-axis labels
  plot <- plot +
    labs(x = paste0(modelType, " Groups"), y = "p value (log10)")
  
  # Add dashed line for significant p values increasing
  plot <- plot + geom_hline(yintercept=hLine1, linetype="dashed", color = "red")
  
  # Add dashed line for significant p values decreasing
  plot <- plot + geom_hline(yintercept=hLine2, linetype="dashed", color = "red")
  
  # Increase text size of labels
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Remove legend
  plot <- plot + theme(legend.position = "none")
  
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Add text indicating the number of significant taxa increasing over time
  plot <- plot + geom_text(data = quants1, aes(x = groupNumbers, y = quant, label = labels1), size = 5)
  
  # Add text indicating the number of significant taxa decreasing over time
  plot <- plot + geom_text(data = quants2, aes(x = groupNumbers, y = quant, label = labels2), size = 5)
  
  # Add up arrow and annotation
  plot2 <- plot + annotate("segment", 
                           x=0.45, xend=0.45,
                           y=0.25, yend=2,
                           col="black", arrow=arrow(length=unit(0.25, "cm")))
  
  plot2 <- plot2 + annotate("text", x=0.5, y=1.5, hjust = 0, label = "Increases in abundance")
  
  # Add down arrow and annotation
  plot2 <- plot2 + annotate("segment", 
                            x=0.45, xend=0.45,
                            y=-0.25, yend=-2,
                            col="black", arrow=arrow(length=unit(0.25, "cm")))
  plot2 <- plot2 + annotate("text", x=0.5, y=-1.5, hjust = 0, label = "Decreases in abundance")
  
  print(plot2)
  
} # for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] )

dev.off()


# knitr::knit_exit()


## ----Figure5-code, fig.height=5, fig.width=7, warning=FALSE, include=FALSE----
message("Figure 5\n")

outputDirNew <- paste0(outputDir, "Figure5/")
dir.create(outputDirNew, showWarnings = FALSE)

viruses <- c("Viruses_noname", "C2likevirus")

qFrame <- read.table(filePath.5, sep="\t",header = TRUE, check.names = FALSE)

qFrame$adjPValues4 <- ifelse(qFrame$adjPValues4 == 0, 2e-16,
                             qFrame$adjPValues4)

qFrame$log <- ifelse(qFrame$direction == "increase", -log10(qFrame$adjPValues4),
                     log10(qFrame$adjPValues4))
qFrame <- na.omit(qFrame)

# Label taxa as significant or not based on adjusted p-value
qFrame$Significance <- ifelse(qFrame$adjPValues4 < 0.05, "Significant",
                              "Not Significant")
qFrame <- qFrame[!( qFrame$qbugName %in% viruses ),]

timeFrames <- c("BLto1M", "BLto6M", "BLto12M", "BLto18M", "BLto24M")
filename.12 <- paste0(level, "_analysis3_pvalue_boxplots.pdf")
filePath <- paste0(outputDirNew, filename.12)
pdf(filePath, width = 7, height = 5)

for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] ) {
  
  # message(timeFrame)
  df <- qFrame[ qFrame$LengthTime == timeFrame, ]
  endPoint <- gsub("BLto", "", timeFrame)
  
  f.test <- fisher.test(df$groupNumbers, df$Significance, workspace = 2e8)
  f.pVal <- paste0("p = ", format(f.test$p.value, digits = 3))
  # f.pVal <- roundP(f.test$p.value)
  
  num_patients <- unique(df$Num_Patients)
  title.lab <- paste0( str_to_title(level),  " level taxonomic changes over time (BL to ", endPoint, ")"); title.lab
  subtitle.lab <- paste0("Fisher's Exact Test (", f.pVal, ")")
  
  hLine1 <- -log10(0.05)
  hLine2 <- log10(0.05)
  
  Q2u <- length(which(df$log > hLine1 & df$groupNumbers == 2))
  Q2d <- length(which(df$log < hLine2 & df$groupNumbers == 2))
  Q1u <- length(which(df$log > hLine1 & df$groupNumbers == 1))
  Q1d <- length(which(df$log < hLine2 & df$groupNumbers == 1))
  
  labels1 <- c(Q1u, Q2u)
  labels2 <- c(Q1d, Q2d)
  
  getWeightGroupLocation1 <- function(df, hLine1) {
    
    df <- na.omit(df)
    quants1 <- data.table(df)[, list(quant = as.numeric(max(log))), by = groupNumbers]
    quants1 <- quants1[order(quants1$groupNumbers),]
    for (k in 1:nrow(quants1)) {
      if (quants1$quant[k] < hLine1) {
        quants1$quant[k] <- hLine1
      }
    }
    quants1$quant <- quants1$quant + 0.5
    
    return(quants1)
  }
  getWeightGroupLocation2 <- function(df, hLine2) {
    
    df <- na.omit(df)
    quants2 <- data.table(df)[, list(quant = as.numeric(min(log))), by = groupNumbers]
    quants2 <- quants2[order(quants2$groupNumbers),]
    for (k in 1:nrow(quants2)) {
      if (quants2$quant[k] > hLine2) {
        quants2$quant[k] <- hLine2
      }
    }
    
    quants2$quant <- quants2$quant - 0.5
    
    return(quants2)
  }
  
  quants1 <- getWeightGroupLocation1(df, hLine1)
  quants2 <- getWeightGroupLocation2(df, hLine2)
  
  ## Generate plot ##
  
  plot <- ggboxplot(
    df, x = "groupNumbers", y = "log",
    add = "jitter", shape = 1, 
    palette = c("orange", "blue", "green3", "black", "red"),
    fill = "groupNumbers",
    # color = "groupNumbers",
    # facet.by = "groupNumbers",
    scales = "free",
    # label="qbugName",
    outlier.shape = 1
  )
  
  # plot <- plot +
  #   xlim(0,4)
  
caption.lab <- "lme(  Taxa ~ Timepoint, method = 'ML', random = ~1 | PatientID )"
# Add title & subtitle
plot <- plot +
  labs(title = title.lab,
       subtitle = subtitle.lab,
       caption = caption.lab)

plot <- plot +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5)
  )
  
  # Edit x- and y-axis labels
  plot <- plot +
    labs(x = paste0(modelType, " Groups"), y = "p value (log10)")
  
  # Add dashed line for significant p values increasing
  plot <- plot + geom_hline(yintercept=hLine1, linetype="dashed", color = "red")
  
  # Add dashed line for significant p values decreasing
  plot <- plot + geom_hline(yintercept=hLine2, linetype="dashed", color = "red")
  
  # Increase text size of labels
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Remove legend
  plot <- plot + theme(legend.position = "none")
  
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Add text indicating the number of significant taxa increasing over time
  plot <- plot + geom_text(data = quants1, aes(x = groupNumbers, y = quant, label = labels1), size = 5)
  
  # Add text indicating the number of significant taxa decreasing over time
  plot <- plot + geom_text(data = quants2, aes(x = groupNumbers, y = quant, label = labels2), size = 5)
  
  # Add up arrow and annotation
  plot2 <- plot + annotate("segment", 
                           x=0.45, xend=0.45,
                           y=0.25, yend=2,
                           col="black", arrow=arrow(length=unit(0.25, "cm")))
  
  plot2 <- plot2 + annotate("text", x=0.5, y=1.5, hjust = 0, label = "Increases in abundance")
  
  # Add down arrow and annotation
  plot2 <- plot2 + annotate("segment", 
                            x=0.45, xend=0.45,
                            y=-0.25, yend=-2,
                            col="black", arrow=arrow(length=unit(0.25, "cm")))
  plot2 <- plot2 + annotate("text", x=0.5, y=-1.5, hjust = 0, label = "Decreases in abundance")
  
  print(plot2)
  
} # for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] )

dev.off()


# knitr::knit_exit()

## ----Figure6-code, fig.height=5, fig.width=7, warning=FALSE, include=FALSE----
message("Figure 6\n")

outputDirNew <- paste0(outputDir, "Figure6/")
dir.create(outputDirNew, showWarnings = FALSE)

viruses <- c("Viruses_noname", "C2likevirus")

qFrame <- read.table(filePath.5, sep="\t",header = TRUE, check.names = FALSE)

qFrame$adjPValues3 <- ifelse(qFrame$adjPValues3 == 0, 2e-16,
                             qFrame$adjPValues3)

qFrame$log <- ifelse(qFrame$direction == "increase", -log10(qFrame$adjPValues3),
                     log10(qFrame$adjPValues3))
qFrame <- na.omit(qFrame)

# Label taxa as significant or not based on adjusted p-value
qFrame$Significance <- ifelse(qFrame$adjPValues3 < 0.05, "Significant",
                              "Not Significant")
qFrame <- qFrame[!( qFrame$qbugName %in% viruses ),]

timeFrames <- c("BLto1M", "BLto6M", "BLto12M", "BLto18M", "BLto24M")
filename.11 <- paste0(level, "_analysis4_pvalue_boxplots.pdf")
filePath <- paste0(outputDirNew, filename.11)
pdf(filePath, width = 7, height = 5)

for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] ) {
  
  # message(timeFrame)
  df <- qFrame[ qFrame$LengthTime == timeFrame, ]
  endPoint <- gsub("BLto", "", timeFrame)
  
  f.test <- fisher.test(df$groupNumbers, df$Significance, workspace = 2e8)
  f.pVal <- paste0("p = ", format(f.test$p.value, digits = 3))
  # f.pVal <- roundP(f.test$p.value)
  
  num_patients <- unique(df$Num_Patients)
  title.lab <- paste0( str_to_title(level),  " level taxonomic changes over time (BL to ", endPoint, ")"); title.lab
  subtitle.lab <- paste0("Fisher's Exact Test (", f.pVal, ")")
  
  hLine1 <- -log10(0.05)
  hLine2 <- log10(0.05)
  
  Q2u <- length(which(df$log > hLine1 & df$groupNumbers == 2))
  Q2d <- length(which(df$log < hLine2 & df$groupNumbers == 2))
  Q1u <- length(which(df$log > hLine1 & df$groupNumbers == 1))
  Q1d <- length(which(df$log < hLine2 & df$groupNumbers == 1))
  
  labels1 <- c(Q1u, Q2u)
  labels2 <- c(Q1d, Q2d)
  
  getWeightGroupLocation1 <- function(df, hLine1) {
    
    df <- na.omit(df)
    quants1 <- data.table(df)[, list(quant = as.numeric(max(log))), by = groupNumbers]
    quants1 <- quants1[order(quants1$groupNumbers),]
    for (k in 1:nrow(quants1)) {
      if (quants1$quant[k] < hLine1) {
        quants1$quant[k] <- hLine1
      }
    }
    quants1$quant <- quants1$quant + 0.5
    
    return(quants1)
  }
  getWeightGroupLocation2 <- function(df, hLine2) {
    
    df <- na.omit(df)
    quants2 <- data.table(df)[, list(quant = as.numeric(min(log))), by = groupNumbers]
    quants2 <- quants2[order(quants2$groupNumbers),]
    for (k in 1:nrow(quants2)) {
      if (quants2$quant[k] > hLine2) {
        quants2$quant[k] <- hLine2
      }
    }
    
    quants2$quant <- quants2$quant - 0.5
    
    return(quants2)
  }
  
  quants1 <- getWeightGroupLocation1(df, hLine1)
  quants2 <- getWeightGroupLocation2(df, hLine2)
  
  ## Generate plot ##
  
  plot <- ggboxplot(
    df, x = "groupNumbers", y = "log",
    add = "jitter", shape = 1, 
    palette = c("orange", "blue", "green3", "black", "red"),
    fill = "groupNumbers",
    # color = "groupNumbers",
    # facet.by = "groupNumbers",
    scales = "free",
    # label="qbugName",
    outlier.shape = 1
  )
  
  # plot <- plot +
  #   xlim(0,4)
  
caption.lab <- "lme(  Taxa ~ Surgery + Timepoint + Surgery:Timepoint, method = 'ML', random = ~1 | PatientID )"
# Add title & subtitle
plot <- plot +
  labs(title = title.lab,
       subtitle = subtitle.lab,
       caption = caption.lab)

plot <- plot +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, size = 8)
  )
  
  # Edit x- and y-axis labels
  plot <- plot +
    labs(x = paste0(modelType, " Groups"), y = "p value (log10)")
  
  # Add dashed line for significant p values increasing
  plot <- plot + geom_hline(yintercept=hLine1, linetype="dashed", color = "red")
  
  # Add dashed line for significant p values decreasing
  plot <- plot + geom_hline(yintercept=hLine2, linetype="dashed", color = "red")
  
  # Increase text size of labels
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Remove legend
  plot <- plot + theme(legend.position = "none")
  
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Add text indicating the number of significant taxa increasing over time
  plot <- plot + geom_text(data = quants1, aes(x = groupNumbers, y = quant, label = labels1), size = 5)
  
  # Add text indicating the number of significant taxa decreasing over time
  plot <- plot + geom_text(data = quants2, aes(x = groupNumbers, y = quant, label = labels2), size = 5)
  
  # Add up arrow and annotation
  plot2 <- plot + annotate("segment", 
                           x=0.45, xend=0.45,
                           y=0.25, yend=2,
                           col="black", arrow=arrow(length=unit(0.25, "cm")))
  
  plot2 <- plot2 + annotate("text", x=0.5, y=1.5, hjust = 0, label = "Increases in abundance")
  
  # Add down arrow and annotation
  plot2 <- plot2 + annotate("segment", 
                            x=0.45, xend=0.45,
                            y=-0.25, yend=-2,
                            col="black", arrow=arrow(length=unit(0.25, "cm")))
  plot2 <- plot2 + annotate("text", x=0.5, y=-1.5, hjust = 0, label = "Decreases in abundance")
  
  print(plot2)
  
} # for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] )

dev.off()


# knitr::knit_exit()

## ----LinearModel-Analyses,  include=FALSE---------------------

outputDirNew <- paste0(outputDir, "LinearModel_Analyses/")
dir.create(outputDirNew, showWarnings = FALSE)

months <- c(0, 1, 6, 12, 18, 24)

outputResults <- paste0(outputDirNew)
# dir.create(outputResults, showWarnings = FALSE)

startAbundanceIndex <- which(colnames(taxa.df)=="ResponderStatus")+1
SampleID <- taxa.df$SampleID
taxa.df2 <- taxa.df[,startAbundanceIndex:ncol(taxa.df)]
taxaNames <- colnames(taxa.df2)
taxa.df2$SampleID <- SampleID

DF <- merge(df.merged, taxa.df2, by = "SampleID")

for (analysis in c("analysis1", "analysis2", "analysis3", "analysis4")) {
  
  groupPValues_Timepoint <- vector()
  groupPValues_Surgery <- vector()
  groupPValues_Interaction <- vector()
  groupNumbers <- vector()
  qbugName <- vector()
  direction <- vector()
  LengthTime <- vector()
  qIndex <- 1

  for( i in 1:length(taxaNames)) {
    BugName <- taxaNames[i]

    for ( t in 2:length(months) ) {

      timeFrame <- months[1:t]
      log.df <- DF[ DF$time %in% timeFrame, ]
      bugCol <- log.df[ , which(colnames(log.df) == BugName) ]
      Month <- log.df$time
      Surgery <- log.df$Surgery
      PatientID <- log.df$PatientID
      Group <- log.df$class
      colEnding <- paste0("BLto", timeFrame[length(timeFrame)],"M")

      qFrameFinal <- data.frame()
      for( j in 1:ng ) {

        lm.df <- data.frame(PatientID, Month, Surgery, bugCol, Group)
        lm.df <- na.omit(lm.df)
        lm.df <- lm.df[lm.df$Group == j,]
        
        tryCatch({
          # Analysis #1
          if (analysis == "analysis1") {
            myLm <- lm(  lm.df$bugCol ~ lm.df$Month) # add surgery type, surgery*time
            p = anova(myLm)$"Pr(>F)"[1]
            groupPValues_Timepoint[qIndex] <- p
          }
          
          # Analysis #2
          if (analysis == "analysis2") {
            myLm <- lm(  lm.df$bugCol ~ lm.df$Surgery + lm.df$Month + lm.df$Surgery:lm.df$Month) # add surgery type, surgery*time
            p = anova(myLm)$"Pr(>F)"[2]
            groupPValues_Timepoint[qIndex] <- p
            p = anova(myLm)$"Pr(>F)"[1]
            groupPValues_Surgery[qIndex] <- p
            p = anova(myLm)$"Pr(>F)"[3]
            groupPValues_Interaction[qIndex] <- p
          }
          
          # Analysis #3
          if (analysis == "analysis3") {
            mlm <- lme( bugCol ~ Month, method = "REML", random = ~1 | PatientID, data = lm.df )
            fit<-anova(mlm)
            p = fit$"p-value"[2]
            groupPValues_Timepoint[qIndex] <- p
          }
          
          # Analysis #4
          if (analysis == "analysis4") {
            mlm <- lme( bugCol ~ Surgery + Month + Surgery:Month, method = "REML", random = ~1 | PatientID, data = lm.df )
            fit<-anova(mlm)
            p = fit$"p-value"[3]
            groupPValues_Timepoint[qIndex] <- p
            p = fit$"p-value"[2]
            groupPValues_Surgery[qIndex] <- p
            p = fit$"p-value"[4]
            groupPValues_Interaction[qIndex] <- p
          }
          
          groupNumbers[qIndex] <- j
          qbugName[qIndex] <- BugName
          
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
        },
        error = function(e){})

      } # for( j in 1:ng )

      qFrame <- data.frame( qbugName, LengthTime, groupNumbers, direction)
      qFrame <- data.frame( qFrame, groupPValues_Timepoint)
      qFrame$adjgroupPValues_Timepoint <- p.adjust(qFrame$groupPValues_Timepoint, method="BH")

      # Analysis #2 & #4
      if (analysis %in% c("analysis2", "analysis4")) {
        qFrame <- data.frame( qFrame, groupPValues_Surgery, groupPValues_Interaction)
        qFrame$adjgroupPValues_Surgery <- p.adjust(qFrame$groupPValues_Surgery, method="BH")
        qFrame$adjgroupPValues_Interaction <- p.adjust(qFrame$groupPValues_Interaction, method="BH")
      }

      qFrameFinal <- rbind(qFrameFinal, qFrame)

    } # for ( t in 2:length(months) )

  } # for( i in 1:length(taxaNames))

  qFrameFinal <- qFrameFinal[ order( qFrameFinal$groupPValues_Timepoint),]

  filename.5 <- paste0(level, "_", analysis, ".tsv")

  filePath <- paste0(outputResults, filename.5)
  write.table(qFrameFinal, file=filePath, sep="\t", row.names=FALSE)
  
} # for (analysis in c(analysis1, analysis2, analysis3, analysis4))


#---- Taxa ~ Timepoint (by GMM-1 group) linear model plots ----

outputDirNew <- paste0(outputDir, "Taxa_by_Timepoint/")
dir.create(outputDirNew, showWarnings = FALSE)

months <- c(0, 1, 6, 12, 18, 24)

outputResults <- paste0(outputDirNew, "ResultsTables/")
# dir.create(outputResults, showWarnings = FALSE)

startAbundanceIndex <- which(colnames(taxa.df)=="ResponderStatus")+1
SampleID <- taxa.df$SampleID
taxa.df2 <- taxa.df[,startAbundanceIndex:ncol(taxa.df)]
taxa.df2$SampleID <- SampleID

DF <- merge(df.merged, taxa.df2, by = "SampleID")


filename.6 <- paste0(level, "_by_Timepoint_groupedby_", modelType, n_class, "_LinearModelPlots.pdf")
filePath <- paste0(outputDirNew, filename.6)
pdf(filePath)

for ( t in 2:length(months) ) {

  timeFrame <- months[1:t]
  colEnding <- paste0("BLto", timeFrame[length(timeFrame)],"M")

  pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Timepoint_", colEnding))]
  n=length(pValues)
  title.lab <- paste0(colEnding, "\nMLM Timepoint p values (n = ", n, ")")
  hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))

  pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Surgery_", colEnding))]
  n=length(pValues)
  title.lab <- paste0(colEnding, "\nMLM Surgery p values (n = ", n, ")")
  hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))

  pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Observed_Weight_", colEnding))]
  n=length(pValues)
  title.lab <- paste0(colEnding, "\nMLM Observed_Weight p values (n = ", n, ")")
  hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))

  pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Modeled_Weight_", colEnding))]
  n=length(pValues)
  title.lab <- paste0(colEnding, "\nMLM Modeled_Weight p values (n = ", n, ")")
  hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
  

} # for ( t in 2:length(months) )

dev.off()


#------ Taxa ~ group (by Timepoint) linear models -----

outputDirNew <- paste0(outputDir, "Taxa_by_Group_and_Timepoint/")
dir.create(outputDirNew, showWarnings = FALSE)

months <- c(0, 1, 6, 12, 18, 24)

outputResults <- paste0(outputDirNew, "ResultsTables/")
dir.create(outputResults, showWarnings = FALSE)

startAbundanceIndex <- which(colnames(taxa.df)=="ResponderStatus")+1
SampleID <- taxa.df$SampleID
taxa.df2 <- taxa.df[,startAbundanceIndex:ncol(taxa.df)]
taxa.df2$SampleID <- SampleID

DF <- merge(df.merged, taxa.df2, by = "SampleID")

dFrame.Final <- data.frame()
for ( t in 1:length(months) ) {

  timeFrame <- months[t]
  log.df <- DF[ DF$time %in% timeFrame, ]
  Surgery <- log.df$Surgery
  Observed_Weight <- log.df$Weight_kg
  Modeled_Weight <- log.df$modeled
  PatientID <- log.df$PatientID
  Group <- log.df$class

  myT<-log.df[,startAbundanceIndex:ncol(log.df)]

  pVal_Group <- vector()
  bugName <- vector()
  index <- 1

  for (i in 1:ncol(myT)){

    bug<-myT[,i]

    if (mean(bug>0, na.rm = TRUE)>0.1){

      df <- data.frame(bug, Group)
      df <- na.omit(df)

      lm <- lm( bug ~ Group, data = df )
      fit<-anova(lm)

      pVal_Group[index] <- fit$`Pr(>F)`[1]

      bugName[index]<-colnames(myT)[i]
      Timepoint[index] <- t
      index<-index+1

    } # if (mean(bug>0, na.rm = TRUE)>0.1)

  } # for (i in 1:ncol(myT))

  dFrame<-data.frame(bugName, pVal_Group)
  dFrame$Adj_pVal_Group <- p.adjust(dFrame$pVal_Group,method = "BH")
  dFrame$Timepoint <- months[t]

  dFrame.Final <- rbind(dFrame.Final, dFrame)
} # for ( t in 1:length(months) )

filename.7 <- paste0(level, "by_", modelType, n_class, "_groupedby_Timepoint_LinearModelResults.tsv")
filePath <- paste0(outputResults, filename.7)
write.table(dFrame.Final, filePath, sep="\t",quote = FALSE, row.names = FALSE)


#----- Taxa ~ group (by Timepoint) linear model plots -----
tukey.df <- data.frame()
plotList <- list()
index <- 1
for (i in unique(dFrame.Final$bugName)) {
  
  dFrame2 <- dFrame.Final[dFrame.Final$bugName == i,]
  # DF2 <- DF[ !is.na( DF[which(colnames(DF) == divMonth)]), ]
  # DF2 <- DF2[ !is.na( DF2[which(colnames(DF2) == i)]), ]
  DF2 <- DF[DF$time %in% months,]
  
  Group <- DF2$class
  Time <- DF2$time
  bug <- DF2[, which(colnames(DF2) == i)]
  df2 <- data.frame(Group, Time, bug)
  df2 <- na.omit(df2)
  x.lab <- paste0(ng, " Class ", modelType, " Groups")
  y.lab <- i
  
  colors <- c("orange", "blue", "green3", "black", "red")
  colorPalette <- colors[1:n_class]
  
  plot <- ggboxplot(
    df2, x = "Group", y = "bug", color = "black",
    fill = "Group", palette = colorPalette,
    # facet.by = "month",
    scales = "free", add = "jitter"
  ); plot
  
  month.labsP <- paste0(month.labs, "\n", roundP(dFrame2$Adj_pVal_Group[1:nrow(dFrame2)]))
  plot <- facet(plot, facet.by = "Time"
                , panel.labs = list(Time = month.labsP[1:length(months)])
  ); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab)+ theme(legend.position = "none"); plot
  
  plotList[[index]] <- plot
  index <- index + 1
}

filename.8 <- paste0(level, "_by_", modelType, n_class, "_groupedby_Timepoint_LinearModel_BoxPlots.pdf")
filePath <- paste0(outputDirNew, filename.8)
pdf(filePath, width = 10, height = 7)
for (x in 1:length(plotList)) {
  print(plotList[[x]])
}
dev.off()


## ----diversity-setup------------------------------------------

outputDirNew <- paste0(outputDir, "Diversity/")
dir.create(outputDirNew, showWarnings = FALSE)

diversity.df <- taxa.df[which(colnames(taxa.df) %in% c("Bray_Uniqueness", "Kendall_Uniqueness", "ShannonIndex", "SimpsonIndex", "Richness", "Evenness", "SampleID"))]
DF <- merge(df.merged, diversity.df, by = "SampleID")
outputDirResults <- paste0(outputDirNew, "ResultsTables/")
dir.create(outputDirResults, showWarnings = FALSE)

##### Bray-Curtis ~ Weight Group & Timepoint PCoA #####

taxaStart <- which(colnames(main.DF) == "class") + 1

taxaTable <- main.DF[,taxaStart:ncol(main.DF)]
taxaTable$class <- main.DF$class
taxaTable$class <- paste0("Group ", taxaTable$class)
taxaTable$class <- as.factor(taxaTable$class)
taxaTable$Timepoint <- main.DF$time
taxaTable$Timepoint <- as.factor(taxaTable$Timepoint)

taxaTable <- na.omit(taxaTable)
end <- which(colnames(taxaTable) == "class") -1

Month <- vector()
p <- vector()
R2 <- vector()
Fval <- vector()

file.path <- paste0(outputDirNew, level, "_BrayCurtis_by_class_by_Timepoint_PCoA.pdf")
pdf(file.path,width = 10, height = 5)
par(mar=c(5, 5, 5, 5), xpd=FALSE, mfrow = c(2,3))

for (month in unique(taxaTable$Timepoint)) {
  
  taxaTable2 <- taxaTable[taxaTable$Timepoint == month,]
  myMDS <- capscale(taxaTable2[,1:end] ~ 1, distance="bray")
  percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
  pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
  
  
  adon.results<-adonis2(taxaTable2[,1:end] ~ taxaTable2$class, method="bray",perm=999)
  # print(adon.results)
  # capture.output(adon.results, file = paste0(outputDir, level, "_PERMANOVA_class.txt"))
  
  Month <- c(Month, month)
  p <- c(p, adon.results$`Pr(>F)`[1])
  R2 <- c(R2, adon.results$R2[1])
  Fval <- c(Fval, adon.results$F[1])
  
  Title <- paste0("PERMANOVA p = ", adon.results$`Pr(>F)`[1])
  Subtitle <- month.labs[which(unique(taxaTable$Timepoint) == month)]
  
  
  all.combindations <- combn(1:2, 2)
  for (combination in 1:ncol(all.combindations)) {
    PCoA_a <- all.combindations[1,combination]
    PCoA_b <- all.combindations[2,combination]
    
    pcoaPlot <- ordiplot(myMDS, choices = c(PCoA_a,PCoA_b), display = "sites", type = "none", cex.lab = 2,
                         xlab = pcoa_p[PCoA_a], ylab = pcoa_p[PCoA_b], main = paste0(Title, "\n", Subtitle))
    col2=adjustcolor(colorPalette[factor(taxaTable2$class)], alpha.f = 1)
    # points(pcoaPlot, "sites", col = taxaTable2$class, pch = 16, cex = 1.5)
    points(pcoaPlot, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
    
    for (n in length(levels(taxaTable2$class)):1) {
      ordiellipse(pcoaPlot, taxaTable2$class, kind = "se", conf = 0.95, lwd = 2, draw = "lines",
                  col = colorPalette[n], show.groups = levels(factor(taxaTable2$class))[n], label = F,
                  font = 2, cex = 1)
    }
    legend("topright", 
           legend = levels(taxaTable2$class),
           col = colorPalette[length(levels(taxaTable2$class)):1],
           cex = 0.75, pch = 16,
           horiz = FALSE)
    
  }
  
  
  
}

dev.off()

dFrame <- data.frame(Month, p, R2, Fval)
file.path <- paste0(outputDirResults, level, "_BrayCurtis_by_class_by_Timepoint_PERMANOVAResults.tsv")
write.table(dFrame, file.path, sep="\t",quote = FALSE, row.names = FALSE)



## Bray-Curtis Uniqueness ~ Latent Class (& Time)-------------------------------------------------------------
x.lab <- paste0(modelType, " Class")
y.lab <- "Bray-Curtis Uniqueness"
title.lab <- str_to_title(level)
subtitle.lab <- ""

DF <- DF[!is.na(DF$class),]
DF$class <- factor(DF$class)
DF$time <- factor(DF$time, levels = months)

stat.test <- DF %>%
  group_by(time) %>%
  wilcox_test(Bray_Uniqueness ~ class) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

file.path <- paste0(outputDirResults, level, "_BrayCurtis_Uniqueness_by_Timepoint_by_", modelType, n_class, "_Wilcox_Results.tsv")
write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)

stat.test <- stat.test %>%
  add_xy_position(x = "class", fun = "max")

plot <- ggboxplot(
  DF, x = "class", y = "Bray_Uniqueness", color = "black",
  fill = "class", palette = Palette,
  # facet.by = "month",
  add = "jitter",
  scales = "free"
)

plot <- facet(plot, facet.by = "time")

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

fileName <-   paste0(level, "_BrayCurtis_Uniqueness_by_Timepoint_by_", modelType, n_class, "_Wilcox_BoxPlot.pdf")
file.path <- paste0(outputDirNew, fileName)
pdf(file.path, width = 10, height = 7)
print(plot)
dev.off()


#' 
## Kendall Uniqueness ~ Latent Class (& Time)-------------------------------------------------------------
x.lab <- paste0(modelType, " Class")
y.lab <- "Kendall Uniqueness"
title.lab <- str_to_title(level)
subtitle.lab <- ""

DF <- DF[!is.na(DF$class),]
DF$class <- factor(DF$class)
DF$time <- factor(DF$time, levels = months)

stat.test <- DF %>%
  group_by(time) %>%
  wilcox_test(Kendall_Uniqueness ~ class) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

file.path <- paste0(outputDirResults, level, "_Kendall_Uniqueness_by_Timepoint_by_", modelType, n_class, "_Wilcox_Results.tsv")
write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)

stat.test <- stat.test %>%
  add_xy_position(x = "class", fun = "max")

plot <- ggboxplot(
  DF, x = "class", y = "Kendall_Uniqueness", color = "black",
  fill = "class", palette = Palette,
  # facet.by = "month",
  add = "jitter",
  scales = "free"
)

plot <- facet(plot, facet.by = "time")

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

fileName <-   paste0(level, "_Kendall_Uniqueness_by_Timepoint_by_", modelType, n_class, "_Wilcox_BoxPlot.pdf")
file.path <- paste0(outputDirNew, fileName)
pdf(file.path, width = 10, height = 7)
print(plot)
dev.off()


## Shannon Index ~ Latent Class (& Time)-------------------------------------------------------------
x.lab <- paste0(modelType, " Class")
y.lab <- "Shannon Index"
title.lab <- str_to_title(level)
subtitle.lab <- ""

DF <- DF[!is.na(DF$class),]
DF$class <- factor(DF$class)
DF$time <- factor(DF$time, levels = months)

stat.test <- DF %>%
  group_by(time) %>%
  wilcox_test(ShannonIndex ~ class) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

file.path <- paste0(outputDirResults, level, "_ShannonIndex_by_Timepoint_by_", modelType, n_class, "_Wilcox_Results.tsv")
write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)

stat.test <- stat.test %>%
  add_xy_position(x = "class", fun = "max")

plot <- ggboxplot(
  DF, x = "class", y = "ShannonIndex", color = "black",
  fill = "class", palette = Palette,
  # facet.by = "month",
  add = "jitter",
  scales = "free"
)

plot <- facet(plot, facet.by = "time")

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

fileName <-   paste0(level, "_ShannonIndex_by_Timepoint_by_", modelType, n_class, "_Wilcox_BoxPlot.pdf")
file.path <- paste0(outputDirNew, fileName)
pdf(file.path, width = 10, height = 7)
print(plot)
dev.off()


## SimpsonIndex ~ Latent Class (& Time)-------------------------------------------------------------
x.lab <- paste0(modelType, " Class")
y.lab <- "SimpsonIndex"
title.lab <- str_to_title(level)
subtitle.lab <- ""

DF <- DF[!is.na(DF$class),]
DF$class <- factor(DF$class)
DF$time <- factor(DF$time, levels = months)

stat.test <- DF %>%
  group_by(time) %>%
  wilcox_test(SimpsonIndex ~ class) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

file.path <- paste0(outputDirResults, level, "_SimpsonIndex_by_Timepoint_by_", modelType, n_class, "_Wilcox_Results.tsv")
write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)

stat.test <- stat.test %>%
  add_xy_position(x = "class", fun = "max")

plot <- ggboxplot(
  DF, x = "class", y = "SimpsonIndex", color = "black",
  fill = "class", palette = Palette,
  # facet.by = "month",
  add = "jitter",
  scales = "free"
)

plot <- facet(plot, facet.by = "time")

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

fileName <-   paste0(level, "_SimpsonIndex_by_Timepoint_by_", modelType, n_class, "_Wilcox_BoxPlot.pdf")
file.path <- paste0(outputDirNew, fileName)
pdf(file.path, width = 10, height = 7)
print(plot)
dev.off()


## Richness ~ Latent Class (& Time)-------------------------------------------------------------
x.lab <- paste0(modelType, " Class")
y.lab <- "Richness"
title.lab <- str_to_title(level)
subtitle.lab <- ""

DF <- DF[!is.na(DF$class),]
DF$class <- factor(DF$class)
DF$time <- factor(DF$time, levels = months)

stat.test <- DF %>%
  group_by(time) %>%
  wilcox_test(Richness ~ class) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

file.path <- paste0(outputDirResults, level, "_Richness_by_Timepoint_by_", modelType, n_class, "_Wilcox_Results.tsv")
write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)

stat.test <- stat.test %>%
  add_xy_position(x = "class", fun = "max")

plot <- ggboxplot(
  DF, x = "class", y = "Richness", color = "black",
  fill = "class", palette = Palette,
  # facet.by = "month",
  add = "jitter",
  scales = "free"
)

plot <- facet(plot, facet.by = "time")

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

fileName <-   paste0(level, "_Richness_by_Timepoint_by_", modelType, n_class, "_Wilcox_BoxPlot.pdf")
file.path <- paste0(outputDirNew, fileName)
pdf(file.path, width = 10, height = 7)
print(plot)
dev.off()


## Evenness ~ Latent Class (& Time)-------------------------------------------------------------
x.lab <- paste0(modelType, " Class")
y.lab <- "Evenness"
title.lab <- str_to_title(level)
subtitle.lab <- ""

DF <- DF[!is.na(DF$class),]
DF$class <- factor(DF$class)
DF$time <- factor(DF$time, levels = months)

stat.test <- DF %>%
  group_by(time) %>%
  wilcox_test(Evenness ~ class) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

file.path <- paste0(outputDirResults, level, "_Evenness_by_Timepoint_by_", modelType, n_class, "_Wilcox_Results.tsv")
write.table(stat.test, file.path, sep="\t",quote = FALSE, row.names = FALSE)

stat.test <- stat.test %>%
  add_xy_position(x = "class", fun = "max")

plot <- ggboxplot(
  DF, x = "class", y = "Evenness", color = "black",
  fill = "class", palette = Palette,
  # facet.by = "month",
  add = "jitter",
  scales = "free"
)

plot <- facet(plot, facet.by = "time")

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

fileName <-   paste0(level, "_Evenness_by_Timepoint_by_", modelType, n_class, "_Wilcox_BoxPlot.pdf")
file.path <- paste0(outputDirNew, fileName)
pdf(file.path, width = 10, height = 7)
print(plot)
dev.off()








