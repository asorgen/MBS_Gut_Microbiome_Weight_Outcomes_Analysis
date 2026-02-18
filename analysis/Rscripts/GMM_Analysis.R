#' ---
#' title: "Growth Mixture Models"
#' author: "Alicia Sorgen"
#' date: December 7, 2023
#' ---

## ----Setup, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo=FALSE)
rm(list=ls())
set.seed(1989)


## ----Library, include=FALSE-----------------------------------
R <- sessionInfo()
message(R$R.version$version.string)

library(MASS); message("MASS: Version ", packageVersion("MASS"))
library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
library(lcmm); message("lcmm: Version ", packageVersion("lcmm"))
library(stringr); message("stringr: Version ", packageVersion("stringr"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
# library(scales); message("scales: Version ", packageVersion("scales"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))
library(dplyr); message("dplyr: Version ", packageVersion("dplyr"))


## ----Script-Edits---------------------------------------------
ANALYSIS <- "MetaPhlAn2_microbiome"

moduleRoot <- paste0("GMM_Analysis_")
params <- vector()
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
params <- c(params, "MetaPhlAn2")
# params <- c(params, "Kraken2")
# params <- c(params, 0)
# params <- c(params, 1)
# params <- c(params, 6)
# params <- c(params, 12)
# params <- c(params, 18)
params <- c(params, 24)
params <- c(params, "LCGA")
params <- c(params, "Multivariate")
# params <- c(params, "GMM-1")
# params <- c(params, "GMM-2")

end_month <- params[length(params)]


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
  pipeline <- paste0(ANALYSIS,"_analysis_")
  
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
  module <- moduleRoot
  module <- paste0(module, args[5], "_",gsub("-", "", args[4]))
  
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
  
  outputDir <- paste0(moduleDir, "output/")
  if (any(dir(moduleDir) == "output") == FALSE) {
    dir.create(outputDir, showWarnings = FALSE)
  }
  
  files <- list.files(outputDir, recursive = TRUE, full.names = TRUE)
  file.remove(files)
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
rm(params, moduleRoot, ANALYSIS)

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
module <- sapply(strsplit(moduleDir, "/"), "[", length(strsplit(moduleDir, "/")[[1]]))
moduleNum <- as.numeric(sapply(strsplit(module, "_"), "[", 1))

#' 
## ----Functions------------------------------------------------
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]
rm(str)

## ----Input----------------------------------------------------
prevModule <- str_subset(dir(pipeRoot), "WeightMetaMerge")
inputDir = paste0(pipeRoot,"/",prevModule,"/output/")
inputFile <- "metadata.tsv"

classifier <- args[2]
endM <- args[3]
modelType <- args[4]
lmeMethod <- args[5]

prevModule2 <- paste0("TaxaMetaMerge")
inputDir2 = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule2),"/output/")
logCountFile <- paste0(level, "_LogNormalizedCounts_", classifier, ".tsv")


## ----Output----------------------------------------------------
outputDir = file.path(moduleDir,"output/")

## ----Script-Variables-----------------------------------------
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

outcomes <- c("Percent_Loss_kg", "Weight_kg")
outcomeLabel <- c("Weight Loss from BL (%)", "Weight_kg")
repeatedMeasure <- c("time")
covariates <- c("Surgery")
subjects <- c("SampleID", "PatientID")
Palette <- c("orange", "blue", "green3", "black", "red")

## ----Read-Table, echo=FALSE-----------------------------------

# Read in table
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

file.path <- paste0(inputDir2, logCountFile)
taxa.df <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

# Trim variables
df.1 <- myTable[ , which(colnames(myTable) %in% c(outcomes, repeatedMeasure, covariates, subjects))]

## ----Data-Summary---------------------------------------------
n <- nrow(df.1)
n.subjects <- length(unique(df.1$PatientID))

# Percent weight
sum.stats <- df.1 %>%
  group_by(time) %>%
  get_summary_stats(Percent_Loss_kg, type = "mean_sd")

# knitr::kable(sum.stats[, c(1,3)], "pipe", caption = "RYGB Patients at each timepoint", booktabs = TRUE)

# outputName <- paste0("Percent_Loss_kg_Summary.tsv")
# outputPath <- file.path(outputDir, outputName)
# write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)


# Weight
sum.stats <- df.1 %>%
  group_by(time, Surgery) %>%
  get_summary_stats(Weight_kg, type = "mean_sd")

# outputName <- paste0("Weight_kg_Summary.tsv")
# outputPath <- file.path(outputDir, outputName)
# write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
weight.table <- sum.stats[,c("time", "Surgery", "n")]

## ----Figure1, echo=FALSE, fig.height=5, fig.width=10, warning=FALSE----

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
Missing_Data_Model_Formula <- paste0(model$call[2])
Missing_Data_Model_Data <- paste0(model$call[3])
Missing_Data_Model_Type <- paste0(model$call[5])

df.2$Baseline_kg <- assignAllbyTerm(df.2, "time", "0", "Weight_kg", "PatientID")

df.2$`Observed - BL` <- df.2$Weight_kg - df.2$Baseline_kg
df.2$`Modeled - BL` <- df.2$modeled - df.2$Baseline_kg

df.2$`(Observed - BL)/BL` <- df.2$`Observed - BL` / df.2$Baseline_kg
df.2$`(Modeled - BL)/BL` <- df.2$`Modeled - BL` / df.2$Baseline_kg

df.2$observed_weight_change <- df.2$`(Observed - BL)/BL` * 100
df.2$modeled_weight_change <- df.2$`(Modeled - BL)/BL` * 100

write.table(df.2, "~/Downloads/output.tsv",  sep="\t", quote = FALSE, row.names = FALSE)
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
  theme_bw()

filename.1 <- "MLM_Modeled_and_Observed_Weight_Change_by_Surgery.pdf"
filePath <- paste0(outputDir, filename.1)

# pdf(filePath, width = 10, height = 5)
# print(plot)
# dev.off()

# df.4 <- df.2[df.2$time != 0,]
# 
# stat.test <- df.4 %>%
#   group_by(time) %>%
#   wilcox_test(observed_weight_change ~ Surgery)
# 
# outputName <- paste0("observed_weight_change_Surgery_Wilcoxon.tsv")
# outputPath <- file.path(outputDir, outputName)
# write.table(stat.test, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
# 
# sum.stats <- df.4 %>%
#   group_by(time, Surgery) %>%
#   get_summary_stats(observed_weight_change, type = "full")
# 
# outputName <- paste0("observed_weight_change_Surgery_desc_statistics.tsv")
# outputPath <- file.path(outputDir, outputName)
# write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)



# stat.test <- df.4 %>%
#   group_by(time) %>%
#   wilcox_test(modeled_weight_change ~ Surgery)
# 
# outputName <- paste0("modeled_weight_change_Surgery_Wilcoxon.tsv")
# outputPath <- file.path(outputDir, outputName)
# write.table(stat.test, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
# 
# sum.stats <- df.4 %>%
#   group_by(time, Surgery) %>%
#   get_summary_stats(modeled_weight_change, type = "full")
# 
# outputName <- paste0("modeled_weight_change_Surgery_desc_statistics.tsv")
# outputPath <- file.path(outputDir, outputName)
# write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)

## ----Modeling, echo=TRUE-----------------------------------------
if (modelType == "LCGA") {
  if (lmeMethod == "Univariate") {
    allModelResults <- data.frame()
  } else {
    prevModule <- str_subset(dir(pipeRoot), paste0(moduleNum-1, "_GMM_Analysis"))
    inputDir = paste0(pipeRoot,"/",prevModule,"/output/")
    inputFile <- "ALL_Summary.tsv"
    file.path <- paste0(inputDir, inputFile)
    allModelResults <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
  }
} 
if (modelType == "GMM-1") {
  prevModule <- str_subset(dir(pipeRoot), paste0(moduleNum-1, "_GMM_Analysis"))
  inputDir = paste0(pipeRoot,"/",prevModule,"/output/")
  inputFile <- "ALL_Summary.tsv"
  file.path <- paste0(inputDir, inputFile)
  allModelResults <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
  
}
if (modelType == "GMM-2") {
  prevModule <- str_subset(dir(pipeRoot), paste0(moduleNum-1, "_GMM_Analysis"))
  inputDir = paste0(pipeRoot,"/",prevModule,"/output/")
  inputFile <- "ALL_Summary.tsv"
  file.path <- paste0(inputDir, inputFile)
  allModelResults <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
  
}

for (model.v in 1:2) {
  
  outputResults <- paste0(outputDir, "Model", model.v, "/")
  dir.create(outputResults, showWarnings = FALSE)
  
  surgTypes <- "all"
  
  plotList <- list()
  gmmSummary <- data.frame()
  index <- 1
  
  for (surg in surgTypes) {
    
    if (surg == "all") {
      df.2a <- df.2
      hlme_Data <- "df.2"
      
    } else {
      df.2a <- df.2[df.2$Surgery == surg,]
      hlme_Data <- paste0("df.", surg)
      
    }
    df.2a$PatientNum <- as.numeric(factor(df.2a$PatientID))
    
    if (model.v %in% c(1)) {
      
      if (modelType == "LCGA") {
        model1 <- hlme( data = df.2a,
                        modeled_weight_change ~ time + Surgery,
                        subject = "PatientNum",
                        ng = 1
        )
      }
      
      if (modelType == "GMM-1") {
        model1 <- hlme( data = df.2a,
                        modeled_weight_change ~ time + Surgery,
                        subject = "PatientNum",
                        random = ~1,  # GMM-1
                        ng = 1
        )
      }
      
      if (modelType == "GMM-2") {
        model1 <- hlme( data = df.2a,
                        modeled_weight_change ~ time + Surgery,
                        subject = "PatientNum",
                        random = ~1 + time, # GMM-2
                        ng = 1
        )
      }
      
      hlme_Formula <- paste0(model1$call[2])
    }
    if (model.v %in% c(2)) {
      if (modelType == "LCGA") {
        model1 <- hlme( data = df.2a,
                        modeled_weight_change ~ time,
                        subject = "PatientNum",
                        ng = 1
        )
      }
      
      if (modelType == "GMM-1") {
        model1 <- hlme( data = df.2a,
                        modeled_weight_change ~ time,
                        subject = "PatientNum",
                        random = ~1,  # GMM-1
                        ng = 1
        )
      }
      
      if (modelType == "GMM-2") {
        model1 <- hlme( data = df.2a,
                        modeled_weight_change ~ time,
                        subject = "PatientNum",
                        random = ~1 + time, # GMM-2
                        ng = 1
        )
      }
      hlme_Formula <- paste0(model1$call[2])
    }
    
    sumTable <- as.data.frame( summarytable( model1))
    rownames(sumTable)[1] <- paste0("model1_", surg)
    
    n_classes <- c(2,3,4,5)
    for (n_class in n_classes) {
      
      outputClass <- paste0(outputResults, "Class_", n_class, "/")
      dir.create(outputClass, showWarnings = FALSE)
      
      if (model.v %in% c(1)) {
        if (modelType == "LCGA") {
          modelx <- gridsearch( rep = 100,
                                maxiter = 10,
                                minit = model1,
                                hlme( data = df.2a,
                                      modeled_weight_change ~ time + Surgery,
                                      var.time = "time",
                                      subject = "PatientNum",
                                      ng = n_class,
                                      mixture = ~ time
                                )
          )
          
        }
        
        if (modelType == "GMM-1") {
          modelx <- gridsearch( rep = 100,
                                maxiter = 10,
                                minit = model1,
                                hlme( data = df.2a,
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
          modelx <- gridsearch( rep = 100,
                                maxiter = 10,
                                minit = model1,
                                hlme( data = df.2a,
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
      if (model.v %in% c(2)) {
        if (modelType == "LCGA") {
          modelx <- gridsearch( rep = 100,
                                maxiter = 10,
                                minit = model1,
                                hlme( data = df.2a,
                                      modeled_weight_change ~ time,
                                      var.time = "time",
                                      subject = "PatientNum",
                                      ng = n_class,
                                      mixture = ~ time
                                )
          )
          
        }
        
        if (modelType == "GMM-1") {
          modelx <- gridsearch( rep = 100,
                                maxiter = 10,
                                minit = model1,
                                hlme( data = df.2a,
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
          modelx <- gridsearch( rep = 100,
                                maxiter = 10,
                                minit = model1,
                                hlme( data = df.2a,
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
      
      sumTable <- bind_rows(sumTable, smry)
      rownames(sumTable)[n_class] <- paste0("gmm", n_class, "_", surg)
      
      ng <- unique(class.df$class)
      
      df.3 <- merge(df.2a, class.df, by = "PatientNum", all.x = TRUE)
      # predic <- modelx[["pred"]]
      # df.3 <- merge(df.3, predic, by = c("PatientNum", "time"), all.x = TRUE)
      
      plot <- ggplot() +
        xlim(0, 36)+
        labs(x = "Follow-up time (months)", y = "Weight Change (%)",
             title = paste0(modelType, " ", n_class, "-class model (BIC=", BIC_value, ") - ", surg, " patients")) +
        theme_bw()
      
      df.merged <- data.frame()
      for (i in ng) {
        
        df.C <- df.3[df.3$class == i,]
        label <- paste0("Group ", i, " (n = ", length(unique(df.C$PatientID)), ", ", round(smry[,x+i],1), "%)")
        
        # # Find baseline weight for each PatientID at time 0
        # df.C$Baseline_kg <- assignAllbyTerm(df.C, "time", "0", "Weight_kg", "PatientID")
        # 
        # df.C$`Observed - BL` <- df.C$Weight_kg - df.C$Baseline_kg
        # df.C$`Modeled - BL` <- df.C$modeled - df.C$Baseline_kg
        # 
        # df.C$`(Observed - BL)/BL` <- df.C$`Observed - BL` / df.C$Baseline_kg
        # df.C$`(Modeled - BL)/BL` <- df.C$`Modeled - BL` / df.C$Baseline_kg
        # 
        # df.C$observed_weight_change <- df.C$`(Observed - BL)/BL` * 100
        # df.C$modeled_weight_change <- df.C$`(Modeled - BL)/BL` * 100
        df.merged <- rbind(df.merged, df.C)
        
        df.C.O <- df.C %>%
          group_by(time) %>%
          get_summary_stats(observed_weight_change, type = "full")
        
        df.C.M <- df.C %>%
          group_by(time) %>%
          get_summary_stats(modeled_weight_change, type = "full")
        
        # Generate the line graph using ggplot2
        plot <- plot + 
          geom_point(data = df.C.O, aes(x = time, y = median), color = Palette[i], shape = 15, size = 3) +
          geom_errorbar(data=df.C.O, mapping=aes(x=time, ymin=q1, ymax=q3), width=0.2, linewidth=0.5, color=Palette[i]) + 
          geom_line(data = df.C.M, aes(x = time, y = mean), color = Palette[i], size=1) +
          annotate("text", x = 25, y = df.C.M$mean[nrow(df.C.M)], hjust = 0, label = label, color = Palette[i])
        
        
      } # for (i in 1:ng) 
      plotList[[index]] <- plot
      index <- index + 1
      
      # df.4 <- df.merged[df.merged$time != 0,]
      # 
      # stat.test <- df.4 %>%
      #   group_by(time) %>%
      #   wilcox_test(observed_weight_change ~ class)
      # 
      # outputName <- paste0("observed_weight_change_", n_class, "_Class_Wilcoxon_", model.v, ".tsv")
      # outputPath <- file.path(outputClass, outputName)
      # write.table(stat.test, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
      # 
      # sum.stats <- df.4 %>%
      #   group_by(time, class) %>%
      #   get_summary_stats(observed_weight_change, type = "full")
      # 
      # outputName <- paste0("observed_weight_change_", n_class, "_Class_desc_statitistics_", model.v, ".tsv")
      # outputPath <- file.path(outputClass, outputName)
      # write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
      
      
      
      # stat.test <- df.4 %>%
      #   group_by(time) %>%
      #   wilcox_test(modeled_weight_change ~ class)
      # 
      # outputName <- paste0("modeled_weight_change_", n_class, "_Class_Wilcoxon_", model.v, ".tsv")
      # outputPath <- file.path(outputClass, outputName)
      # write.table(stat.test, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
      # 
      # sum.stats <- df.4 %>%
      #   group_by(time, class) %>%
      #   get_summary_stats(modeled_weight_change, type = "full")
      # 
      # outputName <- paste0("modeled_weight_change_", n_class, "_Class_desc_statitistics_", model.v, ".tsv")
      # outputPath <- file.path(outputClass, outputName)
      # write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
      
    } # for (n_class in n_classes)
    
    gmmSummary <- rbind(gmmSummary, sumTable)
    
    
  } # for (surg in surgTypes)
  
  gmmSummary$lme_Formula <- Missing_Data_Model_Formula
  gmmSummary$lme_df <- Missing_Data_Model_Data
  gmmSummary$lme_model <- Missing_Data_Model_Type
  gmmSummary$hlme_Formula <- hlme_Formula
  gmmSummary$hlme_df <- hlme_Data
  for (r in 1:nrow(gmmSummary)) {
    gmmSummary$ClassSize_1[r] <- ifelse(min(gmmSummary[r,which(colnames(gmmSummary) %like% "%class")], na.rm = TRUE) > 1, "Pass", "Fail")
    gmmSummary$ClassSize_5[r] <- ifelse(min(gmmSummary[r,which(colnames(gmmSummary) %like% "%class")], na.rm = TRUE) > 5, "Pass", "Fail")
    
  }
  
  
  filename.2 <- paste0(modelType, "_Modeled_and_Observed_Weight_Change_", model.v, ".pdf")
  filePath <- paste0(outputResults, filename.2)
  
  pdf(filePath, width = 10, height = 5)
  for (p in 1:length(plotList)) {
    print(plotList[[p]])
  }
  dev.off()
  
  outputName <- paste0(modelType, "_Summary_", model.v, ".tsv")
  outputPath <- file.path(outputResults, outputName)
  Surgery <- sapply(strsplit(rownames(gmmSummary), "_"), "[", 2)
  gmmSummary <- cbind(Surgery, gmmSummary)
  write.table(gmmSummary, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
  
  gmmSummary$ModelType <- modelType
  allModelResults <- rbind(allModelResults, gmmSummary)
  


}

outputName <- paste0("ALL_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(allModelResults, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
