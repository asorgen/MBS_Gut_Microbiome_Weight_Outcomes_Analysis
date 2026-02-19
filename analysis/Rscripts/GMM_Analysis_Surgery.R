## ----include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo=FALSE)
set.seed(1989)
rm(list=ls())

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
# library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))
library(dplyr); message("dplyr: Version ", packageVersion("dplyr"))



## ----Script-Edits---------------------------------------------
ANALYSIS <- "MetaPhlAn2_microbiome"

moduleRoot <- paste0("GMM_Analysis_Surgery")
included <- c(  0
                , 1
                , 6
                , 12
                , 18
                , 24
)
params <- vector()
params <- c(params, "~/git/gut-microbiome-bariatric-weight-outcomes")
params <- c(params, "RYGB")


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
  module <- gsub("Surgery", args[2], moduleRoot)
  
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

surg <- args[2]
## ----Read-Table, echo=FALSE-----------------------------------

# Read in table
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

# Trim variables
df.1 <- myTable[ , which(colnames(myTable) %in% c(outcomes, repeatedMeasure, covariates, subjects))]

# Weight
sum.stats <- df.1 %>%
  group_by(time, Surgery) %>%
  get_summary_stats(Weight_kg, type = "mean_sd")

outputName <- paste0("Weight_kg_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(sum.stats, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
weight.table <- sum.stats[,c("time", "Surgery", "n")]

# Trim variables
df.1 <- myTable[ , which(colnames(myTable) %in% c( "Weight_kg", repeatedMeasure, covariates, subjects))]

## ----Modeling, echo=TRUE-----------------------------------------
modelTypes <- c("LCGA", "GMM-1", "GMM-2")

allResults <- data.frame()

for (modelType in modelTypes) {
  
  plotList <- list()
  gmmSummary <- data.frame()
  index <- 1
  
  # Filter dataset by surgery type
  df.1a <- df.1[df.1$Surgery == surg,]
  
  # Set patient ID as number b/c that's how hlme takes it
  df.1a$PatientNum <- as.numeric(factor(df.1a$PatientID))
  
  # Perform mixed linear model
  model <- lme( Weight_kg ~ time, method = "ML", random = ~1 | PatientID, data = df.1a, na.action = na.omit )
  # model$sigma # check standard error
  
  ## RYGB ML = 12.21379
  ## SG ML = 11.02674
  
  # # Perform random effects mixed linear model
  # model <- lme( Weight_kg ~ time, method = "REML", random = ~1 | PatientID, data = df.1a, na.action = na.omit )
  # model$sigma # check standard error
  
  ## RYGB REML = 12.22994
  ## SG REML = 11.0592
  
  # Extract predicted values 
  df.1a$predicted_values <- predict(model, newdata = data.frame(time = df.1a$time, PatientID = df.1a$PatientID))
  df.1a$modeled <- ifelse(is.na(df.1a$Weight_kg) == TRUE, df.1a$predicted_values,
                          df.1a$Weight_kg)
  
  # Calculate observed and modeled % weight change
  df.1a$Baseline_kg <- assignAllbyTerm(df.1a, "time", "0", "Weight_kg", "PatientID")
  
  df.1a$`Observed - BL` <- df.1a$Weight_kg - df.1a$Baseline_kg
  df.1a$`Modeled - BL` <- df.1a$modeled - df.1a$Baseline_kg
  
  df.1a$`(Observed - BL)/BL` <- df.1a$`Observed - BL` / df.1a$Baseline_kg
  df.1a$`(Modeled - BL)/BL` <- df.1a$`Modeled - BL` / df.1a$Baseline_kg
  
  df.1a$observed_weight_change <- df.1a$`(Observed - BL)/BL` * 100
  df.1a$modeled_weight_change <- df.1a$`(Modeled - BL)/BL` * 100
  
  
  # Perform 1-class LCGA & GMMs
  if (modelType == "LCGA") {
    model1 <- hlme( data = df.1a,
                    modeled_weight_change ~ time,
                    subject = "PatientNum",
                    ng = 1
    )
  }
  
  if (modelType == "GMM-1") {
    model1 <- hlme( data = df.1a,
                    modeled_weight_change ~ time,
                    subject = "PatientNum",
                    random = ~1,  # GMM-1
                    ng = 1
    )
  }
  
  if (modelType == "GMM-2") {
    model1 <- hlme( data = df.1a,
                    modeled_weight_change ~ time,
                    subject = "PatientNum",
                    random = ~1 + time, # GMM-2
                    ng = 1
    )
  }
  
  sumTable <- as.data.frame( summarytable( model1))
  rownames(sumTable)[1] <- paste0(modelType, "_1")
  
  # Perform 2-class, 3-class, 4-class, 5-class LCGA & GMMs
  n_classes <- c(2,3,4,5)
  for (n_class in n_classes) {
    if (modelType == "LCGA") {
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
    
    smry <- as.data.frame( summarytable( modelx))
    
    class.df <- modelx$pprob
    x <- which(colnames(smry) == "BIC")
    BIC_value <- round(smry[1,x], 2)
    
    sumTable <- bind_rows(sumTable, smry)
    rownames(sumTable)[n_class] <- paste0(modelType, "_", n_class)
    
    ng <- unique(class.df$class)
    
    df.2 <- merge(df.1a, class.df, by = "PatientNum", all.x = TRUE)
    
    plot <- ggplot() +
      # x lim(0, 36)+
      labs(x = "Follow-up time (months)", y = "Weight Change (%)",
           title = paste0(modelType, " ", n_class, "-class model (BIC=", BIC_value, ") - ", surg, " patients")) +
      theme_bw()
    
    df.merged <- data.frame()
    df.C.O.List <- list()
    df.C.M.List <- list()
    ng <- ng[order(ng)]
    for (i in ng) {
      
      df.C <- df.2[df.2$class == i,]
      label <- paste0("Group ", i, " (n = ", length(unique(df.C$PatientID)), ", ", round(smry[,x+i],1), "%)")
      
      df.merged <- rbind(df.merged, df.C)
      
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
    plotList[[index]] <- plot
    index <- index + 1
    
  } # for (n_class in n_classes)
  gmmSummary <- rbind(gmmSummary, sumTable)
  
  filename.2 <- paste0(modelType, "_Modeled_and_Observed_Weight_Change.pdf")
  filePath <- paste0(outputDir, filename.2)
  
  pdf(filePath, width = 10, height = 5)
  for (p in 1:length(plotList)) {
    print(plotList[[p]])
  }
  dev.off()
  
  outputName <- paste0(modelType, "_Summary.tsv")
  outputPath <- file.path(outputDir, outputName)
  Surgery <- surg
  gmmSummary <- cbind(Surgery, gmmSummary)
  GMM_Type <- modelType
  gmmSummary <- cbind(GMM_Type, gmmSummary)
  
  write.table(gmmSummary, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
  
  allResults <- rbind(allResults, gmmSummary)
  
}# for (modelType in modelTypes)

outputName <- paste0("ALL_Summary.tsv")
outputPath <- file.path(outputDir, outputName)
write.table(allResults, outputPath, sep="\t", quote = FALSE, row.names = FALSE)
