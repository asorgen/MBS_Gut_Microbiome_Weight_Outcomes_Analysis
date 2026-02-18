#Author: Alicia Sorgen
#Date: 2022 October 21
#Description: Generate bar plots


##### Edits for script #####
rm(list=ls())
set.seed(1989)
params <- vector()

ANALYSIS <- "microbiome_n124"
# ANALYSIS <- "microbiome_sex"
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
# params <- c(params, "Quintile")
# params <- c(params, "Quartile")
params <- c(params, "Tertile")
# params <- c(params, "Half")

moduleRoot <- paste0("Assign_WLgroups")

##### Libraries #####
R <- sessionInfo()
message("\n************* R Environment *************\n")
message(R$R.version$version.string)
rm(R)

message("\n************* R packages *************\n")
library(stringr); message("stringr: Version ", packageVersion("stringr"))
library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
library(scales); message("scales: Version ", packageVersion("scales"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
# library(data.table); message("data.table: Version ", packageVersion("data.table"))

##### Set up working environment #####
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- params
}
rm(params)

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
  
  if (length(args) > 1) {
    if (args[2] == "MetaPhlAn2" | args[2] == "Kraken2") {
      module <- paste0(args[2], "_", module)
    } 
    
    if (length(args) > 2) {
      if (args[3] %in% c("Quintile", "Quartile", "Thirds", "Half")) {
        module <- paste0(module, "_", args[3])
      } 
      
    }
  }
  if (exists("end_month") == TRUE) {
    module <- paste0(module, "_BLto", end_month, "M")
    rm(end_month)
  }
  
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
rm(moduleRoot, ANALYSIS)

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
prevModule <- "WeightMetaMerge"
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
inputFile <- "metadata.tsv"

month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")
divisions <- c("Quintile", "Quartile", "Tertile", "Half")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Script variables #####
divisions <- args[2:length(args)]

##### Assign weight loss groups #####

# Read in table
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

# Order Timepoint as factor
myTable$Timepoint <- as.factor(myTable$Timepoint)
myTable$Timepoint <- factor(myTable$Timepoint, levels = c("BL", "ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR"))

# Assign SampleID to as dataframe row names
row.names(myTable) = myTable$SampleID

# Convert the timepoint data to numeric values.
myTable$time = as.numeric(myTable$time)

months <- c("12", "18", "24")
included <- c("0", "1", "6")

# Weight loss group assignment
myTable2 <- myTable[ !is.na( myTable$Percent_Loss_kg ), ] # remove rows with percent change NAs
months <- c(1, 6, 12, 18, 24)

for (m in months) {
  
  df.end <- myTable2[myTable2$time == m,]
  divAssignments_endM <- vector()
  
  for (division in divisions) {
    message(paste0("\n************* ", division, " assignment *************\n"))
    
    if (division == "Quintile") {
      div_endM <- quantile(df.end$Percent_Loss_kg , probs=c(0,.2,.4,.6,.8,1), na.rm = TRUE)
      
      for( i in 1:nrow(df.end)) {
        df.end$divAssignments_endM[i] = getQuintileGroup(  df.end$Percent_Loss_kg[i], div_endM )
      }
      
    }
    if (division == "Quartile") {
      div_endM <- quantile(df.end$Percent_Loss_kg , probs=c(0,.25,.5,.75,1), na.rm = TRUE)
      
      for( i in 1:nrow(df.end)) {
        df.end$divAssignments_endM[i] = getQuartileGroup(  df.end$Percent_Loss_kg[i], div_endM )
      }
      
    }
    if (division == "Tertile") {
      div_endM <- quantile(df.end$Percent_Loss_kg , probs=c(0, (1/3), (2/3), 1), na.rm = TRUE)
      
      for( i in 1:nrow(df.end)) {
        df.end$divAssignments_endM[i] = getTertileGroup(  df.end$Percent_Loss_kg[i], div_endM )
      }
      
    }
    if (division == "Half") {
      div_endM <- quantile(df.end$Percent_Loss_kg , probs=c(0,.5,1), na.rm = TRUE)
      
      for( i in 1:nrow(df.end)) {
        df.end$divAssignments_endM[i] = getHalfGroup(  df.end$Percent_Loss_kg[i], div_endM )
      }
      
    }
    
    PatientID <- df.end$PatientID
    Division <- df.end$divAssignments_endM
    PWL_endM <- df.end$Percent_Loss_kg
    df.div <- data.frame(PatientID, Division)
    names(df.div)[names(df.div) == "Division"] <- paste0(division, "Assignment_", m, "M")
    myTable <- merge(df.div, myTable, by = "PatientID", all = TRUE)
    
    
  }
  
  
} # for (m in months)


file.path <- paste0(outputDir, inputFile)
write.table(myTable, file.path, sep="\t",row.names = FALSE,quote = FALSE)
