#Author: Alicia Sorgen
#Date: 10-12-21
#Description: Merge taxa count tables with metadata

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
# library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
# library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))
# library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))


# Script-specific edits -------------------------------------------------
ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, getOption("mbs.pipe_root", default = stop("Set options(mbs.pipe_root = ...) in your local .Rprofile (gitignored)")))
params <- c(params, "MetaPhlAn2")
# params <- c(params, "Kraken2")

module <- "2.0_TaxaMetaMerge"

levels <- vector()
levels <- c(levels, "phylum")
levels <- c(levels, "class")
levels <- c(levels, "order")
levels <- c(levels, "family")
levels <- c(levels, "genus")
levels <- c(levels, "species")


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

# inputDir = file.path(input_root, "Data")
rawDir   = file.path(input_root, "microbiome", "taxa_raw")
normDir  = file.path(input_root, "microbiome", "taxa_normalized")
# sampleID_corrections.R defines sample ID corrections specific to this dataset.
# This file is not included in the public repository; see .gitignore.
source(file.path(input_root, "metadata", "sampleID_corrections.R"))
rawFile <- "_rawCounts"
logFile <- "_LogNormalizedCounts"
relFile <- "_RelativeAbundanceCounts"

metaDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "WeightMetaMerge"),"/output/")
metaFile <- "metadata.tsv"

##### Set up output #####
output = file.path(moduleDir,"output/")

##### Merge data #####
metaTable <- import_file(paste0(metaDir, metaFile))


for (level in levels) {
  
  for (File in c(rawFile, logFile)) {
     if (File == logFile && level == "species") next  # no species-level log-normalized file in Data
     
     if (File == rawFile) {
        Table <- import_file(file.path(rawDir,  paste0("bs_taxa_", level, "_raw.tsv")))
     } else if (File == logFile) {
        Table <- import_file(file.path(normDir, paste0("bs_taxa_", level, "_norm_log10.tsv")))
        } else {Table <- getRelAbun(import_file(file.path(rawDir, paste0("bs_taxa_", level, "_raw.tsv"))))}
    
     if (level == "genus") Table <- Table[,-which(colnames(Table) %like% "virus")]
     
    SampleID <- sapply(strsplit(Table$SampleID, "_"), "[", 1)
    # Strip aliquot suffixes (e.g. "1-241-00-Aliq1" -> "1-241-00") before format conversion
    SampleID <- sub("-Aliq[0-9]+$", "", SampleID)
    # Convert taxa SampleID format "site-patient-timepoint" -> "BIO-site-patient-ZZtimepoint"
    SampleID <- sapply(strsplit(SampleID, "-"), function(x) {
      if (length(x) == 3) sprintf("BIO-%s-%s-%02d", x[1], x[2], as.integer(x[3]))
      else paste(x, collapse="-")
    })
    SampleID <- gsub(pattern = sampleID_mislabeled, replacement = sampleID_corrected, SampleID)
    
    Table <- cbind(SampleID, Table)
    
    if (File == rawFile) {
      dup <- SampleID[which(duplicated(SampleID) == TRUE)]
      dupRows <- vector()
      
      for (i in dup) {
        
        df <- Table[Table$SampleID == i,]
        rowSum <- rowSums(df[3:ncol(df)])
        dupRows <- c(dupRows, which(rownames(Table) %in% rownames(df)[which(rowSum != max(rowSum))]))
        
      } # for (i in dup)
      
    }    
    Table <- Table[-(dupRows),]
    Table <- Table[,-2]
    
    df <- merge(metaTable, Table, by = "SampleID", all = TRUE)
    
    for (i in 1:nrow(df)) {
      
      if (is.na(df$PatientID[i]) && !is.na(df$SampleID[i])) {
        string <- strsplit(as.character(df$SampleID[i]),split = "-")
        temp_string=do.call(rbind,string)
        if (ncol(temp_string) < 4) next  # skip IDs that don't match BIO-site-patient-tp format
        df$PatientID[i] <- paste0(temp_string[,1], "-", temp_string[,2], "-", temp_string[,3])
        df$Timepoint[i] <- temp_string[,4]
        df$Site[i] <- temp_string[,2]
      } # if (is.na(df$PatientID[i]))
      
    } # for (i in 1:nrow(df))
    
    df$Timepoint <- gsub(pattern = "24", replacement = "TWENTY_FOUR", df$Timepoint)
    df$Timepoint <- gsub(pattern = "18", replacement = "EIGHTEEN", df$Timepoint)
    df$Timepoint <- gsub(pattern = "12", replacement = "TWELVE", df$Timepoint)
    df$Timepoint <- gsub(pattern = "06", replacement = "SIX", df$Timepoint)
    df$Timepoint <- gsub(pattern = "01", replacement = "ONE", df$Timepoint)
    df$Timepoint <- gsub(pattern = "00", replacement = "BL", df$Timepoint)
    
    df$Site <- gsub(pattern = "1", replacement = "Fargo", df$Site)
    df$Site <- gsub(pattern = "2", replacement = "Cleveland", df$Site)
    
    for (i in 1:nrow(df)) {
      
      if (is.na(df$Sex[i])) {
        
        ID <- df$PatientID[i]
        df1 <- df[df$PatientID == ID,]
        df1 <- na.omit(df1)
        
        if (nrow(df1) > 0) {
          df$Sex[i] <- df1$Sex
          df$TypeofSurgery[i] <- df1$TypeofSurgery
          df$Ethnicity[i] <- df1$Ethnicity
          
        } # if (nrow(df1) > 0)
        
      } # if (is.na(df$Sex[i]))
      
    } # for (i in 1:nrow(df))
    
    fileOutput <- paste0(output, level, File, "_", classifier, ".tsv"); fileOutput
    write.table(df, fileOutput, sep="\t",quote = FALSE, row.names = FALSE)
    
  } # for (File in c(rawFile, logFile, relFile))  
  
} # for (level in levels)
