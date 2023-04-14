#Author: Alicia Sorgen
#Date: 2023 Feb 23
#Description: Machine learning model for bariatric surgery weight loss outcomes.

#####- Libraries ----------------------------------------------------------#####
R <- sessionInfo()
message(R$R.version$version.string)

library(stringr); message(">>> stringr: Version ", packageVersion("stringr"), "\n")
library(ROCR); message(">>> ROCR: Version ", packageVersion("ROCR"), "\n")
library(randomForest); message(">>> randomForest: Version ", packageVersion("randomForest"), "\n")
library(rstatix); message(">>> rstatix: Version ", packageVersion("rstatix"), "\n")
library(data.table); message(">>> data.table: Version ", packageVersion("data.table"), "\n")

#####- Edits for script ---------------------------------------------------#####
rm(list=ls())

ANALYSIS <- "ML_microbiome"
params <- vector()
params <- c(params, "~/git/ML_MicrobiomeAnalysis_2023")
params <- c(params, "MetaPhlAn2")
# params <- c(params, 2)
params <- c(params, 3)
params <- c(params, "TaxaOnly")
# params <- c(params, "ClinicalOnly")
# params <- c(params, "Taxa_and_Clinical")
# params <- c(params, "Taxa_Clinical_Energy")
params <- c(params, "RYGB")

moduleRoot <- "Classification_randomForest_Prediction"
included <- c(0, 1, 6, 12, 18, 24)

#####- Set up working environment -----------------------------------------#####
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- params
  rm(params)
}

if (args[1] == "BLJ") {
  message("\n************* Running in BioLockJ *************")
} else {
  message("\n************* Running locally *************")
  gitRoot <- args[1]
  gitInput <- file.path(gitRoot, "analysis", "input")
  gitScripts <- file.path(gitRoot, "analysis", "Rscripts"); rm(gitRoot)

  root <- paste0("~/BioLockJ_pipelines/")
  pipeline <- paste0(ANALYSIS,"_analysis_")
  
  if (any(dir(root) %like% pipeline) == TRUE) {
    root <- paste0(root,"/",str_subset(dir(root), pipeline), "/")
  } else {
    today <- as.character(format(Sys.Date(), "%Y%b%d"))
    root <- paste0(root,ANALYSIS,"_analysis_", today, "/")
    dir.create(root, showWarnings = FALSE)
  }; rm(pipeline, ANALYSIS)
  
  if (any(dir(root) == "input") == FALSE) {
    # rootInput <- paste0(root, "input/")
    # dir.create(rootInput, showWarnings = FALSE)
    
    file.copy(gitInput,
              root,
              recursive = TRUE)
    
  }; rm(gitInput)
  
  module <- moduleRoot
  
  if (length(args) > 2) {
    m1 <- sapply(strsplit(module, "_"), "[", 1)
    m2 <- sapply(strsplit(module, "_"), "[", 2)
    m3 <- sapply(strsplit(module, "_"), "[", 3)
    
    if (args[3] == 2 | args[3] == 3) {
      module <- paste0("TertileGroups", args[3], "_", m2, "_", "Prediction")
    }
    module <- paste0(module, "_", args[4])
    rm(m1, m2, m3)
  }
  
  if (length(args) > 4) {
    m1 <- sapply(strsplit(module, "_"), "[", 1)
    m2 <- sapply(strsplit(module, "_"), "[", 2)
    m3 <- sapply(strsplit(module, "_"), "[", 3)
    m4 <- sapply(strsplit(module, "_"), "[", 4)
    
    module <- paste(m1, m2, m3, args[5], m4, sep = "_")
    rm(m1, m2, m3, m4)
  }
  
  if (exists("end_month") == TRUE) {
    mod1 <- paste0(sapply(strsplit(module, "_Weight"), "[", 1), "_"); mod1
    mod2 <- paste0("_Weight", sapply(strsplit(module, "_Weight"), "[", 2)); mod2
    module <- paste0(mod1, end_month, "M", mod2, "_BLto", included[length(included)], "M"); module
    rm(end_month, mod1, mod2)
  }
  
  if (any(dir(root) %like% module) == TRUE) {
    moduleDir <- paste0(root,str_subset(dir(root), module), "/")
  } else {
    moduleDir <- paste0(root, module, "/")
    dir.create(moduleDir, showWarnings = FALSE)
  }; rm(module, root)
  
  scriptDir <- paste0(moduleDir, "script/")
  if (any(dir(moduleDir) == "script") == FALSE) {
    dir.create(scriptDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/", moduleRoot, ".R")
    file.copy(script,
              scriptDir,
              recursive = TRUE)
  }; rm(scriptDir, moduleRoot)
  
  outputDir <- paste0(moduleDir, "output/")
  if (any(dir(moduleDir) == "output") == FALSE) {
    dir.create(outputDir, showWarnings = FALSE)
  }
  
  files <- list.files(outputDir, recursive = TRUE, full.names = TRUE)
  file.remove(files); rm(files, outputDir)
  
  resourcesDir <- paste0(moduleDir, "resources/")
  if (any(dir(moduleDir) == "resources") == FALSE) {
    dir.create(resourcesDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/functions.R")
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
  }; rm(resourcesDir, gitScripts)
  
  setwd(paste0(moduleDir, "script/"))

}

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())

#####- Set up functions file ----------------------------------------------#####
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]
rm(str)

#####- Set up input -------------------------------------------------------#####
classifier <- args[2]

prevModule <- paste0(classifier, "_TaxaMetaMerge")
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
prevModule <- paste0(classifier, "_Taxa_over_time_24M_Weight_Class_BLto24M")
statsDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")

month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")
division <- c("Tertile")

#####- Set up output ------------------------------------------------------#####
outputDir = file.path(moduleDir,"output/")

#####- Set up script variables ------------------------------------------######
groups <- args[3]
trainingInput <- args[4]
level <- "genus"
Level <- str_to_title(level)
clinicalVariables <- c("Age", "Site", "Sex", "Ethnicity", "Loss_from_BL_BMI", "Loss_from_BL_kg", "Percent_Loss_kg", "PEWL_kg")
Timepoints <- c("BL", "ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR")

# ntreeValues <- c(100, 200, 300, 400, 500)
ntreeValues <- c(500)

assignmentMonths <- included[2:length(included)]
# assignmentMonths <- included[length(included)]

endMonths <- included[2:length(included)]
# endMonths <- included[2]

#####- Import data --------------------------------------------------------#####
fileName <- paste0(inputDir, level, "_LogNormalizedCounts_MetaPhlAn2.tsv")
myT <- read.delim(fileName, sep="\t",header = TRUE)


fileName <- paste0(statsDir, level, "/", division, "/", level,  "_by_Timepoint_Tertile_LM_pValues_BLto24M_MetaPhlAn2.tsv")
taxaNames <- read.delim(fileName, sep="\t",header = TRUE); outputName <- "_Non-rare"
taxaNames <- taxaNames[taxaNames$adjPValues < 0.05,]; outputName <- "_Significant"
# taxaNames <- taxaNames[taxaNames$LengthTime == "BLto1M",]; features <- " Significant change from BL to 1M"; outputName <- "_Significant_BLto1M"
taxaNames <- taxaNames[order(taxaNames$divPValues),]
taxaNames <- unique(taxaNames$qbugName)
# taxaNames <- taxaNames[1:10]; features <- " Significant change from BL to 1M (Top 10)"; outputName <- "_Top10_Significant_BLto1M"

features <- paste0(Level)
outputName <- paste0(Level, outputName)

#####- Set random seed to make results reproducible: ----------------------#####
set.seed(200)

#####- Calculate taxa change: ---------------------------------------------#####
taxa.dfList <- list()
listIndex <- 1
for (endM in endMonths) {
  
  taxaList <- list()
  index <- 1
  for (i in taxaNames) {
    
    Timepoint <- myT$Timepoint
    PatientID <- myT$PatientID
    TaxaCol <- myT[,which(colnames(myT) == i)]
    
    df <- data.frame(PatientID, Timepoint, TaxaCol)
    df <- spread(df, Timepoint, TaxaCol)
    taxaList[[index]] <- df[,which(colnames(df) == Timepoints[which(included == endM)])] - df$BL
    names(taxaList)[index] <- i
    index <- index + 1
    
  } # for (i in taxaNames)
  
  PatientID <- df$PatientID
  taxa.df <- data.frame(PatientID)
  for (i in 1:length(taxaList)) {
    
    taxa.df <- data.frame(taxa.df, taxaList[i])
    
  } # for (i in 1:length(taxaList))
  
  taxa.df <- na.omit(taxa.df)
  
  taxa.dfList[[listIndex]] <- taxa.df
  names(taxa.dfList)[listIndex] <- paste0("Taxa_BLto", endM, "M")
  
  listIndex <- listIndex + 1
  
} # for (endM in 2:length(included))

#####- Set up training sets -----------------------------------------------#####
filtered.dfs <- list()
trainingSets <- list()
validationSets <- list()

df.index <- 1
set.index <- 1
for (assignmentMonth in assignmentMonths) {
  
  outputName2 <- paste0(assignmentMonth, "M_WLG")
  # message(paste0(">>>>> ", assignmentMonth, "M weight loss groups"))
  
  ##- Add weight loss assignment ---------------------------------------------##
  
  for (x in 1:length(taxa.dfList)) {
    
    Taxa.DF.Name <- names(taxa.dfList)[x]
    changeName <- sapply(str_split(Taxa.DF.Name, "_"), "[", 2)
    endM <- sapply(str_split(Taxa.DF.Name, "to"), "[", 2)
    endM <- gsub("M", "", endM)
    TimeChange <- paste0("BL - ", endM, "M")
    outputName3 <- paste0(outputName2, "_", changeName)
    # message(paste0("     >>>>> ", TimeChange))
    
    if (trainingInput == "TaxaOnly") {
      meta.df <- myT[which(myT$Timepoint == Timepoints[x+1]), c("PatientID", paste0("TertileAssignment_", assignmentMonth, "M")) ]
      df <- merge(meta.df, taxa.dfList[[x]], by = "PatientID", all = FALSE)
      names(df)[names(df) == paste0("TertileAssignment_", assignmentMonth, "M")] <- "WeightLossGroup"
      dfName <- paste0(Taxa.DF.Name, "_", assignmentMonth, "M")
    } # if (trainingInput == "TaxaOnly")
    
    if (trainingInput == "ClinicalOnly") {
      df <- myT[which(myT$Timepoint == Timepoints[x+1]), c("PatientID", paste0("TertileAssignment_", assignmentMonth, "M"), clinicalVariables) ]
      names(df)[names(df) == paste0("TertileAssignment_", assignmentMonth, "M")] <- "WeightLossGroup"
      dfName <- paste0("Clinical_", changeName, "_", assignmentMonth, "M")
      for (var in clinicalVariables) {
        for (j in 1:nrow(df)) {
          if (is.na(df[j,which(colnames(df) == var)]) == TRUE) {
            df[j,which(colnames(df) == var)] <- mean(na.omit(df[,which(colnames(df) == var)]))
          }
        }
      }
    } # if (trainingInput == "ClinicalOnly")
    
    if (trainingInput == "Taxa_and_Clinical") {
      meta.df <- myT[which(myT$Timepoint == Timepoints[x+1]), c("PatientID", paste0("TertileAssignment_", assignmentMonth, "M"), clinicalVariables) ]
      df <- merge(meta.df, taxa.dfList[[x]], by = "PatientID", all = FALSE)
      names(df)[names(df) == paste0("TertileAssignment_", assignmentMonth, "M")] <- "WeightLossGroup"
      dfName <- paste0("Clinical_", Taxa.DF.Name, "_", assignmentMonth, "M")
      for (var in clinicalVariables) {
        for (j in 1:nrow(df)) {
          if (is.na(df[j,which(colnames(df) == var)]) == TRUE) {
            df[j,which(colnames(df) == var)] <- mean(na.omit(df[,which(colnames(df) == var)]))
          }
        }
      }
      
    } # if (trainingInput == "Taxa_and_Clinical")
    
    if (groups == 2) {
      df$WeightLossGroup <- ifelse(df$WeightLossGroup == 1, "Bottom 33%",
                                   ifelse(df$WeightLossGroup == 2, "Top 67%",
                                          ifelse(df$WeightLossGroup == 3, "Top 67%", NA)))
    } else {
      df$WeightLossGroup <- ifelse(df$WeightLossGroup == 1, "Bottom 33%",
                                   ifelse(df$WeightLossGroup == 2, "Middle 33%",
                                          ifelse(df$WeightLossGroup == 3, "Top 33%", NA)))
      
    } # if (groups == 2)
    
    df <- df[!is.na(df$WeightLossGroup),]
    df$WeightLossGroup <- factor(df$WeightLossGroup)
    factorLength <- length(unique(df$WeightLossGroup))
    filtered.dfs[[df.index]] <- df
    names(filtered.dfs)[df.index] <- dfName
    
    #####- Calculate the size of each of the data sets
    splitSet <- 3
    data_set_size <- floor(nrow(df) * 1/splitSet)
    
    #####- Generate a random sample of "data_set_size" subset
    randomIDs <- list()
    PatientIDs <- df$PatientID
    for (nSet in 1:splitSet) {
      
      if (nSet != splitSet) {
        subset <- sample(PatientIDs, size = data_set_size)
        randomIDs[[nSet]] <- subset
        PatientIDs <- PatientIDs[-which(PatientIDs %in% subset)]
      } else {
        subset <- PatientIDs
        randomIDs[[nSet]] <- subset
        PatientIDs <- PatientIDs[-which(PatientIDs %in% subset)]
        
      }
    }
    
    #####- Assign the data to the correct sets
    for (id in 1:length(randomIDs)) {
      
      IDs <- randomIDs[[id]]
      training <- df[!(df$PatientID %in% IDs),]
      rownames(training) <- training[,1]
      training <- training[,-1]
      trainingSets[[set.index]] <- training
      names(trainingSets)[set.index] <- paste0(dfName, "_Set", id)
      
      validation <- df[df$PatientID %in% IDs,]
      rownames(validation) <- validation[,1]
      validation <- validation[,-1]
      validation_N <- nrow(validation)
      validationSets[[set.index]] <- validation
      names(validationSets)[set.index] <- paste0(dfName, "_Set", id)
      
      set.index <- set.index + 1
      
    } # for (id in 1:length(randomIDs))
    
    df.index <- df.index + 1
  } # for (x in 1:length(taxa.dfList))
  
} # for (a in 2:length(included))



#####- Clean up environment -----------------------------------------------#####
rm(df, meta.df, randomIDs, taxa.df, taxaList, training, validation)

#####- Training -----------------------------------------------------------#####
index <- 1
AssignmentMonth <- vector()
TimeChanges <- vector()
ntree_parameter <- vector()
mtry_parameter <- vector()
OOB_estimate <- vector()
TrainingSet <- vector()
rf.Training <- list()

for (a in 1:length(trainingSets)) {
  
  DF.Name <- names(trainingSets)[a]
  assignmentMonth <- sapply(str_split(DF.Name, "_Set"), "[", 1)
  assignmentMonth <- sapply(str_split(assignmentMonth, "to"), "[", 2)
  assignmentMonth <- gsub("M", "", assignmentMonth)
  endM <- sapply(str_split(assignmentMonth, "_"), "[", 1)
  assignmentMonth <- sapply(str_split(assignmentMonth, "_"), "[", 2)
  set <- sapply(str_split(DF.Name, "_Set"), "[", 2)
  
  outputName2 <- paste0(assignmentMonth, "M_WLG")
  outputDir_WLG <- paste0(outputDir, outputName2, "/")
  dir.create(outputDir_WLG, showWarnings = FALSE)
  
  # message(paste0(">>>>> ", assignmentMonth, " weight loss groups"))
  
  changeName <- paste0("BLto", endM, "M")
  TimeChange <- paste0("BL - ", endM, "M")
  outputName3 <- paste0(outputName2, "_", changeName)
  # message(paste0("     >>>>> ", TimeChange))
  
  AUC.df <- data.frame()

  # message(paste0("          >>>>> Training #", set))
  training <- trainingSets[[a]]
  Nf <- length(2:ncol(training))
  training_N <- nrow(training)
  
  if (set == "1") {
    file.path <- paste0(outputDir_WLG, "VariableImportancePlot_", outputName, "_", outputName3, ".pdf")
    pdf(file.path, width = 15, height = 10)
  }
  
  for (ntreeValue in ntreeValues) {
    
    # mtryValues <- 1:floor(sqrt(Nf))
    mtryValues <- 1
    
    for (mtryValue in mtryValues) {
      
      rf_classifier <- randomForest(WeightLossGroup ~ ., # "WeightLossGroup ~ ." we want to predict WeightLossGroup usings each of the remaining columns of data
                                    data=training,
                                    ntree=ntreeValue, # defines the number of trees to be generated; it is typical to test a range of values for this parameter (i.e. 100, 200, 300, 400, 500) and choose the one that minimises the OOB estimate of error rate
                                    # mtry=mtryValue, # the number of features used in the construction of the tree; these features are selected at random, which is where the "random" in randomForests comes from; the default is sqrt(Nf)
                                    importance=TRUE) # enables the algorithm to calculate variable importance
      
      
      rf.Training[[index]] <-  rf_classifier
      # names(rf.Training)[index] <- paste0(DF.Name, "_", ntreeValue, "_", mtryValue)
      names(rf.Training)[index] <- DF.Name
      
      OOB <- rf_classifier$confusion
      OOB <- OOB[,-ncol(OOB)]
      
      correctPredictions <- 0
      for (i in 1:factorLength) {
        correctPredictions <- correctPredictions + OOB[i,i]
      } # for (i in 1:factorLength)
      
      title.lab <- paste0(features, " (Nf = ", Nf, "); Training Set ", set, " (No = ", training_N, ")\nntree = ", ntreeValue, ", mtry = ", mtryValue, ", OOB estimate = ", round((sum(OOB) - correctPredictions) / sum(OOB) * 100, 2), "%")
      varImpPlot(rf_classifier, main = title.lab)
      
      OOB_estimate[index] <- (sum(OOB) - correctPredictions) / sum(OOB)
      ntree_parameter[index] <- rf_classifier$ntree
      mtry_parameter[index] <- rf_classifier$mtry
      TimeChanges[index] <- TimeChange
      AssignmentMonth[index] <- assignmentMonth
      TrainingSet[index] <- set
      
      index <- index + 1
      
    } # for (mtryValue in mtryValues)
  } # for (ntreeValue in ntreeValues)
  if (set == splitSet) {
    dev.off()
  }
  
} # for (a in 1:length(trainingSets))

trainingParameters <- data.frame(AssignmentMonth, TimeChanges, TrainingSet, ntree_parameter, mtry_parameter, OOB_estimate)
file.path <- paste0(outputDir, outputName, "_MachineLearningTraining.tsv")
write.table(trainingParameters, file=file.path, sep="\t", row.names=FALSE)

#####- Validation ---------------------------------------------------------#####
index <- 1
Accuracy <- vector()
AssignmentMonth <- vector()
TimeChanges <- vector()
ntree_parameter <- vector()
mtry_parameter <- vector()
TrainingSet <- vector()
Pred.Table <- list()
Pred.ROC <- list()
AUC.df <- data.frame()
ntreeIndex <- vector()
mtryIndex <- vector()

for (a in 1:length(rf.Training)) {
  
  DF.Name <- names(rf.Training)[a]
  dfName <- sapply(str_split(DF.Name, "_Set"), "[", 1)
  assignmentMonth <- sapply(str_split(dfName, "to"), "[", 2)
  assignmentMonth <- gsub("M", "", assignmentMonth)
  endM <- sapply(str_split(assignmentMonth, "_"), "[", 1)
  assignmentMonth <- sapply(str_split(assignmentMonth, "_"), "[", 2)
  set <- sapply(str_split(DF.Name, "_Set"), "[", 2)
  
  outputName2 <- paste0(assignmentMonth, "M_WLG")
  # message(paste0(">>>>> ", assignmentMonth, " weight loss groups"))
  outputDir_WLG <- paste0(outputDir, outputName2, "/")
  dir.create(outputDir_WLG, showWarnings = FALSE)
  
  changeName <- paste0("BLto", endM, "M")
  TimeChange <- paste0("BL - ", endM, "M")
  outputName3 <- paste0(outputName2, "_", changeName)
  # message(paste0("     >>>>> ", TimeChange))
  
  validation <- validationSets[[which(names(validationSets) == DF.Name)]]
  validation_N <- nrow(validation)
  
  rf_classifier <- rf.Training[[a]]
  ntreeValue <- rf_classifier$ntree
  mtryValue <- rf_classifier$mtry
  
  # message(paste0("          >>>>> Validation #", set))
  ##- Validation set assessment #1: looking at confusion matrix
  prediction_for_table <- predict(rf_classifier,validation[,-1]) # column 1 is the WeightLossGroup column we are predicting
  
  results <- cbind(prediction_for_table, validation[,1])
  Accuracy[a] <- sum(prediction_for_table==validation[,1]) / nrow(validation)
  AssignmentMonth[a] <- assignmentMonth
  TimeChanges[a] <- TimeChange
  ntree_parameter[a] <- ntreeValue
  mtry_parameter[a] <- mtryValue
  TrainingSet[a] <- set
  
  
  ##- Validation set assessment #2: ROC curves and AUC
  prediction_for_roc_curve <- predict(rf_classifier,validation[,-1],type="prob")
  
  if (set == 1) {
    prediction_for_table_Merged <- prediction_for_table
    prediction_for_roc_curve_Merged <- prediction_for_roc_curve
    file.path <- paste0(outputDir_WLG, "ROC_Curve_", outputName, "_", outputName3, ".pdf")
    pdf(file.path, width = 6, height = 5)
    
  } else {
    prediction_for_table_Merged <- c(prediction_for_table_Merged, prediction_for_table)
    prediction_for_roc_curve_Merged <- rbind(prediction_for_roc_curve_Merged, prediction_for_roc_curve)
  }
  
  ##- Use pretty colours:
  pretty_colours <- c("#F8766D","#00BA38","#619CFF")
  
  ##- Specify the different classes
  classes <- levels(validation$WeightLossGroup)
  
  
  AUC <- vector()
  
  for (i in 1:length(classes)) {
    ##- Define which observations belong to class[i]
    true_values <- ifelse(validation[,1]==classes[i],1,0)
    
    # test <- cbind(validation[,1], true_values)
    ##- Assess the performance of classifier for class[i]
    pred <- prediction(prediction_for_roc_curve[,i],true_values)
    # test <- cbind(prediction_for_roc_curve[,i], true_values)
    
    perf <- performance(pred, "tpr", "fpr")
    ##- Calculate the AUC and print it to screen
    auc.perf <- performance(pred, measure = "auc")
    RF.AUC <- unlist(slot(auc.perf, "y.values"))
    AUC[i] <- RF.AUC
    
    if (i==1)    {
      plot(perf,
           avg = "threshold",
           # main=title.lab,
           lwd = 3,
           col=pretty_colours[i])
    }    else    {
      plot(perf,
           avg = "threshold",
           # main=title.lab,
           lwd = 3,
           col=pretty_colours[i],
           add=TRUE)
    }
    
    # message(paste0(classes[i], " AUC = ", RF.AUC))
    
  } # for (i in 1:length(classes))
  abline(0,1)
  AUC.df <- rbind(AUC.df, AUC)
  for (i in 1:length(classes)) {
    names(AUC.df)[i] <- classes[i]
  }
  mean(AUC)
  
  ##- Add title
  title.lab <- paste0(features, " ", assignmentMonth, "-month WLG BL to ", endM, "M (Nf = ", Nf, ")\n Accuracy = ", round(Accuracy[a] * 100, 2), "%; AUC = ", round(mean(AUC), 3))
  subtitle.lab <- paste0("Training set ", set, " (", training_N, ") Validation set ", set, " (", validation_N, ") ntree = ", ntreeValue, ", mtry = ", mtryValue)
  title(main = title.lab, sub = subtitle.lab)
  
  ##- Add a legend
  legend.labs <- paste0(classes, " (AUC = ", round(AUC, 2), ")")
  legend.labs <- paste0(classes)
  legend(0.5, 0.3, legend=c(legend.labs),
         col=c(pretty_colours), lty = 19, cex=0.9, box.lty=0)
  
  if (set == splitSet) {
    Pred.Table[[index]] <- prediction_for_table_Merged
    names(Pred.Table)[index] <- dfName
    Pred.ROC[[index]] <- prediction_for_roc_curve_Merged
    names(Pred.ROC)[index] <- dfName
    ntreeIndex[index] <- ntreeValue
    mtryIndex[index] <- mtryValue
    
    dev.off()
    index <- index + 1
  }
  
} # for (a in 1:length(trainingSets))

validationResults <- data.frame(AssignmentMonth, TimeChanges, TrainingSet, ntree_parameter, mtry_parameter, Accuracy, AUC.df)
validationResults$Mean.AUC <- rowMeans(AUC.df)
file.path <- paste0(outputDir, outputName, "_MachineLearningValidation.tsv")
write.table(validationResults, file=file.path, sep="\t", row.names=FALSE)

#####- Validation Merge ---------------------------------------------------#####
index <- 1
Accuracy <- vector()
AssignmentMonth <- vector()
TimeChanges <- vector()
ntree_parameter <- vector()
mtry_parameter <- vector()
AUC.df <- data.frame()

for (a in 1:length(Pred.ROC)) {
  
  
  df <- filtered.dfs[[a]]
  dfName <- names(filtered.dfs)[[a]]
  
  assignmentMonth <- sapply(str_split(dfName, "to"), "[", 2)
  assignmentMonth <- gsub("M", "", assignmentMonth)
  endM <- sapply(str_split(assignmentMonth, "_"), "[", 1)
  assignmentMonth <- sapply(str_split(assignmentMonth, "_"), "[", 2)
  
  outputName2 <- paste0(assignmentMonth, "M_WLG")
  outputDir_WLG <- paste0(outputDir, outputName2, "/")
  dir.create(outputDir_WLG, showWarnings = FALSE)
  
  changeName <- paste0("BLto", endM, "M")
  TimeChange <- paste0("BL - ", endM, "M")
  outputName3 <- paste0(outputName2, "_", changeName)
  ntreeValue <- ntreeIndex[a]
  mtryValue <- mtryIndex[a]
  
  df2 <- df[,1:2]
  prediction_for_roc_curve <- Pred.ROC[[a]]
  prediction_for_table <- Pred.Table[[a]]
  prediction_for_table <- as.data.frame(prediction_for_table)
  PatientID <- rownames(prediction_for_table)
  prediction_for_table$PatientID <- PatientID
  df2 <- merge(prediction_for_table, df2, by = "PatientID")
  
  Accuracy[a] <- sum(df2$prediction_for_table==df2$WeightLossGroup) / nrow(df2)
  AssignmentMonth[a] <- assignmentMonth
  TimeChanges[a] <- TimeChange
  ntree_parameter[a] <- ntreeValue
  mtry_parameter[a] <- mtryValue
  TrainingSet[a] <- set
  
  ##- Use pretty colours:
  pretty_colours <- c("#F8766D","#00BA38","#619CFF")
  
  ##- Specify the different classes
  classes <- levels(df$WeightLossGroup)
  
  AUC <- vector()
  
  prediction_for_roc_curve <- as.data.frame(prediction_for_roc_curve)
  PatientID <- rownames(prediction_for_roc_curve)
  prediction_for_roc_curve$PatientID <- PatientID
  df2 <- merge(prediction_for_roc_curve, df2, by = "PatientID")
  
  file.path <- paste0(outputDir_WLG, "Merged_ROC_Curve_", outputName, "_", outputName3, ".pdf")
  pdf(file.path, width = 6, height = 5)
  
  for (i in 1:length(classes)) {
    ##- Define which observations belong to class[i]
    true_values <- ifelse(df2$WeightLossGroup==classes[i],1,0)
    
    ##- Assess the performance of classifier for class[i]
    pred <- prediction(df2[,i+1],true_values)
    perf <- performance(pred, "tpr", "fpr")
    
    ##- Calculate the AUC and print it to screen
    auc.perf <- performance(pred, measure = "auc")
    AUC[i] <- auc.perf@y.values[[1]]
    
    
    if (i==1)    {
      plot(perf,
           # avg = "threshold",
           # main=title.lab,
           lwd = 3,
           col=pretty_colours[i])
    }    else    {
      plot(perf,
           # avg = "threshold",
           # main=title.lab,
           lwd = 3,
           col=pretty_colours[i],
           add=TRUE)
    }
    
  } # for (i in 1:length(classes))
  
  abline(0,1)
  
  AUC.df <- rbind(AUC.df, AUC)
  for (i in 1:length(classes)) {
    names(AUC.df)[i] <- classes[i]
  }
  
  
  ##- Add title
  title.lab <- paste0(features, " ", assignmentMonth, "-month WLG BL to ", endM, "M (Nf = ", Nf, ")\n Accuracy = ", round(Accuracy[a] * 100, 2), "%; AUC = ", round(mean(AUC), 3))
  subtitle.lab <- paste0("Merged validation results; ntree = ", ntreeValue, ", mtry = ", mtryValue)
  title(main = title.lab, sub = subtitle.lab)
  
  ##- Add a legend
  legend.labs <- paste0(classes)
  # legend.labs <- paste0(classes, " (AUC = ", round(AUC, 2), ")")
  legend(0.7, 0.3, legend=c(legend.labs),
         col=c(pretty_colours), lty = 19, cex=0.9, box.lty=0)
  
  
  dev.off()
  
} # for (a in 1:length(trainingSets))

mergedValResults <- data.frame(AssignmentMonth, TimeChanges, ntree_parameter, mtry_parameter, Accuracy, AUC.df)
mergedValResults$Mean.AUC <- rowMeans(AUC.df)
file.path <- paste0(outputDir, outputName, "_MachineLearningValidation_Merged.tsv")
write.table(mergedValResults, file=file.path, sep="\t", row.names=FALSE)




