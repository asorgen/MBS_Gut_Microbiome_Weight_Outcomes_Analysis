#Author: Alicia Sorgen
#Date: 2023 Feb 23
#Description: Machine learning model for bariatric surgery weight loss outcomes.

#####- Libraries ------------------------------------------------#####
library(ROCR)
library(randomForest)

#####- Edits for script -----------------------------------------#####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
params <- c(params, "MetaPhlAn2")
# params <- c(params, 2)
params <- c(params, 3)
params <- c(params, "TaxaOnly")
# params <- c(params, "ClinicalOnly")
# params <- c(params, "Taxa_and_Clinical")
# params <- c(params, 1)
# params <- c(params, 6)
# params <- c(params, 12)
# params <- c(params, 18)
# params <- c(params, 24)
# end_month <- params[length(params)]

moduleRoot <- "Classification_randomForest_Prediction"
included <- c(0, 1, 6, 12, 18, 24)

#####- Set up working environment -------------------------------#####
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
  
  if (length(args) > 2) {
    m1 <- sapply(strsplit(module, "_"), "[", 1)
    m2 <- sapply(strsplit(module, "_"), "[", 2)
    m3 <- sapply(strsplit(module, "_"), "[", 3)
    
    if (args[3] == 2 | args[3] == 3) {
      module <- paste0("TertileGroups", args[3], "_", m2, "_", "Prediction")
    }
    module <- paste0(module, "_", args[4])
    
  }
  
  if (exists("end_month") == TRUE) {
    mod1 <- paste0(sapply(strsplit(module, "_Weight"), "[", 1), "_"); mod1
    mod2 <- paste0("_Weight", sapply(strsplit(module, "_Weight"), "[", 2)); mod2
    module <- paste0(mod1, end_month, "M", mod2, "_BLto", included[length(included)], "M"); module
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

}

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())

#####- Set up functions file ------------------------------------#####
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]
rm(str)

#####- Set up input ---------------------------------------------#####
classifier <- args[2]

prevModule <- paste0(classifier, "_TaxaMetaMerge")
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
prevModule <- paste0(classifier, "_Taxa_over_time_24M_Weight_Class_BLto24M")
statsDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")

month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")
division <- c("Tertile")

#####- Set up output --------------------------------------------#####
outputDir = file.path(moduleDir,"output/")


#####- Import data ----------------------------------------------#####
groups <- args[3]
trainingInput <- args[4]
level <- "genus"

fileName <- paste0(inputDir, level, "_LogNormalizedCounts_MetaPhlAn2.tsv")
myT <- read.delim(fileName, sep="\t",header = TRUE)

fileName <- paste0(statsDir, level, "/", division, "/", level,  "_by_Timepoint_Tertile_LM_pValues_BLto24M_MetaPhlAn2.tsv")
taxaNames <- read.delim(fileName, sep="\t",header = TRUE); features <- " Non-rare (>10%)"; outputName <- "_Non-rare"
taxaNames <- taxaNames[taxaNames$adjPValues < 0.05,]; features <- " Significant change over time"; outputName <- "_Significant"
# taxaNames <- taxaNames[taxaNames$LengthTime == "BLto1M",]; features <- " Significant change from BL to 1M"; outputName <- "_Significant_BLto1M"
taxaNames <- taxaNames[order(taxaNames$divPValues),]
taxaNames <- unique(taxaNames$qbugName)
# taxaNames <- taxaNames[1:10]; features <- " Significant change from BL to 1M (Top 10)"; outputName <- "_Top10_Significant_BLto1M"

features <- paste0(level, features)
outputName <- paste0(level, outputName)

#####- Set random seed to make results reproducible: ------------#####
set.seed(200)

#####- Calculate taxa change (BL-1M): ---------------------------#####

taxaList <- list()
index <- 1
for (i in taxaNames) {
  
  Timepoint <- myT$Timepoint
  PatientID <- myT$PatientID
  TaxaCol <- myT[,which(colnames(myT) == i)]
  
  df <- data.frame(PatientID, Timepoint, TaxaCol)
  df <- spread(df, Timepoint, TaxaCol)
  taxaList[[index]] <- df$ONE - df$BL
  names(taxaList)[index] <- i
  index <- index + 1
  
} # for (i in taxaNames)

PatientID <- df$PatientID
taxa.df <- data.frame(PatientID)
for (i in 1:length(taxaList)) {
  
  taxa.df <- data.frame(taxa.df, taxaList[i])
  
}

taxa.df <- na.omit(taxa.df)

#####- Add weight loss assignment -------------------------------#####
clinicalVariables <- c("Age", "Site", "Sex", "Ethnicity", "Surgery", "Height_m", "Weight_kg", "BMI_kgm2", "Loss_from_BL_BMI", "Loss_from_BL_kg", "Percent_Loss_kg", "PEWL_kg")

if (trainingInput == "TaxaOnly") {
  meta.df <- myT[which(myT$Timepoint == "BL"), c("PatientID", "TertileAssignment_24M") ]
  df <- merge(meta.df, taxa.df, by = "PatientID", all = FALSE)
}

if (trainingInput == "ClinicalOnly") {
  df <- myT[which(myT$Timepoint == "BL"), c("PatientID", "TertileAssignment_24M", clinicalVariables) ]
}

if (trainingInput == "Taxa_and_Clinical") {
  meta.df <- myT[which(myT$Timepoint == "BL"), c("PatientID", "TertileAssignment_24M", clinicalVariables) ]
  df <- merge(meta.df, taxa.df, by = "PatientID", all = FALSE)
}


if (groups == 2) {
  df$TertileAssignment_24M <- ifelse(df$TertileAssignment_24M == 1, "Bottom 33%",
                                     ifelse(df$TertileAssignment_24M == 2, "Top 67%",
                                            ifelse(df$TertileAssignment_24M == 3, "Top 67%", NA)))
  outputName <- paste0(outputName, "_Bottom_3rd_v_Top2_3rds")
} else {
  df$TertileAssignment_24M <- ifelse(df$TertileAssignment_24M == 1, "Bottom 33%",
                                     ifelse(df$TertileAssignment_24M == 2, "Middle 33%",
                                            ifelse(df$TertileAssignment_24M == 3, "Top 33%", NA)))
  
}
df <- df[!is.na(df$TertileAssignment_24M),]
df$TertileAssignment_24M <- factor(df$TertileAssignment_24M)
factorLength <- length(unique(df$TertileAssignment_24M))

#####- Calculate the size of each of the data sets: -------------#####
data_set_size <- floor(nrow(df) * 2/3)

#####- Generate a random sample of "data_set_size" indexes ------#####
indexes <- sample(1:nrow(df), size = data_set_size)
randomIDs <- df$PatientID[indexes]

# indexes <- sample(3:ncol(df), size = Nf)

#####- Assign the data to the correct sets ----------------------#####
training <- df[which(df$PatientID %in% randomIDs),]
print("Training set weight loss breakdown")
print(table(training$TertileAssignment_24M))
rownames(training) <- training[,1]
training <- training[,-1]
Nf <- length(2:ncol(training))
training_N <- nrow(training)

validation1 <- df[-(which(df$PatientID %in% randomIDs)),]
print("Validation set weight loss breakdown")
print(table(validation1$TertileAssignment_24M))
rownames(validation1) <- validation1[,1]
validation1 <- validation1[,-1]
validation1_N <- nrow(validation1)


#####- Perform training: ----------------------------------------#####
# ntreeValues <- c(100, 200, 300, 400, 500)
# mtryValues <- 1:floor(sqrt(Nf))

ntreeValues <- 100
mtryValues <- 1

print(paste0("# of features: ", Nf))
      
ntree_parameter <- vector()
mtry_parameter <- vector()
OOB_estimate <- vector()
index <- 1
file.path <- paste0(outputDir, "VariableImportancePlot_", outputName, ".pdf")
pdf(file.path, width = 15, height = 10)
for (ntree in ntreeValues) {
  
  for (mtry in mtryValues) {
    
    rf_classifier <- randomForest(TertileAssignment_24M ~ ., # "TertileAssignment_24M ~ ." we want to predict TertileAssignment_24M usings each of the remaining columns of data
                                  data=training, 
                                  ntree=ntree, # defines the number of trees to be generated; it is typical to test a range of values for this parameter (i.e. 100, 200, 300, 400, 500) and choose the one that minimises the OOB estimate of error rate
                                  mtry=mtry, # the number of features used in the construction of the tree; these features are selected at random, which is where the "random" in randomForests comes from; the default is sqrt(Nf)
                                  importance=TRUE) # enables the algorithm to calculate variable importance
    
    
    # rf_classifier
    OOB <- rf_classifier$confusion
    OOB <- OOB[,-ncol(OOB)]

    correctPredictions <- 0
    for (i in 1:factorLength) {
      correctPredictions <- correctPredictions + OOB[i,i]
    }
    title.lab <- paste0(features, " (n = ", Nf, ")\nntree = ", ntree, ", mtry = ", mtry, ", OOB estimate = ", round((sum(OOB) - correctPredictions) / sum(OOB) * 100, 2), "%")
    varImpPlot(rf_classifier, main = title.lab)
    
    OOB_estimate[index] <- (sum(OOB) - correctPredictions) / sum(OOB)
    ntree_parameter[index] <- ntree
    mtry_parameter[index] <- mtry
    index <- index + 1
  }
}
dev.off()

trainingParameters <- data.frame(ntree_parameter, mtry_parameter, OOB_estimate)
trainingParameters <- trainingParameters[order(trainingParameters$OOB_estimate),]

print(paste0("Lowest OOB estimate: ", trainingParameters$OOB_estimate[1]))

#####- Validation -----------------------------------------------#####
file.path <- paste0(outputDir, "ROC_Curve_", outputName, ".pdf")
pdf(file.path, width = 6, height = 5)

Accuracy <- vector()
AUC.df <- data.frame()
index <- 1
for (j in 1:nrow(trainingParameters)) {
  
  ntree <- trainingParameters$ntree_parameter[j]
  mtry <- trainingParameters$mtry_parameter[j]
  # print(paste0("ntree = ", ntree))
  # print(paste0("mtry = ", mtry))
  
  rf_classifier <- randomForest(TertileAssignment_24M ~ ., # "TertileAssignment_24M ~ ." we want to predict TertileAssignment_24M usings each of the remaining columns of data
                                data=training, 
                                ntree=ntree, # defines the number of trees to be generated; it is typical to test a range of values for this parameter (i.e. 100, 200, 300, 400, 500) and choose the one that minimises the OOB estimate of error rate
                                mtry=mtry, # the number of features used in the construction of the tree; these features are selected at random, which is where the "random" in randomForests comes from; the default is sqrt(Nf)
                                importance=TRUE) # enables the algorithm to calculate variable importance
  
  ##- Validation set assessment #1: looking at confusion matrix 
  prediction_for_table <- predict(rf_classifier,validation1[,-1]) # column 1 is the TertileAssignment_24M column we are predicting
  
  results <- cbind(prediction_for_table, validation1[,1])
  Accuracy[index] <- sum(prediction_for_table==validation1[,1]) / nrow(validation1)
  
  
  ##- Validation set assessment #2: ROC curves and AUC
  prediction_for_roc_curve <- predict(rf_classifier,validation1[,-1],type="prob")
  
  ##- Use pretty colours:
  pretty_colours <- c("#F8766D","#00BA38","#619CFF")
  
  ##- Specify the different classes 
  classes <- levels(validation1$TertileAssignment_24M)
  title.lab <- paste0(features, " (n = ", Nf, ")\nTraining set (", training_N, ") Validation set (", validation1_N, ")\nntree = ", ntree, ", mtry = ", mtry, ", Accuracy = ", round(Accuracy[index] * 100, 2), "%")
  
  AUC <- vector()
  
  for (i in 1:length(classes)) {
    ##- Define which observations belong to class[i]
    true_values <- ifelse(validation1[,1]==classes[i],1,0)
    
    ##- Assess the performance of classifier for class[i]
    pred <- prediction(prediction_for_roc_curve[,i],true_values)
    perf <- performance(pred, "tpr", "fpr")
    if (i==1)
    {
      plot(perf,main=title.lab,col=pretty_colours[i]) 
    }
    else
    {
      plot(perf,main=title.lab,col=pretty_colours[i],add=TRUE) 
    }
    
    ##- Calculate the AUC and print it to screen
    auc.perf <- performance(pred, measure = "auc")
    AUC[i] <- auc.perf@y.values[[1]]
    # print(paste0(classes[i]))
    # print(paste0("AUC = ", AUC[i]))
  }
  
  AUC.df <- rbind(AUC.df, AUC)
  for (i in 1:length(classes)) {
    names(AUC.df)[i] <- classes[i]
  }
  
  ##- Add a legend
  legend.labs <- paste0(classes, " (AUC = ", round(AUC, 2), ")")
  legend(0.5, 0.3, legend=c(legend.labs),
         col=c(pretty_colours), lty = 19, cex=0.9, box.lty=0)
  index <- index + 1
  
} # for (j in 1:nrow(trainingParameters))


dev.off()

trainingParameters$Accuracy <- Accuracy
trainingParameters <- cbind(trainingParameters, AUC.df)

file.path <- paste0(outputDir, "TrainingParameters_", outputName, ".tsv")
write.table(trainingParameters, file.path,sep="\t",quote = FALSE, row.names = FALSE)
