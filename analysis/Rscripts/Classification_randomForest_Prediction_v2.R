#Author: Alicia Sorgen
#Date: 2023 Feb 23
#Description: Machine learning model for bariatric surgery weight loss outcomes.

#####- Libraries ------------------------------------------------#####
R <- sessionInfo()
message(R$R.version$version.string)

library(stringr); message(">>> stringr: Version ", packageVersion("stringr"), "\n")
library(ROCR); message(">>> ROCR: Version ", packageVersion("ROCR"), "\n")
library(randomForest); message(">>> randomForest: Version ", packageVersion("randomForest"), "\n")
library(rstatix); message(">>> rstatix: Version ", packageVersion("rstatix"), "\n")

#####- Edits for script -----------------------------------------#####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/git/gut-microbiome-bariatric-weight-outcomes")
params <- c(params, "MetaPhlAn2")
params <- c(params, 2)
# params <- c(params, 3)
# params <- c(params, "TaxaOnly")
# params <- c(params, "ClinicalOnly")
# params <- c(params, "Taxa_and_Clinical")
params <- c(params, "Taxa_Clinical_Energy")
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
taxaNames <- read.delim(fileName, sep="\t",header = TRUE); outputName <- "_Non-rare"
taxaNames <- taxaNames[taxaNames$adjPValues < 0.05,]; outputName <- "_Significant"
# taxaNames <- taxaNames[taxaNames$LengthTime == "BLto1M",]; features <- " Significant change from BL to 1M"; outputName <- "_Significant_BLto1M"
taxaNames <- taxaNames[order(taxaNames$divPValues),]
taxaNames <- unique(taxaNames$qbugName)
# taxaNames <- taxaNames[1:10]; features <- " Significant change from BL to 1M (Top 10)"; outputName <- "_Top10_Significant_BLto1M"

features <- paste0(level)
outputName <- paste0(level, outputName)

#####- Set random seed to make results reproducible: ------------#####
set.seed(200)

#####- Calculate taxa change: -----------------------------------#####
Timepoints <- c("BL", "ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR")
taxa.dfList <- list()
listIndex <- 1
for (endM in 2:length(included)) {
  
  taxaList <- list()
  index <- 1
  for (i in taxaNames) {
    
    Timepoint <- myT$Timepoint
    PatientID <- myT$PatientID
    TaxaCol <- myT[,which(colnames(myT) == i)]
    
    df <- data.frame(PatientID, Timepoint, TaxaCol)
    df <- spread(df, Timepoint, TaxaCol)
    taxaList[[index]] <- df[,which(colnames(df) == Timepoints[endM])] - df$BL
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
  listIndex <- listIndex + 1
  
} # for (endM in 2:length(included))

#####- Training -------------------------------------------------#####
main.df <- data.frame()
for (a in 2:length(included)) {
  
  assignmentMonth <- included[a]
  outputName2 <- paste0(assignmentMonth, "M_WLG")
  print(paste0(">>>>> ", assignmentMonth))
  
  #####- Add weight loss assignment
  clinicalVariables <- c("Age", "Site", "Sex", "Ethnicity", "Surgery", "Loss_from_BL_BMI", "Loss_from_BL_kg", "Percent_Loss_kg", "PEWL_kg")
  
  for (x in 1:length(taxa.dfList)) {
    
    endM <- included[x+1]
    TimeChange <- paste0("BL - ", endM, "M")
    outputName3 <- paste0(outputName2, "_BLto", endM, "M")
    print(paste0("     >>>>> ", TimeChange))
    
    if (trainingInput == "TaxaOnly") {
      meta.df <- myT[which(myT$Timepoint == Timepoints[endM]), c("PatientID", paste0("TertileAssignment_", assignmentMonth, "M")) ]
      df <- merge(meta.df, taxa.dfList[[x]], by = "PatientID", all = FALSE)
      names(df)[names(df) == paste0("TertileAssignment_", assignmentMonth, "M")] <- "WeightLossGroup"
    } # if (trainingInput == "TaxaOnly")
    
    if (trainingInput == "ClinicalOnly") {
      df <- myT[which(myT$Timepoint == Timepoints[x+1]), c("PatientID", paste0("TertileAssignment_", assignmentMonth, "M"), clinicalVariables) ]
      names(df)[names(df) == paste0("TertileAssignment_", assignmentMonth, "M")] <- "WeightLossGroup"
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
    
    #####- Calculate the size of each of the data sets
    splitSet <- 3
    data_set_size <- floor(nrow(df) * 1/splitSet)
    
    #####- Generate a random sample of "data_set_size" indexes
    randomIDs <- list()
    PatientIDs <- df$PatientID
    for (nSet in 1:splitSet) {
      
      if (nSet != splitSet) {
        indexes <- sample(PatientIDs, size = data_set_size)
        randomIDs[[nSet]] <- indexes
        PatientIDs <- PatientIDs[-which(PatientIDs %in% indexes)]
      } else {
        indexes <- PatientIDs
        randomIDs[[nSet]] <- indexes
        PatientIDs <- PatientIDs[-which(PatientIDs %in% indexes)]
        
      }
    }
    
    #####- Assign the data to the correct sets
    trainingSets <- list()
    validationSets <- list()
    for (id in 1:length(randomIDs)) {
      
      IDs <- randomIDs[[id]]
      training <- df[!(df$PatientID %in% IDs),]
      rownames(training) <- training[,1]
      training <- training[,-1]
      trainingSets[[id]] <- training
      
      validation <- df[df$PatientID %in% IDs,]
      rownames(validation) <- validation[,1]
      validation <- validation[,-1]
      validation_N <- nrow(validation)
      validationSets[[id]] <- validation
      
    } # for (id in 1:length(randomIDs))
    
    #####- Perform training
    # ntreeValues <- c(100, 200, 300, 400, 500)
    # mtryValues <- 1:floor(sqrt(Nf))
    
    ntreeValue <- 100
    mtryValue <- 1
    AUC.df <- data.frame()
    file.path <- paste0(outputDir, "VariableImportancePlot_", outputName, "_", outputName3, ".pdf")
    pdf(file.path, width = 15, height = 10)
    
    for (set in 1:length(trainingSets)) {
      
      print(paste0("          >>>>> Training #", set))
      # outputName4 <- paste0(outputName3, "_Set", set)
      training <- trainingSets[[set]]
      Nf <- length(2:ncol(training))
      training_N <- nrow(training)
      
      validation <- validationSets[[set]]
      validation_N <- nrow(validation)
      
      rf_classifier <- randomForest(WeightLossGroup ~ ., # "WeightLossGroup ~ ." we want to predict WeightLossGroup usings each of the remaining columns of data
                                    data=training,
                                    ntree=ntreeValue, # defines the number of trees to be generated; it is typical to test a range of values for this parameter (i.e. 100, 200, 300, 400, 500) and choose the one that minimises the OOB estimate of error rate
                                    mtry=mtryValue, # the number of features used in the construction of the tree; these features are selected at random, which is where the "random" in randomForests comes from; the default is sqrt(Nf)
                                    importance=TRUE) # enables the algorithm to calculate variable importance
      
      
      # rf_classifier
      OOB <- rf_classifier$confusion
      OOB <- OOB[,-ncol(OOB)]
      
      correctPredictions <- 0
      for (i in 1:factorLength) {
        correctPredictions <- correctPredictions + OOB[i,i]
      } # for (i in 1:factorLength)
      title.lab <- paste0(features, " (n = ", Nf, ")\nntree = ", ntreeValue, ", mtry = ", mtryValue, ", OOB estimate = ", round((sum(OOB) - correctPredictions) / sum(OOB) * 100, 2), "%")
      varImpPlot(rf_classifier, main = title.lab)
      
      OOB_estimate <- (sum(OOB) - correctPredictions) / sum(OOB)
      ntree_parameter <- ntreeValue
      mtry_parameter <- mtryValue
      
      if (set == 1) {
        trainingParameters <- data.frame(assignmentMonth, TimeChange, ntree_parameter, mtry_parameter, OOB_estimate)
        names(trainingParameters)[names(trainingParameters) == "OOB_estimate"] <- paste0("OOB_estimate_Train", set)
        
        # print(paste0("# of features: ", Nf))
        # print(paste0("Training set #", set, " (n = ", training_N, ")"))
        # print(paste0("Validation set #", set, " (n = ", validation_N, ")"))
        
      } else {
        trainingParameters$OOB_estimate <- OOB_estimate
        names(trainingParameters)[names(trainingParameters) == "OOB_estimate"] <- paste0("OOB_estimate_Train", set)
        
        # print(paste0("Training set #", set, " (n = ", training_N, ")"))
        # print(paste0("Validation set #", set, " (n = ", validation_N, ")"))
        
      } # if (set == 1)
      
      # file.path <- paste0(outputDir, "ROC_Curve_", outputName, "_", outputName3, ".pdf")
      # pdf(file.path, width = 6, height = 5)
      
      print(paste0("          >>>>> Validation #", set))
      ##- Validation set assessment #1: looking at confusion matrix 
      prediction_for_table <- predict(rf_classifier,validation[,-1]) # column 1 is the WeightLossGroup column we are predicting
      
      results <- cbind(prediction_for_table, validation[,1])
      Accuracy <- sum(prediction_for_table==validation[,1]) / nrow(validation)
      trainingParameters$Accuracy <- Accuracy
      names(trainingParameters)[names(trainingParameters) == "Accuracy"] <- paste0("Accuracy_Train", set)
      
      
      ##- Validation set assessment #2: ROC curves and AUC
      prediction_for_roc_curve <- predict(rf_classifier,validation[,-1],type="prob")
      
      if (set == 1) {
        # mergedResults <- cbind(results, prediction_for_roc_curve)
        prediction_for_table_Merged <- prediction_for_table
        prediction_for_roc_curve_Merged <- prediction_for_roc_curve
        
      } else {
        # mergedResultsx <- cbind(results, prediction_for_roc_curve)
        # mergedResults <- rbind(mergedResults, mergedResultsx)
        prediction_for_table_Merged <- c(prediction_for_table_Merged, prediction_for_table)
        prediction_for_roc_curve_Merged <- rbind(prediction_for_roc_curve_Merged, prediction_for_roc_curve)
      }
      
      # ##- Use pretty colours:
      # pretty_colours <- c("#F8766D","#00BA38","#619CFF")
      
      ##- Specify the different classes
      classes <- levels(validation$WeightLossGroup)
      
      # title.lab <- paste0(features, " ", assignmentMonth, "-month WLG BL to ", endM, "M (Nf = ", Nf, ")\nTraining set ", set, " (", training_N, ") Validation set ", set, " (", validation_N, ")\nntree = ", ntreeValue, ", mtry = ", mtryValue, ", Accuracy = ", round(Accuracy * 100, 2), "%")
      
      AUC <- vector()
      
      for (i in 1:length(classes)) {
        ##- Define which observations belong to class[i]
        true_values <- ifelse(validation[,1]==classes[i],1,0)
        
        ##- Assess the performance of classifier for class[i]
        pred <- prediction(prediction_for_roc_curve[,i],true_values)
        perf <- performance(pred, "tpr", "fpr")
        #   if (i==1)
        #   {
        #     plot(perf,main=title.lab,col=pretty_colours[i]) 
        #   }
        #   else
        #   {
        #     plot(perf,main=title.lab,col=pretty_colours[i],add=TRUE) 
        #   }
        
        ##- Calculate the AUC and print it to screen
        auc.perf <- performance(pred, measure = "auc")
        AUC <- auc.perf@y.values[[1]]
        trainingParameters$AUC <- AUC
        names(trainingParameters)[names(trainingParameters) == "AUC"] <- paste0("AUC_Train", set, "_", classes[i])
        
        # print(paste0(classes[i]))
        # print(paste0("AUC = ", AUC[i]))
      } # for (i in 1:length(classes))
      
      # AUC.df <- rbind(AUC.df, AUC)
      # for (i in 1:length(classes)) {
      #   names(AUC.df)[i] <- classes[i]
      # }
      
      # ##- Add a legend
      # legend.labs <- paste0(classes, " (AUC = ", round(AUC, 2), ")")
      # legend(0.5, 0.3, legend=c(legend.labs),
      #        col=c(pretty_colours), lty = 19, cex=0.9, box.lty=0)
      # dev.off()
      
    } # for (set in 1:length(trainingSets))
    dev.off()
    
    file.path <- paste0(outputDir, "ROC_Curve_", outputName, "_", outputName3, ".pdf")
    pdf(file.path, width = 6, height = 5)
    Accuracy <- sum(prediction_for_table_Merged==df[,2]) / nrow(df)
    trainingParameters$Accuracy <- Accuracy
    
    ##- Use pretty colours:
    pretty_colours <- c("#F8766D","#00BA38","#619CFF")
    
    ##- Specify the different classes 
    classes <- levels(df$WeightLossGroup)
    
    
    AUC <- vector()
    
    for (i in 1:length(classes)) {
      ##- Define which observations belong to class[i]
      true_values <- ifelse(df[,2]==classes[i],1,0)
      
      ##- Assess the performance of classifier for class[i]
      pred <- prediction(prediction_for_roc_curve_Merged[,i],true_values)
      perf <- performance(pred, "tpr", "fpr")
      
      ##- Calculate the AUC and print it to screen
      auc.perf <- performance(pred, measure = "auc")
      AUC <- auc.perf@y.values[[1]]
      trainingParameters$AUC <- AUC
      names(trainingParameters)[names(trainingParameters) == "AUC"] <- paste0("AUC_", classes[i])
      
      title.lab <- paste0(features, " ", assignmentMonth, "-month WLG BL to ", endM, "M (Nf = ", Nf, ")\nMerged validation results\nAccuracy = ", round(Accuracy * 100, 2), "%; AUC = ", AUC)
      
      if (i==1)
      {
        plot(perf,main=title.lab,col=pretty_colours[i]) 
      }
      else
      {
        plot(perf,main=title.lab,col=pretty_colours[i],add=TRUE) 
      }
      
      
      # print(paste0(classes[i]))
      # print(paste0("AUC = ", AUC[i]))
    }
    
    # AUC.df <- rbind(AUC.df, AUC)
    # for (i in 1:length(classes)) {
    #   names(AUC.df)[i] <- classes[i]
    # }
    
    ##- Add a legend
    legend.labs <- paste0(classes)
    legend(0.5, 0.3, legend=c(legend.labs),
           col=c(pretty_colours), lty = 19, cex=0.9, box.lty=0)
    
    
    dev.off()
    
    
    main.df <- rbind(main.df, trainingParameters)
    
  } # for (x in 1:length(taxa.dfList))
  
} # for (a in 2:length(included))

file.path <- paste0(outputDir, "MachineLearningResults.tsv")
write.table(main.df, file=file.path, sep="\t", row.names=FALSE)






