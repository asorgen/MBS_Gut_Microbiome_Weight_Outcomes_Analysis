#####- Import data ----------------------------------------------#####
rm(list=ls())
groups <- 2
level <- "genus"

features <- "Clinical Data"
outputName <- "ClinicalData"
fileName <- paste0("~/BioLockJ_pipelines/microbiome_n124_analysis_2023Jan30/07_MetaPhlAn2_TaxaMetaMerge/output/", level, "_LogNormalizedCounts_MetaPhlAn2.tsv"); countNorm <- "LogNorm"
myT <- read.delim(fileName, sep="\t",header = TRUE)

# features <- paste0(countNorm, " ", features)
# outputName <- paste0(countNorm, "_", outputName)
outputDir <- "~/Desktop/MachineLearning/"
dir.create(outputDir, showWarnings = FALSE)

#####- Set random seed to make results reproducible: ------------#####
set.seed(200)

#####- Add weight loss assignment -------------------------------#####
df <- myT[which(myT$Timepoint == "ONE"), c("PatientID", "TertileAssignment_24M", "Age", "Site", "Sex", "Ethnicity", "Surgery", "Loss_from_BL_BMI", "PEWL_kg", "Percent_Loss_kg") ]

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

weightColStart <- which(colnames(df) == "Surgery") + 1

for (i in weightColStart:ncol(df)) {
  for (j in 1:nrow(df)) {
    if (is.na(df[j,i]) == TRUE) {
      df[j,i] <- mean(na.omit(df[,i]))
    }
  }
}

#####- Calculate the size of each of the data sets: -------------#####
data_set_size <- floor(nrow(df) * 2/3)

#####- Generate a random sample of "data_set_size" indexes ------#####
indexes <- sample(1:nrow(df), size = data_set_size)
randomIDs <- df$PatientID[indexes]

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
ntreeValues <- c(100, 200, 300, 400, 500)
mtryValues <- 1:floor(sqrt(Nf))


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
file.path <- paste0(outputDir, "TrainingParameters_", outputName, ".tsv")
write.table(trainingParameters, file.path,sep="\t",quote = FALSE, row.names = FALSE)

print(paste0("Lowest OOB estimate: ", trainingParameters$OOB_estimate[1]))

#####- Validation -----------------------------------------------#####
file.path <- paste0(outputDir, "ROC_Curve_", outputName, ".pdf")
pdf(file.path, width = 6, height = 5)
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
  
  ##- Validation set assessment #1: looking at confusion matrix ------
  prediction_for_table <- predict(rf_classifier,validation1[,-1]) # column 1 is the TertileAssignment_24M column we are predicting
  # table(observed=validation1[,1],predicted=prediction_for_table)
  
  ##- Validation set assessment #2: ROC curves and AUC ---------------
  prediction_for_roc_curve <- predict(rf_classifier,validation1[,-1],type="prob")
  
  # Use pretty colours:
  pretty_colours <- c("#F8766D","#00BA38","#619CFF")
  
  # Specify the different classes 
  classes <- levels(validation1$TertileAssignment_24M)
  title.lab <- paste0(features, " (n = ", Nf, ")\nTraining set (", training_N, ") Validation set (", validation1_N, ")\nntree = ", ntree, ", mtry = ", mtry, ", OOB estimate = ", round(trainingParameters$OOB_estimate[j] * 100, 2), "%")
  
  AUC <- vector()
  
  for (i in 1:length(classes)) {
    ##- Define which observations belong to class[i] -----------------
    true_values <- ifelse(validation1[,1]==classes[i],1,0)
    
    ##- Assess the performance of classifier for class[i] ------------
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
    
    ##- Calculate the AUC and print it to screen ---------------------
    auc.perf <- performance(pred, measure = "auc")
    AUC[i] <- auc.perf@y.values[[1]]
    # print(paste0(classes[i]))
    # print(paste0("AUC = ", AUC[i]))
  }
  
  ##- Add a legend ---------------------------------------------------
  legend.labs <- paste0(classes, " (AUC = ", round(AUC, 2), ")")
  legend(0.5, 0.3, legend=c(legend.labs),
         col=c(pretty_colours), lty = 19, cex=0.9, box.lty=0)
  
}

dev.off()
