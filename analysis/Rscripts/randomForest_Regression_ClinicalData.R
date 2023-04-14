#####- Libraries ------------------------------------------------#####
library(ROCR)
library(randomForest)

#####- Import data ----------------------------------------------#####
rm(list=ls())
groups <- 2
level <- "genus"

features <- "Clinical Data"
outputName <- "ClinicalData"
fileName <- paste0("~/BioLockJ_pipelines/microbiome_n124_analysis_2023Jan30/07_MetaPhlAn2_TaxaMetaMerge/output/", level, "_LogNormalizedCounts_MetaPhlAn2.tsv"); countNorm <- "LogNorm"
myT <- read.delim(fileName, sep="\t",header = TRUE)

# fileName <- paste0("~/BioLockJ_pipelines/microbiome_n124_analysis_2023Jan30/08_MetaPhlAn2_Taxa_over_time_24M_Weight_Class_BLto24M/output/", level, "/Tertile/", level,  "_by_Timepoint_Tertile_LM_pValues_BLto24M_MetaPhlAn2.tsv")
# taxaNames <- read.delim(fileName, sep="\t",header = TRUE); features <- " Non-rare (>10%)"; outputName <- "_Non-rare"
# taxaNames <- taxaNames[taxaNames$adjPValues < 0.05,]; features <- " Significant change over time"; outputName <- "_Significant"
# taxaNames <- taxaNames[taxaNames$LengthTime == "BLto1M",]; features <- " Significant change from BL to 1M"; outputName <- "_Significant_BLto1M"
# taxaNames <- taxaNames[order(taxaNames$divPValues),]
# taxaNames <- unique(taxaNames$qbugName)
# taxaNames <- taxaNames[1:10]; features <- " Significant change from BL to 1M (Top 10)"; outputName <- "_Top10_Significant_BLto1M"

# features <- paste0(countNorm, " ", level, features)
# outputName <- paste0(countNorm, "_", level, outputName)
outputDir <- "~/Desktop/MachineLearning/"
dir.create(outputDir, showWarnings = FALSE)

#####- Set random seed to make results reproducible: ------------#####
set.seed(200)

#####- Add weight loss assignment -------------------------------#####
# trainingVariable <- "PEWL_kg"; varLabel <- "%EWL"
trainingVariable <- "Percent_Loss_kg"; varLabel <- "%WL"
clinicalVars <- c("PatientID", "Age", "Site", "Sex", "Ethnicity", "Surgery", "Loss_from_BL_BMI", "Percent_Loss_kg", "PEWL_kg")
clinicalVars <- clinicalVars[-which(clinicalVars == trainingVariable)]

meta.df <- myT[which(myT$Timepoint == "ONE"), clinicalVars ]
WL_24M <- myT[which(myT$Timepoint == "TWENTY_FOUR"), c("PatientID", trainingVariable) ]

df <- merge(WL_24M, meta.df, by = "PatientID")
# df <- merge(WL_24M, taxa.df, by = "PatientID")

df <- df[!is.na(df[,2]),]

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
rownames(training) <- training[,1]
training <- training[,-1]
colnames(training)[1] <- "trainingVariable"
Nf <- length(2:ncol(training))
training_N <- nrow(training)

validation1 <- df[-(which(df$PatientID %in% randomIDs)),]
rownames(validation1) <- validation1[,1]
validation1 <- validation1[,-1]
validation1_N <- nrow(validation1)
colnames(validation1)[1] <- "trainingVariable"


#####- Perform training: ----------------------------------------#####
ntreeValues <- c(100, 200, 300, 400, 500)
mtryValues <- 1:floor(sqrt(Nf))

ntree_parameter <- vector()
mtry_parameter <- vector()
Adj_r_sq <- vector()
index <- 1
print(paste0("# of features: ", Nf))

file.path <- paste0(outputDir, "Actual_v_Predicted_", trainingVariable, "_", outputName, ".pdf")
pdf(file.path, width = 5, height = 5)
for (ntree in ntreeValues) {
  
  for (mtry in mtryValues) {
    
    rf_regression <- randomForest(trainingVariable ~ ., # "PEWL_kg ~ ." we want to predict PEWL_kg usings each of the remaining columns of data
                                  data=training, 
                                  ntree=ntree, # defines the number of trees to be generated; it is typical to test a range of values for this parameter (i.e. 100, 200, 300, 400, 500) and choose the one that minimises the OOB estimate of error rate
                                  mtry=mtry, # the number of features used in the construction of the tree; these features are selected at random, which is where the "random" in randomForests comes from; the default is sqrt(Nf)
                                  importance=TRUE) # enables the algorithm to calculate variable importance
    
    
    trainPredictions <- rf_regression$predicted
    trainActual <- rf_regression$y
    
    title.lab <- paste0(features, " (Nf = ", Nf, ")\nTraining set (No = ", training_N, ") Validation set (No = ", validation1_N, ")\nntree = ", ntree, ", mtry = ", mtry)
    # varImpPlot(rf_regression, main = title.lab)
    fit <- summary(lm(trainPredictions~trainActual))
    Adj_r_sq <- fit$adj.r.squared
    title.lab <- paste0(features, " (Nf = ", Nf, ")\nTraining set (No = ", training_N, ") adj r^2 = ", round(Adj_r_sq, 3), "\nntree = ", ntree, ", mtry = ", mtry)
    
    plot(trainActual, trainPredictions, main = title.lab,
         xlab = paste0("Actual ", varLabel), ylab = paste0("Predicted ", varLabel))
    abline(lm(trainPredictions~trainActual), col="red", lty = 20) # regression line (y~x)
    
    ##- Validation set assessment #1: looking at confusion matrix ----
    Predicted <- predict(rf_regression, validation1[,-1]) 
    
    PatientID <- rownames(validation1)
    Actual <- validation1[,1]
    
    df2 <- data.frame(PatientID, Actual, Predicted)
    fit <- summary(lm(Predicted~Actual, data = df2))
    Adj_r_sq <- fit$adj.r.squared
    
    title.lab <- paste0(features, " (Nf = ", Nf, ")\nValidation set (No = ", validation1_N, ") adj r^2 = ", round(Adj_r_sq, 3), "\nntree = ", ntree, ", mtry = ", mtry)
    plot(Actual, Predicted, main = title.lab,
         xlab = paste0("Actual ", varLabel), ylab = paste0("Predicted ", varLabel))
    abline(lm(Predicted~Actual), col="red", lty = 20) # regression line (y~x)
    # abline(0,1,col="red", lty = 20)
    
    index <- index + 1
  }
}
dev.off()

