#####- Libraries ------------------------------------------------#####
library(ROCR)
library(randomForest)

#####- Import data ----------------------------------------------#####
rm(list=ls())

data("iris")

#####- Set random seed to make results reproducible: ------------#####
set.seed(17)

# # I want to form predictions on the setosa species so I will extract those observations for validation2
# setosaRows <- which(iris$Species == "setosa")
# validation2 <- iris[setosaRows,]
# iris2 <- iris[-(setosaRows),]

#####- Calculate the size of each of the data sets: -------------#####
data_set_size <- floor(nrow(iris)/2)

#####- Generate a random sample of "data_set_size" indexes ------#####
indexes <- sample(1:nrow(iris), size = data_set_size)

#####- Assign the data to the correct sets ----------------------#####
training <- iris[indexes,]
print(table(training$Species))

validation1 <- iris[-indexes,]
print(table(validation1$Species))

# I try to keep e^Nf < No (Nf = number of features, No = number of observations) to minimise overfitting
# No = data_set_size
Nf <- floor(log(data_set_size))

#####- Perform training: ----------------------------------------#####
rf_classifier <- randomForest(Species ~ ., # "Species ~ ." we want to predict Species usings each of the remaining columns of data
                              data=training, 
                              ntree=500, # defines the number of trees to be generated; it is typical to test a range of values for this parameter (i.e. 100, 200, 300, 400, 500) and choose the one that minimises the OOB estimate of error rate
                              mtry=2, # the number of features used in the construction of the tree; these features are selected at random, which is where the "random" in randomForests comes from; the default is sqrt(Nf)
                              importance=TRUE) # enables the algorithm to calculate variable importance

rf_classifier
##--------------------------------------------------------------------
# Call:
#   randomForest(formula = Species ~ ., data = training, ntree = 100,      mtry = 2, importance = TRUE) 
# Type of random forest: classification
# Number of trees: 100
# No. of variables tried at each split: 2
# 
# OOB estimate of  error rate: 2.67%
# Confusion matrix:
#   setosa versicolor virginica class.error
# setosa         24          0         0  0.00000000
# versicolor      0         24         0  0.00000000
# virginica       0          2        25  0.07407407
##--------------------------------------------------------------------

# OOB estimate is calculated by counting the # of misclassified points in the training set (2 virginica observations = 2) and dividing by the total number of observations (2/75 ~= 2.67%)
# OOB estimate is a useful measure to discriminate b/w different random forest classifiers
# You could vary the # of trees or the # of variables considered and select the combination with the smallest value
# Datasets with higher Nf - good idea to use cross-validation to perform feature selection using OOB error rate


#####- Assess importance ----------------------------------------#####
# Look at the importance that our classifier has assigned to each variable 
varImpPlot(rf_classifier)

# Feature's importance assessed based on:
# MeanDecreaseAccuracy - rough estimate of the loss in prediction performance when that variable is omitted from training
# MeanDecreaseGini - GINI is a measure of node impurity; if you used this feature to split the data, how pure will the nodes be?
## Highest purity = each not contains only elements of a single class
## Assessing the decrease in GINI when that feature is omitted leads to an understanding of how important that feature is

#####- Validation set assessment #1: looking at confusion matrix #####
prediction_for_table <- predict(rf_classifier,validation1[,-5]) # column 5 is the Species column we are predicting

table(observed=validation1[,5],predicted=prediction_for_table)

 ##-------------------------------------------------------------------
# predicted
# observed     setosa versicolor virginica
# setosa         26          0         0
# versicolor      0         23         3
# virginica       0          1        22
##-------------------------------------------------------------------

#####- Validation set assessment #2: ROC curves and AUC ---------#####
# Calculate the probability of new observations belonging to each class
# prediction_for_roc_curve will be a matrix with dimensions data_set_size x number_of_classes
prediction_for_roc_curve <- predict(rf_classifier,validation1[,-5],type="prob")

# Use pretty colours:
pretty_colours <- c("#F8766D","#00BA38","#619CFF")

# Specify the different classes 
classes <- levels(validation1$Species)

# For each class
for (i in 1:length(classes)) {
  # Define which observations belong to class[i]
  true_values <- ifelse(validation1[,5]==classes[i],1,0)
  
  # Assess the performance of classifier for class[i]
  pred <- prediction(prediction_for_roc_curve[,i],true_values)
  perf <- performance(pred, "tpr", "fpr")
  if (i==1)
  {
    plot(perf,main="ROC Curve",col=pretty_colours[i]) 
  }
  else
  {
    plot(perf,main="ROC Curve",col=pretty_colours[i],add=TRUE) 
  }
  
  # Add a legend
  legend(0.75, 0.25, legend=c(classes),
         col=c(pretty_colours), lty = 19, cex=0.8)
  
  # Calculate the AUC and print it to screen
  auc.perf <- performance(pred, measure = "auc")
  AUC <- auc.perf@y.values[[1]]
  print(paste0(classes[i]))
  print(paste0("AUC = ", AUC))
}

##- ntree = 100 ------------------------------------------------------
# [1] "setosa"
# [1] "AUC = 1"
# [1] "versicolor"
# [1] "AUC = 0.991365777080063"
# [1] "virginica"
# [1] "AUC = 0.99247491638796"

##- ntree = 200 ------------------------------------------------------
# [1] "setosa"
# [1] "AUC = 1"
# [1] "versicolor"
# [1] "AUC = 0.992150706436421"      *
# [1] "virginica"
# [1] "AUC = 0.99247491638796"

##- ntree = 300 ------------------------------------------------------
# [1] "setosa"
# [1] "AUC = 1"
# [1] "versicolor"
# [1] "AUC = 0.991365777080063"
# [1] "virginica"
# [1] "AUC = 0.99247491638796"

##- ntree = 400 ------------------------------------------------------
# [1] "setosa"
# [1] "AUC = 1"
# [1] "versicolor"
# [1] "AUC = 0.992150706436421"      *
# [1] "virginica"
# [1] "AUC = 0.99247491638796"

##- ntree = 500 ------------------------------------------------------
# [1] "setosa"
# [1] "AUC = 1"
# [1] "versicolor"
# [1] "AUC = 0.991365777080063"
# [1] "virginica"
# [1] "AUC = 0.99247491638796"