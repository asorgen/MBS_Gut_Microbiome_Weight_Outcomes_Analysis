#' ---
#' title: "Tertile Analysis"
#' author: "Alicia Sorgen"
#' date: "`r Sys.Date()`"
#' output: 
#'   html_document: 
#'     theme: journal
#'     fig_width: 3
#'     fig_caption: yes
#'     fig_height: 3
#' ---
#' 
## ----include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo=FALSE)
rm(list=ls())
set.seed(1989)


#' 
#' 
#' 
#' # Setup
#' 
## ----Library, include=FALSE-----------------------------------
R <- sessionInfo()
message(R$R.version$version.string)

library(MASS); message("MASS: Version ", packageVersion("MASS"))
library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
library(lcmm); message("lcmm: Version ", packageVersion("lcmm"))
library(stringr); message("stringr: Version ", packageVersion("stringr"))
# library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
# library(scales); message("scales: Version ", packageVersion("scales"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))


#' 
#' ``r R$R.version$version.string``
#' 
#' **Packages:**
#' 
#' * `` ggpubr ``: Version ``r packageVersion("ggpubr")``
#' 
#' * `` nlme ``: Version ``r packageVersion("nlme")``
#' 
#' * `` MASS ``: Version ``r packageVersion("MASS")``
#' 
#' * `` ggplot2 ``: Version ``r packageVersion("ggplot2")``
#' 
#' * `` lcmm ``: Version ``r packageVersion("lcmm")``
#' 
#' * `` data.table ``: Version ``r packageVersion("data.table")``
#' 
#' * `` stringr ``: Version ``r packageVersion("stringr")``
#' 
#' * `` rstatix ``: Version ``r packageVersion("rstatix")``
#' 
## ----Script-Edits---------------------------------------------
rm(list=ls())

ANALYSIS <- "MetaPhlAn2_microbiome"

moduleRoot <- paste0("Tertile_Analysis")
params <- vector()
params <- vector()
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
# params <- c(params, "MetaPhlAn2")
# params <- c(params, "Kraken2")
# params <- c(params, 0)
# params <- c(params, 1)
# params <- c(params, 6)
# params <- c(params, 12)
# params <- c(params, 18)
params <- c(params, 24)
# end_month <- params[length(params)]
params <- c(params, "Tertile")


level <- "genus"


#' 
## ----Working-Environment, include=FALSE-----------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- params
}; rm(params)

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
    file.copy(gitInput,
              root,
              recursive = TRUE)
    
  }; rm(gitInput)
  
  module <- moduleRoot
  
  if (any(dir(root) %like% module) == TRUE) {
    moduleDir <- paste0(root,str_subset(dir(root), module), "/")
  } else {
    moduleDir <- paste0(root, module, "/")
    dir.create(moduleDir, showWarnings = FALSE)
  }; rm(root, module)
  
  scriptDir <- paste0(moduleDir, "script/")
  if (any(dir(moduleDir) == "script") == FALSE) {
    dir.create(scriptDir, showWarnings = FALSE)
  }; rm(scriptDir)
  
  outputDir <- paste0(moduleDir, "output/")
  if (any(dir(moduleDir) == "output") == FALSE) {
    dir.create(outputDir, showWarnings = FALSE)
  }
  
  files <- list.files(outputDir, recursive = TRUE, full.names = TRUE); rm(outputDir)
  file.remove(files); rm(files)
  
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
              recursive = TRUE); rm(script)
  }; rm(gitScripts, moduleRoot, resourcesDir)
  
  setwd(paste0(moduleDir, "script/")); rm(moduleDir)
  
}

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())


#' 
## ----Functions------------------------------------------------
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]
rm(str)

#' 
#' ### Input
## -------------------------------------------------------------
prevModule <- str_subset(dir(pipeRoot), "WeightMetaMerge")
inputDir = paste0(pipeRoot,"/",prevModule,"/output/")
inputFile <- "metadata.tsv"

endM <- args[2]
division <- args[length(args)]

prevModule2 <- paste0("TaxaMetaMerge")
inputDir2 = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule2),"/output/")
fileNames <- dir(inputDir2)
logCountFile <- fileNames[which(fileNames %like% paste0(level, "_LogNormalizedCounts"))]


#' 
#' **Input:** 
#' 
#' ``r prevModule``
#' 
#' * ``r inputFile``
#' 
#' ``r prevModule2``
#' 
#' * ``r logCountFile``
#' 
#' 
## -------------------------------------------------------------
outputDir = file.path(moduleDir,"output/")

#' 
#' ### Output
#' 
#' **Output directory:** ``r outputDir``
#' 
## ----Script-Variables-----------------------------------------
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

# outcomes <- c("Weight_kg", "BMI_kgm2", "Loss_from_BL_kg", "Percent_Loss_kg", "PEWL_kg")
# outcomeLabels <- c("Weight (kg)", "BMI (kg/m^2)", "Weight Loss from BL (kg)", "Weight Loss from BL (%)", "Excess Weight Loss from BL (%)")
outcomes <- c("Percent_Loss_kg", "Weight_kg")
outcomeLabel <- c("Weight Loss from BL (%)", "Weight_kg")
repeatedMeasure <- c("time")
covariates <- c("Surgery")
subjects <- c("SampleID", "PatientID")

#' 
#' 
#' # Dataset 
#' 
## ----Read-Table, echo=FALSE-----------------------------------
##### Read in table and data prep #####

# Read in table
filePath <- paste0(inputDir, inputFile)
myTable <- read.table(filePath, sep="\t", header = TRUE, check.names = FALSE)

filePath <- paste0(inputDir2, logCountFile)
taxa.df <- read.table(filePath, sep="\t", header = TRUE, check.names = FALSE)

divMonth <- paste0(division, "Assignment_", endM, "M")

# Trim variables
df.1 <- taxa.df[ , which(colnames(taxa.df) %in% c(outcomes, repeatedMeasure, covariates, subjects, divMonth))]
# file.path <- paste0(outputDir, "test.tsv")
# write.table(df.1, file.path, sep="\t",quote = FALSE, row.names = FALSE)


#' 
#' 
#' 
## ----Figure1-code, fig.height=5, fig.width=10, include=FALSE----

df.2 <- df.1[!is.na(df.1$TertileAssignment_24M),]
df.2$TertileAssignment_24M <- ifelse(df.2$TertileAssignment_24M == 1, 3,
                                     ifelse(df.2$TertileAssignment_24M == 3, 1,
                                            df.2$TertileAssignment_24M))

n <- length(unique(df.2$PatientID))
ng <- length(unique(df.2$TertileAssignment_24M))
message(ng)


Palette <- c("orange", "blue", "green3", "black", "red")

plot <- ggplot() +
  xlim(0, 36)+
  labs(x = "Follow-up time (months)", y = "Weight Change (%)") +
  theme_bw()

df.merged <- data.frame()
for (i in 1:ng) {
  
  df.C <- df.2[df.2$TertileAssignment_24M == i,]
  
  label <- paste0("Group ", i, " (n = ", length(unique(df.C$PatientID)), ", ", round((length(unique(df.C$PatientID))/n)*100,1), "%)")
  message(label)
  model <- lme( Weight_kg ~ time, method = "ML", random = ~1 | PatientID, data = df.C, na.action = na.omit )
  message("\n************* ", i)
  df.C$predicted_values <- predict(model, newdata = data.frame(time = df.C$time, PatientID = df.C$PatientID))
  df.C$modeled <- ifelse(is.na(df.C$Weight_kg) == TRUE, df.C$predicted_values, 
                          df.C$Weight_kg)
  
  # Find baseline weight for each PatientID at time 0
  # baseline_weights <- df.C %>%
  #   filter(time == 0) %>%
  #   select(PatientID, Baseline_kg = Weight_kg)
  # 
  # df.C <- merge(df.C, baseline_weights, by = "PatientID", all.x = TRUE)
  df.C$Baseline_kg <- assignAllbyTerm(df.C, "time", "0", "Weight_kg", "PatientID")
  
  df.C$`Observed - BL` <- df.C$Weight_kg - df.C$Baseline_kg
  df.C$`Modeled - BL` <- df.C$modeled - df.C$Baseline_kg
  
  df.C$`(Observed - BL)/BL` <- df.C$`Observed - BL` / df.C$Baseline_kg
  df.C$`(Modeled - BL)/BL` <- df.C$`Modeled - BL` / df.C$Baseline_kg
  
  df.C$observed_weight_change <- df.C$`(Observed - BL)/BL` * 100
  df.C$modeled_weight_change <- df.C$`(Modeled - BL)/BL` * 100
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
  
  
}


# print(plot)
filename.1 <- "Tertile_Modeled_and_Observed_Weight_Change.pdf"
filePath <- paste0(outputDir, filename.1)

pdf(filePath, width = 10, height = 5)
print(plot)
dev.off()
message("\n************* Setup Complete *************")



## ----Figure1-print,  echo=FALSE, fig.height=5, fig.width=7, warning=FALSE----
# print(plot)


#' 
#' **Figure 1.** Percent weight change trajectories for patients in classified in the top (Group 1), middle (Group 2), and bottom (Group 3) weight loss groups at 24-months (`r filename.1`).
#' 
#' 
#' 
## ----Figure2-code, fig.height=5, fig.width=10, include=FALSE----

df.2 <- df.2[df.2$time == 0,]

f.test <- fisher.test(df.2$Surgery, df.2$TertileAssignment_24M, workspace = 2e8)

freq <- as.data.frame(table(df.2$Surgery, df.2$TertileAssignment_24M))
names(freq)[names(freq) == "Var1"] <- "Surgery"
names(freq)[names(freq) == "Var2"] <- "Group"

n_RYGB <- sum(freq$Freq[which(freq$Surgery == "RYGB")])
n_SG <- sum(freq$Freq[which(freq$Surgery == "SG")])

freq$Percent <- ifelse(freq$Surgery == "RYGB", (freq$Freq / n_RYGB) * 100,
                       (freq$Freq / n_SG) * 100)

plot <- ggplot(freq, aes(fill=Group, y=Freq, x=Surgery)) + 
  geom_bar(position="stack", stat="identity", color = "black")+
  labs(y = "Number of patients")+
  scale_fill_manual(values=Palette[1:ng])

subtitle.lab <- paste0("Fisher test ", round(f.test$p.value, 3))
plot <- plot + labs(subtitle = subtitle.lab)

plot <- plot + geom_text(aes(label = paste0(round(Percent, 1), "%")), size = 3, hjust = 0.5, position = "stack", vjust = 1.5, colour = "black")
# print(plot)
filename.2 <- "Tertile_Groups_Surgery_Fisher.pdf"
filePath <- paste0(outputDir, filename.2)

pdf(filePath, width = 10, height = 5)
print(plot)
dev.off()


## ----Figure2-print, echo=FALSE, fig.height=5, fig.width=7, warning=FALSE----
# print(plot)


#' 
#' **Figure 2.** Weight loss groups per surgery type (`r filename.2`).
#' 
#' 
## ----Overall-MLM, fig.height=5, fig.width=10, warning=FALSE, include=FALSE----
months <- c(0, 1, 6, 12, 18, 24)

outputResults <- paste0(outputDir, "ResultsTables/")
dir.create(outputResults, showWarnings = FALSE)

startAbundanceIndex <- which(colnames(taxa.df)=="ResponderStatus")+1
SampleID <- taxa.df$SampleID
taxa.df2 <- taxa.df[,startAbundanceIndex:ncol(taxa.df)]
taxa.df2$SampleID <- SampleID

DF <- merge(df.merged, taxa.df2, by = "SampleID")

#----- Taxa ~ Timepoint + Surgery + Observed/Modeled weight mixed linear models -----
results.df.list <- list()
df.index <- 1
for ( t in 2:length(months) ) {
  
  timeFrame <- months[1:t]
  log.df <- DF[ DF$time %in% timeFrame, ]
  startAbundanceIndex <- which(colnames(log.df)=="modeled_weight_change")+1

  ## Parse out metadata from counts
  myT<-log.df[,startAbundanceIndex:ncol(log.df)]
  
  PatientID <- log.df$PatientID
  Timepoint <- log.df$time
  Surgery <- log.df$Surgery
  Observed_Weight <- log.df$Weight_kg
  Modeled_Weight <- log.df$modeled
  Group <- log.df$TertileAssignment_24M
  
  #----- Mixed Linear Modeling -----
  bugName <- vector()
  pVal_Timepoint <- vector()
  pVal_Surgery <- vector()
  pVal_Observed_Weight <- vector()
  pVal_Modeled_Weight <- vector()
  
  index <- 1
  
  for (i in 1:ncol(myT)){
    
    bug<-myT[,i]
    
    if (mean(bug>0, na.rm = TRUE)>0.1){
      
      df <- data.frame(bug, PatientID, Timepoint, Surgery, Observed_Weight, Modeled_Weight, Group)
      df <- na.omit(df)
      df$Timepoint <- as.factor(df$Timepoint)
      
      tryCatch({
        mlm <- lme( bug ~ Timepoint + Surgery +  Observed_Weight, method = "REML", random = ~1 | PatientID, data = df )
        fit<-anova(mlm)
        
        pVal_Timepoint[index] <- fit$"p-value"[2]
        pVal_Surgery[index] <- fit$"p-value"[3]
        pVal_Observed_Weight[index] <- fit$"p-value"[4]
        
        mlm <- lme( bug ~ Timepoint + Surgery +  Modeled_Weight, method = "REML", random = ~1 | PatientID, data = df )
        fit<-anova(mlm)
        
        pVal_Modeled_Weight[index] <- fit$"p-value"[4]
        
        bugName[index]<-colnames(myT)[i]
        
        index<-index+1
        
      }, error = function(e){})
    } # if (mean(bug>0, na.rm = TRUE)>0.1)
    
  } # for (i in 1:ncol(myT))

  dFrame<-data.frame(pVal_Timepoint, pVal_Surgery, pVal_Observed_Weight, pVal_Modeled_Weight)
  dFrame$Adj_pVal_Timepoint <- p.adjust(dFrame$pVal_Timepoint,method = "BH")
  dFrame$Adj_pVal_Surgery <- p.adjust(dFrame$pVal_Surgery,method = "BH")
  dFrame$Adj_pVal_Observed_Weight <- p.adjust(dFrame$pVal_Observed_Weight,method = "BH")
  dFrame$Adj_pVal_Modeled_Weight <- p.adjust(dFrame$pVal_Modeled_Weight,method = "BH")
  
  names(dFrame) <- paste0(names(dFrame), "_BLto", timeFrame[length(timeFrame)],"M")
  dFrame <- data.frame(bugName, dFrame)
  #----- Mixed Linear Model END -----
  
  results.df.list[[df.index]] <- dFrame
  df.index <- df.index + 1
  
} # for ( t in 2:length(months) )

merged.results <- results.df.list[[1]]
if (length(results.df.list) > 1) {
  for (l in 2:length(results.df.list)) {

    merged.results <- merge(merged.results, results.df.list[[l]], by = "bugName")

  } # for (l in 2:length(results.df.list))
}

filename.3 <- paste0(level,  "_by_Timepoint+Surgery+Observed_Tertile_Modeled_Weight_MixedLinearModelResults.tsv")
write.table(merged.results, paste0(outputResults, filename.3),sep="\t",row.names = FALSE,quote = FALSE)
#-----  end  ---------


#' 
#' 
## ----LinearModel-Analysis-p-values, fig.height=5, fig.width=10, warning=FALSE, include=FALSE----
months <- c(0, 1, 6, 12, 18, 24)

outputResults <- paste0(outputDir, "ResultsTables/")
dir.create(outputResults, showWarnings = FALSE)

startAbundanceIndex <- which(colnames(taxa.df)=="ResponderStatus")+1
SampleID <- taxa.df$SampleID
taxa.df2 <- taxa.df[,startAbundanceIndex:ncol(taxa.df)]
taxa.df2$SampleID <- SampleID

DF <- merge(df.merged, taxa.df2, by = "SampleID")

#----- Taxa ~ Timepoint + Surgery + Observed/Modeled weight mixed linear model plots -----
filename.4 <- paste0(level,  "_by_Timepoint+Surgery+Observed_Tertile_Modeled_Weight_MixedLinearModelPlots.pdf")
filePath <- paste0(outputDir, filename.4)
pdf(filePath)

par(mfrow=c(2,2))

groupPValues1 <- vector()
groupPValues2 <- vector()
groupPValues3 <- vector()
groupPValues4 <- vector()
groupNumbers <- vector()
qbugName <- vector()
direction <- vector()
LengthTime <- vector()
qIndex <- 1

for( i in 1:nrow(merged.results)) {
  BugName <- merged.results$bugName[i]

  for ( t in 2:length(months) ) {

    timeFrame <- months[1:t]
    log.df <- DF[ DF$time %in% timeFrame, ]
    bugCol <- log.df[ , which(colnames(log.df) == BugName) ]
    Month <- log.df$time
    Surgery <- log.df$Surgery
    Observed_Weight <- log.df$Weight_kg
    Modeled_Weight <- log.df$modeled
    PatientID <- log.df$PatientID
    Group <- log.df$TertileAssignment_24M

    colEnding <- paste0("BLto", timeFrame[length(timeFrame)],"M")
    Timepoint_p <- merged.results[ i , which(colnames(merged.results) == paste0("Adj_pVal_Timepoint_", colEnding))]
    Surgery_p <- merged.results[ i , which(colnames(merged.results) == paste0("Adj_pVal_Surgery_", colEnding))]
    Observed_Weight_p <- merged.results[ i , which(colnames(merged.results) == paste0("Adj_pVal_Observed_Weight_", colEnding))]
    Modeled_Weight_p <- merged.results[ i , which(colnames(merged.results) == paste0("Adj_pVal_Modeled_Weight_", colEnding))]

    df <- data.frame(PatientID, Surgery, Observed_Weight, bugCol)
    df <- na.omit(df)
    patient_n <- length(unique(df$PatientID))
    n <- nrow(df)

    title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Observed_Weight_p, digits = 3), "; n = ", n)
    plot(  df$Observed_Weight, df$bugCol,
           ylab =paste0(BugName, " log10 Abundance") ,
           xlab="Observed Weight (kg)",
           main=title.lab,
           cex.main = 1 ,
           col = ifelse( df$Surgery == "SG", "red", "blue" )
    )
    legend("topright", legend=c("SG", "RYGB"),
           col=c("red", "blue"), pch = 21, cex=0.8)

    df <- data.frame(PatientID, Surgery, Modeled_Weight, bugCol)
    df <- na.omit(df)
    patient_n <- length(unique(df$PatientID))
    n <- nrow(df)
    title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Modeled_Weight_p, digits = 3), "; n = ", n)
    plot(  df$Modeled_Weight, df$bugCol,
           ylab =paste0(BugName, " log10 Abundance") ,
           xlab="Modeled Weight (kg)",
           main=title.lab,
           cex.main = 1 ,
           col = ifelse( df$Surgery == "SG", "red", "blue" )
    )
    legend("topright", legend=c("SG", "RYGB"),
           col=c("red", "blue"), pch = 21, cex=0.8)
    
    df <- data.frame(PatientID, Surgery, Month, bugCol)
    df <- na.omit(df)
    patient_n <- length(unique(df$PatientID))
    n <- nrow(df)
    
    title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Timepoint_p, digits = 3), "; n = ", n)
    plot(  df$Month, df$bugCol, xaxt = 'n',
           ylab =BugName ,
           xlab="Timepoint (months)",
           main=title.lab,
           cex.main = 1 ,
           col = ifelse( df$Surgery == "SG", "red", "blue" )
    )
    axis(1, at = seq(0, max(df$Month), by = 1))
    legend("top", legend=c("SG", "RYGB"),
           col=c("red", "blue"), pch = 21, cex=0.8)

    df <- data.frame(PatientID, Surgery, bugCol)
    df <- na.omit(df)
    patient_n <- length(unique(df$PatientID))
    n <- nrow(df)
    title.lab <- paste0( BugName ," ", colEnding, " - ", patient_n, " patients\nMLM adj p = ", format(Surgery_p, digits = 3), "; n = ", n)
    boxplot( df$bugCol ~  df$Surgery ,
             xlab="Surgery Type",
             ylab =BugName,
             main = title.lab,
             cex.main = 1
    )
    
    qFrameFinal <- data.frame()
    for( j in 1:ng ) {

      lm.df <- data.frame(PatientID, Month, Surgery, bugCol, Group)
      lm.df <- na.omit(lm.df)
      lm.df <- lm.df[lm.df$Group == j,]

      myLm <- lm(  lm.df$bugCol ~ lm.df$Month) # add surgery type, surgery*time
      p = anova(myLm)$"Pr(>F)"[1]
      groupPValues1[qIndex] <- p
      
      title.lab <- paste0( "Group ", j, " (", length(unique(lm.df$PatientID)), " patients)\n", 
                     BugName , 
                     "\nTimepoint (p = ", format(p, digits = 3), ")")
      plot(  lm.df$Month, lm.df$bugCol ,
             ylab =BugName ,
             xlab="Timepoint (months)",
             main=title.lab,
             cex.main = 1,
             col = ifelse( lm.df$Surgery == "SG", "red", "blue" )
      )
      legend("topright", legend=c("SG", "RYGB"),
             col=c("red", "blue"), pch = 21, cex=0.8)

      # myLm <- lm(  lm.df$bugCol ~ lm.df$Month + lm.df$Surgery + lm.df$Month:lm.df$Surgery) # add surgery type, surgery*time
      myLm <- lm(  lm.df$bugCol ~ lm.df$Surgery + lm.df$Month + lm.df$Surgery:lm.df$Month) # add surgery type, surgery*time
      p = anova(myLm)$"Pr(>F)"[2]
      groupPValues2[qIndex] <- p

      # mlm <- lme( bugCol ~ Month + Surgery + Month:Surgery, method = "REML", random = ~1 | PatientID, data = lm.df )
      mlm <- lme( bugCol ~ Surgery + Month + Surgery:Month, method = "REML", random = ~1 | PatientID, data = lm.df )
      fit<-anova(mlm)
      p = fit$"p-value"[3]
      groupPValues3[qIndex] <- p

      mlm <- lme( bugCol ~ Month, method = "REML", random = ~1 | PatientID, data = lm.df )
      fit<-anova(mlm)
      p = fit$"p-value"[2]
      groupPValues4[qIndex] <- p

      
      groupNumbers[qIndex] <- j
      qbugName[qIndex] <- BugName


      averages <- lm.df %>%
        group_by(Month) %>%
        get_summary_stats(bugCol, type = "mean_sd")

      BL_avg <- averages$mean[1]
      end_avg <- averages$mean[nrow(averages)]

      direction[qIndex] <- ifelse(BL_avg > end_avg, "decrease",
                                  ifelse(BL_avg < end_avg, "increase",
                                         "same"))
      LengthTime[qIndex] <- colEnding

      qIndex <- qIndex + 1

    } # for( j in 1:ng )

    qFrame <- data.frame( qbugName, LengthTime, groupNumbers,groupPValues1, groupPValues2, groupPValues3, groupPValues4, direction)
    qFrame$adjPValues1 <- p.adjust(qFrame$groupPValues1, method="BH")
    qFrame$adjPValues2 <- p.adjust(qFrame$groupPValues2, method="BH")
    qFrame$adjPValues3 <- p.adjust(qFrame$groupPValues3, method="BH")
    qFrame$adjPValues4 <- p.adjust(qFrame$groupPValues4, method="BH")

    qFrameFinal <- rbind(qFrameFinal, qFrame)

  } # for ( t in 2:length(months) )

} # for( i in 1:nrow(merged.results))
dev.off()
#-----  end  ---------

#' 
#' 
## ----print-qFrameFinal, echo=FALSE, warning=FALSE, fig.height=5, fig.width=10----
months <- c(0, 1, 6, 12, 18, 24)

outputResults <- paste0(outputDir, "ResultsTables/")
dir.create(outputResults, showWarnings = FALSE)

startAbundanceIndex <- which(colnames(taxa.df)=="ResponderStatus")+1
SampleID <- taxa.df$SampleID
taxa.df2 <- taxa.df[,startAbundanceIndex:ncol(taxa.df)]
taxa.df2$SampleID <- SampleID

DF <- merge(df.merged, taxa.df2, by = "SampleID")

#----- Taxa ~ Timepoint (by lcga2 group) linear models -----
qFrameFinal <- qFrameFinal[ order( qFrameFinal$groupPValues1),]

filename.5 <- paste0(level, "_by_Timepoint_groupedby_lcga2_LinearModelResults.tsv")
filePath <- paste0(outputResults, filename.5)
write.table(qFrameFinal, file=filePath, sep="\t", row.names=FALSE)
#-----  end  ---------


#' 
#' 
## ----Figure3-code, fig.height=5, fig.width=7, warning=FALSE, include=FALSE----
##### p value box plot #####
viruses <- c("Viruses_noname", "C2likevirus")

filePath <- paste0(outputResults, filename.5)
qFrame <- read.table(filePath, sep="\t",header = TRUE, check.names = FALSE)

qFrame$adjPValues1 <- ifelse(qFrame$adjPValues1 == 0, 2e-16,
                             qFrame$adjPValues1)

qFrame$log <- ifelse(qFrame$direction == "increase", -log10(qFrame$adjPValues1),
                     log10(qFrame$adjPValues1))
qFrame <- na.omit(qFrame)

# Label taxa as significant or not based on adjusted p-value
qFrame$Significance <- ifelse(qFrame$adjPValues1 < 0.05, "Significant",
                              "Not Significant")
qFrame <- qFrame[!( qFrame$qbugName %in% viruses ),]

timeFrames <- c("BLto1M", "BLto6M", "BLto12M", "BLto18M", "BLto24M")
filename.9 <- paste0(level, "_Tertile_analysis1_pvalue_boxplots.pdf")
filePath <- paste0(outputDir, filename.9)
pdf(filePath, width = 7, height = 5)

for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] ) {
  
  # message(timeFrame)
  df <- qFrame[ qFrame$LengthTime == timeFrame, ]
  endPoint <- gsub("BLto", "", timeFrame)
  
  f.test <- fisher.test(df$groupNumbers, df$Significance, workspace = 2e8)
  f.pVal <- paste0("p = ", format(f.test$p.value, digits = 3))
  # f.pVal <- roundP(f.test$p.value)
  
  title.lab <- paste0( str_to_title(level),  " level taxonomic changes over time (BL to ", endPoint, ")"); title.lab
  subtitle.lab <- paste0("Fisher's Exact Test (", f.pVal, ")")
  caption.lab <- "lm(  Taxa ~ Timepoint )"
  
  hLine1 <- -log10(0.05)
  hLine2 <- log10(0.05)
  
  Q3u <- length(which(df$log > hLine1 & df$groupNumbers == 3))
  Q3d <- length(which(df$log < hLine2 & df$groupNumbers == 3))
  Q2u <- length(which(df$log > hLine1 & df$groupNumbers == 2))
  Q2d <- length(which(df$log < hLine2 & df$groupNumbers == 2))
  Q1u <- length(which(df$log > hLine1 & df$groupNumbers == 1))
  Q1d <- length(which(df$log < hLine2 & df$groupNumbers == 1))
  
  labels1 <- c(Q1u, Q2u, Q3u)
  labels2 <- c(Q1d, Q2d, Q3d)
  
  getWeightGroupLocation1 <- function(df, hLine1) {
    
    df <- na.omit(df)
    quants1 <- data.table(df)[, list(quant = as.numeric(max(log))), by = groupNumbers]
    quants1 <- quants1[order(quants1$groupNumbers),]
    for (k in 1:nrow(quants1)) {
      if (quants1$quant[k] < hLine1) {
        quants1$quant[k] <- hLine1
      }
    }
    quants1$quant <- quants1$quant + 0.5
    
    return(quants1)
  }
  getWeightGroupLocation2 <- function(df, hLine2) {
    
    df <- na.omit(df)
    quants2 <- data.table(df)[, list(quant = as.numeric(min(log))), by = groupNumbers]
    quants2 <- quants2[order(quants2$groupNumbers),]
    for (k in 1:nrow(quants2)) {
      if (quants2$quant[k] > hLine2) {
        quants2$quant[k] <- hLine2
      }
    }
    
    quants2$quant <- quants2$quant - 0.5
    
    return(quants2)
  }
  
  quants1 <- getWeightGroupLocation1(df, hLine1)
  quants2 <- getWeightGroupLocation2(df, hLine2)
  
  ## Generate plot ##
  
  plot <- ggboxplot(
    df, x = "groupNumbers", y = "log",
    add = "jitter", shape = 1, 
    palette = c("orange", "blue", "green3", "black", "red"),
    fill = "groupNumbers",
    # color = "groupNumbers",
    # facet.by = "groupNumbers",
    scales = "free",
    # label="qbugName",
    outlier.shape = 1
  )
  
  # plot <- plot +
  #   xlim(0,4)
  
  # Add title & subtitle
  plot <- plot +
    labs(title = title.lab,
         subtitle = subtitle.lab,
         caption = caption.lab)
  
  plot <- plot +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5)
    )
  
  # Edit x- and y-axis labels
  plot <- plot +
    labs(x = "Tertile Groups", y = "p value (log10)")
  
  # Add dashed line for significant p values increasing
  plot <- plot + geom_hline(yintercept=hLine1, linetype="dashed", color = "red")
  
  # Add dashed line for significant p values decreasing
  plot <- plot + geom_hline(yintercept=hLine2, linetype="dashed", color = "red")
  
  # Increase text size of labels
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Remove legend
  plot <- plot + theme(legend.position = "none")
  
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Add text indicating the number of significant taxa increasing over time
  plot <- plot + geom_text(data = quants1, aes(x = groupNumbers, y = quant, label = labels1), size = 5)
  
  # Add text indicating the number of significant taxa decreasing over time
  plot <- plot + geom_text(data = quants2, aes(x = groupNumbers, y = quant, label = labels2), size = 5)
  
  # Add up arrow and annotation
  plot2 <- plot + annotate("segment", 
                           x=0.45, xend=0.45,
                           y=0.25, yend=2,
                           col="black", arrow=arrow(length=unit(0.25, "cm")))
  
  plot2 <- plot2 + annotate("text", x=0.5, y=1.5, hjust = 0, label = "Increases in abundance")
  
  # Add down arrow and annotation
  plot2 <- plot2 + annotate("segment", 
                            x=0.45, xend=0.45,
                            y=-0.25, yend=-2,
                            col="black", arrow=arrow(length=unit(0.25, "cm")))
  plot2 <- plot2 + annotate("text", x=0.5, y=-1.5, hjust = 0, label = "Decreases in abundance")
  
  print(plot2)
  
} # for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] )

dev.off()


# knitr::knit_exit()

## ----Figure3-print, echo=FALSE, fig.height=5, fig.width=7, warning=FALSE----
print(plot2)


#' 
#' **Figure 3.** Boxplots of p-values generated by linear modeling (Taxa ~ Timepoint) for each tertile group (`r filename.9`).
#' 
#' 
## ----Figure4-code, fig.height=5, fig.width=7, warning=FALSE, include=FALSE----
##### p value box plot #####
viruses <- c("Viruses_noname", "C2likevirus")

filePath <- paste0(outputResults, filename.5)
qFrame <- read.table(filePath, sep="\t",header = TRUE, check.names = FALSE)

qFrame$adjPValues2 <- ifelse(qFrame$adjPValues2 == 0, 2e-16,
                             qFrame$adjPValues2)

qFrame$log <- ifelse(qFrame$direction == "increase", -log10(qFrame$adjPValues2),
                     log10(qFrame$adjPValues2))
qFrame <- na.omit(qFrame)

# Label taxa as significant or not based on adjusted p-value
qFrame$Significance <- ifelse(qFrame$adjPValues2 < 0.05, "Significant",
                              "Not Significant")
qFrame <- qFrame[!( qFrame$qbugName %in% viruses ),]

timeFrames <- c("BLto1M", "BLto6M", "BLto12M", "BLto18M", "BLto24M")
filename.10 <- paste0(level, "_Tertile_analysis2_pvalue_boxplots.pdf")
filePath <- paste0(outputDir, filename.10)
pdf(filePath, width = 7, height = 5)

for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] ) {
  
  # message(timeFrame)
  df <- qFrame[ qFrame$LengthTime == timeFrame, ]
  endPoint <- gsub("BLto", "", timeFrame)
  
  f.test <- fisher.test(df$groupNumbers, df$Significance, workspace = 2e8)
  f.pVal <- paste0("p = ", format(f.test$p.value, digits = 3))
  # f.pVal <- roundP(f.test$p.value)
  
  title.lab <- paste0( str_to_title(level),  " level taxonomic changes over time (BL to ", endPoint, ")"); title.lab
  subtitle.lab <- paste0("Fisher's Exact Test (", f.pVal, ")")
  
  hLine1 <- -log10(0.05)
  hLine2 <- log10(0.05)
  
  Q3u <- length(which(df$log > hLine1 & df$groupNumbers == 3))
  Q3d <- length(which(df$log < hLine2 & df$groupNumbers == 3))
  Q2u <- length(which(df$log > hLine1 & df$groupNumbers == 2))
  Q2d <- length(which(df$log < hLine2 & df$groupNumbers == 2))
  Q1u <- length(which(df$log > hLine1 & df$groupNumbers == 1))
  Q1d <- length(which(df$log < hLine2 & df$groupNumbers == 1))
  
  labels1 <- c(Q1u, Q2u, Q3u)
  labels2 <- c(Q1d, Q2d, Q3d)
  
  getWeightGroupLocation1 <- function(df, hLine1) {
    
    df <- na.omit(df)
    quants1 <- data.table(df)[, list(quant = as.numeric(max(log))), by = groupNumbers]
    quants1 <- quants1[order(quants1$groupNumbers),]
    for (k in 1:nrow(quants1)) {
      if (quants1$quant[k] < hLine1) {
        quants1$quant[k] <- hLine1
      }
    }
    quants1$quant <- quants1$quant + 0.5
    
    return(quants1)
  }
  getWeightGroupLocation2 <- function(df, hLine2) {
    
    df <- na.omit(df)
    quants2 <- data.table(df)[, list(quant = as.numeric(min(log))), by = groupNumbers]
    quants2 <- quants2[order(quants2$groupNumbers),]
    for (k in 1:nrow(quants2)) {
      if (quants2$quant[k] > hLine2) {
        quants2$quant[k] <- hLine2
      }
    }
    
    quants2$quant <- quants2$quant - 0.5
    
    return(quants2)
  }
  
  quants1 <- getWeightGroupLocation1(df, hLine1)
  quants2 <- getWeightGroupLocation2(df, hLine2)
  
  ## Generate plot ##
  
  plot <- ggboxplot(
    df, x = "groupNumbers", y = "log",
    add = "jitter", shape = 1, 
    palette = c("orange", "blue", "green3", "black", "red"),
    fill = "groupNumbers",
    # color = "groupNumbers",
    # facet.by = "groupNumbers",
    scales = "free",
    # label="qbugName",
    outlier.shape = 1
  )
  
  # plot <- plot +
  #   xlim(0,4)
  caption.lab <- "lm(  Taxa ~ Surgery + Timepoint + Surgery:Timepoint )"

  # Add title & subtitle
  plot <- plot +
    labs(title = title.lab,
         subtitle = subtitle.lab,
         caption = caption.lab)
  
  plot <- plot +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5)
    )
  
  # Edit x- and y-axis labels
  plot <- plot +
    labs(x = "Tertile Groups", y = "p value (log10)")
  
  # Add dashed line for significant p values increasing
  plot <- plot + geom_hline(yintercept=hLine1, linetype="dashed", color = "red")
  
  # Add dashed line for significant p values decreasing
  plot <- plot + geom_hline(yintercept=hLine2, linetype="dashed", color = "red")
  
  # Increase text size of labels
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Remove legend
  plot <- plot + theme(legend.position = "none")
  
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Add text indicating the number of significant taxa increasing over time
  plot <- plot + geom_text(data = quants1, aes(x = groupNumbers, y = quant, label = labels1), size = 5)
  
  # Add text indicating the number of significant taxa decreasing over time
  plot <- plot + geom_text(data = quants2, aes(x = groupNumbers, y = quant, label = labels2), size = 5)
  
  # Add up arrow and annotation
  plot2 <- plot + annotate("segment", 
                           x=0.45, xend=0.45,
                           y=0.25, yend=2,
                           col="black", arrow=arrow(length=unit(0.25, "cm")))
  
  plot2 <- plot2 + annotate("text", x=0.5, y=1.5, hjust = 0, label = "Increases in abundance")
  
  # Add down arrow and annotation
  plot2 <- plot2 + annotate("segment", 
                            x=0.45, xend=0.45,
                            y=-0.25, yend=-2,
                            col="black", arrow=arrow(length=unit(0.25, "cm")))
  plot2 <- plot2 + annotate("text", x=0.5, y=-1.5, hjust = 0, label = "Decreases in abundance")
  
  print(plot2)
  
} # for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] )

dev.off()


# knitr::knit_exit()

## ----Figure4-print, echo=FALSE, fig.height=5, fig.width=7, warning=FALSE----
print(plot2)


#' 
#' **Figure 4.** Boxplots of p-values generated by linear modeling (Taxa ~ Timepoint + Surgery + Timepoint:Surgery) for each tertile group (`r filename.10`).
#' 
#' 
## ----Figure5-code, fig.height=5, fig.width=7, warning=FALSE, include=FALSE----
##### p value box plot #####
viruses <- c("Viruses_noname", "C2likevirus")

filePath <- paste0(outputResults, filename.5)
qFrame <- read.table(filePath, sep="\t",header = TRUE, check.names = FALSE)

qFrame$adjPValues4 <- ifelse(qFrame$adjPValues4 == 0, 2e-16,
                             qFrame$adjPValues4)

qFrame$log <- ifelse(qFrame$direction == "increase", -log10(qFrame$adjPValues4),
                     log10(qFrame$adjPValues4))
qFrame <- na.omit(qFrame)

# Label taxa as significant or not based on adjusted p-value
qFrame$Significance <- ifelse(qFrame$adjPValues4 < 0.05, "Significant",
                              "Not Significant")
qFrame <- qFrame[!( qFrame$qbugName %in% viruses ),]

timeFrames <- c("BLto1M", "BLto6M", "BLto12M", "BLto18M", "BLto24M")
filename.12 <- paste0(level, "_Tertile_analysis3_pvalue_boxplots.pdf")
filePath <- paste0(outputDir, filename.12)
pdf(filePath, width = 7, height = 5)

for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] ) {
  
  # message(timeFrame)
  df <- qFrame[ qFrame$LengthTime == timeFrame, ]
  endPoint <- gsub("BLto", "", timeFrame)
  
  f.test <- fisher.test(df$groupNumbers, df$Significance, workspace = 2e8)
  f.pVal <- paste0("p = ", format(f.test$p.value, digits = 3))
  # f.pVal <- roundP(f.test$p.value)
  
  title.lab <- paste0( str_to_title(level),  " level taxonomic changes over time (BL to ", endPoint, ")"); title.lab
  subtitle.lab <- paste0("Fisher's Exact Test (", f.pVal, ")")
  
  hLine1 <- -log10(0.05)
  hLine2 <- log10(0.05)
  
  Q3u <- length(which(df$log > hLine1 & df$groupNumbers == 3))
  Q3d <- length(which(df$log < hLine2 & df$groupNumbers == 3))
  Q2u <- length(which(df$log > hLine1 & df$groupNumbers == 2))
  Q2d <- length(which(df$log < hLine2 & df$groupNumbers == 2))
  Q1u <- length(which(df$log > hLine1 & df$groupNumbers == 1))
  Q1d <- length(which(df$log < hLine2 & df$groupNumbers == 1))
  
  labels1 <- c(Q1u, Q2u, Q3u)
  labels2 <- c(Q1d, Q2d, Q3d)
  
  getWeightGroupLocation1 <- function(df, hLine1) {
    
    df <- na.omit(df)
    quants1 <- data.table(df)[, list(quant = as.numeric(max(log))), by = groupNumbers]
    quants1 <- quants1[order(quants1$groupNumbers),]
    for (k in 1:nrow(quants1)) {
      if (quants1$quant[k] < hLine1) {
        quants1$quant[k] <- hLine1
      }
    }
    quants1$quant <- quants1$quant + 0.5
    
    return(quants1)
  }
  getWeightGroupLocation2 <- function(df, hLine2) {
    
    df <- na.omit(df)
    quants2 <- data.table(df)[, list(quant = as.numeric(min(log))), by = groupNumbers]
    quants2 <- quants2[order(quants2$groupNumbers),]
    for (k in 1:nrow(quants2)) {
      if (quants2$quant[k] > hLine2) {
        quants2$quant[k] <- hLine2
      }
    }
    
    quants2$quant <- quants2$quant - 0.5
    
    return(quants2)
  }
  
  quants1 <- getWeightGroupLocation1(df, hLine1)
  quants2 <- getWeightGroupLocation2(df, hLine2)
  
  ## Generate plot ##
  
  plot <- ggboxplot(
    df, x = "groupNumbers", y = "log",
    add = "jitter", shape = 1, 
    palette = c("orange", "blue", "green3", "black", "red"),
    fill = "groupNumbers",
    # color = "groupNumbers",
    # facet.by = "groupNumbers",
    scales = "free",
    # label="qbugName",
    outlier.shape = 1
  )
  
  # plot <- plot +
  #   xlim(0,4)
  
caption.lab <- "lme(  Taxa ~ Timepoint, method = 'REML', random = ~1 | PatientID )"
# Add title & subtitle
plot <- plot +
  labs(title = title.lab,
       subtitle = subtitle.lab,
       caption = caption.lab)

plot <- plot +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5)
  )
  
  # Edit x- and y-axis labels
  plot <- plot +
    labs(x = "Tertile Groups", y = "p value (log10)")
  
  # Add dashed line for significant p values increasing
  plot <- plot + geom_hline(yintercept=hLine1, linetype="dashed", color = "red")
  
  # Add dashed line for significant p values decreasing
  plot <- plot + geom_hline(yintercept=hLine2, linetype="dashed", color = "red")
  
  # Increase text size of labels
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Remove legend
  plot <- plot + theme(legend.position = "none")
  
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Add text indicating the number of significant taxa increasing over time
  plot <- plot + geom_text(data = quants1, aes(x = groupNumbers, y = quant, label = labels1), size = 5)
  
  # Add text indicating the number of significant taxa decreasing over time
  plot <- plot + geom_text(data = quants2, aes(x = groupNumbers, y = quant, label = labels2), size = 5)
  
  # Add up arrow and annotation
  plot2 <- plot + annotate("segment", 
                           x=0.45, xend=0.45,
                           y=0.25, yend=2,
                           col="black", arrow=arrow(length=unit(0.25, "cm")))
  
  plot2 <- plot2 + annotate("text", x=0.5, y=1.5, hjust = 0, label = "Increases in abundance")
  
  # Add down arrow and annotation
  plot2 <- plot2 + annotate("segment", 
                            x=0.45, xend=0.45,
                            y=-0.25, yend=-2,
                            col="black", arrow=arrow(length=unit(0.25, "cm")))
  plot2 <- plot2 + annotate("text", x=0.5, y=-1.5, hjust = 0, label = "Decreases in abundance")
  
  print(plot2)
  
} # for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] )

dev.off()


# knitr::knit_exit()

## ----Figure5-print, echo=FALSE, fig.height=5, fig.width=7, warning=FALSE----
print(plot2)


#' 
#' **Figure 5.** Boxplots of p-values generated by random effects mixed linear modeling (Taxa ~ Timepoint + Surgery) for each tertile group (`r filename.12`).
#' 
## ----Figure6-code, fig.height=5, fig.width=7, warning=FALSE, include=FALSE----
##### p value box plot #####
viruses <- c("Viruses_noname", "C2likevirus")

filePath <- paste0(outputResults, filename.5)
qFrame <- read.table(filePath, sep="\t",header = TRUE, check.names = FALSE)

qFrame$adjPValues3 <- ifelse(qFrame$adjPValues3 == 0, 2e-16,
                             qFrame$adjPValues3)

qFrame$log <- ifelse(qFrame$direction == "increase", -log10(qFrame$adjPValues3),
                     log10(qFrame$adjPValues3))
qFrame <- na.omit(qFrame)

# Label taxa as significant or not based on adjusted p-value
qFrame$Significance <- ifelse(qFrame$adjPValues3 < 0.05, "Significant",
                              "Not Significant")
qFrame <- qFrame[!( qFrame$qbugName %in% viruses ),]

timeFrames <- c("BLto1M", "BLto6M", "BLto12M", "BLto18M", "BLto24M")
filename.11 <- paste0(level, "_Tertile_analysis4_pvalue_boxplots.pdf")
filePath <- paste0(outputDir, filename.11)
pdf(filePath, width = 7, height = 5)

for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] ) {
  
  # message(timeFrame)
  df <- qFrame[ qFrame$LengthTime == timeFrame, ]
  endPoint <- gsub("BLto", "", timeFrame)
  
  f.test <- fisher.test(df$groupNumbers, df$Significance, workspace = 2e8)
  f.pVal <- paste0("p = ", format(f.test$p.value, digits = 3))
  # f.pVal <- roundP(f.test$p.value)
  
  title.lab <- paste0( str_to_title(level),  " level taxonomic changes over time (BL to ", endPoint, ")"); title.lab
  subtitle.lab <- paste0("Fisher's Exact Test (", f.pVal, ")")
  
  hLine1 <- -log10(0.05)
  hLine2 <- log10(0.05)
  
  Q3u <- length(which(df$log > hLine1 & df$groupNumbers == 3))
  Q3d <- length(which(df$log < hLine2 & df$groupNumbers == 3))
  Q2u <- length(which(df$log > hLine1 & df$groupNumbers == 2))
  Q2d <- length(which(df$log < hLine2 & df$groupNumbers == 2))
  Q1u <- length(which(df$log > hLine1 & df$groupNumbers == 1))
  Q1d <- length(which(df$log < hLine2 & df$groupNumbers == 1))
  
  labels1 <- c(Q1u, Q2u, Q3u)
  labels2 <- c(Q1d, Q2d, Q3d)
  
  getWeightGroupLocation1 <- function(df, hLine1) {
    
    df <- na.omit(df)
    quants1 <- data.table(df)[, list(quant = as.numeric(max(log))), by = groupNumbers]
    quants1 <- quants1[order(quants1$groupNumbers),]
    for (k in 1:nrow(quants1)) {
      if (quants1$quant[k] < hLine1) {
        quants1$quant[k] <- hLine1
      }
    }
    quants1$quant <- quants1$quant + 0.5
    
    return(quants1)
  }
  getWeightGroupLocation2 <- function(df, hLine2) {
    
    df <- na.omit(df)
    quants2 <- data.table(df)[, list(quant = as.numeric(min(log))), by = groupNumbers]
    quants2 <- quants2[order(quants2$groupNumbers),]
    for (k in 1:nrow(quants2)) {
      if (quants2$quant[k] > hLine2) {
        quants2$quant[k] <- hLine2
      }
    }
    
    quants2$quant <- quants2$quant - 0.5
    
    return(quants2)
  }
  
  quants1 <- getWeightGroupLocation1(df, hLine1)
  quants2 <- getWeightGroupLocation2(df, hLine2)
  
  ## Generate plot ##
  
  plot <- ggboxplot(
    df, x = "groupNumbers", y = "log",
    add = "jitter", shape = 1, 
    palette = c("orange", "blue", "green3", "black", "red"),
    fill = "groupNumbers",
    # color = "groupNumbers",
    # facet.by = "groupNumbers",
    scales = "free",
    # label="qbugName",
    outlier.shape = 1
  )
  
  # plot <- plot +
  #   xlim(0,4)
  
caption.lab <- "lme(  Taxa ~ Surgery + Timepoint + Surgery:Timepoint, method = 'REML', random = ~1 | PatientID )"
# Add title & subtitle
plot <- plot +
  labs(title = title.lab,
       subtitle = subtitle.lab,
       caption = caption.lab)

plot <- plot +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, size = 8)
  )
  
  # Edit x- and y-axis labels
  plot <- plot +
    labs(x = "Tertile Groups", y = "p value (log10)")
  
  # Add dashed line for significant p values increasing
  plot <- plot + geom_hline(yintercept=hLine1, linetype="dashed", color = "red")
  
  # Add dashed line for significant p values decreasing
  plot <- plot + geom_hline(yintercept=hLine2, linetype="dashed", color = "red")
  
  # Increase text size of labels
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Remove legend
  plot <- plot + theme(legend.position = "none")
  
  plot <- plot +
    theme(text = element_text(size = 15))
  
  # Add text indicating the number of significant taxa increasing over time
  plot <- plot + geom_text(data = quants1, aes(x = groupNumbers, y = quant, label = labels1), size = 5)
  
  # Add text indicating the number of significant taxa decreasing over time
  plot <- plot + geom_text(data = quants2, aes(x = groupNumbers, y = quant, label = labels2), size = 5)
  
  # Add up arrow and annotation
  plot2 <- plot + annotate("segment", 
                           x=0.45, xend=0.45,
                           y=0.25, yend=2,
                           col="black", arrow=arrow(length=unit(0.25, "cm")))
  
  plot2 <- plot2 + annotate("text", x=0.5, y=1.5, hjust = 0, label = "Increases in abundance")
  
  # Add down arrow and annotation
  plot2 <- plot2 + annotate("segment", 
                            x=0.45, xend=0.45,
                            y=-0.25, yend=-2,
                            col="black", arrow=arrow(length=unit(0.25, "cm")))
  plot2 <- plot2 + annotate("text", x=0.5, y=-1.5, hjust = 0, label = "Decreases in abundance")
  
  print(plot2)
  
} # for ( timeFrame in timeFrames[1:length(unique(qFrame$LengthTime))] )

dev.off()


# knitr::knit_exit()

## ----Figure6-print, echo=FALSE, fig.height=5, fig.width=7, warning=FALSE----
print(plot2)


#' 
#' **Figure 6.** Boxplots of p-values generated by random effects mixed linear modeling (Taxa ~ Timepoint + Surgery + Timepoint:Surgery) for each tertile group (`r filename.11`).
#' 
#' 
## ----LinearModel-Analyses,  include=FALSE---------------------
months <- c(0, 1, 6, 12, 18, 24)

outputResults <- paste0(outputDir, "ResultsTables/")
dir.create(outputResults, showWarnings = FALSE)

startAbundanceIndex <- which(colnames(taxa.df)=="ResponderStatus")+1
SampleID <- taxa.df$SampleID
taxa.df2 <- taxa.df[,startAbundanceIndex:ncol(taxa.df)]
taxaNames <- colnames(taxa.df2)
taxa.df2$SampleID <- SampleID

DF <- merge(df.merged, taxa.df2, by = "SampleID")

for (analysis in c("analysis1", "analysis2", "analysis3", "analysis4")) {
  
  groupPValues_Timepoint <- vector()
  groupPValues_Surgery <- vector()
  groupPValues_Interaction <- vector()
  groupNumbers <- vector()
  qbugName <- vector()
  direction <- vector()
  LengthTime <- vector()
  qIndex <- 1

  for( i in 1:length(taxaNames)) {
    BugName <- taxaNames[i]

    for ( t in 2:length(months) ) {

      timeFrame <- months[1:t]
      log.df <- DF[ DF$time %in% timeFrame, ]
      bugCol <- log.df[ , which(colnames(log.df) == BugName) ]
      Month <- log.df$time
      Surgery <- log.df$Surgery
      PatientID <- log.df$PatientID
      Group <- log.df$TertileAssignment_24M
      colEnding <- paste0("BLto", timeFrame[length(timeFrame)],"M")

      qFrameFinal <- data.frame()
      for( j in 1:ng ) {

        lm.df <- data.frame(PatientID, Month, Surgery, bugCol, Group)
        lm.df <- na.omit(lm.df)
        lm.df <- lm.df[lm.df$Group == j,]
        
        tryCatch({
          # Analysis #1
          if (analysis == "analysis1") {
            myLm <- lm(  lm.df$bugCol ~ lm.df$Month) # add surgery type, surgery*time
            p = anova(myLm)$"Pr(>F)"[1]
            groupPValues_Timepoint[qIndex] <- p
          }
          
          # Analysis #2
          if (analysis == "analysis2") {
            myLm <- lm(  lm.df$bugCol ~ lm.df$Surgery + lm.df$Month + lm.df$Surgery:lm.df$Month) # add surgery type, surgery*time
            p = anova(myLm)$"Pr(>F)"[2]
            groupPValues_Timepoint[qIndex] <- p
            p = anova(myLm)$"Pr(>F)"[1]
            groupPValues_Surgery[qIndex] <- p
            p = anova(myLm)$"Pr(>F)"[3]
            groupPValues_Interaction[qIndex] <- p
          }
          
          # Analysis #3
          if (analysis == "analysis3") {
            mlm <- lme( bugCol ~ Month, method = "REML", random = ~1 | PatientID, data = lm.df )
            fit<-anova(mlm)
            p = fit$"p-value"[2]
            groupPValues_Timepoint[qIndex] <- p
          }
          
          # Analysis #4
          if (analysis == "analysis4") {
            mlm <- lme( bugCol ~ Surgery + Month + Surgery:Month, method = "REML", random = ~1 | PatientID, data = lm.df )
            fit<-anova(mlm)
            p = fit$"p-value"[3]
            groupPValues_Timepoint[qIndex] <- p
            p = fit$"p-value"[2]
            groupPValues_Surgery[qIndex] <- p
            p = fit$"p-value"[4]
            groupPValues_Interaction[qIndex] <- p
          }
          
          groupNumbers[qIndex] <- j
          qbugName[qIndex] <- BugName
          
          averages <- lm.df %>%
            group_by(Month) %>%
            get_summary_stats(bugCol, type = "mean_sd")
          
          BL_avg <- averages$mean[1]
          end_avg <- averages$mean[nrow(averages)]
          
          direction[qIndex] <- ifelse(BL_avg > end_avg, "decrease",
                                      ifelse(BL_avg < end_avg, "increase",
                                             "same"))
          LengthTime[qIndex] <- colEnding
          
          qIndex <- qIndex + 1
        },
        error = function(e){})

      } # for( j in 1:ng )

      qFrame <- data.frame( qbugName, LengthTime, groupNumbers, direction)
      qFrame <- data.frame( qFrame, groupPValues_Timepoint)
      qFrame$adjgroupPValues_Timepoint <- p.adjust(qFrame$groupPValues_Timepoint, method="BH")

      # Analysis #2 & #4
      if (analysis %in% c("analysis2", "analysis4")) {
        qFrame <- data.frame( qFrame, groupPValues_Surgery, groupPValues_Interaction)
        qFrame$adjgroupPValues_Surgery <- p.adjust(qFrame$groupPValues_Surgery, method="BH")
        qFrame$adjgroupPValues_Interaction <- p.adjust(qFrame$groupPValues_Interaction, method="BH")
      }

      qFrameFinal <- rbind(qFrameFinal, qFrame)

    } # for ( t in 2:length(months) )

  } # for( i in 1:length(taxaNames))

  qFrameFinal <- qFrameFinal[ order( qFrameFinal$groupPValues_Timepoint),]

  filename.5 <- paste0(level, "_Tertile_", analysis, ".tsv")

  filePath <- paste0(outputResults, filename.5)
  write.table(qFrameFinal, file=filePath, sep="\t", row.names=FALSE)
  
} # for (analysis in c(analysis1, analysis2, analysis3, analysis4))


#-----  end  ---------



#' 
#' 
## ----p-value-plots, fig.height=5, fig.width=10, warning=FALSE, include=FALSE----
months <- c(0, 1, 6, 12, 18, 24)

outputResults <- paste0(outputDir, "ResultsTables/")
dir.create(outputResults, showWarnings = FALSE)

startAbundanceIndex <- which(colnames(taxa.df)=="ResponderStatus")+1
SampleID <- taxa.df$SampleID
taxa.df2 <- taxa.df[,startAbundanceIndex:ncol(taxa.df)]
taxa.df2$SampleID <- SampleID

DF <- merge(df.merged, taxa.df2, by = "SampleID")


#---- Taxa ~ Timepoint (by lcga2 group) linear model plots ----
filename.6 <- paste0(level, "_by_Timepoint_groupedby_Tertile_LinearModelPlots.pdf")
filePath <- paste0(outputDir, filename.6)
pdf(filePath)

for ( t in 2:length(months) ) {

  timeFrame <- months[1:t]
  colEnding <- paste0("BLto", timeFrame[length(timeFrame)],"M")

  pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Timepoint_", colEnding))]
  n=length(pValues)
  title.lab <- paste0(colEnding, "\nMLM Timepoint p values (n = ", n, ")")
  hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))

  pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Surgery_", colEnding))]
  n=length(pValues)
  title.lab <- paste0(colEnding, "\nMLM Surgery p values (n = ", n, ")")
  hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))

  pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Observed_Weight_", colEnding))]
  n=length(pValues)
  title.lab <- paste0(colEnding, "\nMLM Observed_Weight p values (n = ", n, ")")
  hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))

  pValues <- merged.results[,which(colnames(merged.results) == paste0("pVal_Modeled_Weight_", colEnding))]
  n=length(pValues)
  title.lab <- paste0(colEnding, "\nMLM Modeled_Weight p values (n = ", n, ")")
  hist(pValues, xlab = "p value bins", main = title.lab, breaks = seq(0, 1, 0.05))
  

} # for ( t in 2:length(months) )

dev.off()
#-----  end  ---------


#' 
## ----Taxa-Group-Boxplots, fig.height=5, fig.width=10, warning=FALSE, include=FALSE----
months <- c(0, 1, 6, 12, 18, 24)

outputResults <- paste0(outputDir, "ResultsTables/")
dir.create(outputResults, showWarnings = FALSE)

startAbundanceIndex <- which(colnames(taxa.df)=="ResponderStatus")+1
SampleID <- taxa.df$SampleID
taxa.df2 <- taxa.df[,startAbundanceIndex:ncol(taxa.df)]
taxa.df2$SampleID <- SampleID

DF <- merge(df.merged, taxa.df2, by = "SampleID")

#------ Taxa ~ lcga2 group (by Timepoint) linear models -----
dFrame.Final <- data.frame()
for ( t in 1:length(months) ) {

  timeFrame <- months[t]
  log.df <- DF[ DF$time %in% timeFrame, ]
  Surgery <- log.df$Surgery
  Observed_Weight <- log.df$Weight_kg
  Modeled_Weight <- log.df$modeled
  PatientID <- log.df$PatientID
  Group <- log.df$TertileAssignment_24M

  myT<-log.df[,startAbundanceIndex:ncol(log.df)]

  pVal_Group <- vector()
  bugName <- vector()
  index <- 1

  for (i in 1:ncol(myT)){

    bug<-myT[,i]

    if (mean(bug>0, na.rm = TRUE)>0.1){

      df <- data.frame(bug, Group)
      df <- na.omit(df)

      lm <- lm( bug ~ Group, data = df )
      fit<-anova(lm)

      pVal_Group[index] <- fit$`Pr(>F)`[1]

      bugName[index]<-colnames(myT)[i]
      Timepoint[index] <- t
      index<-index+1

    } # if (mean(bug>0, na.rm = TRUE)>0.1)

  } # for (i in 1:ncol(myT))

  dFrame<-data.frame(bugName, pVal_Group)
  dFrame$Adj_pVal_Group <- p.adjust(dFrame$pVal_Group,method = "BH")
  dFrame$Timepoint <- months[t]

  dFrame.Final <- rbind(dFrame.Final, dFrame)
} # for ( t in 1:length(months) )

filename.7 <- paste0(level, "by_lcga2_groupedby_Timepoint_LinearModelResults.tsv")
filePath <- paste0(outputResults, filename.7)
write.table(dFrame.Final, filePath, sep="\t",quote = FALSE, row.names = FALSE)
#-----  end  ---------


#----- Taxa ~ lcga2 group (by Timepoint) linear model plots -----
tukey.df <- data.frame()
plotList <- list()
index <- 1
for (i in unique(dFrame.Final$bugName)) {

  dFrame2 <- dFrame.Final[dFrame.Final$bugName == i,]
  # DF2 <- DF[ !is.na( DF[which(colnames(DF) == divMonth)]), ]
  # DF2 <- DF2[ !is.na( DF2[which(colnames(DF2) == i)]), ]
  DF2 <- DF[DF$time %in% months,]

  Group <- DF2$TertileAssignment_24M
  Time <- DF2$time
  bug <- DF2[, which(colnames(DF2) == i)]
  df2 <- data.frame(Group, Time, bug)
  df2 <- na.omit(df2)
  x.lab <- paste0(ng, " Class LCGA Groups")
  y.lab <- i

  colors <- c("orange", "blue", "green3", "black", "red")
  colorPalette <- colors[1:ng]

  plot <- ggboxplot(
    df2, x = "Group", y = "bug", color = "black",
    fill = "Group", palette = colorPalette,
    # facet.by = "month",
    scales = "free", add = "jitter"
  ); plot

  month.labsP <- paste0(month.labs, "\n", roundP(dFrame2$Adj_pVal_Group[1:nrow(dFrame2)]))
  plot <- facet(plot, facet.by = "Time"
                , panel.labs = list(Time = month.labsP[1:length(months)])
  ); plot

  plot <- plot + labs(x=x.lab, y = y.lab)+ theme(legend.position = "none"); plot

  plotList[[index]] <- plot
  index <- index + 1
}

filename.8 <- paste0(level, "_by_lcga2_groupedby_Timepoint_LinearModel_BoxPlots.pdf")
filePath <- paste0(outputDir, filename.8)
pdf(filePath, width = 10, height = 7)
for (x in 1:length(plotList)) {
  print(plotList[[x]])
}
dev.off()
#-----  end  ---------


#' 
## ----KNIT_EXIT, include=FALSE---------------------------------
knitr::knit_exit()


#' 
