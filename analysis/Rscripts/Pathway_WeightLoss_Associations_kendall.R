#Author: Alicia Sorgen
#Date: 2022 August 16
#Description: Kendall correlation modeling weight loss by microbiome

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
params <- c(params, "HumanN2")
params <- c(params, 0)
params <- c(params, 1)
params <- c(params, 6)
params <- c(params, 12)
params <- c(params, 18)
params <- c(params, 24)
end_month <- params[length(params)]


moduleRoot <- "Pathway_WeightLoss_Associations_kendall"

##### Libraries #####
R <- sessionInfo()
message(R$R.version$version.string)

library(stringr); message("stringr: Version ", packageVersion("stringr"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))

##### Set up working environment #####
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
  
  if (args[2] == "MetaPhlAn2" | args[2] == "Kraken2" | args[2] == "HumanN2") {
    module <- paste0(args[2], "_", module)
  } 
  
  if (exists("end_month") == TRUE) {
    module <- paste0(module, "_BLto", end_month, "M")
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
  # file.remove(files)
  
  resourcesDir <- paste0(moduleDir, "resources/")
  if (any(dir(moduleDir) == "resources") == FALSE) {
    dir.create(resourcesDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/functions.R")
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
  }
  setwd(paste0(moduleDir, "script/"))
  # rm(files, gitInput, gitRoot, gitScripts, module, moduleDir, outputDir, params, pipeline, resourcesDir, root, rootInput, script, today)
  
}

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
classifier <- args[2]
included <- args[3:length(args)]

prevModule <- paste0(classifier, "_PathMetaMerge")
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
CountFile <- "PathAbundance_CPM_HumanN2.tsv"

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Prep for analysis #####
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

factorLevels <- c("BL-1M", "BL-6M", "BL-12M", "BL-18M", "BL-24M")
surgeryPalette <- c("steelblue", "gold")

file.path <- paste0(inputDir, CountFile)
logCounts <- read.delim(file.path,sep="\t",header = TRUE,check.names = FALSE)

##### Correlation #####


startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1

# Parse out metadata from counts
myT<-logCounts[,startAbundanceIndex:ncol(logCounts)]

PatientID <- logCounts$PatientID
Timepoint_kg <- paste0("BL-", logCounts$time)
Percent_Loss_kg <- logCounts$Percent_Loss_kg

Path <- vector()
Path_Timepoint <- vector()
Weight_Loss_Timepoint <- vector()
p_Kendall <- vector()
R_Kendall <- vector()
p_Spearman <- vector()
R_Spearman <- vector()
n <- vector()
plotList <- list()
dFrame <- data.frame()
index <- 1

for (i in 1:ncol(myT)) {
  
  pw <- myT[,i]
  pwName <- colnames(myT)[i]
  
  df <- data.frame( pw, PatientID, Timepoint_kg, Percent_Loss_kg )
  
  for (t in 1:length(included)) {
    
    pathMonth <- included[t]
    
    # id rows
    xM <- which(logCounts$time == pathMonth)
    
    # abundance
    pw_xM <- df[xM, "pw"]
    names(pw_xM) <- df[xM, "PatientID"]
    
    # Add pw info to each row
    df$pw_xM <- pw_xM[df$PatientID]
    
    tIndex <- 1
    for (w in 2:length(included)) {
      
      PWL_BL_xM <- paste0("BL-", included[w])
      
      df2 <- df[ df$Timepoint_kg == PWL_BL_xM, ]
      df2 <- na.omit(df2)
      
      if ( mean(df2$pw_xM > 0, na.rm = TRUE) > 0.1 ) {
        
        myKen <- cor.test( df2$pw_xM, df2$Percent_Loss_kg, method = "kendall" )
        mySpear <- cor.test( df2$pw_xM, df2$Percent_Loss_kg, method = "spearman" )
        
        Path[index] <- pwName
        Path_Timepoint[index] <- pathMonth
        Weight_Loss_Timepoint[index] <- paste0(PWL_BL_xM, "M")
        p_Kendall[index] <- myKen$p.value
        R_Kendall[index] <- myKen$estimate
        p_Spearman[index] <- mySpear$p.value
        R_Spearman[index] <- mySpear$estimate
        n[index] <- nrow(df2)
        
        index <- index + 1
        tIndex <- tIndex + 1
        
      } # if ( mean(df2$pw > 0, na.rm = TRUE) > 0.1 )
      
    } # for (w in 2:length(included))
    
  }  # for (m in included)
  
  
} # for (i in 1:ncol(myT))

dFrame <- data.frame( Path, Path_Timepoint, Weight_Loss_Timepoint, n, p_Kendall, R_Kendall, p_Spearman, R_Spearman)

final <- data.frame()
for (x in unique(dFrame$Path_Timepoint)) {
  
  dFrame.Filt <- dFrame[dFrame$Path_Timepoint == x,]
  
  for (y in unique(dFrame$Weight_Loss_Timepoint)) {
    
    dFrame.Filt2 <- dFrame.Filt[dFrame.Filt$Weight_Loss_Timepoint == y,]
    dFrame.Filt2$pAdj_Kendall <- p.adjust(dFrame.Filt2$p_Kendall, method = "BH")
    dFrame.Filt2$pAdj_Spearman <- p.adjust(dFrame.Filt2$p_Spearman, method = "BH")
    final <- rbind(final, dFrame.Filt2)
    
  }
  
}

final <- final[ order(final$p_Kendall), ]

file.path <- paste0(outputDir, "path_by_weight_loss_ken_spear_results_BLto", included[length(included)], "M.tsv")
write.table(final, file.path, sep="\t",quote = FALSE, row.names = FALSE)





##### Plot generation #####
stats.df<-final

startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1

# Parse out metadata from counts
myT<-logCounts[,startAbundanceIndex:ncol(logCounts)]

PatientID <- logCounts$PatientID
Timepoint_kg <- paste0("BL-", logCounts$time, "M")
Percent_Loss_kg <- logCounts$Percent_Loss_kg
SurgeryType <- logCounts$Surgery

plotList <- list()
sig_plotList <- list()
dFrame <- data.frame()
index <- 1
sig_index <- 1

for (i in 1:ncol(myT)) {
  
  pw <- myT[,i]
  pwName <- colnames(myT)[i]
  
  df <- data.frame( pw, PatientID, Timepoint_kg, Percent_Loss_kg, SurgeryType )
  
  for (t in 1:length(included)) {
    
    pathMonth <- included[t]
    # message("pathMonth = ", pathMonth)
    
    # id rows
    xM <- which(logCounts$time == included[t])
    
    # abundance
    pw_xM <- df[xM, "pw"]
    names(pw_xM) <- df[xM, "PatientID"]
    
    # Add pw info to each row
    df$pw_xM <- pw_xM[df$PatientID]
    
    # if (t == 1) {
    for (w in 2:length(included)) {
      
      PWL_BL_xM <- paste0("BL-", included[w], "M")
      
      df2 <- df[ df$Timepoint_kg == PWL_BL_xM, ]
      df2 <- na.omit(df2)
      
      if ( mean(df2$pw_xM > 0, na.rm = TRUE) > 0.1 ) {
        
        stats.df_1 <- stats.df[stats.df$Path == pwName,]
        stats.df_2 <- stats.df_1[stats.df_1$Path_Timepoint == pathMonth,]
        stats.df_3 <- stats.df_2[stats.df_2$Weight_Loss_Timepoint == PWL_BL_xM,]
        pAdj <- roundP(stats.df_3$pAdj_Kendall[1])
        p <- roundP(stats.df_3$p_Kendall[1])
        title.lab <- paste0("Kendall ", p, ", adj ", pAdj)
        subtitle.lab <- paste0("Spearman R = ", round(stats.df_3$R_Spearman[1], 3), "\n", pwName)
        color <- ifelse(stats.df_3$pAdj[1] < 0.05, "red", "black")
        
        plot <- ggscatter(df2, x = "pw_xM", y = "Percent_Loss_kg", color = "SurgeryType", palette = surgeryPalette,
                          shape = 19, size = 2.5 # Points shape and size
                          # add = "reg.line",  # Add regression line
                          # conf.int = TRUE, # Add confidence interval
                          # add.params = list(fill = "lightgray")
        ); plot
        
        x.lab <- paste0("Pathway abundance at ", month.labs[t])
        y.lab <- paste0("Weight loss (%): Baseline - ", month.labs[w])
        plot <- plot + labs(x=x.lab, y = y.lab, title = title.lab, subtitle = subtitle.lab); plot
        # plot <- plot + theme(plot.title = element_text(colour = color)); plot
        
        plotList[[index]] <- plot
        
        if (stats.df_3$pAdj_Kendall[1] < 0.05) {
          sig_plotList[[sig_index]] <- plot
          sig_index <- sig_index + 1
        }
        index <- index + 1
        
      } # if ( mean(df2$pw > 0, na.rm = TRUE) > 0.1 )
    } # for (w in 2:length(included))
  }  # for (t in 1:length(included))
} # for (i in 1:ncol(myT))

# file.path <- paste0(outputDir, "All_path_by_weight_loss_kendall_BLto", included[length(included)], "M.pdf")
# pdf(file.path, width = 10, height = 10)
# mod <- length(plotList) %% 4
# if ( mod == 0 ) {
#   
#   for (i in seq(1, length(plotList), 4)) {
#     grid.arrange(plotList[[i]], plotList[[i+1]], 
#                  plotList[[i+2]],plotList[[i+3]],
#                  ncol = 2, nrow = 2)
#   }
#   
# } else {
#   
#   for (i in seq(1, (length(plotList)-mod), 4)) {
#     grid.arrange(plotList[[i]], plotList[[i+1]], 
#                  plotList[[i+2]],plotList[[i+3]],
#                  ncol = 2, nrow = 2)
#     j <- i
#   }
#   j <- j + 1
#   if ( mod == 3) {
#     grid.arrange(plotList[[j]], plotList[[j+1]], 
#                  plotList[[j+2]],
#                  ncol = 2, nrow = 2)
#   }
#   
#   if ( mod == 2) {
#     grid.arrange(plotList[[j]], plotList[[j+1]], 
#                  ncol = 2, nrow = 2)
#   }
#   
#   if ( mod == 1) {
#     grid.arrange(plotList[[j]],
#                  ncol = 2, nrow = 2)
#   }
#   
# }
# dev.off()


file.path <- paste0(outputDir, "Significant_path_by_weight_loss_kendall_BLto", included[length(included)], "M.pdf")
pdf(file.path, width = 10, height = 10)
mod <- length(sig_plotList) %% 4
if ( mod == 0 ) {
  
  for (i in seq(1, length(sig_plotList), 4)) {
    grid.arrange(sig_plotList[[i]], sig_plotList[[i+1]], 
                 sig_plotList[[i+2]],sig_plotList[[i+3]],
                 ncol = 2, nrow = 2)
  }
  
} else {
  
  for (i in seq(1, (length(sig_plotList)-mod), 4)) {
    grid.arrange(sig_plotList[[i]], sig_plotList[[i+1]], 
                 sig_plotList[[i+2]],sig_plotList[[i+3]],
                 ncol = 2, nrow = 2)
    j <- i
  }
  j <- j + 1
  if ( mod == 3) {
    grid.arrange(sig_plotList[[j]], sig_plotList[[j+1]], 
                 sig_plotList[[j+2]],
                 ncol = 2, nrow = 2)
  }
  
  if ( mod == 2) {
    grid.arrange(sig_plotList[[j]], sig_plotList[[j+1]], 
                 ncol = 2, nrow = 2)
  }
  
  if ( mod == 1) {
    grid.arrange(sig_plotList[[j]],
                 ncol = 2, nrow = 2)
  }
  
}
dev.off()


##### Summary tables (unadjusted p values) #####

PathTime <- vector()
WeightLossTime <- vector()
nSignif <- vector()
index <- 1

for (t in 1:length(included)) {
  
  stats.df_1 <- stats.df[stats.df$Path_Timepoint == included[t], ]
  
  for (w in 2:length(included)) {
    
    stats.df_2 <- stats.df_1[stats.df_1$Weight_Loss_Timepoint == paste0("BL-", included[w]), ]
    nSignif[index] <- ifelse( nrow(stats.df_2) > 0, length(which(stats.df_2$p_Kendall < 0.05)), NA )
    PathTime[index] <- month.labs[t]
    WeightLossTime[index] <- paste0("BL-", included[w], "M")
    index <- index + 1
    
  }
}
dFrame <- data.frame( PathTime, WeightLossTime, nSignif )
dFrame$PathTime <- as.factor(dFrame$PathTime)
dFrame$PathTime <- factor(dFrame$PathTime, levels = month.labs)
dFrame$WeightLossTime <- as.factor(dFrame$WeightLossTime)
dFrame$WeightLossTime <- factor(dFrame$WeightLossTime, levels = factorLevels[1:(length(included)-1)] )

x.lab <- "Weight Loss"
y.lab <- str_to_title(level)
# Heatmap 
plot <- ggplot(dFrame, aes(WeightLossTime, PathTime, fill= nSignif)) + 
  geom_tile()+
  # scale_fill_distiller(palette = "RdPu"); plot
  scale_fill_gradient(low="white", high="blue"); plot

plot <- plot + labs(x=x.lab, y = y.lab); plot

plot <- plot + labs(fill = "Significance\n(p values)"); plot

file.path <- paste0(outputDir, "path_by_weight_loss_ken_spear_results_BLto", included[length(included)], "M_unadjP_heatmap.pdf")
pdf(file.path, width = 7, height = 7)
print(plot)
dev.off()

dFrame2 <- spread( dFrame, PathTime, nSignif )

file.path <- paste0(outputDir, "path_by_weight_loss_ken_spear_results_BLto", included[length(included)], "M_summary_unadjP.tsv")
write.table(dFrame2, file.path, sep="\t",quote = FALSE, row.names = FALSE)


##### Summary tables (adjusted p values) #####

PathTime <- vector()
WeightLossTime <- vector()
nSignif <- vector()
index <- 1

for (t in 1:length(included)) {
  
  stats.df_1 <- stats.df[stats.df$Path_Timepoint == included[t], ]
  
  for (w in 2:length(included)) {
    
    stats.df_2 <- stats.df_1[stats.df_1$Weight_Loss_Timepoint == paste0("BL-", included[w]), ]
    nSignif[index] <- ifelse( nrow(stats.df_2) > 0, length(which(stats.df_2$pAdj_Kendall < 0.05)), NA )
    PathTime[index] <- month.labs[t]
    WeightLossTime[index] <- paste0("BL-", included[w], "M")
    index <- index + 1
    
  }
}
dFrame <- data.frame( PathTime, WeightLossTime, nSignif )
dFrame$PathTime <- as.factor(dFrame$PathTime)
dFrame$PathTime <- factor(dFrame$PathTime, levels = month.labs)
dFrame$WeightLossTime <- as.factor(dFrame$WeightLossTime)
dFrame$WeightLossTime <- factor(dFrame$WeightLossTime, levels = factorLevels[1:(length(included)-1)] )

x.lab <- "Weight Loss"
y.lab <- str_to_title(level)
# Heatmap 
plot <- ggplot(dFrame, aes(WeightLossTime, PathTime, fill= nSignif)) + 
  geom_tile()+
  # scale_fill_distiller(palette = "RdPu"); plot
  scale_fill_gradient(low="white", high="blue"); plot

plot <- plot + labs(x=x.lab, y = y.lab); plot

plot <- plot + labs(fill = "Significance\n(p values)"); plot

file.path <- paste0(outputDir, "path_by_weight_loss_ken_spear_results_BLto", included[length(included)], "M_adjP_heatmap.pdf")
pdf(file.path, width = 7, height = 7)
print(plot)
dev.off()

dFrame2 <- spread( dFrame, PathTime, nSignif )

file.path <- paste0(outputDir, "path_by_weight_loss_ken_spear_results_BLto", included[length(included)], "M_summary_adjP.tsv")
write.table(dFrame2, file.path, sep="\t",quote = FALSE, row.names = FALSE)
