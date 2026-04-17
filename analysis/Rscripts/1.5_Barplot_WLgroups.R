#Author: Alicia Sorgen
#Date: 2022 October 21
#Description: Generate bar plots


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
library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
library(scales); message("scales: Version ", packageVersion("scales"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))


# Script-specific edits -------------------------------------------------
params <- vector()
params <- c(params, getOption("mbs.pipe_root", default = stop("Set options(mbs.pipe_root = ...) in your local .Rprofile (gitignored)")))
params <- c(params, "NA")
params <- c(params, "NA")

module <- "1.5_Barplot_WLgroups"


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
prevModule <- "WeightMetaMerge"
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
inputFile <- "metadata.tsv"

month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")
divisions <- c("Quintile", "Quartile", "Tertile", "Half")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Prep data table for analysis #####
statusColors <- c("tomato", "steelblue")
surgeryColors <- c("#3171BC", "#E8C241")
tags <- c("a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i.", "j.", "k.", "l.")

# Read in table
file.path <- paste0(inputDir, inputFile)
myTable <- read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

file.path <- paste0(inputDir, "metadata_OLD.tsv")
divisions <- c("Quintile", "Quartile", "Tertile", "Half")

if (file.exists(file.path) == FALSE) {
  
  write.table(myTable, file.path, sep="\t",row.names = FALSE,quote = FALSE)
  
  
  # Order Timepoint as factor
  myTable$Timepoint <- as.factor(myTable$Timepoint)
  myTable$Timepoint <- factor(myTable$Timepoint, levels = c("BL", "ONE", "SIX", "TWELVE", "EIGHTEEN", "TWENTY_FOUR"))
  
  # Assign SampleID to as dataframe row names
  row.names(myTable) = myTable$SampleID
  
  # Convert the timepoint data to numeric values.
  myTable$time = as.numeric(myTable$time)
  
  months <- c("12", "18", "24")
  included <- c("0", "1", "6")
  
  ##### Weight loss group assignment #####
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
  
  
  file.path <- paste0(inputDir, inputFile)
  write.table(myTable, file.path, sep="\t",row.names = FALSE,quote = FALSE)
  
  
}


##### Bar plots PEWL between responder status groups at each timepoint #####
message("\n************* Bar plots PEWL between responder status groups at each timepoint *************\n")
months <- c("12", "18", "24")
included <- c("0", "1", "6")
index <- 1
plotList <- list()

for (month in months) {
  
  included <- c(included, month)
  
  PatientID <- myTable$PatientID
  Timepoint <- myTable$time
  PEWL <- myTable$PEWL_kg
  
  # Create smaller data table with relevant columns
  df <- data.frame(PatientID, Timepoint, PEWL)
  df <- na.omit(df)
  
  # Pick out the patients designated at a particular responder status
  rIDs <- myTable$PatientID[which(myTable$ResponderStatus == paste0(month, "-month responder"))]
  nrIDs <- myTable$PatientID[which(myTable$ResponderStatus == paste0(month, "-month non-responder"))]
  
  # Filter table to only contain designated IDs
  df <- df[df$PatientID %in% c(rIDs, nrIDs),]
  
  # Filter table to only contain needed timepoints
  df <- df[df$Timepoint %in% included[2:length(included)],]
  
  df$Status <- ifelse(df$PatientID %in% rIDs, "Responder",
                      "Non-responder")
  stat.test <- df %>%
    group_by(Timepoint) %>%
    wilcox_test(PEWL ~ Status) %>%
    # adjust_pvalue(method = "BH") %>%
    add_significance()
   
  stat.test <- stat.test %>%
    add_xy_position(x = "Timepoint", fun = "mean_se")
  
  plot <- ggbarplot(df, "Timepoint", "PEWL",
                    fill = "Status",
                    color = "black", 
                    palette = statusColors,
                    add = "mean_se",
                    label = FALSE, position = position_dodge())
  
  plot <- plot +
    labs(x = "Time (months)", y = "Excess weight loss (%)")
  
  plot <- plot +
    stat_pvalue_manual(
      stat.test,
      bracket.nudge.y = 0.5,
      # bracket.shorten = 1,
      bracket.size = 0.5,
      tip.length = 0.01,
      # remove.bracket = TRUE,
      size = 6,
      hide.ns = TRUE,
      # label = "p.adj.signif"
      label = "p.signif"
    )
  plot
  plotList[[index]] <- plot
  
  index <- index + 1
  
} # for (month in months)

plotFileName <- "response_PEWL_barplots.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 15, height = 5)
for (i in c(1)) {
  grid.arrange(plotList[[i]], plotList[[i+1]],
               plotList[[i+2]], 
               ncol = 3, nrow = 1)
}
dev.off()

##### Bar plots PEWL between surgery types at each timepoint #####
message("\n************* Bar plots PEWL between surgery types at each timepoint *************\n")
index <- 1
plotList <- list()

PatientID <- myTable$PatientID
Timepoint <- myTable$time
PEWL <- myTable$PEWL_kg
Surgery <- myTable$Surgery

# Create smaller data table with relevant columns
df <- data.frame(PatientID, Timepoint, PEWL, Surgery)
df <- na.omit(df)
df <- df[!(df$Timepoint == 0),]

stat.test <- df %>%
  group_by(Timepoint) %>%
  wilcox_test(PEWL ~ Surgery) %>%
  # adjust_pvalue(method = "BH") %>%
  add_significance()

stat.test <- stat.test %>%
  add_xy_position(x = "Timepoint", fun = "mean_se")

plot <- ggbarplot(df, "Timepoint", "PEWL",
                  fill = "Surgery",
                  color = "black", 
                  palette = surgeryColors,
                  add = "mean_se",
                  label = FALSE, position = position_dodge())


plot <- plot +
  labs(x = "Time (months)", y = "Excess weight loss (%)")

plot <- plot +
  stat_pvalue_manual(
    stat.test,
    bracket.nudge.y = 0.5,
    # bracket.shorten = 1,
    bracket.size = 0.5,
    tip.length = 0.01,
    # remove.bracket = TRUE,
    size = 6,
    hide.ns = TRUE,
    # label = "p.adj.signif"
    label = "p.signif"
  )
plot
plotList[[index]] <- plot

plotFileName <- "surgery_PEWL_barplots.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 5, height = 5)
for (i in c(1)) {
  grid.arrange(plotList[[i]],
               ncol = 1, nrow = 1)
}
dev.off()

##### Bar plots Percent_Loss (BL-12M) between surgery types at each timepoint #####
index <- 1
plotList <- list()

PatientID <- myTable$PatientID
Timepoint <- myTable$time
Percent_Loss <- myTable$Percent_Loss_kg
Surgery <- myTable$Surgery
Weight <- myTable$Weight_kg


# Create smaller data table with relevant columns
df <- data.frame(PatientID, Timepoint, Percent_Loss, Surgery, Weight)
df <- na.omit(df)
df <- df[!(df$Timepoint %in% c(18, 24)),]

stat.test <- df %>%
  group_by(Timepoint) %>%
  wilcox_test(Weight ~ Surgery) %>%
  # adjust_pvalue(method = "BH") %>%
  add_significance()

stat.test <- stat.test %>%
  add_xy_position(x = "Timepoint", fun = "mean_se")

plot <- ggbarplot(df, "Timepoint", "Weight",
                  fill = "Surgery",
                  color = "black", 
                  palette = surgeryColors,
                  add = "mean_se",
                  label = FALSE, position = position_dodge()); plot

plot <- plot +
  labs(x = "Time (months)", y = "Weight (kg)")

plot <- plot +
  stat_pvalue_manual(
    stat.test,
    bracket.nudge.y = 0.5,
    # bracket.shorten = 1,
    bracket.size = 0.5,
    tip.length = 0.01,
    # remove.bracket = TRUE,
    size = 6,
    hide.ns = TRUE,
    # label = "p.adj.signif"
    label = "p.signif"
  )
plot
plotList[[index]] <- plot

plotFileName <- "surgery_Weight_barplots_BLto12M.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 5, height = 5)
for (i in c(1)) {
  grid.arrange(plotList[[i]],
               ncol = 1, nrow = 1)
}
dev.off()


df <- df[!(df$Timepoint %in% c(0)),]
stat.test <- df %>%
  group_by(Timepoint) %>%
  wilcox_test(Percent_Loss ~ Surgery) %>%
  # adjust_pvalue(method = "BH") %>%
  add_significance()

stat.test <- stat.test %>%
  add_xy_position(x = "Timepoint", fun = "mean_se")

plot <- ggbarplot(df, "Timepoint", "Percent_Loss",
                  fill = "Surgery",
                  color = "black", 
                  palette = surgeryColors,
                  add = "mean_se",
                  label = FALSE, position = position_dodge())

plot <- plot +
  labs(x = "Time (months)", y = "Weight loss (%)")

plot <- plot +
  stat_pvalue_manual(
    stat.test,
    bracket.nudge.y = 0.5,
    # bracket.shorten = 1,
    bracket.size = 0.5,
    tip.length = 0.01,
    # remove.bracket = TRUE,
    size = 6,
    hide.ns = TRUE,
    # label = "p.adj.signif"
    label = "p.signif"
  )
plot
plotList[[2]] <- plot

plotFileName <- "surgery_Percent_Loss_barplots_BLto12M.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 5, height = 5)
for (i in c(2)) {
  grid.arrange(plotList[[i]],
               ncol = 1, nrow = 1)
}
dev.off()




##### Bar plots Percent_Loss (BL-12M) between weight loss groups at each timepoint #####
for (division in divisions) {
  plotList <- list()
  PatientID <- myTable$PatientID
  Timepoint <- myTable$time
  Surgery <- myTable$Surgery
  # Quintile <- myTable$QuintileAssignment_12M
  Division <- myTable[,which(colnames(myTable) == paste0(division, "Assignment_12M"))]
  
  # Create smaller data table with relevant columns
  df <- data.frame(PatientID, Timepoint, Surgery, Division)
  df <- na.omit(df)
  df <- df[(df$Timepoint %in% c(12)),]
  
  df2 <- as.data.frame(table(df$Surgery, df$Division))
  names(df2)[names(df2) == "Var1"] <- "Surgery"
  names(df2)[names(df2) == "Var2"] <- "Division"
  names(df2)[names(df2) == "Freq"] <- "Frequency"
  
  # Chi-squared test
  stats.summ <- chisq_test(df$Surgery, df$Division)
  p_value <- roundP(stats.summ$p)
  
  subtitle.lab <- paste0("Chi-squared ", p_value)
  
  plot <- ggplot(df2, aes(fill=Surgery, y=Frequency, x=Division)) + 
    geom_bar(position="stack", stat="identity", color = "black")+
    labs(x = division, y = "Number of patients")+
    scale_fill_manual(values=surgeryColors); plot
  
  plot <- plot + labs(subtitle = subtitle.lab); plot
  
  plot <- plot + geom_text(aes(label = paste0(Frequency)), size = 3, hjust = 0.5, position = "stack", vjust = 1.5, colour = "white"); plot
  
  plotList[[1]] <- plot
  
  plotFileName <- paste0("surgery_", division, "_stackedbarplots_BLto12M.pdf")
  file.path <- paste0(outputDir, plotFileName)
  pdf(file.path, width = 5, height = 5)
  for (i in c(1)) {
    grid.arrange(plotList[[i]],
                 ncol = 1, nrow = 1)
  }
  dev.off()
  
}
##### Bar plots Percent_Loss (BL-24M) between surgery types at each timepoint #####
index <- 1
plotList <- list()

PatientID <- myTable$PatientID
Timepoint <- myTable$time
Percent_Loss <- myTable$Percent_Loss_kg
Surgery <- myTable$Surgery
Weight <- myTable$Weight_kg


# Create smaller data table with relevant columns
df <- data.frame(PatientID, Timepoint, Percent_Loss, Surgery, Weight)
df <- na.omit(df)

stat.test <- df %>%
  group_by(Timepoint) %>%
  wilcox_test(Weight ~ Surgery) %>%
  # adjust_pvalue(method = "BH") %>%
  add_significance()

stat.test <- stat.test %>%
  add_xy_position(x = "Timepoint", fun = "mean_se")

plot <- ggbarplot(df, "Timepoint", "Weight",
                  fill = "Surgery",
                  color = "black", 
                  palette = surgeryColors,
                  add = "mean_se",
                  label = FALSE, position = position_dodge()); plot

plot <- plot +
  labs(x = "Time (months)", y = "Weight (kg)")

plot <- plot +
  stat_pvalue_manual(
    stat.test,
    bracket.nudge.y = 0.5,
    # bracket.shorten = 1,
    bracket.size = 0.5,
    tip.length = 0.01,
    # remove.bracket = TRUE,
    size = 6,
    hide.ns = TRUE,
    # label = "p.adj.signif"
    label = "p.signif"
  )
plot
plotList[[index]] <- plot

plotFileName <- "surgery_Weight_barplots_BLto24M.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 5, height = 5)
for (i in c(1)) {
  grid.arrange(plotList[[i]],
               ncol = 1, nrow = 1)
}
dev.off()


df <- df[!(df$Timepoint %in% c(0)),]
stat.test <- df %>%
  group_by(Timepoint) %>%
  wilcox_test(Percent_Loss ~ Surgery) %>%
  # adjust_pvalue(method = "BH") %>%
  add_significance()

stat.test <- stat.test %>%
  add_xy_position(x = "Timepoint", fun = "mean_se")

plot <- ggbarplot(df, "Timepoint", "Percent_Loss",
                  fill = "Surgery",
                  color = "black", 
                  palette = surgeryColors,
                  add = "mean_se",
                  label = FALSE, position = position_dodge())

plot <- plot +
  labs(x = "Time (months)", y = "Weight loss (%)")

plot <- plot +
  stat_pvalue_manual(
    stat.test,
    bracket.nudge.y = 0.5,
    # bracket.shorten = 1,
    bracket.size = 0.5,
    tip.length = 0.01,
    # remove.bracket = TRUE,
    size = 6,
    hide.ns = TRUE,
    # label = "p.adj.signif"
    label = "p.signif"
  )
plot
plotList[[2]] <- plot

plotFileName <- "surgery_Percent_Loss_barplots_BLto24M.pdf"
file.path <- paste0(outputDir, plotFileName)
pdf(file.path, width = 5, height = 5)
for (i in c(2)) {
  grid.arrange(plotList[[i]],
               ncol = 1, nrow = 1)
}
dev.off()





##### Bar plots Percent_Loss (BL-24M) between weight loss groups at each timepoint #####
for (division in divisions) {
  plotList <- list()
  PatientID <- myTable$PatientID
  Timepoint <- myTable$time
  Surgery <- myTable$Surgery
  # Quintile <- myTable$QuintileAssignment_24M
  Division <- myTable[,which(colnames(myTable) == paste0(division, "Assignment_24M"))]
  
  # Create smaller data table with relevant columns
  df <- data.frame(PatientID, Timepoint, Surgery, Division)
  df <- na.omit(df)
  df <- df[(df$Timepoint %in% c(24)),]
  
  df2 <- as.data.frame(table(df$Surgery, df$Division))
  names(df2)[names(df2) == "Var1"] <- "Surgery"
  names(df2)[names(df2) == "Var2"] <- "Division"
  names(df2)[names(df2) == "Freq"] <- "Frequency"
  
  # Chi-squared test
  stats.summ <- chisq_test(df$Surgery, df$Division)
  p_value <- roundP(stats.summ$p)
  
  subtitle.lab <- paste0("Chi-squared ", p_value)
  
  plot <- ggplot(df2, aes(fill=Surgery, y=Frequency, x=Division)) + 
    geom_bar(position="stack", stat="identity", color = "black")+
    labs(x = division, y = "Number of patients")+
    scale_fill_manual(values=surgeryColors); plot
  
  plot <- plot + labs(subtitle = subtitle.lab); plot
  
  plot <- plot + geom_text(aes(label = paste0(Frequency)), size = 3, hjust = 0.5, position = "stack", vjust = 1.5, colour = "white"); plot
  
  plotList[[1]] <- plot
  
  plotFileName <- paste0("surgery_", division, "_stackedbarplots_BLto24M.pdf")
  file.path <- paste0(outputDir, plotFileName)
  pdf(file.path, width = 5, height = 5)
  for (i in c(1)) {
    grid.arrange(plotList[[i]],
                 ncol = 1, nrow = 1)
  }
  dev.off()
  
}

##### Plots for %EWL correlations (Kendall) - 12, 18, and 24 months only #####
PatientID <- myTable$PatientID
Timepoint <- myTable$Timepoint
PEWL <- myTable$PEWL_kg

df <- data.frame(PatientID, Timepoint, PEWL)
df <- na.omit(df)

df <- spread(df, Timepoint, PEWL)

MONTHS <- c("TWELVE", "EIGHTEEN", "TWENTY_FOUR")
mo <- c(12,18,24)
y <- 2
all_combinations <- combn(MONTHS,y)
all_combinations <- all_combinations[,-1]

plotList <- list()
index <- 1
tagIndex <- 1

for (i in 1:ncol(all_combinations)) {
  
  m1 <- all_combinations[1,i]
  m2 <- all_combinations[2,i]
  
  col1 <- which(MONTHS == m1)
  col2 <- which(MONTHS == m2)
  
  plot <- ggscatter(df, x = m1, y = m2,
                    shape = 21, size = 2.5, # Points shape and size
                    add = "reg.line",  # Add regression line
                    conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "kendall"),
                    add.params = list(fill = "lightgray")
  )
  tag <- tags[tagIndex]
  x.lab <- paste0(mo[col1], "-month excess weight loss (%)")
  y.lab <- paste0(mo[col2], "-month excess weight loss (%)")
  plot <- plot + labs(x=x.lab, y = y.lab, tag = tag)
  plot <- plot + geom_abline(intercept = 0, slope = 1, color = "red", linetype="dashed")
  # plot <- plot + stat_regline_equation()
  plot
  plotList[[i]] <- plot
  tagIndex <- tagIndex + 1
  
}


file.path <- paste0(outputDir, "excess_weight_loss_kendall_plots_12_18_24_only.pdf")
pdf(file.path, width = 7, height = 3.5)
for (i in seq(1, length(plotList), 2)) {
  # message(i)
  grid.arrange(plotList[[i]], plotList[[i+1]], 
               ncol = 2, nrow = 1)
}
dev.off()

