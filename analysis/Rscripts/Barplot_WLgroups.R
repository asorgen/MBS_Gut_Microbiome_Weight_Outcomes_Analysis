#Author: Alicia Sorgen
#Date: 2022 October 21
#Description: Generate bar plots


##### Edits for script #####
rm(list=ls())
params <- vector()

ANALYSIS <- "microbiome_n124"
params <- c(params, "~/git/gut-microbiome-bariatric-weight-outcomes")
params <- c(params, "NA")
params <- c(params, "NA")

moduleRoot <- paste0("Barplot_WLgroups")

##### Libraries #####
R <- sessionInfo()
message("\n************* R Environment *************\n")
message(R$R.version$version.string)
rm(R)

message("\n************* R packages *************\n")
library(stringr); message("stringr: Version ", packageVersion("stringr"))
library(ggplot2); message("ggplot2: Version ", packageVersion("ggplot2"))
library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
library(scales); message("scales: Version ", packageVersion("scales"))
library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
library(nlme); message("nlme: Version ", packageVersion("nlme"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
# library(data.table); message("data.table: Version ", packageVersion("data.table"))

##### Set up working environment #####
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args <- params
}
rm(params)

if (args[1] == "BLJ") {
  message("\n************* Running in BioLockJ *************")
} else {
  message("\n************* Running locally *************")
  gitRoot <- args[1]
  gitInput <- file.path(gitRoot, "analysis", "input")
  gitScripts <- file.path(gitRoot, "analysis", "Rscripts")
  rm(gitRoot)
  
  root <- paste0("~/BioLockJ_pipelines/")
  pipeline <- paste0(ANALYSIS,"_analysis_")
  
  if (any(dir(root) %like% pipeline) == TRUE) {
    root <- paste0(root,"/",str_subset(dir(root), pipeline), "/")
  } else {
    today <- as.character(format(Sys.Date(), "%Y%b%d"))
    root <- paste0(root,ANALYSIS,"_analysis_", today, "/")
    dir.create(root, showWarnings = FALSE)
  }
  rm(pipeline)
  
  if (any(dir(root) == "input") == FALSE) {
    file.copy(gitInput,
              root,
              recursive = TRUE)
    
  }
  rm(gitInput)
  module <- moduleRoot
  
  if (args[2] == "MetaPhlAn2" | args[2] == "Kraken2") {
    module <- paste0(args[2], "_", module)
  } 
  
  if (args[3] %in% c("Quintile", "Quartile", "Thirds", "Half")) {
    module <- paste0(module, "_", args[3])
  } 
  
  if (exists("end_month") == TRUE) {
    module <- paste0(module, "_BLto", end_month, "M")
    rm(end_month)
  }
  
  if (any(dir(root) %like% module) == TRUE) {
    moduleDir <- paste0(root,str_subset(dir(root), module), "/")
  } else {
    moduleDir <- paste0(root, module, "/")
    dir.create(moduleDir, showWarnings = FALSE)
  }
  rm(module, root)
  
  scriptDir <- paste0(moduleDir, "script/")
  if (any(dir(moduleDir) == "script") == FALSE) {
    dir.create(scriptDir, showWarnings = FALSE)
  }
  rm(scriptDir)
  
  outputDir <- paste0(moduleDir, "output/")
  if (any(dir(moduleDir) == "output") == FALSE) {
    dir.create(outputDir, showWarnings = FALSE)
  }
  
  files <- list.files(outputDir, recursive = TRUE, full.names = TRUE)
  file.remove(files)
  rm(outputDir, files)
  
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
              recursive = TRUE)
    rm(script)
    
  }
  rm(resourcesDir, gitScripts)
  
  setwd(paste0(moduleDir, "script/"))
  rm(moduleDir)
  
}
rm(moduleRoot, ANALYSIS)

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

