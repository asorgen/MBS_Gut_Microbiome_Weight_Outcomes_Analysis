#Author: Alicia Sorgen
#Date: 2022 August 16
#Description: Kendall correlation modeling weight loss by microbiome

##### Edits for script #####
rm(list=ls())

ANALYSIS <- "microbiome_n124"
params <- vector()
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
params <- c(params, "MetaPhlAn2")
params <- c(params, 0)
params <- c(params, 1)
params <- c(params, 6)
params <- c(params, 12)
params <- c(params, 18)
params <- c(params, 24)
end_month <- params[length(params)]

levels <- vector()
levels <- c(levels, "phylum")
# levels <- c(levels, "class")
# levels <- c(levels, "order")
# levels <- c(levels, "family")
levels <- c(levels, "genus")
# levels <- c(levels, "species")

moduleRoot <- "Microbiome_WeightLoss_Associations_kendall"

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
  
  if (args[2] == "MetaPhlAn2" | args[2] == "Kraken2") {
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

prevModule <- paste0(classifier, "_TaxaMetaMerge")
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")
CountFile <- paste0("_LogNormalizedCounts_", classifier, ".tsv")

##### Set up output #####
outputDir = file.path(moduleDir,"output/")

##### Prep for analysis #####
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")

factorLevels <- c("BL-1M", "BL-6M", "BL-12M", "BL-18M", "BL-24M")
surgeryPalette <- c("steelblue", "gold")

##### Kendall correlation #####
for (level in levels) {

  message(paste0("\n***** Starting ", level, " *****\n"))

  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)

  file.path <- paste0(inputDir,level, CountFile)
  logCounts<-read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)

  startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1

  # Parse out metadata from counts
  myT<-logCounts[,startAbundanceIndex:ncol(logCounts)]

  PatientID <- logCounts$PatientID
  Timepoint_kg <- paste0("BL-", logCounts$time)
  Percent_Loss_kg <- logCounts$Percent_Loss_kg

  Taxa <- vector()
  Taxa_Timepoint <- vector()
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

    bug <- myT[,i]
    bugName <- colnames(myT)[i]

    if ( mean(bug > 0, na.rm = TRUE) > 0.1 ) {

      df <- data.frame( bug, PatientID, Timepoint_kg, Percent_Loss_kg )

      for (t in 1:length(included)) {

        taxaMonth <- included[t]

        # id rows
        xM <- which(logCounts$time == taxaMonth)

        # abundance
        bug_xM <- df[xM, "bug"]
        names(bug_xM) <- df[xM, "PatientID"]

        # Add bug info to each row
        df$bug_xM <- bug_xM[df$PatientID]

        tIndex <- 1
        for (w in 2:length(included)) {

          PWL_BL_xM <- paste0("BL-", included[w])

          df2 <- df[ df$Timepoint_kg == PWL_BL_xM, ]
          df2 <- na.omit(df2)

          myKen <- cor.test( df2$bug_xM, df2$Percent_Loss_kg, method = "kendall" )
          mySpear <- cor.test( df2$bug_xM, df2$Percent_Loss_kg, method = "spearman" )

          Taxa[index] <- bugName
          Taxa_Timepoint[index] <- taxaMonth
          Weight_Loss_Timepoint[index] <- PWL_BL_xM
          p_Kendall[index] <- myKen$p.value
          R_Kendall[index] <- myKen$estimate
          p_Spearman[index] <- mySpear$p.value
          R_Spearman[index] <- mySpear$estimate
          n[index] <- nrow(df2)

          index <- index + 1
          tIndex <- tIndex + 1

        } # for (w in 2:length(included))

      }  # for (m in included)

    } # if ( mean(bug > 0, na.rm = TRUE) > 0.1 )

  } # for (i in 1:ncol(myT))

  dFrame <- data.frame( Taxa, Taxa_Timepoint, Weight_Loss_Timepoint, n, p_Kendall, R_Kendall, p_Spearman, R_Spearman)

  final <- data.frame()
  for (x in unique(dFrame$Taxa_Timepoint)) {

    dFrame.Filt <- dFrame[dFrame$Taxa_Timepoint == x,]

    for (y in unique(dFrame$Weight_Loss_Timepoint)) {

      dFrame.Filt2 <- dFrame.Filt[dFrame.Filt$Weight_Loss_Timepoint == y,]
      dFrame.Filt2$pAdj_Kendall <- p.adjust(dFrame.Filt2$p_Kendall, method = "BH")
      dFrame.Filt2$pAdj_Spearman <- p.adjust(dFrame.Filt2$p_Spearman, method = "BH")
      final <- rbind(final, dFrame.Filt2)

    }

  }

  final <- final[ order(final$p_Kendall), ]

  file.path <- paste0(outputLevel, level, "_taxa_by_weight_loss_ken_spear_results_BLto", included[length(included)], "M.tsv")
  write.table(final, file.path, sep="\t",quote = FALSE, row.names = FALSE)


} # for (level in levels)



##### Plot generation #####
for (level in levels) {
  
  message(paste0("\n***** Starting ", level, " *****\n"))
  
  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  
  file.path <- paste0(inputDir,level, CountFile)
  logCounts<-read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
  
  file.path <- paste0(outputLevel, level, "_taxa_by_weight_loss_ken_spear_results_BLto", included[length(included)], "M.tsv")
  stats.df<-read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
  
  startAbundanceIndex <- which(colnames(logCounts)=="ResponderStatus")+1
  
  # Parse out metadata from counts
  myT<-logCounts[,startAbundanceIndex:ncol(logCounts)]
  
  PatientID <- logCounts$PatientID
  Timepoint_kg <- paste0("BL-", logCounts$time)
  Percent_Loss_kg <- logCounts$Percent_Loss_kg
  SurgeryType <- logCounts$Surgery
  
  plotList <- list()
  dFrame <- data.frame()
  index <- 1
  
  for (i in 1:ncol(myT)) {
    
    bug <- myT[,i]
    bugName <- colnames(myT)[i]
    
    stats.df_1 <- stats.df[stats.df$Taxa == bugName,]
    
    if ( mean(bug > 0, na.rm = TRUE) > 0.1 ) {
      
      df <- data.frame( bug, PatientID, Timepoint_kg, Percent_Loss_kg, SurgeryType )
      
      for (t in 1:length(included)) {
        
        taxaMonth <- included[t]
        # message("taxaMonth = ", taxaMonth)
        
        stats.df_2 <- stats.df_1[stats.df_1$Taxa_Timepoint == taxaMonth,]
        
        # id rows
        xM <- which(logCounts$time == included[t])
        
        # abundance
        bug_xM <- df[xM, "bug"]
        names(bug_xM) <- df[xM, "PatientID"]
        
        # Add bug info to each row
        df$bug_xM <- bug_xM[df$PatientID]
        
        # if (t == 1) {
          for (w in 2:length(included)) {

            PWL_BL_xM <- paste0("BL-", included[w])
            
            stats.df_3 <- stats.df_2[stats.df_2$Weight_Loss_Timepoint == PWL_BL_xM,]
            pAdj <- roundP(stats.df_3$pAdj_Kendall[1])
            p <- roundP(stats.df_3$p_Kendall[1])
            title.lab <- paste0("Kendall ", p, ", adj ", pAdj)
            subtitle.lab <- paste0("Spearman R = ", round(stats.df_3$R_Spearman[1], 3))
            color <- ifelse(stats.df_3$pAdj[1] < 0.05, "red", "black")
            
            df2 <- df[ df$Timepoint_kg == PWL_BL_xM, ]
            df2 <- na.omit(df2)
            
            plot <- ggscatter(df2, x = "bug_xM", y = "Percent_Loss_kg", color = "SurgeryType", palette = surgeryPalette,
                              shape = 19, size = 2.5 # Points shape and size
                              # add = "reg.line",  # Add regression line
                              # conf.int = TRUE, # Add confidence interval
                              # add.params = list(fill = "lightgray")
            ); plot
            
            x.lab <- paste0(bugName, " relative abundance at ", month.labs[t])
            y.lab <- paste0("Weight loss (%): Baseline - ", month.labs[w])
            plot <- plot + labs(x=x.lab, y = y.lab, title = title.lab, subtitle = subtitle.lab); plot
            # plot <- plot + theme(plot.title = element_text(colour = color)); plot
            
            plotList[[index]] <- plot
            index <- index + 1
          }
        # } else {
        #   for (w in t:length(included)) {
        # 
        #     PWL_BL_xM <- paste0("BL-", included[w])
        #     
        #     stats.df_3 <- stats.df_2[stats.df_2$Weight_Loss_Timepoint == PWL_BL_xM,]
        #     p <- roundP(stats.df_3$pAdj[1])
        #     title.lab <- paste0("R^2 = ", round(stats.df_3$R_Spearman[1], 3), ", adj ", p)
        #     color <- ifelse(stats.df_3$pAdj[1] < 0.05, "red", "black")
        # 
        #     df2 <- df[ df$Timepoint_kg == PWL_BL_xM, ]
        #     
        #     plot <- ggscatter(df2, x = "bug_xM", y = "Percent_Loss_kg",
        #                       shape = 21, size = 2.5, # Points shape and size
        #                       add = "reg.line",  # Add regression line
        #                       conf.int = TRUE, # Add confidence interval
        #                       add.params = list(fill = "lightgray")
        #     ); plot
        #     
        #     x.lab <- paste0(bugName, " relative abundance at ", month.labs[t])
        #     y.lab <- paste0("Weight loss (%): Baseline - ", month.labs[w])
        #     plot <- plot + labs(x=x.lab, y = y.lab, title = title.lab); plot
        #     plot <- plot + theme(plot.title = element_text(colour = color))
        # 
        #     plotList[[index]] <- plot
        #     index <- index + 1
        #   }
        # }
        
      }  # for (m in included)
      
    } # if ( mean(bug > 0, na.rm = TRUE) > 0.1 )
    
  } # for (i in 1:ncol(myT))
  
  file.path <- paste0(outputLevel, level, "_taxa_by_weight_loss_kendall_BLto", included[length(included)], "M.pdf")
  pdf(file.path, width = 10, height = 10)
  mod <- length(plotList) %% 4
  if ( mod == 0 ) {
    
    for (i in seq(1, length(plotList), 4)) {
      grid.arrange(plotList[[i]], plotList[[i+1]], 
                   plotList[[i+2]],plotList[[i+3]],
                   ncol = 2, nrow = 2)
    }
    
  } else {
    
    for (i in seq(1, (length(plotList)-mod), 4)) {
      grid.arrange(plotList[[i]], plotList[[i+1]], 
                   plotList[[i+2]],plotList[[i+3]],
                   ncol = 2, nrow = 2)
      j <- i
    }
    j <- j + 1
    if ( mod == 3) {
      grid.arrange(plotList[[j]], plotList[[j+1]], 
                   plotList[[j+2]],
                   ncol = 2, nrow = 2)
    }
    
    if ( mod == 2) {
      grid.arrange(plotList[[j]], plotList[[j+1]], 
                   ncol = 2, nrow = 2)
    }
    
    if ( mod == 1) {
      grid.arrange(plotList[[j]],
                   ncol = 2, nrow = 2)
    }
    
  }
  dev.off()
  

} # for (level in levels)

##### Summary tables (unadjusted p values) #####
for (level in levels) {
  
  message(paste0("\n***** Starting ", level, " *****\n"))
  
  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  
  file.path <- paste0(outputLevel, level, "_taxa_by_weight_loss_ken_spear_results_BLto", included[length(included)], "M.tsv")
  stats.df<-read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
  
  TaxaTime <- vector()
  WeightLossTime <- vector()
  nSignif <- vector()
  index <- 1
  
  for (t in 1:length(included)) {
    
    stats.df_1 <- stats.df[stats.df$Taxa_Timepoint == included[t], ]
    
    for (w in 2:length(included)) {
      
      stats.df_2 <- stats.df_1[stats.df_1$Weight_Loss_Timepoint == paste0("BL-", included[w]), ]
      nSignif[index] <- ifelse( nrow(stats.df_2) > 0, length(which(stats.df_2$p_Kendall < 0.05)), NA )
      TaxaTime[index] <- month.labs[t]
      WeightLossTime[index] <- paste0("BL-", included[w], "M")
      index <- index + 1
      
    }
  }
  dFrame <- data.frame( TaxaTime, WeightLossTime, nSignif )
  dFrame$TaxaTime <- as.factor(dFrame$TaxaTime)
  dFrame$TaxaTime <- factor(dFrame$TaxaTime, levels = month.labs)
  dFrame$WeightLossTime <- as.factor(dFrame$WeightLossTime)
  dFrame$WeightLossTime <- factor(dFrame$WeightLossTime, levels = factorLevels[1:(length(included)-1)] )

  x.lab <- "Weight Loss"
  y.lab <- str_to_title(level)
  # Heatmap 
  plot <- ggplot(dFrame, aes(WeightLossTime, TaxaTime, fill= nSignif)) + 
    geom_tile()+
    # scale_fill_distiller(palette = "RdPu"); plot
    scale_fill_gradient(low="white", high="blue"); plot
  
  plot <- plot + labs(x=x.lab, y = y.lab); plot
  
  plot <- plot + labs(fill = "Significance\n(p values)"); plot
  
  file.path <- paste0(outputLevel, level, "_taxa_by_weight_loss_ken_spear_results_BLto", included[length(included)], "M_unadjP_heatmap.pdf")
  pdf(file.path, width = 7, height = 7)
  print(plot)
  dev.off()
  
  dFrame2 <- spread( dFrame, TaxaTime, nSignif )
  
  file.path <- paste0(outputLevel, level, "_taxa_by_weight_loss_ken_spear_results_BLto", included[length(included)], "M_summary_unadjP.tsv")
  write.table(dFrame2, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  
} # for (level in levels)

##### Summary tables (adjusted p values) #####
for (level in levels) {
  
  message(paste0("\n***** Starting ", level, " *****\n"))
  
  outputLevel <- paste0(outputDir, level, "/")
  dir.create(outputLevel, showWarnings = FALSE)
  
  file.path <- paste0(outputLevel, level, "_taxa_by_weight_loss_ken_spear_results_BLto", included[length(included)], "M.tsv")
  stats.df<-read.table(file.path, sep="\t", header = TRUE, check.names = FALSE)
  
  TaxaTime <- vector()
  WeightLossTime <- vector()
  nSignif <- vector()
  index <- 1
  
  for (t in 1:length(included)) {
    
    stats.df_1 <- stats.df[stats.df$Taxa_Timepoint == included[t], ]
    
    for (w in 2:length(included)) {
      
      stats.df_2 <- stats.df_1[stats.df_1$Weight_Loss_Timepoint == paste0("BL-", included[w]), ]
      nSignif[index] <- ifelse( nrow(stats.df_2) > 0, length(which(stats.df_2$pAdj_Kendall < 0.05)), NA )
      TaxaTime[index] <- month.labs[t]
      WeightLossTime[index] <- paste0("BL-", included[w], "M")
      index <- index + 1
      
    }
  }
  dFrame <- data.frame( TaxaTime, WeightLossTime, nSignif )
  dFrame$TaxaTime <- as.factor(dFrame$TaxaTime)
  dFrame$TaxaTime <- factor(dFrame$TaxaTime, levels = month.labs)
  dFrame$WeightLossTime <- as.factor(dFrame$WeightLossTime)
  dFrame$WeightLossTime <- factor(dFrame$WeightLossTime, levels = factorLevels[1:(length(included)-1)] )
  
  x.lab <- "Weight Loss"
  y.lab <- str_to_title(level)
  # Heatmap 
  plot <- ggplot(dFrame, aes(WeightLossTime, TaxaTime, fill= nSignif)) + 
    geom_tile()+
    # scale_fill_distiller(palette = "RdPu"); plot
    scale_fill_gradient(low="white", high="blue"); plot

  plot <- plot + labs(x=x.lab, y = y.lab); plot
  
  plot <- plot + labs(fill = "Significance\n(p values)"); plot
  
  file.path <- paste0(outputLevel, level, "_taxa_by_weight_loss_ken_spear_results_BLto", included[length(included)], "M_adjP_heatmap.pdf")
  pdf(file.path, width = 7, height = 7)
  print(plot)
  dev.off()
  
  dFrame2 <- spread( dFrame, TaxaTime, nSignif )
  
  file.path <- paste0(outputLevel, level, "_taxa_by_weight_loss_ken_spear_results_BLto", included[length(included)], "M_summary_adjP.tsv")
  write.table(dFrame2, file.path, sep="\t",quote = FALSE, row.names = FALSE)
  
  
} # for (level in levels)

