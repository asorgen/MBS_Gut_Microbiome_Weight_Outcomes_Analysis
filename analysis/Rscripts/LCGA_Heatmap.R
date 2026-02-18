#Author: Alicia Sorgen
#Date: 2023 Nov 29
#Description: Linear modeling for taxa by time

##----- Edits for script --------------------------------------------------------------------------
rm(list=ls())
set.seed(1989)
date = "2024Mar05"
ANALYSIS <- "MetaPhlAn2_microbiome"
params <- vector()
params <- c(params, "~/git/BS_MicrobiomeAnalysis_2022")
params <- c(params, "MetaPhlAn2")
# params <- c(params, 0)
# params <- c(params, 1)
# params <- c(params, 6)
# params <- c(params, 12)
# params <- c(params, 18)
params <- c(params, 24)


levels <- vector()
levels <- c(levels, "phylum")
# levels <- c(levels, "class")
# levels <- c(levels, "order")
# levels <- c(levels, "family")
levels <- c(levels, "genus")
# levels <- c(levels, "species")

moduleRoot <- "LCGA_Heatmap"

##----- Libraries ---------------------------------------------------------------------------------
R <- sessionInfo()
message(R$R.version$version.string); rm(R)

library(stringr); message("stringr: Version ", packageVersion("stringr"))
# library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
# library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
# library(nlme); message("nlme: Version ", packageVersion("nlme"))
# library(data.table); message("data.table: Version ", packageVersion("data.table"))
# library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
# library(ggrepel); message("ggrepel: Version ", packageVersion("ggrepel"))
library(ComplexHeatmap); message("ComplexHeatmap: Version ", packageVersion("ComplexHeatmap"))

##----- Set up working environment ----------------------------------------------------------------
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
  pipeline <- paste0(ANALYSIS,"_analysis_", date)
  
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
  
  if (exists("end_month") == TRUE) {
    mod1 <- paste0(sapply(strsplit(module, "_Weight"), "[", 1), "_"); mod1
    mod2 <- paste0("_Weight", sapply(strsplit(module, "_Weight"), "[", 2)); mod2
    module <- paste0(mod1, end_month, "M", mod2, "_BLto24M"); module
  }
  
  if (any(dir(root) %like% module) == TRUE) {
    moduleDir <- paste0(root,str_subset(dir(root), module), "/")
  } else {
    moduleDir <- paste0(root, module, "/")
    dir.create(moduleDir, showWarnings = FALSE)
  }; rm(root, module)
  
  scriptDir <- paste0(moduleDir, "script/")
  if (any(dir(moduleDir) == "script") == FALSE) {
    dir.create(scriptDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/", moduleRoot, ".R"); script
    file.copy(script,
              scriptDir,
              recursive = TRUE)
  }; rm(scriptDir, moduleRoot)
  
  outputDir <- paste0(moduleDir, "output/")
  if (any(dir(moduleDir) == "output") == FALSE) {
    dir.create(outputDir, showWarnings = FALSE)
  }
  
  files <- list.files(outputDir, recursive = TRUE, full.names = TRUE); rm(outputDir)
  # file.remove(files)
  rm(files)
  
  resourcesDir <- paste0(moduleDir, "resources/")
  if (any(dir(moduleDir) == "resources") == FALSE) {
    dir.create(resourcesDir, showWarnings = FALSE)
    
    script = paste0(gitScripts,"/functions.R")
    file.copy(script,
              resourcesDir,
              recursive = TRUE)
  }; rm(gitScripts, resourcesDir)
  setwd(paste0(moduleDir, "script/")); rm(moduleDir)
  
}

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())

##----- Set up functions file ---------------------------------------------------------------------
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

str <- sapply(strsplit(pipeRoot, "_analysis"), "[", 1)
str <- str_split(str, "/")
str <- do.call(rbind, str)
ANALYSIS <- str[length(str)]
rm(str)

##----- Set up input ------------------------------------------------------------------------------

inputDir = file.path(pipeRoot, "input")
count.InputDir = file.path(inputDir, paste0("MetaPhlan2TaxaTables/"))
taxonomyFile <- "Full_Taxonomy.tsv"

prevModule <- str_subset(dir(pipeRoot), "Univariate_RYGB_LCGA_Model2_2Class_Analysis"); message(prevModule)
LCGADir = paste0(pipeRoot,"/",prevModule,"/output"); message(LCGADir)
pvalue.File <- list.files(file.path(LCGADir, "LinearModel_Analysis_pvalues", "ResultsTables"), full.names = TRUE); message(pvalue.File)

##----- Set up output -----------------------------------------------------------------------------
outputDir = file.path(moduleDir,"output/")

##----- Script variables --------------------------------------------------------------------------
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")
new_getlog10p<-function(pval,coefficient){
  
  log10p<-sapply(1:length(pval), function(x){
    if (coefficient[x] == "increase") return(-log10(pval[x]))
    else return(log10(pval[x]))
  })
}

##----- Read table ----------------------------------------------------

taxonomyTable <- read.delim(paste0(count.InputDir, taxonomyFile), sep="\t",header = TRUE)
message(">>> taxa table")

level = "genus"
myT <- read.table(pvalue.File, sep="\t",header = TRUE, check.names = FALSE)
message(">>> p-value table")
##---- Set up for heatmap -----
myT$groupPValues4 <- ifelse(myT$groupPValues4 == 0, 2e-16, myT$groupPValues4)

log10p <- new_getlog10p(myT$groupPValues4, myT$direction)
comparisons <- paste0(myT$LengthTime)
taxaNames <- paste0(myT$qbugName, "&", myT$groupNumbers)
adjp <- myT$adjPValues4

df <- data.frame(taxaNames, comparisons, log10p)
adjP.df <- data.frame(taxaNames, comparisons, adjp)


df <- spread(df, comparisons, log10p)
adjP.df <- spread(adjP.df, comparisons, adjp)

df <- df[!(df$taxaNames %like% "Other"),]
adjP.df <- adjP.df[!(adjP.df$taxaNames %like% "Other"),]

taxaNames <- df$taxaNames
BLto1M <- df$BLto1M
BLto6M <- df$BLto6M
BLto12M <- df$BLto12M
BLto18M <- df$BLto18M
BLto24M <- df$BLto24M
df <- data.frame(taxaNames, BLto1M, BLto6M, BLto12M, BLto18M, BLto24M)
BLto1M <- adjP.df$BLto1M
BLto6M <- adjP.df$BLto6M
BLto12M <- adjP.df$BLto12M
BLto18M <- adjP.df$BLto18M
BLto24M <- adjP.df$BLto24M
adjP.df <- data.frame(taxaNames, BLto1M, BLto6M, BLto12M, BLto18M, BLto24M)

rownames(df) <- df$taxaNames
rownames(adjP.df) <- adjP.df$taxaNames

df <- df[,-1]
adjP.df <- adjP.df[,-1]

df <- df[order(rownames(df)),]
adjP.df <- adjP.df[order(rownames(adjP.df)),]

if (level == "genus") {
  
  # Find rows where all values are between the specified range
  df <- df[!rowSums(adjP.df >= 0.05, na.rm = TRUE) == ncol(adjP.df), ]
  adjP.df <- adjP.df[!rowSums(adjP.df >= 0.05, na.rm = TRUE) == ncol(adjP.df), ]
}

df.1 <- df[which(rownames(df) %like% "&1"),]
df.2 <- df[which(rownames(df) %like% "&2"),]

# rownames(df.1) <- gsub("&1", "", rownames(df.1))
group1_Taxa <- nrow(df.1)
# rownames(df.2) <- gsub("&2", "", rownames(df.2))
group2_Taxa <- nrow(df.2)
df <- rbind(df.1, df.2)
df <- as.matrix(df)

adjP.df.1 <- adjP.df[which(rownames(adjP.df) %like% "&1"),]
# rownames(adjP.df.1) <- gsub("&1", "", rownames(adjP.df.1))
group1_Taxa <- nrow(adjP.df.1)
adjP.df.2 <- adjP.df[which(rownames(adjP.df) %like% "&2"),]
# rownames(adjP.df.2) <- gsub("&2", "", rownames(adjP.df.2))
group2_Taxa <- nrow(adjP.df.2)
adjP.df <- rbind(adjP.df.1, adjP.df.2)
adjP.df <- as.matrix(adjP.df)



##---- Heatmap -----
# df<-t(df)
# df.adj<-t(df.adj)

# Specify the range
lower_limit <- log10(0.05)
upper_limit <- -log10(0.05)


ra = rowAnnotation(
  empty = anno_empty(border = FALSE),
  foo2 = anno_block(gp = gpar(fill = 2:3), labels = c("Group 1" , "Group 2"))
)

splitRow = c(rep(1, each = group1_Taxa), rep(2, each = group2_Taxa))

genera <- sapply(strsplit(rownames(df.1), "&"), "[", 1)
taxonomyTable1 <- taxonomyTable[,c("Phylum", "Genus")]
taxonomyTable1 <- taxonomyTable1[!duplicated(taxonomyTable1$Genus),]
taxonomyTable1 <- taxonomyTable1[taxonomyTable1$Genus %in% genera,]
taxonomyTable1 <- taxonomyTable1[order(taxonomyTable1$Genus), ]
phyla1 <- taxonomyTable1$Phylum

genera <- sapply(strsplit(rownames(df.2), "&"), "[", 1)
taxonomyTable2 <- taxonomyTable[,c("Phylum", "Genus")]
taxonomyTable2 <- taxonomyTable2[!duplicated(taxonomyTable2$Genus),]
taxonomyTable2 <- taxonomyTable2[taxonomyTable2$Genus %in% genera,]
taxonomyTable2 <- taxonomyTable2[order(taxonomyTable2$Genus), ]
phyla2 <- taxonomyTable2$Phylum
phyla <- c(phyla1, phyla2)

col2 = list(Phylum = c("Actinobacteria" = "yellow2", 
                       "Bacteroidetes" = "limegreen",
                       "Firmicutes" = "skyblue1",
                       "Fusobacteria" = "orange",
                       "Proteobacteria" = "purple1",
                       "Verrucomicrobia" = "hotpink"))

right_ha <- rowAnnotation(
  Phylum = phyla, col = col2
)

rownames(df) <- sapply(strsplit(rownames(df), "&"), "[", 1)
# Replace all NAs with 1
adjP.df <- replace(adjP.df, is.na(adjP.df), 1)

pdf(file.path(outputDir, paste0("genus_HeatmapAndCluster_Classified.pdf")),height = 8)
Heatmap(df,
        # top_annotation = ha,
        left_annotation = ra,
        right_annotation = right_ha,
        rect_gp = gpar(col = "white", lwd = 1),
        name = "log10 p-value",
        row_split = splitRow,
        column_title = NULL,
        row_title = NULL,
        border = TRUE,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 6),
        # heatmap_height = unit(40, "cm"),
        layer_fun = function(j, i, x, y, w, h, fill) {
          ind_mat = restore_matrix(j, i, x, y)
          for(ir in seq_len(nrow(ind_mat))) {
            for(ic in seq_len(ncol(ind_mat))) {
              ind1 = ind_mat[ir, ic] # previous column
              v1 = adjP.df[i[ind1], j[ind1]]
              # print(v1)
              if(v1 < 0.05) { # if they have the same sign
                grid.points(x[c(ind1)], y[c(ind1)],
                            pch = 8, size = unit(0.5, "mm"), 
                            gp = gpar(col = "black"))
              }
            } # for(ic in seq_len(ncol(ind_mat)))
          } # for(ir in seq_len(nrow(ind_mat)))
        }
)
dev.off()
##---- end -----

