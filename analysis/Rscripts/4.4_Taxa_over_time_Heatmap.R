#Author: Alicia Sorgen
#Date: 2023 Nov 29
#Description: Linear modeling for taxa by time

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
# library(ggpubr); message("ggpubr: Version ", packageVersion("ggpubr"))
# library(tidyr); message("tidyr: Version ", packageVersion("tidyr"))
# library(rstatix); message("rstatix: Version ", packageVersion("rstatix"))
# library(nlme); message("nlme: Version ", packageVersion("nlme"))
library(data.table); message("data.table: Version ", packageVersion("data.table"))
# library(gridExtra); message("gridExtra: Version ", packageVersion("gridExtra"))
# library(ggrepel); message("ggrepel: Version ", packageVersion("ggrepel"))
library(ComplexHeatmap); message("ComplexHeatmap: Version ", packageVersion("ComplexHeatmap"))


# Script-specific edits -------------------------------------------------
params <- vector()
params <- c(params, getOption("mbs.pipe_root", default = stop("Set options(mbs.pipe_root = ...) in your local .Rprofile (gitignored)")))
params <- c(params, "MetaPhlAn2")
# params <- c(params, 0)
# params <- c(params, 1)
# params <- c(params, 6)
# params <- c(params, 12)
# params <- c(params, 18)
params <- c(params, 24)

module <- "4.4_Taxa_over_time_Heatmap"

levels <- vector()
levels <- c(levels, "phylum")
# levels <- c(levels, "class")
# levels <- c(levels, "order")
# levels <- c(levels, "family")
levels <- c(levels, "genus")
# levels <- c(levels, "species")


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


##----- Set up input ------------------------------------------------------------------------------

classifier <- args[2]
inputDir = file.path(pipeRoot, "input")
taxonomyFile <- "Full_Taxonomy.tsv"

prevModule <- paste0("Taxa_over_time_MLM")
inputDir2 = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), prevModule),"/output/")

##----- Set up output -----------------------------------------------------------------------------
outputDir = file.path(moduleDir,"output/")

##----- Script variables --------------------------------------------------------------------------
month.labs <- c("Baseline", "1 month post-op", "6 months post-op", "12 months post-op",
                "18 months post-op", "24 months post-op")
endM <- args[length(args)]
comparisons <- c("1_v_0", "6_v_0", "12_v_0", "18_v_0", "24_v_0",
                 "6_v_1", "12_v_1", "18_v_1", "24_v_1",
                 "12_v_6", "18_v_6", "24_v_6",
                 "18_v_12", "24_v_12",
                 "24_v_18")
##----- Generate heatmaps (unadjusted p-values) ----------------------------------------------------

taxonomyTable <- read.delim(file.path(inputDir, "microbiome", taxonomyFile), sep="\t", header = TRUE)

for (level in levels) {
  ##---- Read in data -----
  # outputLevel <- paste0(outputDir, level, "/")
  # dir.create(outputLevel, showWarnings = FALSE)
  
  inputLevel <- paste0(inputDir2, level, "/")
  file.path <- paste0(inputLevel, level, "~Timepoint+Surgery+Weight_kg_MixedLinearModelResults_Heatmap.tsv")
  myT <- read.table(file.path, sep="\t",header = TRUE, check.names = FALSE)
  ##---- end -----
  
  ##---- Set up for heatmap -----
  taxaNames <- myT$bugName
  taxaNames <-  gsub("_", " ", taxaNames)
  if (level == "genus") {
    taxaNames <-  gsub("noname", "(unnamed genus)", taxaNames)
    taxaNames <-  gsub("unclassified", "(unclassified genus)", taxaNames)
  }
  
  adjp_List <- list()
  pval_List <- list()
  index <- 1
  for (comparison in comparisons) {
    
    t1 <- sapply(strsplit(comparison, "_v_"), "[", 2)
    t1 <- ifelse(t1 == "0", "BL", paste0(t1, "M"))
    t2 <- sapply(strsplit(comparison, "_v_"), "[", 1)
    t2 <- ifelse(t2 == "0", "BL", paste0(t2, "M"))
    
    adjp <- myT[,which(colnames(myT) == paste0("Adj_p_Timepoint", comparison))]
    p <- myT[,which(colnames(myT) == paste0("p_Timepoint", comparison))]
    s <- myT[,which(colnames(myT) == paste0("s_Timepoint", comparison))]
    log10p <- getlog10p(p, s)
    
    pval_List[[index]] <- log10p
    adjp_List[[index]] <- adjp
    names(pval_List)[index] <- paste0(t1, "-", t2)
    names(adjp_List)[index] <- paste0(t1, "-", t2)
    index <- index + 1
  } # for (comparison in comparisons)
  
  df <- do.call(rbind, pval_List)
  df.adj <- do.call(rbind, adjp_List)
  colnames(df)<- taxaNames
  colnames(df.adj)<- taxaNames
  ##---- end -----
  
  
  ##---- Heatmap -----
  df<-t(df)
  df.adj<-t(df.adj)
  
  # Specify the range
  lower_limit <- log10(0.05)
  upper_limit <- -log10(0.05)
  
  if (level == "genus") {

    # Find rows where all values are between the specified range
    df <- df[!rowSums(df.adj >= 0.05) == ncol(df.adj), ]
    df.adj <- df.adj[!rowSums(df.adj >= 0.05) == ncol(df.adj), ]
  }
  
  times<-sapply(colnames(df), function(x){strsplit(x,"-")[[1]][1]})
  times <- as.factor(times)
  times <- factor(times, levels = c("BL", "1M", "6M", "12M", "18M"))
  
  col = list(Time = c("BL" = "lightpink", "1M" = "yellow","6M" = "blue","12M" = "grey","18M" = "purple"))
  
  ha <- HeatmapAnnotation(
    Time = times, col = col
  )
  
  if (level == "phylum") {
    phylumplot <- Heatmap(df,
                          top_annotation = ha,
                          rect_gp = gpar(col = "white", lwd = 1),
                          row_names_gp = gpar(fontsize = 10),
                          heatmap_height = unit(10, "cm"),
                          name = "log10 p-value",
                          layer_fun = function(j, i, x, y, w, h, fill) {
                            ind_mat = restore_matrix(j, i, x, y)
                            for(ir in seq_len(nrow(ind_mat))) {
                              for(ic in seq_len(ncol(ind_mat))) {
                                ind1 = ind_mat[ir, ic] # previous column
                                v1 = df.adj[i[ind1], j[ind1]]
                                # print(v1)
                                if(v1 < 0.05) { # if they have the same sign
                                  grid.points(x[c(ind1)], y[c(ind1)],
                                              pch = 8, size = unit(2, "mm"))
                                }
                              } # for(ic in seq_len(ncol(ind_mat)))
                            } # for(ir in seq_len(nrow(ind_mat)))
                          }
    )
  }
  if (level == "genus") {
    genusplot <- Heatmap(df,
                          top_annotation = ha,
                          rect_gp = gpar(col = "white", lwd = 1),
                          row_names_gp = gpar(fontsize = 6),
                          heatmap_height = unit(20, "cm"),
                          name = "log10 p-value",
                          layer_fun = function(j, i, x, y, w, h, fill) {
                            ind_mat = restore_matrix(j, i, x, y)
                            for(ir in seq_len(nrow(ind_mat))) {
                              for(ic in seq_len(ncol(ind_mat))) {
                                ind1 = ind_mat[ir, ic] # previous column
                                v1 = df.adj[i[ind1], j[ind1]]
                                # print(df.adj[53,1])
                                if(v1 < 0.05) { # if they have the same sign
                                  grid.points(x[c(ind1)], y[c(ind1)],
                                              pch = 8, size = unit(0.5, "mm"))
                                }
                              } # for(ic in seq_len(ncol(ind_mat)))
                            } # for(ir in seq_len(nrow(ind_mat)))
                          }
    )
  }
  
  ##---- end -----
  
} # for (level in levels)

pdf(file.path(outputDir, paste0("phylum_HeatmapAndCluster.pdf")),height = 8)
print(phylumplot)
dev.off()

pdf(file.path(outputDir, paste0("genus_HeatmapAndCluster.pdf")),height = 8)
print(genusplot)
dev.off()



##----- classified genera -----
level = "genus"
inputLevel <- paste0(inputDir2, level, "/")
file.path <- paste0(inputLevel, level, "~Timepoint+Surgery+Weight_kg_MixedLinearModelResults_Heatmap_Classified.tsv")
myT <- read.table(file.path, sep="\t",header = TRUE, check.names = FALSE)
##---- end -----

##---- Set up for heatmap -----
taxaNames <- myT$bugName
taxaNames <-  gsub("_", " ", taxaNames)

adjp_List <- list()
pval_List <- list()
index <- 1
for (comparison in comparisons) {
  
  t1 <- sapply(strsplit(comparison, "_v_"), "[", 2)
  t1 <- ifelse(t1 == "0", "BL", paste0(t1, "M"))
  t2 <- sapply(strsplit(comparison, "_v_"), "[", 1)
  t2 <- ifelse(t2 == "0", "BL", paste0(t2, "M"))
  
  adjp <- myT[,which(colnames(myT) == paste0("Adj_p_Timepoint", comparison))]
  p <- myT[,which(colnames(myT) == paste0("p_Timepoint", comparison))]
  s <- myT[,which(colnames(myT) == paste0("s_Timepoint", comparison))]
  log10p <- getlog10p(p, s)
  
  pval_List[[index]] <- log10p
  adjp_List[[index]] <- adjp
  names(pval_List)[index] <- paste0(t1, "-", t2)
  names(adjp_List)[index] <- paste0(t1, "-", t2)
  index <- index + 1
} # for (comparison in comparisons)

df <- do.call(rbind, pval_List)
df.adj <- do.call(rbind, adjp_List)
colnames(df)<- taxaNames
colnames(df.adj)<- taxaNames
##---- end -----


##---- Heatmap -----
df<-t(df)
df.adj<-t(df.adj)

# Specify the range
lower_limit <- log10(0.05)
upper_limit <- -log10(0.05)

if (level == "genus") {
  
  # Find rows where all values are between the specified range
  df <- df[!rowSums(df.adj >= 0.05) == ncol(df.adj), ]
  df.adj <- df.adj[!rowSums(df.adj >= 0.05) == ncol(df.adj), ]
}

times<-sapply(colnames(df), function(x){strsplit(x,"-")[[1]][1]})
times <- as.factor(times)
times <- factor(times, levels = c("BL", "1M", "6M", "12M", "18M"))

col = list(Time = c("BL" = "lightpink", "1M" = "yellow","6M" = "blue","12M" = "grey","18M" = "purple"))

genera <- rownames(df)
taxonomyTable <- taxonomyTable[,c("Phylum", "Genus")]
taxonomyTable <- taxonomyTable[!duplicated(taxonomyTable$Genus),]
taxonomyTable <- taxonomyTable[taxonomyTable$Genus %in% genera,]
taxonomyTable <- taxonomyTable[order(taxonomyTable$Genus), ]

# Keep only genera present in taxonomy table; align row order to match
df     <- df[rownames(df) %in% taxonomyTable$Genus, , drop = FALSE]
df     <- df[order(rownames(df)), , drop = FALSE]
df.adj <- df.adj[rownames(df.adj) %in% taxonomyTable$Genus, , drop = FALSE]
df.adj <- df.adj[order(rownames(df.adj)), , drop = FALSE]

phyla <- taxonomyTable$Phylum

col2 = list(Phylum = c("Actinobacteria" = "yellow2",
                    "Ascomycota" = "chocolate",
                    "Bacteroidetes" = "limegreen",
                    "Cyanobacteria" = "darkgreen",
                    "Deinococcus_Thermus" = "salmon",
                    "Euryarchaeota" = "gray60",
                    "Firmicutes" = "skyblue1",
                    "Fusobacteria" = "orange",
                    "Proteobacteria" = "purple1",
                    "Synergistetes" = "darkorange3",
                    "Tenericutes" = "mediumpurple",
                    "Verrucomicrobia" = "hotpink"))

top_ha <- HeatmapAnnotation(
  Time = times, col = col
)

right_ha <- rowAnnotation(
  Phylum = phyla, col = col2
)

pdf(file.path(outputDir, paste0("genus_HeatmapAndCluster_Classified.pdf")),height = 8)
Heatmap(df,
        top_annotation = ha,
        left_annotation = right_ha,
        rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize = 6),
        heatmap_height = unit(20, "cm"),
        name = "log10 p-value",
        layer_fun = function(j, i, x, y, w, h, fill) {
          ind_mat = restore_matrix(j, i, x, y)
          for(ir in seq_len(nrow(ind_mat))) {
            for(ic in seq_len(ncol(ind_mat))) {
              ind1 = ind_mat[ir, ic] # previous column
              v1 = df.adj[i[ind1], j[ind1]]
              # print(df.adj[53,1])
              if(v1 < 0.05) { # if they have the same sign
                grid.points(x[c(ind1)], y[c(ind1)],
                            pch = 8, size = unit(0.5, "mm"))
              }
            } # for(ic in seq_len(ncol(ind_mat)))
          } # for(ir in seq_len(nrow(ind_mat)))
        }
)
dev.off()
##---- end -----
