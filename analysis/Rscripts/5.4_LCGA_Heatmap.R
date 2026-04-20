# Author: Alicia Sorgen
# Description: Heatmap of taxa trajectories by LCGA group (RYGB, 2-class, Model2)

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

library(stringr);       message("stringr: Version ",       packageVersion("stringr"))
library(data.table);    message("data.table: Version ",    packageVersion("data.table"))
library(tidyr);         message("tidyr: Version ",         packageVersion("tidyr"))
library(ComplexHeatmap); message("ComplexHeatmap: Version ", packageVersion("ComplexHeatmap"))

# Script-specific edits -------------------------------------------------
params <- vector()
params <- c(params, getOption("mbs.pipe_root", default = stop("Set options(mbs.pipe_root = ...) in your local .Rprofile (gitignored)")))

# Set working environment ---------------------------------------------------------
H1("Working Environment")
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
   args <- params
}

message("\n************* Running locally *************")
proj_root <- args[1]
message("Project root directory: ", proj_root, "\n")
inputEnv  <- Sys.getenv("INPUT_ROOT");   input_root <- if (nchar(inputEnv)  > 0) inputEnv  else file.path(dirname(proj_root), "Data")
script_root <- file.path(proj_root, basename(proj_root), "analysis", "Rscripts")

resultsEnv <- Sys.getenv("RESULTS_ROOT"); pipeRoot <- if (nchar(resultsEnv) > 0) resultsEnv else file.path(proj_root, gsub("Analysis", "Results", basename(proj_root)))

module    <- "5.4_LCGA_Heatmap"
moduleDir <- file.path(pipeRoot, module)

outputDir <- file.path(moduleDir, "output/")

rm(proj_root, inputEnv, resultsEnv, module, params)

# Set functions ------------------------------------------------
H1("Functions")
message("Loading functions from: ", script_root)
funcScript <- file.path(script_root, "functions.R")
source(funcScript); rm(funcScript)

# Input ----------------------------------------------------------------
H1("Input")
taxonomyFile  <- "Full_Taxonomy.tsv"
taxonomyTable <- read.delim(file.path(input_root, "microbiome", taxonomyFile), sep = "\t", header = TRUE)
message(">>> taxa table")

prevModule  <- str_subset(dir(pipeRoot), "5.2_LCGA2_RYGB_2Class")
LCGADir     <- file.path(pipeRoot, prevModule, "output")
message("LCGA module dir: ", LCGADir)
pvalue.File <- list.files(file.path(LCGADir, "08_LinearModel_Analysis_pvalues", "ResultsTables"), full.names = TRUE)
message("p-value file: ", pvalue.File)

# Output ----------------------------------------------------------------
H1("Output")

# Script variables ----------------------------------------------------------------
H1("Variables")
level <- "genus"

new_getlog10p <- function(pval, coefficient) {
   sapply(seq_along(pval), function(x) {
      if (coefficient[x] == "increase") -log10(pval[x]) else log10(pval[x])
   })
}

# Build heatmap matrices -------------------------------------------------------
H1("Build heatmap matrices")

myT <- read.table(pvalue.File, sep = "\t", header = TRUE, check.names = FALSE)
message(">>> p-value table")

myT$groupPValues4 <- ifelse(myT$groupPValues4 == 0, 2e-16, myT$groupPValues4)

log10p      <- new_getlog10p(myT$groupPValues4, myT$direction)
comparisons <- paste0(myT$LengthTime)
taxaNames   <- paste0(myT$qbugName, "&", myT$groupNumbers)
adjp        <- myT$adjPValues4

df      <- spread(data.frame(taxaNames, comparisons, log10p), comparisons, log10p)
adjP.df <- spread(data.frame(taxaNames, comparisons, adjp),   comparisons, adjp)

df      <- df[!(df$taxaNames      %like% "Other"), ]
adjP.df <- adjP.df[!(adjP.df$taxaNames %like% "Other"), ]

taxaNames <- df$taxaNames
df      <- data.frame(taxaNames, BLto1M = df$BLto1M, BLto6M = df$BLto6M,
                      BLto12M = df$BLto12M, BLto18M = df$BLto18M, BLto24M = df$BLto24M)
adjP.df <- data.frame(taxaNames, BLto1M = adjP.df$BLto1M, BLto6M = adjP.df$BLto6M,
                      BLto12M = adjP.df$BLto12M, BLto18M = adjP.df$BLto18M, BLto24M = adjP.df$BLto24M)

rownames(df)      <- df$taxaNames;      df      <- df[, -1]
rownames(adjP.df) <- adjP.df$taxaNames; adjP.df <- adjP.df[, -1]

df      <- df[order(rownames(df)), ]
adjP.df <- adjP.df[order(rownames(adjP.df)), ]

if (level == "genus") {
   keep    <- !rowSums(adjP.df >= 0.05, na.rm = TRUE) == ncol(adjP.df)
   df      <- df[keep, ]
   adjP.df <- adjP.df[keep, ]
}

df.1 <- df[which(rownames(df) %like% "&1"), ]
df.2 <- df[which(rownames(df) %like% "&2"), ]
group1_Taxa <- nrow(df.1)
group2_Taxa <- nrow(df.2)
df <- as.matrix(rbind(df.1, df.2))

adjP.df.1 <- adjP.df[which(rownames(adjP.df) %like% "&1"), ]
adjP.df.2 <- adjP.df[which(rownames(adjP.df) %like% "&2"), ]
adjP.df   <- as.matrix(rbind(adjP.df.1, adjP.df.2))

# Phylum annotations -------------------------------------------------------
H1("Phylum annotations")

col2 <- list(Phylum = c(
   "Actinobacteria"  = "yellow2",
   "Bacteroidetes"   = "limegreen",
   "Euryarchaeota"   = "coral",
   "Firmicutes"      = "skyblue1",
   "Fusobacteria"    = "orange",
   "Proteobacteria"  = "purple1",
   "Verrucomicrobia" = "hotpink"
))

get_phyla <- function(mat, taxonomy) {
   genera <- sapply(strsplit(rownames(mat), "&"), "[", 1)
   t2 <- taxonomy[, c("Phylum", "Genus")]
   t2 <- t2[!duplicated(t2$Genus), ]
   t2$Phylum[match(genera, t2$Genus)]
}

phyla <- c(get_phyla(df.1, taxonomyTable), get_phyla(df.2, taxonomyTable))

rownames(df) <- sapply(strsplit(rownames(df), "&"), "[", 1)
adjP.df      <- replace(adjP.df, is.na(adjP.df), 1)

ra <- rowAnnotation(
   empty = anno_empty(border = FALSE),
   foo2  = anno_block(gp = gpar(fill = 2:3), labels = c("Group 1", "Group 2"))
)
right_ha <- rowAnnotation(Phylum = phyla, col = col2)
splitRow  <- c(rep(1, group1_Taxa), rep(2, group2_Taxa))

# Heatmap -------------------------------------------------------
H1("Heatmap")

pdf(file.path(outputDir, "genus_HeatmapAndCluster_Classified.pdf"), height = 8)
Heatmap(df,
        left_annotation  = ra,
        right_annotation = right_ha,
        rect_gp          = gpar(col = "white", lwd = 1),
        name             = "log10 p-value",
        row_split        = splitRow,
        column_title     = NULL,
        row_title        = NULL,
        border           = TRUE,
        cluster_columns  = FALSE,
        cluster_rows     = FALSE,
        row_names_gp     = gpar(fontsize = 6),
        layer_fun = function(j, i, x, y, w, h, fill) {
           ind_mat <- restore_matrix(j, i, x, y)
           for (ir in seq_len(nrow(ind_mat))) {
              for (ic in seq_len(ncol(ind_mat))) {
                 ind1 <- ind_mat[ir, ic]
                 if (adjP.df[i[ind1], j[ind1]] < 0.05) {
                    grid.points(x[ind1], y[ind1], pch = 8,
                                size = unit(0.5, "mm"), gp = gpar(col = "black"))
                 }
              }
           }
        }
)
dev.off()

message("Output written to: ", outputDir)
