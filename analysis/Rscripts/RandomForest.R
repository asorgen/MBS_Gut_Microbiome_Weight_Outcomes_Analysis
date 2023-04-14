#Author: Alicia Sorgen
#Date: 08-16-21
#Description: Merge taxa count tables with metadata

## Libraries
library("tidyverse")
library("SIAMCAT")
library("yaml")
library("stringr")

rm(list=ls())

date <- Sys.Date()
date <- format(date, "%Y%b%d")
date <- "2021Aug31"

# root <- paste0("~/BioLockJ_pipelines/actigraph_analysis_", date)
root <- paste0("~/BioLockJ_pipelines/actigraph_update_", date)
root <- dir(root, pattern=paste0("RandomForest"), full.names=TRUE)

if ((dirname(dirname(dirname(getwd()))) == "/mnt/efs/pipelines")) {
  message("************* Running in BioLockJ *************")
  # studies <- commandArgs(trailingOnly = TRUE)
} else {
  setwd(paste0(root, "/script"))
}


##### Prep #####
pipeRoot = dirname(dirname(getwd()))
inputDir = paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "TaxaMetaRelAb"),"/output/")
moduleDir <- dirname(getwd())
output = file.path(moduleDir,"output/")

funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

inputFile <- "_metaCounts_relabundance.tsv"

# levels <- c("phylum", "class", "order", "family", "genus", "species")
levels <- c("genus", "species")

for (level in levels) {
  
  #Loading tables
  # level <- "phylum"
  metaCounts <- read.table(paste0(inputDir, level, inputFile), header = TRUE, quote='')
  rownames(metaCounts) <- metaCounts$SampleID
  taxaStart <- which(colnames(metaCounts) == "OpTime")+1
  counts <- metaCounts[,taxaStart:ncol(metaCounts)]
  feat.all <- t(counts)
  meta <- metaCounts[,1:taxaStart-1]
  stopifnot(all(meta$SampleID %in% colnames(feat.all)))

  #Set parameters
  norm.method="log.std"
  n.p<-list(log.n0=1e-05,sd.min.q=0.1,n.p=2,norm.margin=1)
  num.folds=8
  num.resample=10
  ml.method="randomForest"
  modsel.crit=list("pr")
  min.nonzero.coeff=1
  param.fs.ss=list(thres.fs=800,method.fs="AUC")
  
  #Model building
  feat.train <- feat.all[,as.character(meta %>% pull(SampleID))]
  
  siamcat <- siamcat(feat=feat.train, meta=meta,
                     label = 'Assessment', case=1)
  siamcat <- normalize.features(siamcat, norm.method = norm.method,
                                norm.param = n.p, feature.type = 'original',
                                verbose=3)
  siamcat <- create.data.split(siamcat, num.folds = num.folds,
                               num.resample = num.resample,stratify=TRUE)
  siamcat <- train.model(siamcat,
                         method = ml.method,
                         modsel.crit=modsel.crit,
                         min.nonzero.coeff = min.nonzero.coeff,
                         perform.fs = FALSE)
  siamcat <- make.predictions(siamcat)
  siamcat <- evaluate.predictions(siamcat)
  save(siamcat, file=paste0(output, level, "_", ml.method,'_model.RData'))
  
  write.table(siamcat@pred_matrix,paste0(output, level, "_", ml.method,'_model.tsv'),sep="\t")
  
  
}