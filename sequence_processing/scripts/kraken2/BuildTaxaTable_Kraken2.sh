#!/bin/bash

# Shell script to automate Kraken2 
# Author = Alicia Sorgen
# Date = 2021 Oct 11

################## Edit ##############################
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
inputDir=${Dir}/03_kraken2
outputDir=${Dir}/06_kraken2TaxaTables
scriptPath=${ScriptDir}/kraken2
funcScript=${ScriptDir}/functions.R
logDir=${Dir}/outputLogs
######################################################

export inputDir
export outputDir
export funcScript
export scriptPath
export logDir

# Submit job
sbatch \
--job-name=k2Tables \
${scriptPath}/BuildTaxaTable_Kraken2.slurm 

