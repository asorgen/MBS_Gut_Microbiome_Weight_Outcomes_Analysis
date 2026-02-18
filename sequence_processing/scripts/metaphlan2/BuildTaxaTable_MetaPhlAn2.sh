#!/bin/bash

# Shell script to automate MetaPhlAn2 
# Author = Alicia Sorgen
# Date = 2021 Oct 11

################## Edit ##############################
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
inputDir=${Dir}/04_metaphlan2
outputDir=${Dir}/07_metaphlan2TaxaTables
scriptPath=${ScriptDir}/metaphlan2
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
--job-name=m2test \
--time=02:00:00 \
--output=${slurm_out}/metaphlan2/%x.%j.out \
--error=${slurm_out}/metaphlan2/%x.%j.err \
${scriptPath}/BuildTaxaTable_MetaPhlAn2.slurm

