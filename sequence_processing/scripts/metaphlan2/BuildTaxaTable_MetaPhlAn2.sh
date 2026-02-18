#!/bin/bash

# Shell script to automate MetaPhlAn2 
# Author = Alicia Sorgen
# Date = 2021 Oct 11

################## Edit ##############################
Dir=/nobackup/afodor_research/asorgen/BariatricSurgery
inputDir=${Dir}/04_metaphlan2
outputDir=${Dir}/07_metaphlan2TaxaTables
scriptPath=${Dir}/scripts/metaphlan2
funcScript=${Dir}/scripts/functions.R
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
${scriptPath}/BuildTaxaTable_MetaPhlAn2.slurm 

