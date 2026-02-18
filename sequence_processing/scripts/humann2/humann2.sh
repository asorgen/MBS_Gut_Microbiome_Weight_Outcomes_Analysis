#!/bin/bash

# Shell script to automate HumanN2 
# Author = Alicia Sorgen
# Date = 2021 Oct 11

################## Edit ##############################
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
inputDir=${Dir}/02_pearMerge
outputDir=${Dir}/05_humann2
scriptPath=${ScriptDir}/humann2
######################################################

# mkdir ${outputDir} # create new directory

# mkdir ${outputDir}/pathCoverage
# mkdir ${outputDir}/pathAbundance
# mkdir ${outputDir}/geneFamilies

ID=BIO-2-010-01_CAACACCT-TGGTAGCT_S072_L006

# for ID in $(cat $inputList); do
    
    # Assign variable names
    export ID
    export INPUT=${inputDir}/${ID}_merge.assembled.fastq
    export outputDir

    # Submit job
    sbatch \
    --job-name=humann_${ID} \
    --time=10:00:00 \
    ${scriptPath}/humann2.slurm 

# done

