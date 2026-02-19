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

# for ID in $(cat $inputList); do
    
    # Assign variable names
    export ID
    export INPUT=${inputDir}/${ID}_merge.assembled.fastq
    export outputDir

    # Submit job
    sbatch \
    --job-name=humann_${ID} \
    --time=10:00:00 \
    --output=${slurm_out}/humann2/%x.%j.out \
    --error=${slurm_out}/humann2/%x.%j.err \
    ${scriptPath}/humann2.slurm

# done

