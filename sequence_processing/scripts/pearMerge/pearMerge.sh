#!/bin/bash

# Shell script to automate PEAR merge
# Author = Alicia Sorgen
# Date = 2021 Oct 8

########### Edit directory paths ###########
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
inputDir=${Dir}/01_kneaddata_hg19
outputDir=${Dir}/02_pearMerge
R1File=R1_001_kneaddata_paired_1.fastq
R2File=R1_001_kneaddata_paired_2.fastq
scriptPath=${ScriptDir}/pearMerge
############################################


# rm -r ${Dir}/pearMerge # remove old directory
# mkdir ${outputDir} # create new directory


for ID in $(cat $inputList); do
    
    # Assign variable names
    export ID
    export INPUT=${inputDir}/${ID}
    export R1=${INPUT}/${ID}_${R1File}
    export R2=${INPUT}/${ID}_${R2File}
    export merge=${outputDir}/${ID}_merge

    # Submit job
    sbatch \
    --job-name=pear_${ID} \
    --output=${slurm_out}/pearMerge/%x.%j.out \
    --error=${slurm_out}/pearMerge/%x.%j.err \
    ${scriptPath}/pearMerge.slurm 

done


