#!/bin/bash

# Shell script to automate Kneaddata 
# Author = Alicia Sorgen
# Date = 2021 Oct 7

source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
inputDir=${rawSeqDir}
scriptPath=${ScriptDir}/kneaddata_hg19

# rm -r ${Dir}/kneaddata_hg19 # remove old directory
mkdir ${Dir}/01_kneaddata_hg19 # create new directory

outputDir=${Dir}/01_kneaddata_hg19

for ID in $(cat $inputList); do
    
    # Assign variable names
    export ID
    export INPUT=${inputDir}
    export R1=${INPUT}/${ID}_R1_001.fastq.gz
    export R2=${INPUT}/${ID}_R2_001.fastq.gz
    export PEout=${outputDir}/${ID}

    # Submit job
    sbatch \
    --job-name=hg19_${ID} \
    ${scriptPath}/kneaddata_hg19.slurm 

done