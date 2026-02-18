#!/bin/bash

# Shell script to automate MetaPhlAn2 
# Author = Alicia Sorgen
# Date = 2021 Oct 11

################## Edit ##############################
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
inputDir=${Dir}/01_kneaddata_hg19
outputDir=${Dir}/04_metaphlan2
R1File=R1_001_kneaddata_paired_1.fastq
R2File=R1_001_kneaddata_paired_2.fastq
scriptPath=${ScriptDir}/metaphlan2
######################################################

# rm -r ${outputDir} # remove old directory
# mkdir ${outputDir} # create new directory


for ID in $(cat $inputList); do
    
    # Assign variable names
    export ID
    export database=${metaphlan2DB}
    export INPUT=${inputDir}/${ID}
    export R1=${INPUT}/${ID}_${R1File}
    export R2=${INPUT}/${ID}_${R2File}
    export bowtieOut=${outputDir}/${ID}.bowtie2.bz2
    export profileOut=${outputDir}/${ID}_profile.txt

    # Submit job
    sbatch \
    --job-name=metaphlan_${ID} \
    --time=5:00:00 \
    ${scriptPath}/metaphlan2.slurm 

done

