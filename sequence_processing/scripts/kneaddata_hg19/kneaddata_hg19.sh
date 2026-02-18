#!/bin/bash

# Shell script to automate Kneaddata 
# Author = Alicia Sorgen
# Date = 2021 Oct 7

Dir=/nobackup/afodor_research/asorgen/BariatricSurgery/
inputDir=/scratch/asorgen/metagenome/BariatricSurgery/rawSeqs
scriptPath=${Dir}/scripts/kneaddata_hg19

# rm -r ${Dir}/kneaddata_hg19 # remove old directory
mkdir ${Dir}/01_kneaddata_hg19 # create new directory

outputDir=${Dir}/01_kneaddata_hg19
# inputList=${Dir}/BariatricSurgery_FileList.txt
# inputList=${Dir}/BariatricSurgery_FileList_2021Dec08.txt
inputList=${Dir}/BariatricSurgery_FileList_2022Nov18.txt

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