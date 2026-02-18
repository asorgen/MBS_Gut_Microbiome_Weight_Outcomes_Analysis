#!/bin/bash

# Shell script to automate Kraken2 
# Author = Alicia Sorgen
# Date = 2021 Oct 11

################## Edit ##############################
Dir=/nobackup/afodor_research/asorgen/BariatricSurgery
inputDir=${Dir}/01_kneaddata_hg19
outputDir=${Dir}/03_kraken2
# inputList=${Dir}/BariatricSurgery_FileList.txt
inputList=${Dir}/BariatricSurgery_FileList_2022Nov18.txt
R1File=R1_001_kneaddata_paired_1.fastq
R2File=R1_001_kneaddata_paired_2.fastq
scriptPath=${Dir}/scripts/kraken2
database=/nobackup/afodor_research/databases/kraken2_standardDB
kraken2Path=/users/asorgen/PROGRAMS/kraken2
######################################################

# rm -r ${Dir}/kraken2 # remove old directory
# mkdir ${outputDir} # create new directory


for ID in $(cat $inputList); do
    
    # Assign variable names
    export ID
    export database
    export kraken2Path
    export INPUT=${inputDir}/${ID}
    export R1=${INPUT}/${ID}_${R1File}
    export R2=${INPUT}/${ID}_${R2File}
    export report=${outputDir}/${ID}_k2_REPORT.txt
    export K2out=${outputDir}/${ID}_k2_OUTPUT.txt

    # Submit job
    sbatch \
    --job-name=k2_${ID} \
    --time=10:00:00 \
    ${scriptPath}/kraken2.slurm 

done