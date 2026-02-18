#!/bin/bash

# Shell script to automate PEAR merge
# Author = Alicia Sorgen
# Date = 2021 Oct 8

########### Edit directory paths ###########
Dir=/nobackup/afodor_research/asorgen/BariatricSurgery
inputDir=${Dir}/01_kneaddata_hg19
# inputDir=/scratch/asorgen/metagenome/BariatricSurgery/01_kneaddata_hg19
# inputList=${Dir}/BariatricSurgery_FileList.txt
inputList=${Dir}/BariatricSurgery_FileList_2022Nov18.txt
outputDir=${Dir}/02_pearMerge
#outputDir=/scratch/asorgen/metagenome/BariatricSurgery/02_pearMerge
R1File=R1_001_kneaddata_paired_1.fastq
R2File=R1_001_kneaddata_paired_2.fastq
scriptPath=${Dir}/scripts/pearMerge
############################################


# rm -r ${Dir}/pearMerge # remove old directory
# mkdir ${outputDir} # create new directory
# ID=BIO-2-031-12_ATCGTGGA-GGCCGTTG_S022_L002


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
    ${scriptPath}/pearMerge.slurm 

done


