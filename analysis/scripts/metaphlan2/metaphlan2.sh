#!/bin/bash

# Shell script to automate MetaPhlAn2 
# Author = Alicia Sorgen
# Date = 2021 Oct 11

################## Edit ##############################
Dir=/nobackup/afodor_research/asorgen/BariatricSurgery
inputDir=${Dir}/01_kneaddata_hg19
outputDir=${Dir}/04_metaphlan2
# inputList=${Dir}/BariatricSurgery_FileList.txt
inputList=${Dir}/BariatricSurgery_FileList_2022Nov18.txt
R1File=R1_001_kneaddata_paired_1.fastq
R2File=R1_001_kneaddata_paired_2.fastq
scriptPath=${Dir}/scripts/metaphlan2
database=/nobackup/afodor_research/databases/metaphlan2_db/mpa_v20_m200
######################################################

# rm -r ${outputDir} # remove old directory
# mkdir ${outputDir} # create new directory


for ID in $(cat $inputList); do
    
    # Assign variable names
    export ID
    export database
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

