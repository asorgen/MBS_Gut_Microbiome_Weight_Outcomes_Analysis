#!/bin/bash

# Shell script to merge MetaPhlAn2 data 
# Author = Alicia Sorgen
# Date = 2021 Oct 12

################# EDIT ####################
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
outputDir=${Dir}/04_metaphlan2
# inputList=${Dir}/BariatricSurgery_FileList.txt    # uses original file list, not 2022 update
###########################################

# mkdir ${outputDir}/bowtie2Output
mv ${outputDir}/*.bowtie2.bz2 ${outputDir}/bowtie2Output

# mkdir ${outputDir}/profileOutput
mv ${outputDir}/*_profile.txt ${outputDir}/profileOutput

# module load metaphlan/3.0.7
module load metaphlan/2.8.1

# rm ${outputDir}/merged_abundance_table.txt
merge_metaphlan_tables.py ${outputDir}/profileOutput/*_profile.txt > ${outputDir}/merged_abundance_table.txt