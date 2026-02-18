#!/bin/bash

# Shell script to automate PEAR merge
# Author = Alicia Sorgen
# Date = 2021 Oct 8

########### Edit directory paths ###########
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
inputDir=${Dir}/02_pearMerge
# inputList=${Dir}/BariatricSurgery_FileList.txt    # uses original file list, not 2022 update
scriptPath=${ScriptDir}/pearMerge
############################################


echo -e "SampleID" > ${Dir}/sampleid.tsv
echo -e "ReadCount" > ${Dir}/readcount.tsv


for ID in $(cat $inputList); do
    
    num=`echo $(cat ${inputDir}/${ID}_merge.assembled.fastq|wc -l)/4|bc`
    echo -e ${ID} >> ${Dir}/sampleid.tsv
    echo -e ${num} >> ${Dir}/readcount.tsv

done

paste ${Dir}/sampleid.tsv ${Dir}/readcount.tsv > ${Dir}/seqCounts.tsv

rm ${Dir}/sampleid.tsv
rm ${Dir}/readcount.tsv

