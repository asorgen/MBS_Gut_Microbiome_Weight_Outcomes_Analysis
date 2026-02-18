#!/bin/bash

#SBATCH --partition=Orion											    # Partition/queue requested on server
#SBATCH --time=04:00:00												    # Time limit (hrs:min:sec)
#SBATCH --nodes=1													    # Number of nodes requested
#SBATCH --ntasks-per-node=16										    # Number of CPUs (processor cores/tasks)
#SBATCH --mem=125gb													    # Memory limit

source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
logDir=${Dir}/outputLogs
jobName=$SLURM_JOB_NAME
outFile=$SLURM_JOB_NAME.$SLURM_JOB_ID

echo Output file: $outFile

mv ${logDir}/${jobName}.* ${logDir}/oldLogs
mv ${logDir}/oldLogs/${outFile}.* ${logDir}

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

### Load modules
# module load anaconda3/2020.11
module load kneaddata/0.7.2 

echo -e "\n${yellow}[Running Kneaddata for sample ${ID}]\n"
/bin/date
SECONDS=0


kneaddata \
--input1 $R1 \
--input2 $R2 \
--reference-db ${kneaddataDB} \
-t 4 \
-p 2 \
--remove-intermediate-output \
--output $PEout \
--trimmomatic-options ILLUMINACLIP:${trimmomaticAdapters}:2:30:10:8:true \
--trimmomatic-options LEADING:3 \
--trimmomatic-options TRAILING:3 \
--trimmomatic-options SLIDINGWINDOW:4:15 \
--trimmomatic-options MINLEN:50

# newlogDir=${logDir}/knead_hg19_2021Dec08
# mkdir ${newlogDir}/
# mv ${outFile}.* ${newlogDir}

duration=$SECONDS

hour=`echo "$(( ${duration} / 3600))"`
minute=`echo "$(( (${duration} / 60) - (hour * 60) ))"`
second=`echo "$((${duration} % 60))"`


echo "Job time: ${hour} hour(s) $minute minute(s) and ${second} seconds"
