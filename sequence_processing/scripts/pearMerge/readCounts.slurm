#!/bin/bash

#SBATCH --partition=Orion											    # Partition/queue requested on server
#SBATCH --time=05:00:00												    # Time limit (hrs:min:sec)
#SBATCH --nodes=1													    # Number of nodes requested
#SBATCH --ntasks-per-node=16										    # Number of CPUs (processor cores/tasks)
#SBATCH --mem=125gb													    # Memory limit
#SBATCH --job-name=readCounts										    # Job name

source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
logDir=${Dir}/outputLogs
jobName=$SLURM_JOB_NAME
outFile=$SLURM_JOB_NAME.$SLURM_JOB_ID

echo Output file: $outFile
echo time=5:00:00

mv ${logDir}/${jobName}.* ${logDir}/oldLogs
mv ${logDir}/oldLogs/${outFile}.* ${logDir}

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

### Load modules
# module load metaphlan/3.0.7

/bin/date

SECONDS=0

${ScriptDir}/pearMerge/readCounts.sh


duration=$SECONDS

hour=`echo "$(( ${duration} / 3600))"`
minute=`echo "$(( (${duration} / 60) - (hour * 60) ))"`
second=`echo "$((${duration} % 60))"`


echo "Job time: ${hour} hour(s) $minute minute(s) and ${second} seconds"
