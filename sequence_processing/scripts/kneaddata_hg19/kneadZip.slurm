#!/bin/bash

#SBATCH --job-name=KneadData_zip 																# Job name
#SBATCH --partition=Orion 																		# Partition/queue requested on server
#SBATCH --time=100:00:00																		# Time limit (hrs:min:sec)
#SBATCH --nodes=1 																				# Number of nodes requested
#SBATCH --ntasks-per-node=16 																	# Number of CPUs (processor cores/tasks)
#SBATCH --mem=125gb																				# Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL 																# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=${email}

SECONDS=0

source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
logDir=${Dir}/outputLogs
jobName=$SLURM_JOB_NAME
outFile=$SLURM_JOB_NAME.$SLURM_JOB_ID

echo Output file: $outFile
echo time=100:00:00

mv ${logDir}/${jobName}.* ${logDir}/oldLogs
mv ${logDir}/oldLogs/${outFile}.* ${logDir}

### Display the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host: `hostname`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

### Load modules

# zip -r ${Dir}/01_kneaddata_hg19.zip ${Dir}/01_kneaddata_hg19/
unzip ${Dir}/01_kneaddata_hg19.zip -d ${Dir}_GMHI/BS_data -x "*/*_unmatched_*" "*/*_contam_*" "*/*.log"

duration=$SECONDS

hour=`echo "$(( ${duration} / 3600))"`
minute=`echo "$(( (${duration} / 60) - (hour * 60) ))"`
second=`echo "$((${duration} % 60))"`


echo "Job time: ${hour} hour(s) $minute minute(s) and ${second} seconds"
