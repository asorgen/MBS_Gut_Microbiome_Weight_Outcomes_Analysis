#!/bin/bash

#SBATCH --job-name=KneadData_zip 																# Job name
#SBATCH --partition=Orion 																		# Partition/queue requested on server
#SBATCH --time=100:00:00																		# Time limit (hrs:min:sec)
#SBATCH --nodes=1 																				# Number of nodes requested
#SBATCH --ntasks-per-node=16 																	# Number of CPUs (processor cores/tasks)
#SBATCH --mem=125gb																				# Memory limit
#SBATCH --output=/nobackup/afodor_research/asorgen/BariatricSurgery/outputLogs/%x.%j.out 		# Standard output
#SBATCH --error=/nobackup/afodor_research/asorgen/BariatricSurgery/outputLogs/%x.%j.err 		# Error log
#SBATCH --mail-type=BEGIN,END,FAIL 																# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=asorgen@uncc.edu 															# Specified email address

SECONDS=0

logDir=/nobackup/afodor_research/asorgen/BariatricSurgery/outputLogs
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

# zip -r /nobackup/afodor_research/asorgen/BariatricSurgery/01_kneaddata_hg19.zip /nobackup/afodor_research/asorgen/BariatricSurgery/01_kneaddata_hg19/
unzip /nobackup/afodor_research/asorgen/BariatricSurgery/01_kneaddata_hg19.zip -d /nobackup/afodor_research/asorgen/BariatricSurgery_GMHI/BS_data -x "*/*_unmatched_*" "*/*_contam_*" "*/*.log"

duration=$SECONDS

hour=`echo "$(( ${duration} / 3600))"`
minute=`echo "$(( (${duration} / 60) - (hour * 60) ))"`
second=`echo "$((${duration} % 60))"`


echo "Job time: ${hour} hour(s) $minute minute(s) and ${second} seconds"
