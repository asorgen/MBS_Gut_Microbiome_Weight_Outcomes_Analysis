#!/bin/bash

#SBATCH --partition=Orion											    # Partition/queue requested on server
#SBATCH --time=05:00:00												    # Time limit (hrs:min:sec)
#SBATCH --nodes=1													    # Number of nodes requested
#SBATCH --ntasks-per-node=16										    # Number of CPUs (processor cores/tasks)
#SBATCH --mem=125gb													    # Memory limit
#SBATCH --output=/nobackup/afodor_research/asorgen/BariatricSurgery/outputLogs/%x.%j.out 	# Standard output
#SBATCH --error=/nobackup/afodor_research/asorgen/BariatricSurgery/outputLogs/%x.%j.err 	# Error log
#SBATCH --job-name=readCounts										    # Job name
##SBATCH --mail-type=BEGIN,END,FAIL									    # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=asorgen@uncc.edu								    # Specified email address

logDir=/nobackup/afodor_research/asorgen/BariatricSurgery/outputLogs
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

/nobackup/afodor_research/asorgen/BariatricSurgery/scripts/pearMerge/readCounts.sh


duration=$SECONDS

hour=`echo "$(( ${duration} / 3600))"`
minute=`echo "$(( (${duration} / 60) - (hour * 60) ))"`
second=`echo "$((${duration} % 60))"`


echo "Job time: ${hour} hour(s) $minute minute(s) and ${second} seconds"
