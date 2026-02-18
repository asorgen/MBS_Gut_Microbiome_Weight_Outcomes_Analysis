#!/bin/bash

#SBATCH --partition=Orion											    # Partition/queue requested on server
#SBATCH --job-name=copy										    # Job name
#SBATCH --time=08:00:00												    # Time limit (hrs:min:sec)
#SBATCH --nodes=1													    # Number of nodes requested
#SBATCH --ntasks-per-node=16										    # Number of CPUs (processor cores/tasks)
#SBATCH --mem=125gb													    # Memory limit
#SBATCH --mail-type=BEGIN,END,FAIL									    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=asorgen@uncc.edu								    # Specified email address
#SBATCH --output=/nobackup/afodor_research/asorgen/BariatricSurgery/outputLogs/%x.%j.out 	# Standard output
#SBATCH --error=/nobackup/afodor_research/asorgen/BariatricSurgery/outputLogs/%x.%j.err 	# Error log

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

### Working script
${ScriptDir}/copy/copy.sh
