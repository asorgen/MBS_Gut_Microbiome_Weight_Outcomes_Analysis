#!/bin/bash

#SBATCH --partition=Orion											    # Partition/queue requested on server
##SBATCH --job-name=k2Tables
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
module load R/4.0.3

echo -e "\n${yellow}[Running Kraken2 taxa tables]\n"
/bin/date

Rscript \
${scriptPath}/BuildTaxaTable_Kraken2.R \
${inputDir} \
${outputDir} \
${funcScript}

Rscript \
${scriptPath}/NormalizeCounts_Kraken2.R \
${outputDir} \
${outputDir} \
${funcScript}


echo -e "\nCOMPLETE"
/bin/date
