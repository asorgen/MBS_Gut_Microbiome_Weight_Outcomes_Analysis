#!/bin/bash

#SBATCH --partition=Orion											    # Partition/queue requested on server
##SBATCH --time=10:00:00												    # Time limit (hrs:min:sec)
#SBATCH --nodes=1													    # Number of nodes requested
#SBATCH --ntasks-per-node=16										    # Number of CPUs (processor cores/tasks)
#SBATCH --mem=125gb													    # Memory limit
#SBATCH --output=/nobackup/afodor_research/asorgen/BariatricSurgery/outputLogs/%x.%j.out 	# Standard output
#SBATCH --error=/nobackup/afodor_research/asorgen/BariatricSurgery/outputLogs/%x.%j.err 	# Error log

logDir=/nobackup/afodor_research/asorgen/BariatricSurgery/outputLogs
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

echo -e "\n${yellow}[Running Kraken2 for sample ${ID}]\n"
/bin/date


${kraken2Path}/kraken2 \
--db ${database} \
--threads 20 \
--report ${report} \
--use-mpa-style \
--paired ${R1} ${R2} > ${K2out}


echo -e "\nCOMPLETED"
/bin/date
