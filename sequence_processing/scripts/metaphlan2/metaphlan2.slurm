#!/bin/bash

#SBATCH --partition=Orion											    # Partition/queue requested on server
##SBATCH --time=2:00:00												    # Time limit (hrs:min:sec)
#SBATCH --nodes=1													    # Number of nodes requested
#SBATCH --ntasks-per-node=16										    # Number of CPUs (processor cores/tasks)
#SBATCH --mem=125gb													    # Memory limit
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
# module load anaconda3/2020.11
module load metaphlan/2.8.1

echo -e "\n${yellow}[Running MetaPhlAn2 for sample ${ID}]\n"
/bin/date


metaphlan2.py \
${R1},${R2} \
--bowtie2db ${database} \
--bowtie2out ${bowtieOut} \
--nproc 8 \
-t rel_ab_w_read_stats \
--input_type fastq \
-o ${profileOut}

echo -e "\n${ID} COMPLETE"
/bin/date
