#!/bin/bash

#SBATCH --partition=Orion											    # Partition/queue requested on server
##SBATCH --time=10:00:00												    # Time limit (hrs:min:sec)
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
module load humann2/0.11.2 

echo -e "\n${yellow}[Running HumanN2 for sample ${ID}]\n"
/bin/date


humann2 \
--input $INPUT \
--nucleotide-database ${humann2NuclDB} \
--protein-database ${humann2ProtDB} \
--metaphlan-options "--bowtie2db ${metaphlan2DB}" \
--remove-temp-output \
--threads 8 \
--output $outputDir 

mv ${outputDir}/${ID}_merge.assembled_pathcoverage.tsv ${outputDir}/pathCoverage
mv ${outputDir}/${ID}_merge.assembled_pathabundance.tsv ${outputDir}/pathAbundance
mv ${outputDir}/${ID}_merge.assembled_genefamilies.tsv ${outputDir}/geneFamilies

echo -e "\n${ID} COMPLETE"
/bin/date

