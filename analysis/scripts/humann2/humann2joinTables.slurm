#!/bin/bash

#SBATCH --partition=Orion											    # Partition/queue requested on server
#SBATCH --job-name=h2Join
#SBATCH --time=10:00:00												    # Time limit (hrs:min:sec)
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
module load humann2/0.11.2 

ROOT=/nobackup/afodor_research/asorgen/BariatricSurgery
genefamiliesINPUT=${ROOT}/05_humann2/geneFamilies
pathabundanceINPUT=${ROOT}/05_humann2/pathAbundance
pathcoverageINPUT=${ROOT}/05_humann2/pathCoverage
output=${ROOT}/08_humann2JoinTables

# mkdir ${output}



# echo -e "\n[ Join genefamilies Tables for HumanN2 ]"
# /bin/date


# humann2_join_tables \
# -i ${genefamiliesINPUT} \
# -o ${output}/humanN2_genefamilies.tsv \
# --file_name genefamilies

# echo -e "\n[ Join genefamilies Tables COMPLETE! ]"
# /bin/date


# echo -e "\n[ Renormalize genefamilies Table for HumanN2 ]"
# /bin/date


# humann2_renorm_table \
# -i ${output}/humanN2_genefamilies.tsv \
# -o ${output}/humanN2_genefamilies_cpm.tsv \
# --units cpm

# echo -e "\n[ Renormalize genefamilies Tables COMPLETE! ]"
# /bin/date





# echo -e "\n[ Join pathabundance Tables for HumanN2 ]"
# /bin/date


# humann2_join_tables \
# -i ${pathabundanceINPUT} \
# -o ${output}/humanN2_pathabundance.tsv \
# --file_name pathabundance

# echo -e "\n[ Join pathabundance Tables COMPLETE! ]"
# /bin/date


echo -e "\n[ Renormalize pathabundance Table for HumanN2 ]"
/bin/date


humann2_renorm_table \
-i ${output}/humanN2_pathabundance.tsv \
-o ${output}/humanN2_pathabundance_cpm.tsv \
--units cpm

echo -e "\n[ Renormalize pathabundance Tables COMPLETE! ]"
/bin/date







# echo -e "\n[ Join pathcoverage Tables for HumanN2 ]"
# /bin/date


# humann2_join_tables \
# -i ${pathcoverageINPUT} \
# -o ${output}/humanN2_pathcoverage.tsv \
# --file_name pathcoverage

# echo -e "\n[ Join pathcoverage Tables COMPLETE! ]"
# /bin/date


# echo -e "\n[ Renormalize pathcoverage Table for HumanN2 ]"
# /bin/date


# humann2_renorm_table \
# -i ${output}/humanN2_pathcoverage.tsv \
# -o ${output}/humanN2_pathcoverage_cpm.tsv \
# --units cpm

# echo -e "\n[ Renormalize pathcoverage Tables COMPLETE! ]"
# /bin/date

