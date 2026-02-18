
# Bariatric Surgery MetagenomicSequencing Processing

## Copy sequence files from all source directories into one scratch folder

`./scripts/copy.sh`
Working script: copy.slurm
Approx run time = 00:38:40


## Run kneaddata on sequences to filter out human DNA

`./scripts/kneaddata/kneaddata_hg19.sh`
Working script: kneaddata_hg19.slurm
Requested walltime = 02:00:00
Database: /nobackup/afodor_research/databases/kneaddata_database/Homo_sapiens
Output: 01_kneaddata_hg19
Logs: outputLogs/knead_hg19

**Trimmomatic parameters**
ILLUMINACLIP:/users/asorgen/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:8:true \
LEADING:3
TRAILING:3
SLIDINGWINDOW:4:15
MINLEN:50

	
## Merge forward and reverse reads using PEAR

`./scripts/pearMerge/pearMerge.sh`
Working script: pearMerge.slurm
Requested walltime = 10:00:00
Output: 02_pearMerge
Logs: outputLogs/pearMerge



## Taxonomic classification using Kraken2

`./scripts/kraken2/kraken2.sh`
Working script: kraken2.slurm
Requested walltime = 01:00:00
Database: /scratch/asorgen/kraken2_standardDB
Output: 03_kraken2
Logs: outputLogs/kraken2


## Taxonomic classification using MetaPhlAn2

`./scripts/metaphlan2/metaphlan2.sh`
Working script: metaphlan2.slurm
Requested walltime = 05:00:00
Database: /scratch/asorgen/mpa_v20_m200
Output: 04_metaphlan2
Logs: outputLogs/metaphlan2


## Merge MetaPhlAn2 output into single file

`./scripts/metaphlan2/metaphlanMerge.sh`
Output: 04_metaphlan2


## Pathway classification using HumanN2

`./scripts/humann2/humann2.sh`
Working script: humann2.slurm
Requested walltime = 10:00:00
Nucleotide database: /nobackup/afodor_research/databases/humanN2/chocophlan
Protein database: /nobackup/afodor_research/databases/humanN2/uniref
Output: 05_humann2
Logs: outputLogs/humann2



## Build Kraken2 taxa tables

Working script: BuildTaxaTable_Kraken2.R
`module load R/4.0.3`
`Rscript /nobackup/afodor_research/asorgen/BariatricSurgery/scripts/kraken2/BuildTaxaTable_Kraken2.R`
Output: 06_kraken2TaxaTables


## Normalize Kraken2 taxa tables

`sbatch scripts/kraken2/NormalizeCounts_Kraken2.slurm`
Working script: NormalizeCounts_Kraken2.R
Output: 06_kraken2TaxaTables


## Build MetaPhlAn2 taxa tables

`./scripts/metaphlan2/BuildTaxaTable_MetaPhlAn2.sh`
Working script: BuildTaxaTable_MetaPhlAn2.slurm
	Calls working script: BuildTaxaTable_MetaPhlAn2.R and 
	Calls working script: NormalizeCounts_MetaPhlAn2.R
Output: 07_metaphlan2TaxaTables
_MetaPhlAn2 default output as relative abundance and raw counts are estimated_



## Join and renormalize HumanN2 tables

Working script: humann2joinTables.slurm
`sbatch scripts/humann2/humann2joinTables.slurm`
Output: 08_humann2JoinTables
