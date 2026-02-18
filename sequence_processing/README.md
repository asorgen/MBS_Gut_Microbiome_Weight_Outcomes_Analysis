# Bariatric Surgery Metagenomic Sequence Processing

Pipeline for quality filtering, host decontamination, taxonomic classification, and functional profiling of paired-end whole metagenome shotgun sequencing data from bariatric surgery patients.

---

## Pipeline Overview

```
Raw paired-end FASTQ reads
        │
        ▼
1. Copy sequences to scratch
        │
        ▼
2. KneadData — quality trim & remove human reads (hg19)
        │
        ▼
3. PEAR — merge paired-end reads
        │
        ├──────────────────────┬──────────────────────┐
        ▼                      ▼                      ▼
4. Kraken2               5. MetaPhlAn2           7. HUMAnN2
   Taxonomic                Taxonomic              Functional
   Classification           Classification         Profiling
        │                      │                      │
        ▼                      ▼                      ▼
6. Build & normalize     8. Build taxa           9. Join & renormalize
   Kraken2 taxa tables      tables                  pathway tables
```

---

## Directory Structure

```
sequence_processing/
├── scripts/
│   ├── copy/               # File staging scripts
│   ├── kneaddata_hg19/     # Host decontamination scripts
│   ├── pearMerge/          # Read merging scripts
│   ├── kraken2/            # Kraken2 classification & table-building scripts
│   ├── metaphlan2/         # MetaPhlAn2 classification & table-building scripts
│   ├── humann2/            # HUMAnN2 profiling scripts
│   └── functions.R         # Shared R utility functions
├── 01_kneaddata_hg19/      # KneadData output (not tracked)
├── 02_pearMerge/           # PEAR merged reads (not tracked)
├── 03_kraken2/             # Kraken2 output (not tracked)
├── 04_metaphlan2/          # MetaPhlAn2 output (not tracked)
├── 05_humann2/             # HUMAnN2 output (not tracked)
├── 06_kraken2TaxaTables/   # Kraken2 taxa tables (not tracked)
├── 07_metaphlan2TaxaTables/# MetaPhlAn2 taxa tables (not tracked)
├── 08_humann2JoinTables/   # HUMAnN2 joined tables (not tracked)
└── outputLogs/             # SLURM job logs (not tracked)
```

---

## Configuration

All personal paths, database locations, tool paths, and SLURM settings are stored in a single configuration file that is **not tracked by git** to keep sensitive paths private.

### Setup

1. Copy the example config to create your own:
   ```bash
   cp scripts/config.sh.example scripts/config.sh
   ```

2. Edit `scripts/config.sh` and fill in your values:

   | Variable | Description |
   |---|---|
   | `Dir` | Base working directory for all pipeline data |
   | `rawSeqDir` | Location of raw input FASTQ files |
   | `ScriptDir` | Path to this `scripts/` directory in your clone |
   | `inputList` | Sample file list (one ID per line) |
   | `kneaddataDB` | KneadData human reference database |
   | `kraken2DB` | Kraken2 standard database |
   | `metaphlan2DB` | MetaPhlAn2 bowtie2 database |
   | `humann2NuclDB` | HUMAnN2 ChocoPhlAn nucleotide database |
   | `humann2ProtDB` | HUMAnN2 UniRef protein database |
   | `kraken2Path` | Path to Kraken2 executable |
   | `pearPath` | Path to PEAR executable |
   | `trimmomaticAdapters` | Path to Trimmomatic adapter FASTA |
   | `email` | Email address for SLURM job notifications |
   | `slurm_out` | Directory for SLURM output/error logs |

All scripts source `config.sh` automatically at runtime — no other path changes are needed.

---

## Steps

### 1. Copy Sequences to Scratch

Copies raw FASTQ files from all source directories into a single scratch folder for processing.

```bash
./scripts/copy/copy.sh
```

- **SLURM script:** `scripts/copy/copy.slurm`
- **Approx. runtime:** 00:38:40

---

### 2. Quality Filtering & Host Decontamination (KneadData)

Trims adapter sequences and low-quality bases, then removes reads mapping to the human reference genome (hg19).

```bash
./scripts/kneaddata_hg19/kneaddata_hg19.sh
```

- **SLURM script:** `scripts/kneaddata_hg19/kneaddata_hg19.slurm`
- **Requested walltime:** 02:00:00
- **Reference database:** `/nobackup/afodor_research/databases/kneaddata_database/Homo_sapiens`
- **Output:** `01_kneaddata_hg19/`
- **Logs:** `outputLogs/knead_hg19/`

**Trimmomatic parameters:**

| Parameter | Value |
|---|---|
| ILLUMINACLIP | `TruSeq3-PE-2.fa:2:30:10:8:true` |
| LEADING | 3 |
| TRAILING | 3 |
| SLIDINGWINDOW | 4:15 |
| MINLEN | 50 |

---

### 3. Merge Paired-End Reads (PEAR)

Merges overlapping forward and reverse reads into single assembled reads.

```bash
./scripts/pearMerge/pearMerge.sh
```

- **SLURM script:** `scripts/pearMerge/pearMerge.slurm`
- **Requested walltime:** 10:00:00
- **Output:** `02_pearMerge/`
- **Logs:** `outputLogs/pearMerge/`

---

### 4. Taxonomic Classification (Kraken2)

Classifies reads against the Kraken2 standard database.

```bash
./scripts/kraken2/kraken2.sh
```

- **SLURM script:** `scripts/kraken2/kraken2.slurm`
- **Requested walltime:** 01:00:00
- **Database:** `/scratch/asorgen/kraken2_standardDB`
- **Output:** `03_kraken2/`
- **Logs:** `outputLogs/kraken2/`

---

### 5. Taxonomic Classification (MetaPhlAn2)

Profiles microbial community composition using MetaPhlAn2 marker genes.

```bash
./scripts/metaphlan2/metaphlan2.sh
```

- **SLURM script:** `scripts/metaphlan2/metaphlan2.slurm`
- **Requested walltime:** 05:00:00
- **Database:** `/scratch/asorgen/mpa_v20_m200`
- **Output:** `04_metaphlan2/`
- **Logs:** `outputLogs/metaphlan2/`

#### Merge MetaPhlAn2 profiles

Merges per-sample MetaPhlAn2 output into a single abundance table.

```bash
./scripts/metaphlan2/metaphlanMerge.sh
```

- **Output:** `04_metaphlan2/`

---

### 6. Build & Normalize Kraken2 Taxa Tables

Constructs and normalizes taxa count tables from Kraken2 output.

```bash
# Build taxa tables
sbatch scripts/kraken2/BuildTaxaTable_Kraken2.slurm

# Normalize counts
sbatch scripts/kraken2/NormalizeCounts_Kraken2.slurm
```

- **R scripts:** `scripts/kraken2/BuildTaxaTable_Kraken2.R`, `scripts/kraken2/NormalizeCounts_Kraken2.R`
- **Output:** `06_kraken2TaxaTables/`

---

### 7. Functional Pathway Profiling (HUMAnN2)

Profiles functional pathways and gene families from merged reads.

```bash
./scripts/humann2/humann2.sh
```

- **SLURM script:** `scripts/humann2/humann2.slurm`
- **Requested walltime:** 10:00:00
- **Nucleotide database:** `/nobackup/afodor_research/databases/humanN2/chocophlan`
- **Protein database:** `/nobackup/afodor_research/databases/humanN2/uniref`
- **Output:** `05_humann2/` (subdirs: `pathCoverage/`, `pathAbundance/`, `geneFamilies/`)
- **Logs:** `outputLogs/humann2/`

---

### 8. Build MetaPhlAn2 Taxa Tables

Constructs and normalizes taxa tables from MetaPhlAn2 profiles. Note: MetaPhlAn2 default output is relative abundance; raw counts are estimated.

```bash
./scripts/metaphlan2/BuildTaxaTable_MetaPhlAn2.sh
```

- **SLURM script:** `scripts/metaphlan2/BuildTaxaTable_MetaPhlAn2.slurm`
- **R scripts:** `scripts/metaphlan2/BuildTaxaTable_MetaPhlAn2.R`, `scripts/metaphlan2/NormalizeCounts_MetaPhlAn2.R`
- **Output:** `07_metaphlan2TaxaTables/`

---

### 9. Join & Renormalize HUMAnN2 Tables

Joins per-sample HUMAnN2 output files and renormalizes for downstream analysis.

```bash
sbatch scripts/humann2/humann2joinTables.slurm
```

- **Output:** `08_humann2JoinTables/`
