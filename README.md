# Gut Microbiome and Bariatric Surgery Weight Outcomes

Analysis code for the following publication:

> Steffen, K.J., Sorgen, A.A., Fodor, A.A., Carroll, I.M., et al. (2024). Early changes in the gut microbiota are associated with weight outcomes over 2 years following metabolic and bariatric surgery. *Obesity* 32(11):1985–1997.

---

## Overview

This repository contains the full analysis pipeline examining whether early post-surgical gut microbiome changes predict long-term weight loss outcomes (up to 2 years) in patients undergoing Roux-en-Y gastric bypass (RYGB) or sleeve gastrectomy (SG). Microbiome profiling was performed using MetaPhlAn2 and Kraken2.

---

## Repository Structure

```
.
├── analysis/
│   ├── BLJ_config_files/   # archived BioLockJ pipeline configuration files
│   ├── Rscripts/           # analysis R scripts
│   └── run_analysis.sh     # pipeline orchestration script
├── Final_Tables_Figures/   # publication-ready outputs
└── Results_Presentations/  # interim results and presentations
```

Results are written to a sibling directory `../MBS_Gut_Microbiome_Weight_Outcomes_Results/` outside this repo. Input data is read from `../../../Data/` and is not included (see [Data Availability](#data-availability)).

---

## Sequence Processing

Raw paired-end whole metagenome shotgun sequencing reads were processed prior to this analysis using the pipeline in the [Bariatric_Metagenomic_Sequence_Processing](https://github.com/asorgen/Bariatric_Metagenomic_Sequence_Processing) repository.

| Step | Tool | Description |
|---|---|---|
| 1 | — | Copy raw FASTQ files to scratch |
| 2 | KneadData | Quality trimming and human read removal (hg19) |
| 3 | PEAR | Paired-end read merging |
| 4 | Kraken2 | Taxonomic classification |
| 5 | MetaPhlAn2 | Taxonomic classification and abundance profiling |
| 6 | HUMAnN2 | Functional pathway and gene family profiling |
| 7 | R | Build and normalize Kraken2 and MetaPhlAn2 taxa tables |
| 8 | HUMAnN2 | Join and renormalize pathway tables |

---

## Statistical Analysis (`analysis/`)

Downstream statistical analysis is implemented in R and orchestrated by [`analysis/run_analysis.sh`](analysis/run_analysis.sh). See [`analysis/README.md`](analysis/README.md) for full usage details.

### Pipeline Sections

| Section | Description |
|---|---|
| 1 | Patient & weight characteristics — BMI, excess weight loss, responder classification, tertile group assignment |
| 2 | Taxa data preparation — merge MetaPhlAn2/Kraken2 tables with patient metadata at all taxonomic levels |
| 3 | Taxa–weight associations — univariate linear models, Kendall rank correlations, surgery-type comparisons |
| 4 | Longitudinal analyses — alpha/beta diversity, microbiome uniqueness, mixed linear models for taxa trajectories |
| 5 | Group analyses & prediction — tertile comparisons, GMM trajectory clustering, LCGA, prediction models |

### Quick Start

```bash
# Run all incomplete modules
bash analysis/run_analysis.sh

# Re-run everything from scratch
bash analysis/run_analysis.sh --clean

# Run a specific section range
bash analysis/run_analysis.sh --from 2 --to 4

# Run a single module
bash analysis/run_analysis.sh 3.1_Microbiome_WeightLoss_Associations_uLM.R
```

### First-Time Setup

`analysis/.Rprofile` must exist and set the `mbs.pipe_root` option to the local repo path. This file is not tracked in git. Create it once:

```r
# analysis/.Rprofile
options(mbs.pipe_root = "/path/to/this/repo")
```

---

## Dependencies

- **R** (≥ 4.0)
- Key R packages: `nlme`, `ggplot2`, `ggpubr`, `vegan`, `rstatix`, `randomForest`, `data.table`, `stringr`, `tidyr`, `gridExtra`, `scales`

---

## Data Availability

Raw sequencing data are available through NCBI. Patient metadata and weight data are not included in this repository to protect participant privacy.
