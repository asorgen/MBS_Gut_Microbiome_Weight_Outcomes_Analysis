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
├── sequence_processing/      # Upstream sequence QC and profiling pipeline
│   └── scripts/              # Shell, SLURM, and R scripts (see sequence_processing/README.md)
├── analysis/
│   ├── BLJ_config_files/     # BioLockJ pipeline configuration files
│   ├── input/
│   │   ├── metadataTables/   # Patient weight and clinical metadata (not included)
│   │   ├── MetaPhlAn2TaxaTables/  # Taxonomic abundance tables (not included)
│   │   └── Kraken2TaxaTables/     # Taxonomic abundance tables (not included)
│   └── Rscripts/             # All analysis R scripts
├── Final_Tables_Figures/     # Publication-ready outputs
└── Results_Presentations/    # Interim results and presentations
```

---

## Sequence Processing (`sequence_processing/`)

Raw paired-end whole metagenome shotgun sequencing reads were processed through the following steps before downstream analysis. See [`sequence_processing/README.md`](sequence_processing/README.md) for full details, script locations, databases, and runtimes.

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

Downstream statistical analysis is implemented in R and managed with [BioLockJ](https://github.com/BioLockJ-Dev-Team/BioLockJ). The pipeline configuration used for the publication is [`analysis/BLJ_config_files/MetaPhlAn2_microbiome_analysis.properties`](analysis/BLJ_config_files/MetaPhlAn2_microbiome_analysis.properties). See [`analysis/BLJ_config_files/README.md`](analysis/BLJ_config_files/README.md) for descriptions of all pipeline configurations.

---

## Dependencies

- **R** (≥ 4.0)
- Key R packages: `nlme`, `ggplot2`, `ggpubr`, `vegan`, `rstatix`, `randomForest`, `data.table`, `stringr`, `tidyr`, `gridExtra`
- **BioLockJ** for pipeline execution
- **KneadData** for quality trimming and host decontamination
- **PEAR** for paired-end read merging
- **Kraken2** and **MetaPhlAn2** for taxonomic profiling
- **HUMAnN2** for functional pathway profiling

---

## Data Availability

Raw sequencing data are available through NCBI. Patient metadata and weight data are not included in this repository to protect participant privacy.
