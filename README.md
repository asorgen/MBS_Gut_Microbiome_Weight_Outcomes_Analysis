# gut-microbiome-bariatric-weight-outcomes

Analysis code for the following publication:

> Steffen, K.J., **Sorgen, A.A.**, Fodor, A.A., Carroll, I.M., et al. (2024). Early changes in the gut microbiota are associated with weight outcomes over 2 years following metabolic and bariatric surgery. *Obesity* 32(11):1985–1997.

---

## Overview

This repository contains the full analysis pipeline examining whether early post-surgical gut microbiome changes predict long-term weight loss outcomes (up to 2 years) in patients undergoing Roux-en-Y gastric bypass (RYGB) or sleeve gastrectomy (SG). Microbiome profiling was performed using MetaPhlAn2 and Kraken2.

---

## Repository Structure

```
.
├── analysis/
│   ├── BLJ_config_files/     # BioLockJ pipeline configuration files
│   ├── input/
│   │   ├── metadataTables/   # Patient weight and clinical metadata (not included)
│   │   ├── MetaPhlAn2TaxaTables/  # Taxonomic abundance tables (not included)
│   │   └── Kraken2TaxaTables/     # Taxonomic abundance tables (not included)
│   ├── Rscripts/             # All analysis R scripts
│   └── scripts/              # Shell scripts for bioinformatics processing
├── Final_Tables_Figures/     # Publication-ready outputs
└── Results_Presentations/    # Interim results and presentations
```

---

## Analysis Scripts (`analysis/Rscripts/`)

| Script | Description |
|---|---|
| `TaxaMetaMerge.R` | Merges MetaPhlAn2/Kraken2 taxonomic tables with patient metadata |
| `WeightMetaMerge.R` | Prepares and merges weight outcome data |
| `PatientCharacteristics.R` | Summarizes patient demographics |
| `ExcessWeightLoss.R` | Calculates excess weight loss (EWL) and assigns responder/non-responder groups |
| `General_Weight_Stats.R` | Descriptive weight outcome statistics |
| `BMI_Results.R` | BMI trajectory summaries |
| `Diversity_Metrics.R` | Computes alpha and beta diversity (Shannon, Bray-Curtis) |
| `Taxa_Changes.R` | Tests for changes in taxa abundance over time |
| `Taxa_over_time_MLM.R` | Mixed linear models for taxa trajectories over time |
| `Taxa_over_time_by_Surgery_MLM.R` | Surgery-type-stratified taxa trajectory models |
| `Taxa_over_time_Heatmap.R` | Heatmap visualization of taxa trajectories |
| `Taxa_over_time_Weight_Class.R` | Taxa trajectories stratified by weight loss class |
| `Taxa_by_Surgery.R` | Taxonomic composition differences by surgery type |
| `Taxa_pValuePlots_M.R` | p-value summary plots for taxa associations |
| `Microbiome_WeightLoss_Associations_kendall.R` | Kendall rank correlations between taxa and weight loss |
| `Microbiome_WeightLoss_Associations_uLM.R` | Univariate linear models for taxa-weight associations |
| `Assign_WLgroups.R` | Assigns patients to weight loss trajectory groups |
| `Barplot_WLgroups.R` | Bar plots of taxonomic composition by weight loss group |
| `Tertile_Ranks.R` / `Tertile_Analysis.R` | Tertile-based stratification and ranked analyses |
| `GMM_Analysis.R` / `GMM_Analysis_Surgery.R` / `GMM1_2_Analysis.R` | Gaussian mixture model clustering of weight trajectories |
| `LCGA2_Analysis.R` | Latent class growth analysis (LCGA) for weight trajectory groups |
| `LCGA_Heatmap.R` / `Prep_Heatmap.R` | Heatmap visualizations of taxa by LCGA group |
| `PredictionModels.R` | Logistic and linear prediction models for weight outcomes |
| `MissingData_Analysis.R` | Evaluates missingness patterns in the dataset |
| `Uniqueness.R` | Microbiome uniqueness analysis |
| `functions.R` | Shared utility functions used across scripts |

Scripts prefixed with `x_` are exploratory or developmental versions not used in the final published analysis.

---

## Pipeline Management

Analysis pipelines are managed with [BioLockJ](https://github.com/BioLockJ-Dev-Team/BioLockJ). Configuration files in `analysis/BLJ_config_files/` define the full reproducible pipeline.

---

## Dependencies

- **R** (≥ 4.0)
- Key R packages: `nlme`, `ggplot2`, `ggpubr`, `vegan`, `rstatix`, `randomForest`, `data.table`, `stringr`, `tidyr`, `gridExtra`
- **BioLockJ** for pipeline execution
- **MetaPhlAn2** and **Kraken2** for upstream taxonomic profiling

---

## Data Availability

Raw sequencing data are available through NCBI. Patient metadata and weight data are not included in this repository to protect participant privacy.
