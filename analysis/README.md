# Analysis

Statistical analysis code for examining gut microbiome changes and weight loss outcomes following bariatric surgery. See the [top-level README](../README.md) for a full project overview.

---

## Directory Structure

```
analysis/
├── BLJ_config_files/     # BioLockJ pipeline configuration files
├── input/
│   ├── metadataTables/   # Patient weight and clinical metadata (not included)
│   ├── MetaPhlAn2TaxaTables/  # Taxonomic abundance tables (not included)
│   └── Kraken2TaxaTables/     # Taxonomic abundance tables (not included)
└── Rscripts/             # All analysis R scripts
```

---

## Pipeline Management

Analysis pipelines are managed with [BioLockJ](https://github.com/BioLockJ-Dev-Team/BioLockJ). Configuration files in `BLJ_config_files/` define the full reproducible pipeline.

---

## R Scripts (`Rscripts/`)

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

