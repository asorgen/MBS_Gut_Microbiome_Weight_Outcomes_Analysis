# Analysis

Statistical analysis code for examining gut microbiome changes and weight loss outcomes following bariatric surgery. See the [top-level README](../README.md) for a full project overview.

---

## Directory Structure

```
analysis/
├── BLJ_config_files/   # archived BioLockJ pipeline configuration files
├── Rscripts/           # analysis R scripts
│   ├── 1.0_–5.3_*.R    # numbered pipeline modules (run in order)
│   ├── functions.R     # shared utility library
│   └── x_*             # exploratory/archived scripts (excluded from pipeline)
└── run_analysis.sh     # pipeline orchestration script
```

---

## Running the Pipeline

All analysis is orchestrated via `run_analysis.sh`. Results are written to `../../MBS_Gut_Microbiome_Weight_Outcomes_Results/` (sibling to the repo root). Input data is read from `../../../Data/`.

```bash
# Run all incomplete modules
bash run_analysis.sh

# Re-run everything from scratch
bash run_analysis.sh --clean

# Run a specific section range (1–5)
bash run_analysis.sh --from 2 --to 4

# Run a single module with default parameters
bash run_analysis.sh 3.1_Microbiome_WeightLoss_Associations_uLM.R

# Run a single module with custom parameters
bash run_analysis.sh 2.0_TaxaMetaMerge.R MetaPhlAn2
```

Each module writes output to `<results_root>/<ModuleName>/output/`, logs to `<results_root>/<ModuleName>/log/<ModuleName>.log`, and creates `STARTED`, `COMPLETE`, or `FAILED` flag files. When a module fails, check the log file first.

### First-Time Setup

`analysis/.Rprofile` must exist and set the `mbs.pipe_root` option to the local repo path. This file is not tracked in git. Create it once:

```r
# analysis/.Rprofile
options(mbs.pipe_root = "/path/to/this/repo")
```

---

## Pipeline Sections

| Section | Scripts | Description |
|---|---|---|
| 1 | `1.0_`–`1.7_` | Patient & weight characteristics — BMI, EWL, responder classification, group assignment |
| 2 | `2.0_` | Taxa data preparation — merge MetaPhlAn2/Kraken2 tables with metadata |
| 3 | `3.0_`–`3.3_` | Taxa–weight associations — linear models, Kendall correlations, surgery comparisons |
| 4 | `4.0_`–`4.5_` | Longitudinal analyses — diversity, uniqueness, MLM trajectories, missing data |
| 5 | `5.0_`–`5.3_` | Group analyses & prediction — tertile comparisons, GMM clustering, LCGA, prediction models |

---

## R Scripts (`Rscripts/`)

### Pipeline Modules

| Script | Description |
|---|---|
| `1.0_PatientCharacteristics.R` | Patient demographics summary |
| `1.1_BMI_Results.R` | BMI trajectory summaries |
| `1.2_ExcessWeightLoss.R` | Calculates EWL; assigns responder/non-responder groups |
| `1.3_WeightMetaMerge.R` | Merges weight outcome data with metadata |
| `1.4_Assign_WLgroups.R` | Assigns patients to weight loss trajectory groups |
| `1.5_Barplot_WLgroups.R` | Bar plots of weight loss group composition |
| `1.6_Tertile_Ranks.R` | Tertile-based group assignment and ranking |
| `1.7_General_Weight_Stats.R` | Descriptive weight outcome statistics |
| `2.0_TaxaMetaMerge.R` | Merges MetaPhlAn2/Kraken2 taxonomic tables with patient metadata |
| `3.0_Taxa_Changes.R` | Tests for changes in taxa abundance over time |
| `3.1_Microbiome_WeightLoss_Associations_uLM.R` | Univariate linear models for taxa–weight associations |
| `3.2_Microbiome_WeightLoss_Associations_kendall.R` | Kendall rank correlations between taxa and weight loss |
| `3.3_Taxa_by_Surgery.R` | Taxonomic composition differences by surgery type |
| `4.0_Diversity_Metrics.R` | Alpha and beta diversity (Shannon, Bray-Curtis) |
| `4.1_Uniqueness.R` | Microbiome uniqueness analysis |
| `4.2_Taxa_over_time_MLM.R` | Mixed linear models for taxa trajectories over time |
| `4.3_Taxa_over_time_by_Surgery_MLM.R` | Surgery-type-stratified taxa trajectory models |
| `4.4_Taxa_over_time_Heatmap.R` | Heatmap visualization of taxa trajectories |
| `4.5_MissingData_Analysis.R` | Evaluates missingness patterns in the dataset |
| `5.0_Tertile_Analysis.R` | Tertile-based group comparisons |
| `5.1_GMM_Analysis_Surgery.R` | Gaussian mixture model clustering of weight trajectories by surgery type |
| `5.2_LCGA2_Analysis.R` | Latent class growth analysis (LCGA) for weight trajectory groups |
| `5.3_PredictionModels.R` | Logistic and linear prediction models for weight outcomes |
| `functions.R` | Shared utility functions sourced by all numbered scripts |
