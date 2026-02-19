# BioLockJ Pipeline Configuration Files

This directory contains [BioLockJ](https://github.com/BioLockJ-Dev-Team/BioLockJ) `.properties` configuration files used to orchestrate the downstream statistical analysis. Each file defines an ordered sequence of R script modules, their input data, and runtime parameters.

---

## Configuration Files

### `MetaPhlAn2_microbiome_analysis.properties` *(used for publication)*

The primary analysis pipeline used for the published results (Steffen, Sorgen et al. 2024, *Obesity*). Uses MetaPhlAn2-classified taxonomic abundance tables as input.

Modules include:
- Patient characteristics and weight loss summary statistics
- Microbiome taxa merging and normalization across taxonomic levels
- Taxa abundance changes at 1 month and 6 months post-surgery
- Mixed linear models (MLMs) of taxa over time
- Kendall correlations between taxa abundance and weight loss outcomes
- Alpha/beta diversity analysis
- Microbiome uniqueness over time
- Latent class growth analysis (LCGA) and Gaussian mixture model (GMM) trajectory modeling
- Weight loss prediction models using microbiome features

---

### `Kraken2_microbiome_analysis.properties`

Mirrors the MetaPhlAn2 primary pipeline but uses Kraken2-classified taxonomic abundance tables as input. Used as a secondary classifier comparison. Focuses on taxa abundance changes over time stratified by weight loss class across multiple timepoints.

---

### `microbiome_sex_analysis.properties`

Extends the primary analysis with sex-stratified comparisons using both MetaPhlAn2 and Kraken2 classifiers. Adds `Taxa_by_Sex.R` to test for sex-based differences in microbiome composition and weight outcomes.

---

### `quintile_analysis.properties`

Investigates alternative weight loss rank-grouping strategies (quintiles, quartiles, halves) alongside the standard tertile groupings used in the primary analysis. Covers taxa abundance over time from baseline to 24 months and baseline to 12 months.

---

### `microbiome_n124_analysis.properties`

The most comprehensive pipeline, integrating MetaPhlAn2, Kraken2, and HUMAnN2 functional pathway data. Includes random forest classification models using:
- Taxa features only
- Clinical features only
- Combined taxa and clinical features

Models are run for both 2-group and 3-group (tertile) weight outcome classifications.
