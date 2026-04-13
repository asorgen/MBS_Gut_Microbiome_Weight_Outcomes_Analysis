#!/bin/bash
# Pipeline orchestration script for MBS_Gut_Microbiome_Weight_Outcomes_Analysis
# Replaces BioLockJ (.properties) pipeline execution — no Java/Docker required.
#
# Usage:
#   bash run_analysis.sh                          # run all incomplete modules
#   bash run_analysis.sh --clean                  # wipe flags and re-run all modules
#   bash run_analysis.sh --from <N>               # start at section N (1–5)
#   bash run_analysis.sh --to <N>                 # stop after section N (1–5)
#   bash run_analysis.sh --from <N> --to <M>      # run sections N through M
#   bash run_analysis.sh <Script.R>               # run a single module (uses default params)
#   bash run_analysis.sh <Script.R> [args]        # run a single module with custom params
#
# Sections:
#   1 = patient & weight characteristics
#   2 = taxa data preparation
#   3 = taxa-weight associations
#   4 = longitudinal, diversity, and uniqueness
#   5 = group analyses and prediction models

set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)                         # absolute path to repo root (one level up from analysis/)
results_root="/Users/aliciasorgen/UNCC/Projects/Bariatric_Surgery/Git_Repositories/MBS_Gut_Microbiome_Weight_Outcomes_Analysis/MBS_Gut_Microbiome_Weight_Outcomes_Results"  # sibling results directory
classifier="MetaPhlAn2"


# ─── Functions ─────────────────────────────────────────────────────────────────

    error() { echo "ERROR: $1" >&2; exit 1; }

    # run_module <ModuleName> <Script.R> [arg1 arg2 ...]
    #   ModuleName : label used for Results/ subdirectory, log, and script snapshot
    #   Script.R   : file under analysis/Rscripts/
    #   args       : params passed to the R script AFTER the repo root path
    run_module() {
        local name=$1
        local script=$2
        shift 2
        local scriptPath="${root}/analysis/Rscripts/${script}"
        [[ ! -f "${scriptPath}" ]] && error "Script not found: analysis/Rscripts/${script}"

        local moduleDir="${results_root}/${name}"
        mkdir -p "${moduleDir}/input" "${moduleDir}/output" "${moduleDir}/log" "${moduleDir}/script"

        # In clean mode, wipe everything and re-run; otherwise skip completed modules
        if [[ "${clean_mode}" == true ]]; then
            rm -f "${moduleDir}/script"/*.R
            rm -f "${moduleDir}/log"/*.log
            rm -f "${moduleDir}/STARTED" "${moduleDir}/COMPLETE" "${moduleDir}/FAILED"
        elif [[ -f "${moduleDir}/COMPLETE" ]]; then
            echo "Skipping ${name} (COMPLETE)"
            return 0
        else
            rm -f "${moduleDir}/STARTED" "${moduleDir}/FAILED"
        fi

        # Copy script and shared functions file into module script dir
        cp "${scriptPath}" "${root}/analysis/Rscripts/functions.R" "${moduleDir}/script/"

        local logFile="${moduleDir}/log/${name}.log"
        local elapsed hh mm ss

        echo -e "Running ${name} (${script})"
        touch "${moduleDir}/STARTED"
        local start=$SECONDS

        if RESULTS_ROOT="${results_root}" INPUT_ROOT="${root}/Results/input" Rscript "${scriptPath}" "${root}" "$@" >> "${logFile}" 2>&1; then
            elapsed=$(( SECONDS - start ))
            hh=$(( elapsed / 3600 )); mm=$(( (elapsed % 3600) / 60 )); ss=$(( elapsed % 60 ))
            rm -f "${moduleDir}/STARTED"
            touch "${moduleDir}/COMPLETE"
            printf "  %-50s [%02d:%02d:%02d]\n" "${name}" "$hh" "$mm" "$ss"
        else
            elapsed=$(( SECONDS - start ))
            hh=$(( elapsed / 3600 )); mm=$(( (elapsed % 3600) / 60 )); ss=$(( elapsed % 60 ))
            rm -f "${moduleDir}/STARTED"
            touch "${moduleDir}/FAILED"
            printf "  %-50s [%02d:%02d:%02d] FAILED\n" "${name}" "$hh" "$mm" "$ss"
            error "Module ${name} failed. Check ${logFile}"
        fi
    }

# ─── Defaults ─────────────────────────────────────────────────────────────────

    from_section=1
    to_section=5
    clean_mode=false

# ─── Single-module mode ────────────────────────────────────────────────────────

    if [[ "${1:-}" == *.R ]]; then
        script=$1; shift
        [[ ! -f "${root}/analysis/Rscripts/${script}" ]] \
            && error "Module not found: analysis/Rscripts/${script}"
        run_module "${script%.R}" "${script}" "$@"
        exit 0
    fi

# ─── Parse --from / --to flags ────────────────────────────────────────────────

    while [[ "${1:-}" == --* ]]; do
        case "$1" in
            --from)  from_section=$2; shift 2 ;;
            --to)    to_section=$2;   shift 2 ;;
            --clean) clean_mode=true; shift ;;
            *) error "Unknown flag: $1" ;;
        esac
    done

    skip() { [[ $1 -lt $from_section ]]; }
    stop() { [[ "${to_section}" != "5" ]] && [[ $1 -ge $to_section ]]; }

# ─── 1. Patient & weight characteristics ──────────────────────────────────────

    skip 1 || run_module 1.0_PatientCharacteristics  1.0_PatientCharacteristics.R
    skip 1 || run_module 1.1_BMI_Results             1.1_BMI_Results.R
    skip 1 || run_module 1.2_ExcessWeightLoss        1.2_ExcessWeightLoss.R
    skip 1 || run_module 1.3_WeightMetaMerge         1.3_WeightMetaMerge.R
    skip 1 || run_module 1.4_Assign_WLgroups         1.4_Assign_WLgroups.R         Tertile
    skip 1 || run_module 1.5_Barplot_WLgroups        1.5_Barplot_WLgroups.R        Tertile
    skip 1 || run_module 1.6_Tertile_Ranks           1.6_Tertile_Ranks.R
    skip 1 || run_module 1.7_General_Weight_Stats    1.7_General_Weight_Stats.R    0 1 6 12 18 24

    stop 1 && exit 0

# ─── 2. Taxa data preparation ──────────────────────────────────────────────────

    skip 2 || run_module 2.0_TaxaMetaMerge    2.0_TaxaMetaMerge.R    ${classifier}
    stop 2 && exit 0

# ─── 3. Taxa-weight associations ──────────────────────────────────────────────

    skip 3 || run_module 3.0_Taxa_Changes_1M                            3.0_Taxa_Changes.R                            ${classifier} 1
    skip 3 || run_module 3.0_Taxa_Changes_6M                            3.0_Taxa_Changes.R                            ${classifier} 6
    skip 3 || run_module 3.1_Microbiome_WeightLoss_Associations_uLM     3.1_Microbiome_WeightLoss_Associations_uLM.R     ${classifier} 0 1 6 12 18 24
    skip 3 || run_module 3.2_Microbiome_WeightLoss_Associations_kendall 3.2_Microbiome_WeightLoss_Associations_kendall.R ${classifier} 0 1 6 12 18 24
    skip 3 || run_module 3.3_Taxa_by_Surgery                            3.3_Taxa_by_Surgery.R                            ${classifier} 0 1 6 12 18 24

    stop 3 && exit 0

# ─── 4. Longitudinal, diversity & uniqueness ──────────────────────────────────

    skip 4 || run_module 4.0_Diversity_Metrics             4.0_Diversity_Metrics.R             ${classifier}
    skip 4 || run_module 4.1_Uniqueness                    4.1_Uniqueness.R                    ${classifier}
    skip 4 || run_module 4.2_Taxa_over_time_MLM            4.2_Taxa_over_time_MLM.R            24
    skip 4 || run_module 4.3_Taxa_over_time_by_Surgery_MLM 4.3_Taxa_over_time_by_Surgery_MLM.R 24
    skip 4 || run_module 4.4_Taxa_over_time_Heatmap        4.4_Taxa_over_time_Heatmap.R        ${classifier} 24
    skip 4 || run_module 4.5_MissingData_Analysis          4.5_MissingData_Analysis.R          ${classifier} 24

    stop 4 && exit 0

# ─── 5. Group analyses & prediction models ────────────────────────────────────

    skip 5 || run_module 5.0_Tertile_Analysis    5.0_Tertile_Analysis.R    24 Tertile
    skip 5 || run_module 5.1_GMM_Analysis_RYGB   5.1_GMM_Analysis_Surgery.R RYGB
    skip 5 || run_module 5.1_GMM_Analysis_SG     5.1_GMM_Analysis_Surgery.R SG
    skip 5 || run_module 5.2_LCGA2_RYGB_2Class   5.2_LCGA2_Analysis.R       LCGA 2 RYGB Univariate Model2
    skip 5 || run_module 5.3_PredictionModels    5.3_PredictionModels.R     ${classifier} 24

    stop 5 && exit 0

echo -e "\nPipeline complete."
