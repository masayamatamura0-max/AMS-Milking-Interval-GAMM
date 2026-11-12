# AMS milking interval GAMM: supplementary material and analysis code

This repository accompanies a manuscript examining early-lactation milking interval and milking interval variability clusters in automated milking systems (AMS). The supplied supplementary document reports **basis-dimension validation (`mgcv::k.check`)** for 11 outcomes, 6 clustering analyses, and parity-specific DIM smooths.

## Contents

- [`supplementary/Supplemental_File_S1_GAMM_model_validation.docx`](supplementary/Supplemental_File_S1_GAMM_model_validation.docx): Supplemental File S1, Tables S1–S12.
- [`analysis/run_all_outcomes_GAMM.R`](analysis/run_all_outcomes_GAMM.R): extended R analysis code based on the original ECM/SCC script.
- [`config/outcome_columns.csv`](config/outcome_columns.csv): optional mapping from outcome identifiers to processed-data CSV column names.

## Outcomes and smoothing basis dimensions

| Outcome | k |
|---|---:|
| Mean milking interval (h) | 25 |
| Milking frequency (no./d) | 25 |
| Milk yield (kg/d) | 25 |
| Energy-corrected milk (kg/d) | 100 |
| AMS concentrate delivered (kg/d) | 25 |
| Milk fat (%) | 25 |
| Milk protein (%) | 50 |
| Milk lactose (%) | 25 |
| log10 somatic cell count | 100 |
| Milk temperature (°C) | 25 |
| Rumination time (h/d) | 25 |

## Running the analysis

Install R packages `mgcv`, `emmeans`, `dplyr`, `tidyr`, `readr`, `stringr`, and `tibble`.

From the repository root:

```bash
Rscript analysis/run_all_outcomes_GAMM.R data/processed results/GAMM_all_outcomes
```

The script expects six processed daily CSV files, with filenames:

```
MI_1_21_daily_for_GAMM.csv
MI_1_42_daily_for_GAMM.csv
MI_1_63_daily_for_GAMM.csv
MIV_1_21_daily_for_GAMM.csv
MIV_1_42_daily_for_GAMM.csv
MIV_1_63_daily_for_GAMM.csv
```

Each CSV must include `Cluster`, `Parity_group`, `DIM_corrected`, `herd_id`, `lactation_segment_id`, and outcome-specific source columns. The original script derives ECM from `Sum_Cond_MilkYield`, `Mean_Daily_MdpFatPercentage`, and `Mean_Daily_MdpProteinPercentage`; this is preserved.

**Important:** The raw input CSVs were not supplied with this repository package. Several outcome source column names are therefore provisional. Compare the actual CSV headers with the `Source_candidates` table in the script. To override, add rows to `config/outcome_columns.csv`, e.g. `Rumination_time,Actual_Rumination_Column`. Do not publish identifiable farm or animal-level data without appropriate permission.

The script fits a preliminary GAMM to estimate residual lag-1 correlation and then a Gaussian `mgcv::bam` model with AR(1), parity-specific smooths, herd and lactation random effects. It exports model information, EMMs, Tukey comparisons, smooth diagnostics, and k-check results.

## Reproducibility limitations

- The supplied source R script only covered ECM and log10 SCC; the additional outcome mappings and script expansion have **not been validated against the six processed input datasets**.
- R is not installed in the preparation environment, so no models have been executed or results regenerated. Existing tables in Supplemental File S1 are retained unchanged; they are **not claimed to have been reproduced** by this package.
- The `k.check` p-values are simulation-based and may vary across runs and package versions.
- The original R script's statistical model, including its EMM and AR(1) handling, is retained as a starting point and should be independently reviewed before publication.

## Citation and access

Please cite the associated manuscript once published. Link directly to this repository's supplementary document in the manuscript or data-availability statement, subject to journal requirements. For a permanent versioned DOI, consider archiving a GitHub release on Zenodo.
