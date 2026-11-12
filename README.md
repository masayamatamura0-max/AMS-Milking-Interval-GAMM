# AMS milking interval GAMM: supplementary material and analysis code

This repository accompanies a manuscript examining early-lactation milking interval and milking interval variability clusters in automated milking systems (AMS). The supplied supplementary document reports **basis-dimension validation (`mgcv::k.check`)** for 11 outcomes, 6 clustering analyses, and parity-specific DIM smooths.

## Contents

- [`supplementary/Supplemental_File_S1_GAMM_model_validation.docx`](supplementary/Supplemental_File_S1_GAMM_model_validation.docx): Supplemental File S1, Tables S1–S12.
- [`analysis/run_all_outcomes_GAMM.R`](analysis/run_all_outcomes_GAMM.R): extended R analysis code based on the original ECM/SCC script.

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

The script expects six processed daily CSV files, with filenames:

```
MI_1_21_daily_for_GAMM.csv
MI_1_42_daily_for_GAMM.csv
MI_1_63_daily_for_GAMM.csv
MIV_1_21_daily_for_GAMM.csv
MIV_1_42_daily_for_GAMM.csv
MIV_1_63_daily_for_GAMM.csv
```

Each CSV must include `Cluster`, `Parity_group`, `DIM_corrected`, `herd_id`, `lactation_segment_id`, and outcome-specific source columns.

The script fits a preliminary GAMM to estimate residual lag-1 correlation and then a Gaussian `mgcv::bam` model with AR(1), parity-specific smooths, herd and lactation random effects. It exports model information, EMMs, Tukey comparisons, smooth diagnostics, and k-check results.

## Citation and access

