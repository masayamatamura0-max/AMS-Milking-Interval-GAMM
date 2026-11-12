# ============================================================
# GAMM REANALYSIS
# ALL 11 SUPPLEMENTARY OUTCOMES
# ALL 6 CLUSTERING ANALYSES
# outcome-specific k from Supplemental File S1
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

library(mgcv)
library(emmeans)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)


# ============================================================
# 1. PATHS
# ============================================================

# Run from the repository root, or pass an explicit input/output directory.
args <- commandArgs(trailingOnly = TRUE)
input_folder <- if (length(args) >= 1) args[1] else file.path("data", "processed")
output_folder <- if (length(args) >= 2) args[2] else file.path("results", "GAMM_all_outcomes")
model_folder <- file.path(output_folder, "Models")
if (!dir.exists(input_folder)) stop("Input directory not found: ", input_folder,
  "\nSee README.md for required processed CSV inputs.")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(model_folder, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. GENERAL SETTINGS
# ============================================================

minimum_DIM <- 1
maximum_DIM <- 305

smooth_k <- NA_integer_ # Set individually for each outcome in analysis_targets

n_threads <- max(
  1,
  parallel::detectCores() - 1
)


# ============================================================
# 3. ECM SOURCE VARIABLES
# ============================================================

milk_variable <-
  "Sum_Cond_MilkYield"

fat_percent_variable <-
  "Mean_Daily_MdpFatPercentage"

protein_percent_variable <-
  "Mean_Daily_MdpProteinPercentage"


# ============================================================
# 4. CLUSTERING ANALYSES
# ============================================================

cluster_targets <- tibble::tribble(
  
  ~Cluster_analysis,
  ~File_name,
  ~Cluster_type,
  ~Window_end,
  
  "MI_1_21",
  "MI_1_21_daily_for_GAMM.csv",
  "Milking interval",
  21,
  
  "MI_1_42",
  "MI_1_42_daily_for_GAMM.csv",
  "Milking interval",
  42,
  
  "MI_1_63",
  "MI_1_63_daily_for_GAMM.csv",
  "Milking interval",
  63,
  
  "MIV_1_21",
  "MIV_1_21_daily_for_GAMM.csv",
  "Milking interval variability",
  21,
  
  "MIV_1_42",
  "MIV_1_42_daily_for_GAMM.csv",
  "Milking interval variability",
  42,
  
  "MIV_1_63",
  "MIV_1_63_daily_for_GAMM.csv",
  "Milking interval variability",
  63
)


# ============================================================
# 5. OUTCOMES
# ============================================================

# The input column names below are candidates, not verified raw-data schema.
# Override using config/outcome_columns.csv (Outcome,Source_column).
# The first three observed column names are grounded in the original ECM script.
# Other candidate names MUST be verified against the processed CSV headers.
outcome_targets <- tibble::tribble(
  ~Outcome, ~Outcome_label, ~Smooth_k, ~Source_candidates,
  "Milking_interval", "Mean milking interval (h)", 25L, "Mean_Daily_MilkInterval_hour",
  "Milking_frequency", "Milking frequency (no./d)", 25L, "Count_Milking",
  "Milk_yield", "Milk yield (kg/d)", 25L, "Sum_Cond_MilkYield",
  "ECM_kg_d", "Energy-corrected milk (kg/d)", 100L, "ECM_kg_d",
  "AMS_concentrate", "AMS concentrate delivered (kg/d)", 25L, "AMS_concentrate_delivered",
  "Milk_fat", "Milk fat (%)", 25L, "Mean_Daily_MdpFatPercentage",
  "Milk_protein", "Milk protein (%)", 50L, "Mean_Daily_MdpProteinPercentage",
  "Milk_lactose", "Milk lactose (%)", 25L, "Mean_Daily_MdpLactosePercentage",
  "Log10_SCC", "Log10 somatic cell count", 100L, "log10_Mean_Daily_MdpSCC",
  "Milk_temperature", "Milk temperature (C)", 25L, "Mean_Daily_MilkTemperature",
  "Rumination_time", "Rumination time (h/d)", 25L, "Mean_Daily_RuminationTime"
)

column_map_path <- file.path("config", "outcome_columns.csv")
column_map <- if (file.exists(column_map_path)) readr::read_csv(column_map_path,
  show_col_types = FALSE) else tibble::tibble(Outcome=character(),Source_column=character())
stopifnot(all(c("Outcome", "Source_column") %in% names(column_map)))

resolve_outcome <- function(dat, outcome_id, default_column) {
  custom <- column_map$Source_column[column_map$Outcome == outcome_id]
  if (length(custom) > 1) stop("Duplicate mapping for ", outcome_id)
  chosen <- if (length(custom) == 1 && nzchar(custom)) custom else default_column
  if (!(chosen %in% names(dat))) stop("Missing outcome column for ", outcome_id,
    ": ", chosen, ". Edit config/outcome_columns.csv to map your actual CSV header.")
  chosen
}

analysis_targets <- tidyr::crossing(cluster_targets, outcome_targets)

# ============================================================
# 6. HELPER FUNCTIONS
# ============================================================


# ------------------------------------------------------------
# AR.start
# ------------------------------------------------------------

prepare_AR_data <- function(dat){
  
  dat %>%
    
    arrange(
      lactation_segment_id,
      DIM_corrected
    ) %>%
    
    group_by(
      lactation_segment_id
    ) %>%
    
    mutate(
      
      AR_start =
        row_number() == 1 |
        is.na(
          lag(DIM_corrected)
        ) |
        (
          DIM_corrected -
            lag(DIM_corrected) != 1
        )
    ) %>%
    
    ungroup()
}


# ------------------------------------------------------------
# Estimate rho from preliminary residuals
# ------------------------------------------------------------

estimate_rho <- function(
    dat,
    residual_vector
){
  
  temp <- dat %>%
    
    mutate(
      residual =
        as.numeric(
          residual_vector
        )
    ) %>%
    
    arrange(
      lactation_segment_id,
      DIM_corrected
    ) %>%
    
    group_by(
      lactation_segment_id
    ) %>%
    
    mutate(
      
      previous_residual =
        lag(residual),
      
      previous_DIM =
        lag(DIM_corrected),
      
      consecutive =
        DIM_corrected -
        previous_DIM == 1
    ) %>%
    
    ungroup()
  
  
  keep <- (
    temp$consecutive &
      is.finite(temp$residual) &
      is.finite(temp$previous_residual)
  )
  
  
  if(
    sum(
      keep,
      na.rm = TRUE
    ) < 3
  ){
    
    return(
      NA_real_
    )
  }
  
  
  rho <- cor(
    
    temp$residual[keep],
    
    temp$previous_residual[keep],
    
    use = "complete.obs"
  )
  
  
  rho <- max(
    min(
      rho,
      0.95
    ),
    -0.95
  )
  
  
  rho
}


# ------------------------------------------------------------
# Extract P-value
# ------------------------------------------------------------

extract_parametric_p <- function(
    anova_object,
    term_name
){
  
  p_table <-
    anova_object$pTerms.table
  
  
  if(
    is.null(p_table) ||
    nrow(p_table) == 0
  ){
    
    return(
      NA_real_
    )
  }
  
  
  row_id <- which(
    rownames(p_table) ==
      term_name
  )
  
  
  if(length(row_id) == 0){
    
    return(
      NA_real_
    )
  }
  
  
  p_column <- grep(
    "p|Pr",
    colnames(p_table),
    ignore.case = TRUE,
    value = TRUE
  )
  
  
  if(length(p_column) == 0){
    
    return(
      NA_real_
    )
  }
  
  
  as.numeric(
    p_table[
      row_id[1],
      p_column[1]
    ]
  )
}


# ------------------------------------------------------------
# Post-AR1 lag-1 correlation by lactation
# ------------------------------------------------------------

calculate_lag1_by_lactation <- function(
    dat,
    residual_vector
){
  
  dat %>%
    
    mutate(
      .residual =
        as.numeric(
          residual_vector
        )
    ) %>%
    
    arrange(
      lactation_segment_id,
      DIM_corrected
    ) %>%
    
    group_by(
      lactation_segment_id
    ) %>%
    
    mutate(
      
      .previous_residual =
        lag(.residual),
      
      .previous_DIM =
        lag(DIM_corrected),
      
      .consecutive =
        DIM_corrected -
        .previous_DIM == 1
    ) %>%
    
    summarise(
      
      n_pairs = sum(
        
        .consecutive &
          is.finite(.residual) &
          is.finite(.previous_residual),
        
        na.rm = TRUE
      ),
      
      lag1_r = {
        
        keep <- (
          .consecutive &
            is.finite(.residual) &
            is.finite(.previous_residual)
        )
        
        if(
          sum(
            keep,
            na.rm = TRUE
          ) >= 3
        ){
          
          cor(
            
            .residual[keep],
            
            .previous_residual[keep],
            
            use = "complete.obs"
          )
          
        } else {
          
          NA_real_
        }
      },
      
      .groups = "drop"
    )
}


# ============================================================
# 7. RESULT OBJECTS
# ============================================================

all_main_results <-
  tibble()

all_emm_results <-
  tibble()

all_pairwise_results <-
  tibble()

all_model_information <-
  tibble()

all_k_check <-
  tibble()

all_smooth_terms <-
  tibble()

all_post_ar1 <-
  tibble()


# ============================================================
# 8. RUN 66 MODELS (6 clustering analyses x 11 outcomes)
# ============================================================

for(
  i in seq_len(
    nrow(
      analysis_targets
    )
  )
){
  
  
  cluster_analysis <-
    analysis_targets$Cluster_analysis[i]
  
  file_name <-
    analysis_targets$File_name[i]
  
  cluster_type <-
    analysis_targets$Cluster_type[i]
  
  window_end <-
    analysis_targets$Window_end[i]
  
  outcome <-
    analysis_targets$Outcome[i]
  
  outcome_label <-
    analysis_targets$Outcome_label[i]

  smooth_k <- analysis_targets$Smooth_k[i]
  source_candidates <- analysis_targets$Source_candidates[i]
  
  
  input_file <- file.path(
    input_folder,
    file_name
  )
  
  
  cat(
    
    "\n\n============================================================\n",
    
    "Analysis: ",
    cluster_analysis,
    "\n",
    
    "Outcome: ",
    outcome_label,
    "\n",
    
    "k: ",
    smooth_k,
    "\n",
    
    "============================================================\n",
    
    sep = ""
  )
  
  
  # ==========================================================
  # READ DATA
  # ==========================================================
  
  dat <- read_csv(
    input_file,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  
  # ==========================================================
  # CALCULATE ECM
  # ==========================================================
  
  ecm_variables <- c(
    
    milk_variable,
    
    fat_percent_variable,
    
    protein_percent_variable
  )
  
  
  missing_ecm_variables <- setdiff(
    
    ecm_variables,
    
    names(dat)
  )
  
  
  if(
    length(
      missing_ecm_variables
    ) > 0
  ){
    
    stop(
      
      paste0(
        
        "Cannot calculate ECM. Missing variables: ",
        
        paste(
          missing_ecm_variables,
          collapse = ", "
        )
      )
    )
  }
  
  
  dat[[milk_variable]] <-
    as.numeric(
      dat[[milk_variable]]
    )
  
  dat[[fat_percent_variable]] <-
    as.numeric(
      dat[[fat_percent_variable]]
    )
  
  dat[[protein_percent_variable]] <-
    as.numeric(
      dat[[protein_percent_variable]]
    )
  
  
  dat <- dat %>%
    
    mutate(
      
      Milk_kg_d =
        .data[[milk_variable]],
      
      Fat_kg_d =
        .data[[milk_variable]] *
        .data[[fat_percent_variable]] /
        100,
      
      Protein_kg_d =
        .data[[milk_variable]] *
        .data[[protein_percent_variable]] /
        100,
      
      ECM_kg_d =
        0.327 * Milk_kg_d +
        12.95 * Fat_kg_d +
        7.2 * Protein_kg_d
    )
  
  
  # Resolve all non-derived outcomes using a configurable raw-data column name.
  if (outcome != "ECM_kg_d") {
    source_column <- resolve_outcome(dat, outcome, source_candidates)
    dat[[outcome]] <- suppressWarnings(as.numeric(dat[[source_column]]))
  }

  # ==========================================================
  # CHECK VARIABLES
  # ==========================================================
  
  required_variables <- c(
    
    outcome,
    
    "Cluster",
    
    "Parity_group",
    
    "DIM_corrected",
    
    "herd_id",
    
    "lactation_segment_id"
  )
  
  
  missing_variables <- setdiff(
    
    required_variables,
    
    names(dat)
  )
  
  
  if(
    length(
      missing_variables
    ) > 0
  ){
    
    warning(
      
      "Missing variables in ",
      file_name,
      ": ",
      paste(
        missing_variables,
        collapse = ", "
      )
    )
    
    next
  }
  
  
  # ==========================================================
  # CLUSTER LEVELS
  # ==========================================================
  
  if(
    cluster_type ==
    "Milking interval"
  ){
    
    cluster_levels <- c(
      
      "Short interval",
      
      "Intermediate interval",
      
      "Long interval"
    )
    
  } else {
    
    cluster_levels <- c(
      
      "Consistent interval",
      
      "Intermediate variability",
      
      "Variable interval"
    )
  }
  
  
  dat$Cluster <- factor(
    dat$Cluster,
    levels = cluster_levels
  )
  
  
  dat$Parity_group <- factor(
    
    dat$Parity_group,
    
    levels = c(
      "Parity 1",
      "Parity 2+"
    )
  )
  
  
  dat$herd_id <-
    factor(
      dat$herd_id
    )
  
  
  dat$lactation_segment_id <-
    factor(
      dat$lactation_segment_id
    )
  
  
  dat$DIM_corrected <-
    as.numeric(
      dat$DIM_corrected
    )
  
  
  dat[[outcome]] <-
    as.numeric(
      dat[[outcome]]
    )
  
  
  # ==========================================================
  # ANALYSIS DATASET
  # ==========================================================
  
  dat2 <- dat %>%
    
    filter(
      
      DIM_corrected >=
        minimum_DIM,
      
      DIM_corrected <=
        maximum_DIM,
      
      is.finite(
        DIM_corrected
      ),
      
      is.finite(
        .data[[outcome]]
      ),
      
      !is.na(
        Cluster
      ),
      
      !is.na(
        Parity_group
      ),
      
      !is.na(
        herd_id
      ),
      
      !is.na(
        lactation_segment_id
      )
    ) %>%
    
    droplevels() %>%
    
    arrange(
      lactation_segment_id,
      DIM_corrected
    )
  
  
  cat(
    
    "Records: ",
    format(
      nrow(dat2),
      big.mark = ","
    ),
    "\n",
    
    "Lactations: ",
    format(
      n_distinct(
        dat2$lactation_segment_id
      ),
      big.mark = ","
    ),
    "\n",
    
    "Herds: ",
    format(
      n_distinct(
        dat2$herd_id
      ),
      big.mark = ","
    ),
    "\n",
    
    sep = ""
  )
  
  
  # ==========================================================
  # MODEL FORMULA
  # ==========================================================
  
  model_formula <- as.formula(
    
    paste0(
      
      outcome,
      
      " ~ Cluster * Parity_group + ",
      
      "s(DIM_corrected, by = Parity_group, ",
      
      "k = smooth_k, bs = 'cr', id = 1) + ",
      
      "s(herd_id, bs = 're') + ",
      
      "s(lactation_segment_id, bs = 're')"
    )
  )
  
  
  # ==========================================================
  # PRELIMINARY MODEL
  # ==========================================================
  
  cat(
    "Fitting preliminary model...\n"
  )
  
  
  preliminary_model <- bam(
    
    formula =
      model_formula,
    
    data =
      dat2,
    
    family =
      gaussian(
        link = "identity"
      ),
    
    method =
      "fREML",
    
    discrete =
      TRUE,
    
    nthreads =
      n_threads,
    
    drop.unused.levels =
      TRUE,
    
    gc.level =
      1
  )
  
  
  # ==========================================================
  # ESTIMATE RHO
  # ==========================================================
  
  preliminary_residuals <- residuals(
    
    preliminary_model,
    
    type = "pearson"
  )
  
  
  rho_est <- estimate_rho(
    
    dat2,
    
    preliminary_residuals
  )
  
  
  cat(
    "Estimated rho: ",
    round(
      rho_est,
      4
    ),
    "\n"
  )
  
  
  rm(
    preliminary_model,
    preliminary_residuals
  )
  
  invisible(
    gc()
  )
  
  
  # ==========================================================
  # AR.start
  # ==========================================================
  
  dat2 <-
    prepare_AR_data(
      dat2
    )
  
  
  # ==========================================================
  # FINAL AR(1) MODEL
  # ==========================================================
  
  cat(
    "Fitting final AR(1) model...\n"
  )
  
  
  final_model <- bam(
    
    formula =
      model_formula,
    
    data =
      dat2,
    
    family =
      gaussian(
        link = "identity"
      ),
    
    method =
      "fREML",
    
    discrete =
      TRUE,
    
    rho =
      rho_est,
    
    AR.start =
      dat2$AR_start,
    
    nthreads =
      n_threads,
    
    drop.unused.levels =
      TRUE,
    
    gc.level =
      1
  )
  
  
  # ==========================================================
  # P-VALUES
  # ==========================================================
  
  model_anova <- anova(
    final_model,
    freq = TRUE
  )
  
  
  cluster_p <-
    extract_parametric_p(
      model_anova,
      "Cluster"
    )
  
  
  parity_p <-
    extract_parametric_p(
      model_anova,
      "Parity_group"
    )
  
  
  interaction_p <-
    extract_parametric_p(
      model_anova,
      "Cluster:Parity_group"
    )
  
  
  # ==========================================================
  # EMM
  # ==========================================================
  
  emm <- emmeans(
    
    final_model,
    
    ~ Cluster |
      Parity_group,
    
    at = list(
      
      DIM_corrected =
        seq(
          1,
          305,
          by = 1
        )
    ),
    
    weights =
      "equal"
  )
  
  
  emm_df <- as.data.frame(
    emm
  ) %>%
    
    mutate(
      
      Cluster_analysis =
        cluster_analysis,
      
      Cluster_type =
        cluster_type,
      
      Clustering_window =
        paste0(
          "DIM 1-",
          window_end
        ),
      
      Outcome =
        outcome,
      
      Outcome_label =
        outcome_label,
      
      Smooth_k =
        smooth_k
    )
  
  
  all_emm_results <-
    bind_rows(
      all_emm_results,
      emm_df
    )
  
  
  # ==========================================================
  # TUKEY
  # ==========================================================
  
  pairwise_df <- as.data.frame(
    
    pairs(
      emm,
      adjust = "tukey"
    )
  ) %>%
    
    mutate(
      
      Cluster_analysis =
        cluster_analysis,
      
      Cluster_type =
        cluster_type,
      
      Clustering_window =
        paste0(
          "DIM 1-",
          window_end
        ),
      
      Outcome =
        outcome,
      
      Outcome_label =
        outcome_label,
      
      Smooth_k =
        smooth_k
    )
  
  
  all_pairwise_results <-
    bind_rows(
      all_pairwise_results,
      pairwise_df
    )
  
  
  # ==========================================================
  # MAIN RESULTS
  # ==========================================================
  
  main_df <- emm_df %>%
    
    mutate(
      
      Cluster_p =
        cluster_p,
      
      Parity_p =
        parity_p,
      
      Interaction_p =
        interaction_p,
      
      Rho_used =
        rho_est,
      
      N_records =
        nrow(
          dat2
        ),
      
      N_lactations =
        n_distinct(
          dat2$lactation_segment_id
        ),
      
      N_herds =
        n_distinct(
          dat2$herd_id
        )
    )
  
  
  all_main_results <-
    bind_rows(
      all_main_results,
      main_df
    )
  
  
  # ==========================================================
  # MODEL INFORMATION
  # ==========================================================
  
  model_summary <-
    summary(
      final_model
    )
  
  
  model_info <- tibble(
    
    Cluster_analysis =
      cluster_analysis,
    
    Cluster_type =
      cluster_type,
    
    Clustering_window =
      paste0(
        "DIM 1-",
        window_end
      ),
    
    Outcome =
      outcome,
    
    Outcome_label =
      outcome_label,
    
    Smooth_k =
      smooth_k,
    
    N_records =
      nrow(
        dat2
      ),
    
    N_lactations =
      n_distinct(
        dat2$lactation_segment_id
      ),
    
    N_herds =
      n_distinct(
        dat2$herd_id
      ),
    
    Rho_used =
      rho_est,
    
    AIC =
      AIC(
        final_model
      ),
    
    BIC =
      BIC(
        final_model
      ),
    
    Deviance_explained =
      model_summary$dev.expl,
    
    Adjusted_R2 =
      model_summary$r.sq,
    
    Cluster_p =
      cluster_p,
    
    Parity_p =
      parity_p,
    
    Interaction_p =
      interaction_p
  )
  
  
  all_model_information <-
    bind_rows(
      all_model_information,
      model_info
    )
  
  
  # ==========================================================
  # SMOOTH TERMS
  # ==========================================================
  
  smooth_df <- as.data.frame(
    model_summary$s.table
  ) %>%
    
    rownames_to_column(
      "Smooth_term"
    ) %>%
    
    mutate(
      
      Cluster_analysis =
        cluster_analysis,
      
      Cluster_type =
        cluster_type,
      
      Clustering_window =
        paste0(
          "DIM 1-",
          window_end
        ),
      
      Outcome =
        outcome,
      
      Outcome_label =
        outcome_label,
      
      Smooth_k =
        smooth_k
    )
  
  
  all_smooth_terms <-
    bind_rows(
      all_smooth_terms,
      smooth_df
    )
  
  
  set.seed(2026 + i)

  # ==========================================================
  # k.check
  # ==========================================================
  
  k_df <- as.data.frame(
    
    mgcv::k.check(
      final_model
    )
  ) %>%
    
    rownames_to_column(
      "Smooth_term"
    ) %>%
    
    filter(
      
      grepl(
        "DIM_corrected",
        Smooth_term,
        fixed = TRUE
      )
    ) %>%
    
    mutate(
      
      Cluster_analysis =
        cluster_analysis,
      
      Cluster_type =
        cluster_type,
      
      Clustering_window =
        paste0(
          "DIM 1-",
          window_end
        ),
      
      Outcome =
        outcome,
      
      Outcome_label =
        outcome_label,
      
      Smooth_k =
        smooth_k
    )
  
  
  all_k_check <-
    bind_rows(
      all_k_check,
      k_df
    )
  
  
  # ==========================================================
  # POST-AR1 AUTOCORRELATION
  # ==========================================================
  
  if(
    !is.null(
      final_model$std.rsd
    )
  ){
    
    final_residuals <-
      as.numeric(
        final_model$std.rsd
      )
    
  } else {
    
    final_residuals <- residuals(
      
      final_model,
      
      type = "pearson"
    )
  }
  
  
  post_pooled_r <- estimate_rho(
    
    dat2,
    
    final_residuals
  )
  
  
  lag1_by_lactation <-
    calculate_lag1_by_lactation(
      
      dat2,
      
      final_residuals
    )
  
  
  valid_lactations <-
    lag1_by_lactation %>%
    
    filter(
      is.finite(
        lag1_r
      )
    )
  
  
  post_ar1_df <- tibble(
    
    Cluster_analysis =
      cluster_analysis,
    
    Cluster_type =
      cluster_type,
    
    Clustering_window =
      paste0(
        "DIM 1-",
        window_end
      ),
    
    Outcome =
      outcome,
    
    Outcome_label =
      outcome_label,
    
    Smooth_k =
      smooth_k,
    
    Rho_used =
      rho_est,
    
    Post_AR1_pooled_r =
      post_pooled_r,
    
    Post_AR1_mean_lactation_r =
      mean(
        valid_lactations$lag1_r,
        na.rm = TRUE
      ),
    
    Post_AR1_median_lactation_r =
      median(
        valid_lactations$lag1_r,
        na.rm = TRUE
      ),
    
    N_lactations_with_lag1 =
      nrow(
        valid_lactations
      )
  )
  
  
  all_post_ar1 <-
    bind_rows(
      all_post_ar1,
      post_ar1_df
    )
  
  
  # ==========================================================
  # PRINT k.check
  # ==========================================================
  
  cat(
    "\nk.check results:\n"
  )
  
  print(
    k_df
  )
  
  
  # ==========================================================
  # CLEAN MEMORY
  # ==========================================================
  
  rm(
    dat,
    dat2,
    final_model,
    final_residuals,
    model_summary
  )
  
  invisible(
    gc()
  )
}


# ============================================================
# 9. SAVE COMBINED RESULTS
# ============================================================

write_csv(
  all_main_results,
  file.path(
    output_folder,
    "ALL_OUTCOMES_main_results_by_outcome.csv"
  )
)

write_csv(
  all_emm_results,
  file.path(
    output_folder,
    "ALL_OUTCOMES_EMM_by_outcome.csv"
  )
)

write_csv(
  all_pairwise_results,
  file.path(
    output_folder,
    "ALL_OUTCOMES_pairwise_Tukey_by_outcome.csv"
  )
)

write_csv(
  all_model_information,
  file.path(
    output_folder,
    "ALL_OUTCOMES_model_information_by_outcome.csv"
  )
)

write_csv(
  all_smooth_terms,
  file.path(
    output_folder,
    "ALL_OUTCOMES_smooth_terms_by_outcome.csv"
  )
)

write_csv(
  all_k_check,
  file.path(
    output_folder,
    "ALL_OUTCOMES_k_check_by_outcome.csv"
  )
)

write_csv(
  all_post_ar1,
  file.path(
    output_folder,
    "ALL_OUTCOMES_post_AR1_autocorrelation_by_outcome.csv"
  )
)


# ============================================================
# 10. FINISHED
# ============================================================
cat("\nAll outcome models completed. Results saved in: ", output_folder, "\n")
