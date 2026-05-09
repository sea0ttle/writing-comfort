# ============================================================
# 03_glmm_review_and_validation.R
#
# Purpose:
#   - Combine cleaned item-level features, respondent covariates,
#     and HB / cluster outputs
#   - Produce descriptive statistics and sample diagnostics
#   - Fit and compare the core, extended, pen-only, random-slope,
#     threshold, and validation GLMM families
#   - Bridge MaxDiff stated importance to item-level evidence
#   - Save a single review bundle for downstream review and plotting
#
# Note: coef_glmm_simple, summarize_glmm_simple, and
#       extract_glmm_simple_coefs live in pipeline_utils.R.
# ============================================================

project_dir <- if (exists("project_dir")) {
  normalizePath(project_dir, mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(project_dir, "config.R"))
source(file.path(project_dir, "pipeline_utils.R"))

ensure_packages(c("tidyverse", "lme4", "performance", "pROC"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(performance)
})

set.seed(config$seed)

paths <- make_pipeline_paths(project_dir)
make_output_dirs(paths)
clean_bundle_path <- file.path(paths$out$clean, "01_clean_data_bundle.rds")
maxdiff_bundle_path <- file.path(paths$out$maxdiff, "02_maxdiff_bundle.rds")

if (!file.exists(clean_bundle_path)) {
  stop("Could not find 01_clean_data_bundle.rds. Run 01_clean_data_and_features.R first.")
}
if (!file.exists(maxdiff_bundle_path)) {
  stop("Could not find 02_maxdiff_bundle.rds. Run 02_maxdiff_hb_and_clusters.R first.")
}

clean_bundle <- readRDS(clean_bundle_path)
hb_bundle <- readRDS(maxdiff_bundle_path)

ctrl <- glmerControl(
  optimizer = config$glmm$optimizer,
  optCtrl   = list(maxfun = config$glmm$max_fun)
)

responses2 <- prepare_respondent_covariates(
  clean_bundle$responses_enriched,
  hb_bundle$cluster_assignments
)

features <- clean_bundle$features_2_final

df <- features %>%
  left_join(responses2, by = "survey_id") %>%
  mutate(
    comfort_rating = safe_num(comfort_rating),
    high = if_else(!is.na(comfort_rating) & comfort_rating >= 6, 1L, 0L, missing = NA_integer_),
    high5 = if_else(!is.na(comfort_rating) & comfort_rating >= 5, 1L, 0L, missing = NA_integer_),
    high6 = if_else(!is.na(comfort_rating) & comfort_rating >= 6, 1L, 0L, missing = NA_integer_),
    high7 = if_else(!is.na(comfort_rating) & comfort_rating >= 7, 1L, 0L, missing = NA_integer_),
    # Keep the current is_fountain logic to match the validated working pipeline.
    is_fountain = if_else(
      !is.na(instrument) & grepl("Fountain", instrument, ignore.case = TRUE),
      1L,
      0L,
      missing = 0L
    ),
    diameter = safe_num(med_diameter_grip),
    weight = safe_num(w_writing),
    length = safe_num(l_writing),
    cg_abs_used = safe_num(cg_est_used_raw),
    cg_rel_used = safe_num(cg_est_used_rel),
    torque_used = safe_num(used_25mm_torque_abs_n_m),
    is_pen_only = product_type %in% c(
      "Ballpoint Pens",
      "Gel Pens",
      "Rollerball Pens",
      "Fountain Pens"
    )
  ) %>%
  group_by(survey_id) %>%
  mutate(n_rated_by_respondent = n()) %>%
  ungroup()

df_glmm <- df %>%
  select(
    survey_id,
    cluster,
    comfort_rating,
    high,
    high5,
    high6,
    high7,
    press_force,
    hand_size,
    handwriting_size,
    balance_preference,
    grip_posture,
    longest_session,
    activity_cluster,
    types_cluster,
    instrument,
    is_fountain,
    diameter,
    weight,
    length,
    med_diameter,
    med_diameter_grip,
    med_diameter_grip_min,
    med_diameter_grip_max,
    l_writing,
    w_writing,
    cgu_est,
    cgp_est,
    cgpn_est,
    cgps_est,
    cg_abs_used,
    cg_rel_used,
    torque_used,
    starts_with("used_25mm_"),
    starts_with("uncapped_25mm_"),
    starts_with("posted_25mm_")
  )

overall_summary <- tibble(
  n_rows = nrow(df),
  n_respondents = n_distinct(df$survey_id),
  mean_items_per_respondent = mean(df$n_rated_by_respondent),
  median_items_per_respondent = median(df$n_rated_by_respondent),
  pct_high5 = mean(df$high5, na.rm = TRUE),
  pct_high6 = mean(df$high6, na.rm = TRUE),
  pct_high7 = mean(df$high7, na.rm = TRUE),
  pct_fountain = mean(df$is_fountain, na.rm = TRUE),
  pct_pen_only = mean(df$is_pen_only, na.rm = TRUE)
)

feature_coverage <- tibble(
  variable = c(
    "product_type", "total_records",
    "diameter", "weight", "length",
    "cg_abs_used", "cg_rel_used", "torque_used",
    "hand_size", "press_force", "handwriting_size",
    "balance_preference", "grip_posture", "longest_session"
  ),
  nonmissing = c(
    sum(!is.na(df$product_type)),
    sum(!is.na(df$total_records)),
    sum(!is.na(df$diameter)),
    sum(!is.na(df$weight)),
    sum(!is.na(df$length)),
    sum(!is.na(df$cg_abs_used)),
    sum(!is.na(df$cg_rel_used)),
    sum(!is.na(df$torque_used)),
    sum(!is.na(df$hand_size)),
    sum(!is.na(df$press_force)),
    sum(!is.na(df$handwriting_size)),
    sum(!is.na(df$balance_preference)),
    sum(!is.na(df$grip_posture)),
    sum(!is.na(df$longest_session))
  )
) %>%
  mutate(
    total_rows = nrow(df),
    pct_nonmissing = nonmissing / total_rows
  )

match_quality_summary <- df %>%
  summarise(
    pct_with_product_type = mean(!is.na(product_type)),
    pct_with_total_records = mean(!is.na(total_records)),
    median_total_records = median(total_records, na.rm = TRUE),
    p25_total_records = quantile(total_records, 0.25, na.rm = TRUE),
    p75_total_records = quantile(total_records, 0.75, na.rm = TRUE)
  )

product_type_counts <- df %>% count(product_type, sort = TRUE)

missing_by_product_type <- df %>%
  group_by(product_type) %>%
  summarise(
    n = n(),
    pct_missing_diameter = mean(is.na(diameter)),
    pct_missing_weight = mean(is.na(weight)),
    pct_missing_length = mean(is.na(length)),
    pct_missing_cg_rel = mean(is.na(cg_rel_used)),
    pct_missing_torque = mean(is.na(torque_used)),
    .groups = "drop"
  ) %>%
  arrange(desc(n))

respondent_item_count <- df %>%
  distinct(survey_id, n_rated_by_respondent) %>%
  count(n_rated_by_respondent, name = "n_respondents")

core_vars <- c(
  "survey_id", "high5", "high6", "high7",
  "diameter", "weight", "length",
  "cg_abs_used", "cg_rel_used", "torque_used",
  "hand_size", "is_fountain"
)

ext_vars <- c(
  core_vars,
  "press_force", "handwriting_size",
  "balance_preference", "grip_posture", "longest_session"
)

df_core <- df %>%
  filter(if_all(all_of(core_vars), ~ !is.na(.x))) %>%
  prep_model_data()

df_ext <- df %>%
  filter(if_all(all_of(ext_vars), ~ !is.na(.x))) %>%
  prep_model_data()

df_pen_core <- df %>%
  filter(is_pen_only) %>%
  filter(if_all(all_of(core_vars), ~ !is.na(.x))) %>%
  prep_model_data()

attrition_summary <- tibble(
  sample_name = c("all_rated_items", "core_same_sample", "extended_same_sample", "pen_only_core_same_sample"),
  n_rows = c(nrow(df), nrow(df_core), nrow(df_ext), nrow(df_pen_core)),
  n_respondents = c(
    n_distinct(df$survey_id),
    n_distinct(df_core$survey_id),
    n_distinct(df_ext$survey_id),
    n_distinct(df_pen_core$survey_id)
  ),
  pct_of_all_rows = c(
    1,
    nrow(df_core) / nrow(df),
    nrow(df_ext) / nrow(df),
    nrow(df_pen_core) / nrow(df)
  )
)

within_person_summary <- df_core %>%
  group_by(survey_id) %>%
  summarise(
    n_items = n(),
    has_within_person_outcome_variation = n_distinct(high6) > 1,
    has_within_person_diameter_variation = n_distinct(diameter) > 1,
    has_within_person_cg_rel_variation = n_distinct(cg_rel_used) > 1,
    has_within_person_fountain_variation = n_distinct(is_fountain) > 1,
    .groups = "drop"
  ) %>%
  summarise(
    respondents = n(),
    pct_gt_1_item = mean(n_items > 1),
    pct_with_outcome_variation = mean(has_within_person_outcome_variation),
    pct_with_diameter_variation = mean(has_within_person_diameter_variation),
    pct_with_cg_rel_variation = mean(has_within_person_cg_rel_variation),
    pct_with_fountain_variation = mean(has_within_person_fountain_variation)
  )

numeric_corr_core <- df_core %>%
  select(diameter, weight, length, cg_abs_used, cg_rel_used, torque_used, is_fountain) %>%
  cor(use = "pairwise.complete.obs")

numeric_corr_core_long <- as.data.frame(as.table(round(numeric_corr_core, 3)), stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rename(var1 = Var1, var2 = Var2, correlation = Freq)

core_rhs <- c(
  null = "1 + (1 | survey_id)",
  fountain_only = "is_fountain + (1 | survey_id)",
  diameter_only = "diameter_z + (1 | survey_id)",
  diameter_fountain = "diameter_z + is_fountain + (1 | survey_id)",
  diameter_hand_fountain = "diameter_z + hand_size + is_fountain + (1 | survey_id)",
  cg_rel_no_interaction = "diameter_z + hand_size + cg_rel_used_z + is_fountain + (1 | survey_id)",
  cg_rel_interaction = "diameter_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id)",
  cg_abs_interaction = "diameter_z + hand_size * cg_abs_used_z + is_fountain + (1 | survey_id)",
  torque_interaction = "diameter_z + hand_size * torque_used_z + is_fountain + (1 | survey_id)",
  full_geom_rel = "diameter_z + weight_z + length_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id)"
)

ext_rhs <- c(
  current_main = "diameter_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id)",
  plus_press_force = "diameter_z + hand_size * cg_rel_used_z + press_force + is_fountain + (1 | survey_id)",
  plus_handwriting = "diameter_z + hand_size * cg_rel_used_z + handwriting_size + is_fountain + (1 | survey_id)",
  plus_balance_pref = "diameter_z + hand_size * cg_rel_used_z + balance_preference + is_fountain + (1 | survey_id)",
  balance_pref_interaction = "diameter_z + hand_size + balance_preference * cg_rel_used_z + is_fountain + (1 | survey_id)",
  plus_grip_posture = "diameter_z + hand_size * cg_rel_used_z + grip_posture + is_fountain + (1 | survey_id)",
  plus_session = "diameter_z + hand_size * cg_rel_used_z + longest_session + is_fountain + (1 | survey_id)",
  original_press_model = "diameter_z * press_force + weight_z * press_force + is_fountain + (1 | survey_id)",
  all_controls = "diameter_z + weight_z + length_z + hand_size * cg_rel_used_z + press_force + handwriting_size + balance_preference + grip_posture + longest_session + is_fountain + (1 | survey_id)"
)

random_slope_rhs <- c(
  intercept_only = "diameter_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id)",
  diameter_slope = "diameter_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id) + (0 + diameter_z || survey_id)",
  cg_rel_slope = "diameter_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id) + (0 + cg_rel_used_z || survey_id)",
  both_slopes = "diameter_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id) + (0 + diameter_z || survey_id) + (0 + cg_rel_used_z || survey_id)"
)

threshold_rhs <- c(
  current_main = "diameter_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id)"
)

core_results <- fit_model_set(df_core, "high6", core_rhs, "core_same_sample", ctrl)
ext_results <- fit_model_set(df_ext, "high6", ext_rhs, "extended_same_sample", ctrl)
pen_results <- fit_model_set(df_pen_core, "high6", core_rhs, "pen_only_core_same_sample", ctrl)
random_slope_results <- fit_model_set(df_core, "high6", random_slope_rhs, "core_random_slopes", ctrl)

thresh5_results <- fit_model_set(df_core, "high5", threshold_rhs, "core_thresholds", ctrl)
thresh6_results <- fit_model_set(df_core, "high6", threshold_rhs, "core_thresholds", ctrl)
thresh7_results <- fit_model_set(df_core, "high7", threshold_rhs, "core_thresholds", ctrl)

threshold_comparison <- bind_rows(
  thresh5_results$comparison,
  thresh6_results$comparison,
  thresh7_results$comparison
)

threshold_coefficients <- bind_rows(
  thresh5_results$coefficients,
  thresh6_results$coefficients,
  thresh7_results$coefficients
)

core_lrt <- make_lrt_table(
  fits = core_results$fits,
  model_names = c(
    "null",
    "fountain_only",
    "diameter_fountain",
    "diameter_hand_fountain",
    "cg_rel_no_interaction",
    "cg_rel_interaction",
    "full_geom_rel"
  ),
  dataset_name = "core_same_sample",
  outcome_name = "high6",
  comparison_name = "nested_core_chain"
)

selected_collinearity <- bind_rows(
  extract_collinearity_tbl(core_results$fits$cg_rel_interaction$fit, "cg_rel_interaction", "core_same_sample", "high6"),
  extract_collinearity_tbl(core_results$fits$full_geom_rel$fit, "full_geom_rel", "core_same_sample", "high6"),
  extract_collinearity_tbl(ext_results$fits$all_controls$fit, "all_controls", "extended_same_sample", "high6")
)

best_models_by_family <- bind_rows(
  core_results$comparison %>% mutate(family = "core"),
  ext_results$comparison %>% mutate(family = "extended"),
  pen_results$comparison %>% mutate(family = "pen_only"),
  random_slope_results$comparison %>% mutate(family = "random_slopes")
) %>%
  filter(is.na(error)) %>%
  group_by(family) %>%
  slice_min(order_by = AIC, n = 3, with_ties = FALSE) %>%
  ungroup()

# ------------------------------------------------------------
# Fair matched-sample comparison: balance vs press vs weight/length variants
# ------------------------------------------------------------
cmp_vars <- c(
  "survey_id", "high6", "diameter", "weight", "length",
  "hand_size", "press_force", "cg_rel_used", "is_fountain"
)

df_cmp <- df %>%
  filter(if_all(all_of(cmp_vars), ~ !is.na(.x))) %>%
  prep_model_data()

m_cmp_balance <- glmer(
  high6 ~ diameter_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id),
  data = df_cmp,
  family = binomial,
  control = ctrl
)

m_cmp_press <- glmer(
  high6 ~ diameter_z * press_force + weight_z * press_force + is_fountain + (1 | survey_id),
  data = df_cmp,
  family = binomial,
  control = ctrl
)

m_cmp_weight_hand <- glmer(
  high6 ~ diameter_z + hand_size * weight_z + is_fountain + (1 | survey_id),
  data = df_cmp,
  family = binomial,
  control = ctrl
)

m_cmp_length_hand <- glmer(
  high6 ~ diameter_z + hand_size * length_z + is_fountain + (1 | survey_id),
  data = df_cmp,
  family = binomial,
  control = ctrl
)

m_cmp_balance_metrics <- binary_glmm_metrics(
  fit = m_cmp_balance,
  data = df_cmp,
  outcome_name = "high6",
  threshold = 0.5,
  re_form = NA
)

m_cmp_press_metrics <- binary_glmm_metrics(
  fit = m_cmp_press,
  data = df_cmp,
  outcome_name = "high6",
  threshold = 0.5,
  re_form = NA
)

m_cmp_weight_hand_metrics <- binary_glmm_metrics(
  fit = m_cmp_weight_hand,
  data = df_cmp,
  outcome_name = "high6",
  threshold = 0.5,
  re_form = NA
)

m_cmp_length_hand_metrics <- binary_glmm_metrics(
  fit = m_cmp_length_hand,
  data = df_cmp,
  outcome_name = "high6",
  threshold = 0.5,
  re_form = NA
)

matched_model_comparison <- bind_rows(
  tibble(
    model = "matched_current_main",
    n = nobs(m_cmp_balance),
    respondents = nlevels(getME(m_cmp_balance, "flist")[[1]]),
    AIC = AIC(m_cmp_balance),
    BIC = BIC(m_cmp_balance),
    logLik = as.numeric(logLik(m_cmp_balance)),
    auc = m_cmp_balance_metrics$auc,
    accuracy = m_cmp_balance_metrics$accuracy,
    misclassification_rate = m_cmp_balance_metrics$misclassification_rate,
    R2_marginal = safe_scalar(performance::r2_nakagawa(m_cmp_balance)$R2_marginal),
    R2_conditional = safe_scalar(performance::r2_nakagawa(m_cmp_balance)$R2_conditional)
  ),
  tibble(
    model = "matched_original_press_model",
    n = nobs(m_cmp_press),
    respondents = nlevels(getME(m_cmp_press, "flist")[[1]]),
    AIC = AIC(m_cmp_press),
    BIC = BIC(m_cmp_press),
    logLik = as.numeric(logLik(m_cmp_press)),
    auc = m_cmp_press_metrics$auc,
    accuracy = m_cmp_press_metrics$accuracy,
    misclassification_rate = m_cmp_press_metrics$misclassification_rate,
    R2_marginal = safe_scalar(performance::r2_nakagawa(m_cmp_press)$R2_marginal),
    R2_conditional = safe_scalar(performance::r2_nakagawa(m_cmp_press)$R2_conditional)
  ),
  tibble(
    model = "matched_weight_hand_interaction",
    n = nobs(m_cmp_weight_hand),
    respondents = nlevels(getME(m_cmp_weight_hand, "flist")[[1]]),
    AIC = AIC(m_cmp_weight_hand),
    BIC = BIC(m_cmp_weight_hand),
    logLik = as.numeric(logLik(m_cmp_weight_hand)),
    auc = m_cmp_weight_hand_metrics$auc,
    accuracy = m_cmp_weight_hand_metrics$accuracy,
    misclassification_rate = m_cmp_weight_hand_metrics$misclassification_rate,
    R2_marginal = safe_scalar(performance::r2_nakagawa(m_cmp_weight_hand)$R2_marginal),
    R2_conditional = safe_scalar(performance::r2_nakagawa(m_cmp_weight_hand)$R2_conditional)
  ),
  tibble(
    model = "matched_length_hand_interaction",
    n = nobs(m_cmp_length_hand),
    respondents = nlevels(getME(m_cmp_length_hand, "flist")[[1]]),
    AIC = AIC(m_cmp_length_hand),
    BIC = BIC(m_cmp_length_hand),
    logLik = as.numeric(logLik(m_cmp_length_hand)),
    auc = m_cmp_length_hand_metrics$auc,
    accuracy = m_cmp_length_hand_metrics$accuracy,
    misclassification_rate = m_cmp_length_hand_metrics$misclassification_rate,
    R2_marginal = safe_scalar(performance::r2_nakagawa(m_cmp_length_hand)$R2_marginal),
    R2_conditional = safe_scalar(performance::r2_nakagawa(m_cmp_length_hand)$R2_conditional)
  )
) %>%
  arrange(AIC)

matched_model_coefficients <- bind_rows(
  coef_glmm_simple(m_cmp_balance, "matched_current_main"),
  coef_glmm_simple(m_cmp_press, "matched_original_press_model")
)

# ------------------------------------------------------------
# Within-person / between-person decomposition
# ------------------------------------------------------------
wb_vars <- c("survey_id", "high6", "diameter", "cg_rel_used", "is_fountain", "hand_size")

df_wb <- df %>%
  filter(if_all(all_of(wb_vars), ~ !is.na(.x))) %>%
  mutate(hand_size = factor(hand_size, levels = c("Small", "Medium", "Large"))) %>%
  group_by(survey_id) %>%
  mutate(
    diameter_mean = mean(diameter, na.rm = TRUE),
    cg_rel_mean = mean(cg_rel_used, na.rm = TRUE),
    fountain_mean = mean(is_fountain, na.rm = TRUE),
    diameter_within = diameter - diameter_mean,
    cg_rel_within = cg_rel_used - cg_rel_mean,
    fountain_within = is_fountain - fountain_mean
  ) %>%
  ungroup() %>%
  mutate(
    diameter_within_z = as.numeric(scale(diameter_within)),
    diameter_between_z = as.numeric(scale(diameter_mean)),
    cg_rel_within_z = as.numeric(scale(cg_rel_within)),
    cg_rel_between_z = as.numeric(scale(cg_rel_mean)),
    fountain_within_z = as.numeric(scale(fountain_within)),
    fountain_between_z = as.numeric(scale(fountain_mean))
  )

m_within_between <- glmer(
  high6 ~
    diameter_within_z + diameter_between_z +
    hand_size * cg_rel_within_z + cg_rel_between_z +
    fountain_within_z + fountain_between_z +
    (1 | survey_id),
  data = df_wb,
  family = binomial,
  control = ctrl
)

m_within_between_metrics <- binary_glmm_metrics(
  fit = m_within_between,
  data = df_wb,
  outcome_name = "high6",
  threshold = 0.5,
  re_form = NA
)

within_between_summary <- tibble(
  model = "within_between_main",
  n = nobs(m_within_between),
  respondents = nlevels(getME(m_within_between, "flist")[[1]]),
  AIC = AIC(m_within_between),
  BIC = BIC(m_within_between),
  logLik = as.numeric(logLik(m_within_between)),
  auc = m_within_between_metrics$auc,
  accuracy = m_within_between_metrics$accuracy,
  misclassification_rate = m_within_between_metrics$misclassification_rate,
  R2_marginal = safe_scalar(performance::r2_nakagawa(m_within_between)$R2_marginal),
  R2_conditional = safe_scalar(performance::r2_nakagawa(m_within_between)$R2_conditional)
)

within_between_coefficients <- coef_glmm_simple(m_within_between, "within_between_main")

# ------------------------------------------------------------
# Product-type factor comparison
# ------------------------------------------------------------
pen_types_keep <- c("Ballpoint Pens", "Gel Pens", "Rollerball Pens", "Fountain Pens")

type_vars <- c("survey_id", "high6", "diameter", "cg_rel_used", "hand_size", "product_type")

df_type <- df %>%
  filter(product_type %in% pen_types_keep) %>%
  filter(if_all(all_of(type_vars), ~ !is.na(.x))) %>%
  mutate(
    is_fountain = if_else(product_type == "Fountain Pens", 1L, 0L),
    product_type = factor(product_type),
    product_type = relevel(product_type, ref = "Ballpoint Pens")
  ) %>%
  prep_model_data()

m_type_binary <- glmer(
  high6 ~ diameter_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id),
  data = df_type,
  family = binomial,
  control = ctrl
)

m_type_factor <- glmer(
  high6 ~ diameter_z + hand_size * cg_rel_used_z + product_type + (1 | survey_id),
  data = df_type,
  family = binomial,
  control = ctrl
)

m_type_binary_metrics <- binary_glmm_metrics(
  fit = m_type_binary,
  data = df_type,
  outcome_name = "high6",
  threshold = 0.5,
  re_form = NA
)

m_type_factor_metrics <- binary_glmm_metrics(
  fit = m_type_factor,
  data = df_type,
  outcome_name = "high6",
  threshold = 0.5,
  re_form = NA
)

product_type_comparison <- bind_rows(
  tibble(
    model = "binary_fountain_same_pen_sample",
    n = nobs(m_type_binary),
    respondents = nlevels(getME(m_type_binary, "flist")[[1]]),
    AIC = AIC(m_type_binary),
    BIC = BIC(m_type_binary),
    logLik = as.numeric(logLik(m_type_binary)),
    auc = m_type_binary_metrics$auc,
    accuracy = m_type_binary_metrics$accuracy,
    misclassification_rate = m_type_binary_metrics$misclassification_rate,
    R2_marginal = safe_scalar(performance::r2_nakagawa(m_type_binary)$R2_marginal),
    R2_conditional = safe_scalar(performance::r2_nakagawa(m_type_binary)$R2_conditional)
  ),
  tibble(
    model = "product_type_factor",
    n = nobs(m_type_factor),
    respondents = nlevels(getME(m_type_factor, "flist")[[1]]),
    AIC = AIC(m_type_factor),
    BIC = BIC(m_type_factor),
    logLik = as.numeric(logLik(m_type_factor)),
    auc = m_type_factor_metrics$auc,
    accuracy = m_type_factor_metrics$accuracy,
    misclassification_rate = m_type_factor_metrics$misclassification_rate,
    R2_marginal = safe_scalar(performance::r2_nakagawa(m_type_factor)$R2_marginal),
    R2_conditional = safe_scalar(performance::r2_nakagawa(m_type_factor)$R2_conditional)
  )
) %>%
  arrange(AIC)

product_type_coefficients <- bind_rows(
  coef_glmm_simple(m_type_binary, "binary_fountain_same_pen_sample"),
  coef_glmm_simple(m_type_factor, "product_type_factor")
)

# ------------------------------------------------------------
# Length / weight validation and stated-vs-revealed bridge
# ------------------------------------------------------------
df_lw <- df %>%
  left_join(hb_bundle$respondent_hb_shares, by = "survey_id") %>%
  filter(
    !is.na(survey_id),
    !is.na(high),
    !is.na(hand_size),
    !is.na(diameter),
    !is.na(weight),
    !is.na(length),
    !is.na(cg_rel_used),
    !is.na(is_fountain),
    !is.na(hb_share_weight),
    !is.na(hb_share_length),
    !is.na(hb_share_diameter),
    !is.na(hb_share_balance)
  ) %>%
  mutate(
    hand_size = factor(hand_size, levels = c("Small", "Medium", "Large")),
    diameter_z = as.numeric(scale(diameter)),
    weight_z = as.numeric(scale(weight)),
    length_z = as.numeric(scale(length)),
    cg_rel_used_z = as.numeric(scale(cg_rel_used)),
    hb_weight_z = as.numeric(scale(hb_share_weight)),
    hb_length_z = as.numeric(scale(hb_share_length)),
    hb_diameter_z = as.numeric(scale(hb_share_diameter)),
    hb_balance_z = as.numeric(scale(hb_share_balance))
  )

matched_sample_summary_lw <- tibble(
  n_rows = nrow(df_lw),
  n_respondents = n_distinct(df_lw$survey_id),
  pct_high = mean(df_lw$high, na.rm = TRUE),
  pct_fountain = mean(df_lw$is_fountain, na.rm = TRUE),
  mean_items_per_respondent = mean(table(df_lw$survey_id))
)

stated_importance_physical <- hb_bundle$stated_importance_summary %>%
  filter(attribute %in% c("Grip Diameter", "Overall Weight", "Length", "Balance")) %>%
  arrange(desc(mean_hb_share))

m_lw_base <- glmer(
  high ~ diameter_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

m_lw_weight <- glmer(
  high ~ diameter_z + hand_size * cg_rel_used_z + weight_z + is_fountain + (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

m_lw_length <- glmer(
  high ~ diameter_z + hand_size * cg_rel_used_z + length_z + is_fountain + (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

m_lw_weight_length <- glmer(
  high ~ diameter_z + hand_size * cg_rel_used_z + weight_z + length_z + is_fountain + (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

m_lw_weight_hand <- glmer(
  high ~ diameter_z + hand_size * cg_rel_used_z + hand_size * weight_z + is_fountain + (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

m_lw_length_hand <- glmer(
  high ~ diameter_z + hand_size * cg_rel_used_z + hand_size * length_z + is_fountain + (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

m_lw_weight_length_hand <- glmer(
  high ~ diameter_z + hand_size * cg_rel_used_z +
    hand_size * weight_z + hand_size * length_z +
    is_fountain + (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

lw_model_comparison <- bind_rows(
  summarize_glmm_simple(m_lw_base, df_lw, "high", "base_current_main", "incremental"),
  summarize_glmm_simple(m_lw_weight, df_lw, "high", "plus_weight", "incremental"),
  summarize_glmm_simple(m_lw_length, df_lw, "high", "plus_length", "incremental"),
  summarize_glmm_simple(m_lw_weight_length, df_lw, "high", "plus_weight_length", "incremental"),
  summarize_glmm_simple(m_lw_weight_hand, df_lw, "high", "plus_weight_hand_interaction", "incremental"),
  summarize_glmm_simple(m_lw_length_hand, df_lw, "high", "plus_length_hand_interaction", "incremental"),
  summarize_glmm_simple(m_lw_weight_length_hand, df_lw, "high", "plus_weight_length_hand_interaction", "incremental")
) %>%
  arrange(AIC)

lw_coefficients <- bind_rows(
  extract_glmm_simple_coefs(m_lw_base, "base_current_main", "incremental"),
  extract_glmm_simple_coefs(m_lw_weight, "plus_weight", "incremental"),
  extract_glmm_simple_coefs(m_lw_length, "plus_length", "incremental"),
  extract_glmm_simple_coefs(m_lw_weight_length, "plus_weight_length", "incremental"),
  extract_glmm_simple_coefs(m_lw_weight_hand, "plus_weight_hand_interaction", "incremental"),
  extract_glmm_simple_coefs(m_lw_length_hand, "plus_length_hand_interaction", "incremental"),
  extract_glmm_simple_coefs(m_lw_weight_length_hand, "plus_weight_length_hand_interaction", "incremental")
)

lw_collinearity <- bind_rows(
  extract_collinearity_tbl(m_lw_weight_length, "plus_weight_length", "incremental", "high6"),
  extract_collinearity_tbl(m_lw_weight_length_hand, "plus_weight_length_hand_interaction", "incremental", "high6")
)

lw_lrt_weight_chain <- tidy_lrt(
  c("base_current_main", "plus_weight", "plus_weight_length"),
  m_lw_base, m_lw_weight, m_lw_weight_length
) %>%
  mutate(chain = "weight_then_both", .before = 1)

lw_lrt_length_chain <- tidy_lrt(
  c("base_current_main", "plus_length", "plus_weight_length"),
  m_lw_base, m_lw_length, m_lw_weight_length
) %>%
  mutate(chain = "length_then_both", .before = 1)

df_lw_wb <- df_lw %>%
  group_by(survey_id) %>%
  mutate(
    diameter_mean = mean(diameter, na.rm = TRUE),
    weight_mean = mean(weight, na.rm = TRUE),
    length_mean = mean(length, na.rm = TRUE),
    cg_rel_mean = mean(cg_rel_used, na.rm = TRUE),
    fountain_mean = mean(is_fountain, na.rm = TRUE),
    diameter_within = diameter - diameter_mean,
    weight_within = weight - weight_mean,
    length_within = length - length_mean,
    cg_rel_within = cg_rel_used - cg_rel_mean,
    fountain_within = is_fountain - fountain_mean
  ) %>%
  ungroup() %>%
  mutate(
    diameter_within_z = as.numeric(scale(diameter_within)),
    diameter_between_z = as.numeric(scale(diameter_mean)),
    weight_within_z = as.numeric(scale(weight_within)),
    weight_between_z = as.numeric(scale(weight_mean)),
    length_within_z = as.numeric(scale(length_within)),
    length_between_z = as.numeric(scale(length_mean)),
    cg_rel_within_z = as.numeric(scale(cg_rel_within)),
    cg_rel_between_z = as.numeric(scale(cg_rel_mean)),
    fountain_within_z = as.numeric(scale(fountain_within)),
    fountain_between_z = as.numeric(scale(fountain_mean))
  )

m_wb_base <- glmer(
  high ~
    diameter_within_z + diameter_between_z +
    hand_size * cg_rel_within_z + cg_rel_between_z +
    fountain_within_z + fountain_between_z +
    (1 | survey_id),
  data = df_lw_wb,
  family = binomial,
  control = ctrl
)

m_wb_weight <- glmer(
  high ~
    diameter_within_z + diameter_between_z +
    weight_within_z + weight_between_z +
    hand_size * cg_rel_within_z + cg_rel_between_z +
    fountain_within_z + fountain_between_z +
    (1 | survey_id),
  data = df_lw_wb,
  family = binomial,
  control = ctrl
)

m_wb_length <- glmer(
  high ~
    diameter_within_z + diameter_between_z +
    length_within_z + length_between_z +
    hand_size * cg_rel_within_z + cg_rel_between_z +
    fountain_within_z + fountain_between_z +
    (1 | survey_id),
  data = df_lw_wb,
  family = binomial,
  control = ctrl
)

m_wb_weight_length <- glmer(
  high ~
    diameter_within_z + diameter_between_z +
    weight_within_z + weight_between_z +
    length_within_z + length_between_z +
    hand_size * cg_rel_within_z + cg_rel_between_z +
    fountain_within_z + fountain_between_z +
    (1 | survey_id),
  data = df_lw_wb,
  family = binomial,
  control = ctrl
)

lw_wb_model_comparison <- bind_rows(
  summarize_glmm_simple(m_wb_base, df_lw_wb, "high", "wb_base", "within_between"),
  summarize_glmm_simple(m_wb_weight, df_lw_wb, "high", "wb_plus_weight", "within_between"),
  summarize_glmm_simple(m_wb_length, df_lw_wb, "high", "wb_plus_length", "within_between"),
  summarize_glmm_simple(m_wb_weight_length, df_lw_wb, "high", "wb_plus_weight_length", "within_between")
) %>%
  arrange(AIC)

lw_wb_coefficients <- bind_rows(
  extract_glmm_simple_coefs(m_wb_base, "wb_base", "within_between"),
  extract_glmm_simple_coefs(m_wb_weight, "wb_plus_weight", "within_between"),
  extract_glmm_simple_coefs(m_wb_length, "wb_plus_length", "within_between"),
  extract_glmm_simple_coefs(m_wb_weight_length, "wb_plus_weight_length", "within_between")
)

lw_wb_lrt_weight_chain <- tidy_lrt(
  c("wb_base", "wb_plus_weight", "wb_plus_weight_length"),
  m_wb_base, m_wb_weight, m_wb_weight_length
) %>%
  mutate(chain = "wb_weight_then_both", .before = 1)

lw_wb_lrt_length_chain <- tidy_lrt(
  c("wb_base", "wb_plus_length", "wb_plus_weight_length"),
  m_wb_base, m_wb_length, m_wb_weight_length
) %>%
  mutate(chain = "wb_length_then_both", .before = 1)

m_bridge_base <- glmer(
  high ~ diameter_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

m_bridge_weight <- glmer(
  high ~ diameter_z + hand_size * cg_rel_used_z + is_fountain +
    weight_z + hb_weight_z + weight_z:hb_weight_z +
    (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

m_bridge_length <- glmer(
  high ~ diameter_z + hand_size * cg_rel_used_z + is_fountain +
    length_z + hb_length_z + length_z:hb_length_z +
    (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

m_bridge_weight_length <- glmer(
  high ~ diameter_z + hand_size * cg_rel_used_z + is_fountain +
    weight_z + hb_weight_z + weight_z:hb_weight_z +
    length_z + hb_length_z + length_z:hb_length_z +
    (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

hb_bridge_model_comparison <- bind_rows(
  summarize_glmm_simple(m_bridge_base, df_lw, "high", "bridge_base", "hb_bridge"),
  summarize_glmm_simple(m_bridge_weight, df_lw, "high", "bridge_weight", "hb_bridge"),
  summarize_glmm_simple(m_bridge_length, df_lw, "high", "bridge_length", "hb_bridge"),
  summarize_glmm_simple(m_bridge_weight_length, df_lw, "high", "bridge_weight_length", "hb_bridge")
) %>%
  arrange(AIC)

hb_bridge_coefficients <- bind_rows(
  extract_glmm_simple_coefs(m_bridge_base, "bridge_base", "hb_bridge"),
  extract_glmm_simple_coefs(m_bridge_weight, "bridge_weight", "hb_bridge"),
  extract_glmm_simple_coefs(m_bridge_length, "bridge_length", "hb_bridge"),
  extract_glmm_simple_coefs(m_bridge_weight_length, "bridge_weight_length", "hb_bridge")
)

hb_bridge_lrt_weight_chain <- tidy_lrt(
  c("bridge_base", "bridge_weight", "bridge_weight_length"),
  m_bridge_base, m_bridge_weight, m_bridge_weight_length
) %>%
  mutate(chain = "bridge_weight_then_both", .before = 1)

hb_bridge_lrt_length_chain <- tidy_lrt(
  c("bridge_base", "bridge_length", "bridge_weight_length"),
  m_bridge_base, m_bridge_length, m_bridge_weight_length
) %>%
  mutate(chain = "bridge_length_then_both", .before = 1)

m_attr_base <- glmer(
  high ~ hand_size + is_fountain + (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

m_attr_diameter <- glmer(
  high ~ hand_size + is_fountain + diameter_z + (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

m_attr_weight <- glmer(
  high ~ hand_size + is_fountain + weight_z + (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

m_attr_length <- glmer(
  high ~ hand_size + is_fountain + length_z + (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

m_attr_balance <- glmer(
  high ~ hand_size + is_fountain + cg_rel_used_z + hand_size:cg_rel_used_z + (1 | survey_id),
  data = df_lw,
  family = binomial,
  control = ctrl
)

get_attr_support <- function(attribute, fit, base_fit, coef_term) {
  an <- as.data.frame(anova(base_fit, fit, test = "Chisq"))
  sm <- as.data.frame(coef(summary(fit)))
  sm$term <- rownames(sm)
  rownames(sm) <- NULL
  this_coef <- sm %>% filter(term == coef_term)

  tibble(
    attribute = attribute,
    revealed_model = paste(deparse(formula(fit)), collapse = " "),
    delta_AIC_vs_base = AIC(fit) - AIC(base_fit),
    delta_logLik_vs_base = as.numeric(logLik(fit) - logLik(base_fit)),
    lrt_p_value = if (nrow(an) >= 2 && "Pr(>Chisq)" %in% names(an)) an$`Pr(>Chisq)`[2] else NA_real_,
    revealed_estimate = if (nrow(this_coef) == 1) this_coef$Estimate else NA_real_,
    revealed_p_value = if (nrow(this_coef) == 1) this_coef$`Pr(>|z|)` else NA_real_,
    revealed_odds_ratio = if (nrow(this_coef) == 1) exp(this_coef$Estimate) else NA_real_
  )
}

revealed_attribute_support <- bind_rows(
  get_attr_support("Grip Diameter", m_attr_diameter, m_attr_base, "diameter_z"),
  get_attr_support("Overall Weight", m_attr_weight, m_attr_base, "weight_z"),
  get_attr_support("Length", m_attr_length, m_attr_base, "length_z"),
  get_attr_support("Balance", m_attr_balance, m_attr_base, "cg_rel_used_z")
)

stated_vs_revealed_summary <- stated_importance_physical %>%
  left_join(revealed_attribute_support, by = "attribute") %>%
  arrange(desc(mean_hb_share))

key_terms_incremental <- lw_coefficients %>%
  filter(
    (model %in% c("plus_weight", "plus_weight_length", "plus_weight_hand_interaction") &
       term %in% c("weight_z", "hand_sizeMedium:weight_z", "hand_sizeLarge:weight_z")) |
      (model %in% c("plus_length", "plus_weight_length", "plus_length_hand_interaction") &
         term %in% c("length_z", "hand_sizeMedium:length_z", "hand_sizeLarge:length_z"))
  )

key_terms_within_between <- lw_wb_coefficients %>%
  filter(term %in% c(
    "diameter_within_z", "diameter_between_z",
    "weight_within_z", "weight_between_z",
    "length_within_z", "length_between_z",
    "cg_rel_within_z", "cg_rel_between_z",
    "hand_sizeMedium:cg_rel_within_z",
    "hand_sizeLarge:cg_rel_within_z"
  ))

key_terms_bridge <- hb_bridge_coefficients %>%
  filter(term %in% c(
    "weight_z", "hb_weight_z", "weight_z:hb_weight_z",
    "length_z", "hb_length_z", "length_z:hb_length_z"
  ))

selected_coefficients <- bind_rows(
  core_results$coefficients %>%
    filter(model %in% c("cg_rel_interaction", "cg_abs_interaction", "torque_interaction", "full_geom_rel")),
  ext_results$coefficients %>%
    filter(model %in% c("current_main", "balance_pref_interaction", "original_press_model", "all_controls")),
  pen_results$coefficients %>%
    filter(model == "cg_rel_interaction"),
  threshold_coefficients
)

final_model_name <- "cg_rel_interaction"
final_model <- core_results$fits[[final_model_name]]$fit
final_model_data <- df_core
model_selection_note <- paste(
  "Selected final model:",
  "high6 ~ diameter_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id).",
  "Rationale: best AIC within the core family, interpretable,",
  "and it outperformed the matched press-force alternative on the exact same sample."
)

review_bundle <- list(
  metadata = list(
    created_at = Sys.time(),
    project_dir = project_dir,
    seed = 812,
    model_selection_note = model_selection_note
  ),
  data = list(
    df = df,
    df_glmm = df_glmm,
    df_core = df_core,
    df_ext = df_ext,
    df_pen_core = df_pen_core,
    df_cmp = df_cmp,
    df_wb = df_wb,
    df_type = df_type,
    df_lw = df_lw,
    df_lw_wb = df_lw_wb
  ),
  descriptive = list(
    overall_summary = overall_summary,
    feature_coverage = feature_coverage,
    match_quality_summary = match_quality_summary,
    product_type_counts = product_type_counts,
    missing_by_product_type = missing_by_product_type,
    respondent_item_count = respondent_item_count,
    attrition_summary = attrition_summary,
    within_person_summary = within_person_summary,
    numeric_corr_core = numeric_corr_core,
    numeric_corr_core_long = numeric_corr_core_long,
    activity_cluster_counts = clean_bundle$activity_cluster_counts,
    types_cluster_counts = clean_bundle$types_cluster_counts,
    hb_cluster_sizes = hb_bundle$cluster_sizes,
    stated_importance_summary = hb_bundle$stated_importance_summary,
    proxy_cv_summary = clean_bundle$proxy_cv_summary
  ),
  model_sets = list(
    core = core_results,
    extended = ext_results,
    pen_only = pen_results,
    random_slopes = random_slope_results,
    thresholds = list(
      high5 = thresh5_results,
      high6 = thresh6_results,
      high7 = thresh7_results,
      comparison = threshold_comparison,
      coefficients = threshold_coefficients
    ),
    core_lrt = core_lrt,
    selected_collinearity = selected_collinearity,
    best_models_by_family = best_models_by_family,
    selected_coefficients = selected_coefficients
  ),
  validation = list(
    matched_model_comparison = matched_model_comparison,
    matched_model_coefficients = matched_model_coefficients,
    within_between_summary = within_between_summary,
    within_between_coefficients = within_between_coefficients,
    product_type_comparison = product_type_comparison,
    product_type_coefficients = product_type_coefficients,
    matched_sample_summary_lw = matched_sample_summary_lw,
    stated_importance_physical = stated_importance_physical,
    lw_model_comparison = lw_model_comparison,
    lw_coefficients = lw_coefficients,
    lw_collinearity = lw_collinearity,
    lw_lrt_weight_chain = lw_lrt_weight_chain,
    lw_lrt_length_chain = lw_lrt_length_chain,
    lw_wb_model_comparison = lw_wb_model_comparison,
    lw_wb_coefficients = lw_wb_coefficients,
    lw_wb_lrt_weight_chain = lw_wb_lrt_weight_chain,
    lw_wb_lrt_length_chain = lw_wb_lrt_length_chain,
    hb_bridge_model_comparison = hb_bridge_model_comparison,
    hb_bridge_coefficients = hb_bridge_coefficients,
    hb_bridge_lrt_weight_chain = hb_bridge_lrt_weight_chain,
    hb_bridge_lrt_length_chain = hb_bridge_lrt_length_chain,
    revealed_attribute_support = revealed_attribute_support,
    stated_vs_revealed_summary = stated_vs_revealed_summary,
    key_terms_incremental = key_terms_incremental,
    key_terms_within_between = key_terms_within_between,
    key_terms_bridge = key_terms_bridge
  ),
  selected_models = list(
    final_model_name = final_model_name,
    final_model = final_model,
    final_model_data = final_model_data,
    matched_current_main = m_cmp_balance,
    matched_original_press_model = m_cmp_press,
    within_between_main = m_within_between,
    product_type_binary = m_type_binary,
    product_type_factor = m_type_factor
  )
)

write_csv(overall_summary, file.path(paths$out$review, "overall_summary.csv"))
write_csv(feature_coverage, file.path(paths$out$review, "feature_coverage.csv"))
write_csv(match_quality_summary, file.path(paths$out$review, "match_quality_summary.csv"))
write_csv(attrition_summary, file.path(paths$out$review, "attrition_summary.csv"))
write_csv(core_results$comparison, file.path(paths$out$review, "core_model_comparison.csv"))
write_csv(core_results$coefficients, file.path(paths$out$review, "core_model_coefficients.csv"))
write_csv(ext_results$comparison, file.path(paths$out$review, "extended_model_comparison.csv"))
write_csv(pen_results$comparison, file.path(paths$out$review, "pen_only_model_comparison.csv"))
write_csv(random_slope_results$comparison, file.path(paths$out$review, "random_slope_model_comparison.csv"))
write_csv(threshold_comparison, file.path(paths$out$review, "threshold_model_comparison.csv"))
write_csv(matched_model_comparison, file.path(paths$out$review, "matched_model_comparison.csv"))
write_csv(within_between_summary, file.path(paths$out$review, "within_between_summary.csv"))
write_csv(product_type_comparison, file.path(paths$out$review, "product_type_comparison.csv"))
write_csv(lw_model_comparison, file.path(paths$out$review, "lw_model_comparison.csv"))
write_csv(lw_wb_model_comparison, file.path(paths$out$review, "lw_wb_model_comparison.csv"))
write_csv(hb_bridge_model_comparison, file.path(paths$out$review, "hb_bridge_model_comparison.csv"))
write_csv(stated_vs_revealed_summary, file.path(paths$out$review, "stated_vs_revealed_summary.csv"))
write_csv(selected_coefficients, file.path(paths$out$review, "selected_coefficients.csv"))
write_csv(selected_collinearity, file.path(paths$out$review, "selected_collinearity.csv"))
saveRDS(review_bundle, file.path(paths$out$review, "03_glmm_review_bundle.rds"))

cat("\n============================================================\n")
cat("GLMM REVIEW + VALIDATION COMPLETE\n")
cat("============================================================\n")
cat(model_selection_note, "\n\n")
cat("Overall analysis sample summary:\n")
print(as_tibble(overall_summary), n = Inf)
cat("\nTop core models by AIC:\n")
print(as_tibble(core_results$comparison %>% arrange(AIC)), n = nrow(core_results$comparison))
cat("\nMatched balance vs press comparison:\n")
print(as_tibble(matched_model_comparison), n = nrow(matched_model_comparison))
cat("\nStated vs revealed summary:\n")
print(as_tibble(stated_vs_revealed_summary), n = nrow(stated_vs_revealed_summary))
cat("\nOutputs written to:\n", paths$out$review, "\n", sep = "")
