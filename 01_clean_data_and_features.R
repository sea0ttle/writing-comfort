# ============================================================
# 01_clean_data_and_features.R
#
# Purpose:
#   - Read the raw product, comfort-item, JRG, and response files
#   - Build the product-side feature lookup
#   - Fit and cross-validate the CG proxy models
#   - Construct features_2_final with CG / torque estimates
#   - Build respondent-side activity and type clusters
#   - Save a single clean-data bundle for downstream scripts
#
# Expected raw files in project_dir:
#   - combined_data_manual_clean_5.csv
#   - Pen Survey Responses - comfort_items.csv
#   - jrg.csv
#   - Pen Survey Responses - responses (6).csv
# ============================================================

project_dir <- if (exists("project_dir")) {
  normalizePath(project_dir, mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(project_dir, "config.R"))
source(file.path(project_dir, "pipeline_utils.R"))

ensure_packages(c("tidyverse", "readr", "glmnet", "cluster", "proxy"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(glmnet)
  library(cluster)
  library(proxy)
})

set.seed(config$seed)

paths <- make_pipeline_paths(project_dir)
make_output_dirs(paths)
check_required_files(paths)

products_raw <- read.csv(paths$raw$combined_data, na.strings = c("", "NA"))
comfort_items_raw <- read.csv(paths$raw$comfort_items, na.strings = c("", "NA"))
responses_raw <- read_csv(paths$raw$responses, show_col_types = FALSE)
jrg <- read_csv(paths$raw$jrg, show_col_types = FALSE) %>%
  mutate(
    across(
      all_of(c("LPn", "LCp", "LPs", "DBa", "D25", "WPn", "WCp", "CGu", "CGp", "CgPn", "CgPs")),
      as.numeric
    )
  )

product_lookup <- build_product_feature_lookup(products_raw)
proxy_models <- fit_jrg_proxy_models(jrg)
proxy_cv_summary <- summarize_jrg_proxy_cv(jrg)
feature_build <- build_features_2_final(comfort_items_raw, product_lookup, proxy_models)
responses_build <- build_responses_enriched(responses_raw)

responses_enriched <- responses_build$responses_enriched
features_2_final <- feature_build$features_2_final

activity_cluster_counts <- responses_enriched %>%
  count(activity_cluster, activity_cluster_label, activity_cluster_name, sort = TRUE)

types_cluster_counts <- responses_enriched %>%
  count(types_cluster, types_cluster_label, types_cluster_name, sort = TRUE)

feature_missing_summary <- tibble(
  variable = c(
    "product_type", "total_records", "l_writing", "w_writing",
    "med_diameter_grip", "med_diameter", "cg_est_used_raw", "cg_est_used_rel",
    "used_25mm_torque_abs_n_m"
  ),
  nonmissing = c(
    sum(!is.na(features_2_final$product_type)),
    sum(!is.na(features_2_final$total_records)),
    sum(!is.na(features_2_final$l_writing)),
    sum(!is.na(features_2_final$w_writing)),
    sum(!is.na(features_2_final$med_diameter_grip)),
    sum(!is.na(features_2_final$med_diameter)),
    sum(!is.na(features_2_final$cg_est_used_raw)),
    sum(!is.na(features_2_final$cg_est_used_rel)),
    sum(!is.na(features_2_final$used_25mm_torque_abs_n_m))
  )
) %>%
  mutate(
    total_rows = nrow(features_2_final),
    pct_nonmissing = nonmissing / total_rows
  )

clean_bundle <- list(
  metadata = list(
    created_at = Sys.time(),
    project_dir = project_dir,
    seed = 812
  ),
  inputs = list(
    products_raw = products_raw,
    comfort_items_raw = comfort_items_raw,
    responses_raw = responses_raw,
    jrg = jrg
  ),
  product_lookup = product_lookup,
  proxy_models = proxy_models,
  proxy_cv_summary = proxy_cv_summary,
  joined_ratings = feature_build$joined_ratings,
  features_2 = feature_build$features_2,
  features_2_pred = feature_build$features_2_pred,
  features_2_final = features_2_final,
  responses_enriched = responses_enriched,
  activities_fit = responses_build$activities_fit,
  types_fit = responses_build$types_fit,
  activity_cluster_counts = activity_cluster_counts,
  types_cluster_counts = types_cluster_counts,
  feature_missing_summary = feature_missing_summary
)

write_csv(features_2_final, file.path(paths$out$clean, "features_2_final.csv"))
write_csv(responses_enriched, file.path(paths$out$clean, "responses_enriched.csv"))
write_csv(proxy_cv_summary, file.path(paths$out$clean, "proxy_cv_summary.csv"))
write_csv(activity_cluster_counts, file.path(paths$out$clean, "activity_cluster_counts.csv"))
write_csv(types_cluster_counts, file.path(paths$out$clean, "types_cluster_counts.csv"))
write_csv(feature_missing_summary, file.path(paths$out$clean, "feature_missing_summary.csv"))

saveRDS(proxy_models$fit_cgu, file.path(paths$out$clean, "model_cgu_linear.rds"))
saveRDS(proxy_models$fit_cgp, file.path(paths$out$clean, "model_cgp_ridge.rds"))
saveRDS(proxy_models$fit_cgpn, file.path(paths$out$clean, "model_cgpn_linear_no_wcp.rds"))
saveRDS(proxy_models$fit_cgps, file.path(paths$out$clean, "model_cgps_linear_with_wcp.rds"))
saveRDS(clean_bundle, file.path(paths$out$clean, "01_clean_data_bundle.rds"))

cat("\n============================================================\n")
cat("CLEAN DATA + FEATURE BUILD COMPLETE\n")
cat("============================================================\n")
cat("Project directory:\n", project_dir, "\n", sep = "")
cat("\nProduct lookup rows:", nrow(product_lookup), "\n")
cat("Features rows:", nrow(features_2_final), "\n")
cat("Respondents:", nrow(responses_enriched), "\n")
cat("\nCG proxy cross-validation summary:\n")
print(proxy_cv_summary, n = nrow(proxy_cv_summary))
cat("\nFeature missingness summary:\n")
print(feature_missing_summary, n = nrow(feature_missing_summary))
cat("\nOutputs written to:\n", paths$out$clean, "\n", sep = "")
