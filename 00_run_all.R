# ============================================================
# 00_run_all.R
#
# Run the full writing-comfort analysis pipeline end-to-end.
# Place this file, the four numbered scripts, pipeline_utils.R,
# config.R, the R/ utility modules, the stan/ model, and the
# data/raw/ CSVs in the same project directory, then run:
#   source("00_run_all.R")
# ============================================================

project_dir <- if (exists("project_dir")) {
  normalizePath(project_dir, mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(project_dir, "config.R"))

run_step <- function(script_name) {
  step_env <- new.env(parent = globalenv())
  step_env$project_dir <- project_dir
  step_env$config      <- config
  sys.source(file.path(project_dir, script_name), envir = step_env)
}

run_step("01_clean_data_and_features.R")
run_step("02_maxdiff_hb_and_clusters.R")
run_step("03_glmm_review_and_validation.R")
run_step("04_final_model_visuals.R")
