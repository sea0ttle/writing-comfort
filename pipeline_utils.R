# ============================================================
# pipeline_utils.R
#
# Thin loader: sources all modular utility files from R/.
# Sourced by every numbered script and by the QMD report.
# To add or modify helpers, edit the relevant R/utils_*.R file.
# ============================================================

.pipeline_utils_dir <- file.path(
  normalizePath(
    if (exists("project_dir")) project_dir else getwd(),
    mustWork = FALSE
  ),
  "R"
)

source(file.path(.pipeline_utils_dir, "utils_io.R"))
source(file.path(.pipeline_utils_dir, "utils_data.R"))
source(file.path(.pipeline_utils_dir, "utils_clustering.R"))
source(file.path(.pipeline_utils_dir, "utils_models.R"))
source(file.path(.pipeline_utils_dir, "utils_catalog.R"))

rm(.pipeline_utils_dir)
