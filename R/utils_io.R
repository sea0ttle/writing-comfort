# ============================================================
# R/utils_io.R
#
# Package management, file-path helpers, and I/O utilities.
# ============================================================

ensure_packages <- function(pkgs) {
  to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install)) install.packages(to_install)
  invisible(pkgs)
}

# ------------------------------------------------------------
# Path helpers
# ------------------------------------------------------------

make_pipeline_paths <- function(project_dir = getwd()) {
  project_dir  <- normalizePath(project_dir, mustWork = TRUE)
  raw_dir      <- file.path(project_dir, "data", "raw")
  output_root  <- file.path(project_dir, "analysis_pipeline_output")

  list(
    project_dir = project_dir,
    raw = list(
      combined_data = file.path(raw_dir, "combined_data_manual_clean_5.csv"),
      comfort_items = file.path(raw_dir, "Pen Survey Responses - comfort_items.csv"),
      jrg           = file.path(raw_dir, "jrg.csv"),
      responses     = file.path(raw_dir, "Pen Survey Responses - responses (6).csv")
    ),
    out = list(
      root    = output_root,
      clean   = file.path(output_root, "01_clean"),
      maxdiff = file.path(output_root, "02_maxdiff"),
      review  = file.path(output_root, "03_review"),
      visuals = file.path(output_root, "04_visuals")
    )
  )
}

make_output_dirs <- function(paths) {
  dirs <- unlist(paths$out, use.names = FALSE)
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

check_required_files <- function(paths) {
  raw_paths <- unlist(paths$raw, use.names = TRUE)
  missing   <- raw_paths[!file.exists(raw_paths)]
  if (length(missing)) {
    stop(
      "Missing raw input files:\n",
      paste(sprintf("  - %s", missing), collapse = "\n"),
      "\n\nExpected location: data/raw/",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
