# ============================================================
# config.R
#
# Central configuration for the writing-comfort pipeline.
# Sourced by 00_run_all.R and all numbered scripts.
# Change parameters here rather than in individual scripts.
# ============================================================

config <- list(

  # Reproducibility seed — used by all scripts that call set.seed()
  seed = 812,

  # Stan / HB MaxDiff sampling settings (script 02)
  stan = list(
    chains          = 4,
    iter            = 2000,
    warmup          = 1000,
    refresh         = 0,
    reuse_existing  = TRUE   # if TRUE, skip re-fitting if hb_fit.rds already exists
  ),

  # GLMM optimizer settings (script 03)
  glmm = list(
    optimizer = "bobyqa",
    max_fun   = 5e5
  )

)
