# ============================================================
# 02_maxdiff_hb_and_clusters.R
#
# Purpose:
#   - Load the cleaned respondent data bundle
#   - Rebuild the MaxDiff tasks from the survey responses
#   - Fit the hierarchical Bayes MaxDiff model in Stan
#   - Cluster respondents from HB utilities
#   - Save respondent HB shares, HB clusters, and share summaries
# ============================================================

project_dir <- if (exists("project_dir")) {
  normalizePath(project_dir, mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(project_dir, "config.R"))
source(file.path(project_dir, "pipeline_utils.R"))

ensure_packages(c("rstan", "tidyverse", "cluster", "proxy", "scales"))

suppressPackageStartupMessages({
  library(rstan)
  library(tidyverse)
  library(cluster)
  library(proxy)
  library(scales)
})

set.seed(config$seed)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

paths <- make_pipeline_paths(project_dir)
make_output_dirs(paths)
clean_bundle_path <- file.path(paths$out$clean, "01_clean_data_bundle.rds")
if (!file.exists(clean_bundle_path)) {
  stop("Could not find clean data bundle. Run 01_clean_data_and_features.R first.")
}

clean_bundle <- readRDS(clean_bundle_path)
responses <- clean_bundle$responses_enriched

bw_design <- default_bw_design()
responses_bw <- derive_bw_mid_options(responses, bw_design = bw_design)
bw_build <- build_bw_tasks_table(responses_bw, bw_design = bw_design)
attrs <- bw_build$attrs

stan_data <- list(
  J = nrow(responses_bw),
  K = bw_build$K_attr,
  T = nrow(bw_build$tasks),
  best = bw_build$tasks$best,
  worst = bw_build$tasks$worst,
  mid = bw_build$tasks$mid,
  id = bw_build$tasks$id
)

hb_fit_path <- file.path(paths$out$maxdiff, "hb_fit.rds")
fit_hb <- NULL

if (config$stan$reuse_existing && file.exists(hb_fit_path)) {
  fit_hb <- tryCatch(readRDS(hb_fit_path), error = function(e) NULL)
}

if (is.null(fit_hb)) {
  hb_model <- stan_model(file = file.path(project_dir, "stan", "hb_maxdiff.stan"))
  fit_hb <- sampling(
    hb_model,
    data    = stan_data,
    chains  = config$stan$chains,
    iter    = config$stan$iter,
    warmup  = config$stan$warmup,
    refresh = config$stan$refresh
  )
  saveRDS(fit_hb, hb_fit_path)
}

hb_fit_summary <- as.data.frame(rstan::summary(fit_hb)$summary) %>%
  tibble::rownames_to_column("parameter") %>%
  as_tibble()

draws_u <- rstan::extract(fit_hb, "u")$u
u_mean <- apply(draws_u, c(2, 3), mean)
colnames(u_mean) <- attrs

hb_cluster_build <- make_hb_cluster_assignments(responses_bw, u_mean, seed = 812)
cluster_assignments <- hb_cluster_build$cluster_assignments
cluster_sizes <- cluster_assignments %>%
  distinct(respondent_id, cluster) %>%
  count(cluster, name = "n")

ci_probs <- c(0.025, 0.975)
thin <- 2
iter_idx <- seq(1, dim(draws_u)[1], by = thin)
S <- length(iter_idx)
C <- length(levels(cluster_assignments$cluster))
K <- dim(draws_u)[3]
cl_vec <- cluster_assignments %>%
  arrange(respondent_id) %>%
  pull(cluster) %>%
  as.factor()
clusters <- levels(cl_vec)
cluster_rows <- lapply(clusters, function(cc) which(cl_vec == cc))

centroid_draw <- array(NA_real_, dim = c(S, C, K))
for (si in seq_along(iter_idx)) {
  s <- iter_idx[si]
  U <- draws_u[s, , ]
  for (ci in seq_along(clusters)) {
    rows <- cluster_rows[[ci]]
    centroid_draw[si, ci, ] <- colMeans(U[rows, , drop = FALSE])
  }
}

share_draw <- array(NA_real_, dim = c(S, C, K))
for (si in seq_len(S)) {
  for (ci in seq_len(C)) {
    share_draw[si, ci, ] <- 100 * softmax(centroid_draw[si, ci, ])
  }
}

share_sum <- bind_rows(lapply(seq_len(C), function(ci) {
  mat <- share_draw[, ci, ]
  tibble(
    cluster = factor(clusters[ci], levels = clusters),
    Attribute = attrs,
    mean = colMeans(mat),
    q025 = apply(mat, 2, quantile, probs = ci_probs[1]),
    q975 = apply(mat, 2, quantile, probs = ci_probs[2])
  )
})) %>%
  standardize_ci_cols()

respondent_share_mat <- t(apply(u_mean, 1, function(x) 100 * softmax(x)))
respondent_hb_shares <- tibble(
  respondent_id = responses_bw$respondent_id,
  survey_id = responses_bw$survey_id,
  hb_share_balance = respondent_share_mat[, "Balance"],
  hb_share_material = respondent_share_mat[, "Barrel Material"],
  hb_share_diameter = respondent_share_mat[, "Grip Diameter"],
  hb_share_texture = respondent_share_mat[, "Grip Texture"],
  hb_share_length = respondent_share_mat[, "Length"],
  hb_share_smoothness = respondent_share_mat[, "Nib/tip Smoothness"],
  hb_share_weight = respondent_share_mat[, "Overall Weight"]
)

stated_importance_summary <- respondent_hb_shares %>%
  summarise(
    Balance = mean(hb_share_balance, na.rm = TRUE),
    `Barrel Material` = mean(hb_share_material, na.rm = TRUE),
    `Grip Diameter` = mean(hb_share_diameter, na.rm = TRUE),
    `Grip Texture` = mean(hb_share_texture, na.rm = TRUE),
    Length = mean(hb_share_length, na.rm = TRUE),
    `Nib/tip Smoothness` = mean(hb_share_smoothness, na.rm = TRUE),
    `Overall Weight` = mean(hb_share_weight, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "attribute", values_to = "mean_hb_share") %>%
  arrange(desc(mean_hb_share))

hb_fit_diagnostics <- tibble(
  n_respondents = nrow(responses_bw),
  n_tasks_used = nrow(bw_build$tasks),
  max_rhat = suppressWarnings(max(hb_fit_summary$Rhat, na.rm = TRUE)),
  min_n_eff = suppressWarnings(min(hb_fit_summary$n_eff, na.rm = TRUE))
)

maxdiff_bundle <- list(
  metadata = list(
    created_at         = Sys.time(),
    seed               = config$seed,
    reuse_existing_fit = config$stan$reuse_existing
  ),
  bw_design = bw_design,
  responses_bw = responses_bw,
  tasks = bw_build$tasks,
  stan_data = stan_data,
  hb_fit_diagnostics = hb_fit_diagnostics,
  hb_fit_summary = hb_fit_summary,
  u_mean = u_mean,
  cluster_assignments = cluster_assignments,
  cluster_sizes = cluster_sizes,
  km_final = hb_cluster_build$km_final,
  smooth_by_cluster = hb_cluster_build$smooth_by_cluster,
  share_sum = share_sum,
  respondent_hb_shares = respondent_hb_shares,
  stated_importance_summary = stated_importance_summary
)

write_csv(cluster_assignments, file.path(paths$out$maxdiff, "hb_cluster_assignments.csv"))
write_csv(cluster_sizes, file.path(paths$out$maxdiff, "hb_cluster_sizes.csv"))
write_csv(share_sum, file.path(paths$out$maxdiff, "hb_cluster_share_summary.csv"))
write_csv(respondent_hb_shares, file.path(paths$out$maxdiff, "respondent_hb_shares.csv"))
write_csv(stated_importance_summary, file.path(paths$out$maxdiff, "stated_importance_summary.csv"))
write_csv(hb_fit_diagnostics, file.path(paths$out$maxdiff, "hb_fit_diagnostics.csv"))
write_csv(hb_fit_summary, file.path(paths$out$maxdiff, "hb_fit_summary.csv"))
saveRDS(maxdiff_bundle, file.path(paths$out$maxdiff, "02_maxdiff_bundle.rds"))

cat("\n============================================================\n")
cat("MAXDIFF HB + RESPONDENT CLUSTERING COMPLETE\n")
cat("============================================================\n")
cat("HB diagnostics:\n")
print(hb_fit_diagnostics, n = nrow(hb_fit_diagnostics))
cat("\nStated importance summary:\n")
print(stated_importance_summary, n = nrow(stated_importance_summary))
cat("\nHB cluster sizes:\n")
print(cluster_sizes, n = nrow(cluster_sizes))
cat("\nOutputs written to:\n", paths$out$maxdiff, "\n", sep = "")
