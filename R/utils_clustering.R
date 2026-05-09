# ============================================================
# R/utils_clustering.R
#
# Clustering and MaxDiff helpers for the writing-comfort
# analysis pipeline.
#
# Sections:
#   1. ggplot2 / tidytext helpers
#   2. MaxDiff / HB helpers
#   3. MaxDiff task construction
#   4. HB cluster assignments
#
# NOTE: The Stan model for HB MaxDiff has been extracted to
#   stan/hb_maxdiff.stan
# and should NOT be inlined here.  Load it in the calling
# script with:
#   stan_model(file = file.path(project_dir, "stan", "hb_maxdiff.stan"))
# ============================================================


# ------------------------------------------------------------
# ggplot2 / tidytext helpers
# ------------------------------------------------------------

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(labels = function(x) gsub(paste0(sep, ".*$"), "", x), ...)
}

standardize_ci_cols <- function(df) {
  nms <- names(df)
  mean_col <- intersect(nms, c("Mean", "mean", "avg", "estimate"))[1]
  lo_col <- intersect(nms, c("CI_low", "ci_low", "low", "lwr", "q025", "p025", "lower"))[1]
  hi_col <- intersect(nms, c("CI_high", "ci_high", "high", "upr", "q975", "p975", "upper"))[1]

  if (any(is.na(c(mean_col, lo_col, hi_col)))) {
    stop(
      "Could not find mean/CI columns. Found columns: ",
      paste(nms, collapse = ", ")
    )
  }

  df %>%
    dplyr::rename(
      Mean = dplyr::all_of(mean_col),
      CI_low = dplyr::all_of(lo_col),
      CI_high = dplyr::all_of(hi_col)
    )
}


# ------------------------------------------------------------
# MaxDiff / HB helpers
# ------------------------------------------------------------

softmax <- function(x) {
  ex <- exp(x - max(x))
  ex / sum(ex)
}

default_bw_design <- function() {
  list(
    A = c("Balance", "Barrel Material", "Nib/tip Smoothness"),
    B = c("Grip Diameter", "Overall Weight", "Balance"),
    C = c("Grip Diameter", "Nib/tip Smoothness", "Grip Texture"),
    D = c("Overall Weight", "Length", "Nib/tip Smoothness"),
    E = c("Grip Diameter", "Length", "Barrel Material"),
    F = c("Balance", "Length", "Grip Texture"),
    G = c("Overall Weight", "Barrel Material", "Grip Texture")
  )
}

build_multiselect_matrix <- function(data, id_col, value_col, sep = ";", min_prevalence = 2) {
  id_sym <- rlang::sym(id_col)
  value_sym <- rlang::sym(value_col)

  long_df <- data %>%
    dplyr::select(!!id_sym, !!value_sym) %>%
    dplyr::mutate(
      !!value_sym := dplyr::coalesce(!!value_sym, ""),
      !!value_sym := stringr::str_replace_all(!!value_sym, "\\s*;\\s*", ";")
    ) %>%
    tidyr::separate_rows(!!value_sym, sep = sep) %>%
    dplyr::mutate(
      !!value_sym := stringr::str_trim(!!value_sym),
      !!value_sym := dplyr::na_if(!!value_sym, ""),
      !!value_sym := dplyr::na_if(!!value_sym, "NA")
    ) %>%
    dplyr::filter(!is.na(!!value_sym)) %>%
    dplyr::distinct(!!id_sym, !!value_sym)

  matrix_df <- long_df %>%
    dplyr::mutate(value = 1L) %>%
    tidyr::pivot_wider(
      names_from = !!value_sym,
      values_from = value,
      values_fill = 0
    ) %>%
    dplyr::arrange(!!id_sym)

  respondent_ids <- matrix_df %>% dplyr::pull(!!id_sym)

  X <- matrix_df %>%
    dplyr::select(-!!id_sym) %>%
    as.matrix()

  storage.mode(X) <- "numeric"

  keep_cols <- colSums(X) >= min_prevalence
  X <- X[, keep_cols, drop = FALSE]

  list(
    long_df = long_df,
    matrix_df = matrix_df,
    respondent_ids = respondent_ids,
    X = X
  )
}

evaluate_pam_binary <- function(dmat, k_values = 2:8, min_cluster_size = 8) {
  purrr::map_dfr(k_values, function(k) {
    fit <- cluster::pam(dmat, k = k, diss = TRUE)
    sizes <- table(fit$clustering)

    tibble::tibble(
      k = k,
      avg_silhouette = fit$silinfo$avg.width,
      objective = fit$objective[2],
      min_cluster_size = min(sizes),
      max_cluster_size = max(sizes),
      n_small_clusters = sum(sizes < min_cluster_size),
      model = list(fit)
    )
  })
}

choose_k <- function(results, silhouette_tolerance = 0.02, min_cluster_size = 8) {
  max_sil <- max(results$avg_silhouette, na.rm = TRUE)

  candidates <- results %>%
    dplyr::filter(
      avg_silhouette >= (max_sil - silhouette_tolerance),
      min_cluster_size >= min_cluster_size
    ) %>%
    dplyr::arrange(k)

  if (nrow(candidates) == 0) {
    candidates <- results %>%
      dplyr::arrange(dplyr::desc(avg_silhouette), dplyr::desc(min_cluster_size), k)
  }

  candidates %>% dplyr::slice(1)
}

fit_multiselect_clusters <- function(data,
                                     id_col,
                                     value_col,
                                     cluster_col,
                                     label_col,
                                     min_prevalence = 2,
                                     k_values = 2:8,
                                     min_cluster_size = 8,
                                     silhouette_tolerance = 0.02,
                                     forced_k = NULL,
                                     exclude_from_labels = "Other",
                                     min_label_cluster_pct = 0.15) {
  built <- build_multiselect_matrix(
    data = data,
    id_col = id_col,
    value_col = value_col,
    min_prevalence = min_prevalence
  )

  dmat <- proxy::dist(built$X, method = "Jaccard")

  pam_results <- evaluate_pam_binary(
    dmat = dmat,
    k_values = k_values,
    min_cluster_size = min_cluster_size
  )

  best_k <- if (is.null(forced_k)) {
    choose_k(
      pam_results,
      silhouette_tolerance = silhouette_tolerance,
      min_cluster_size = min_cluster_size
    ) %>% dplyr::pull(k)
  } else {
    forced_k
  }

  best_pam <- pam_results %>%
    dplyr::filter(k == best_k) %>%
    dplyr::pull(model) %>%
    .[[1]]

  clusters <- tibble::tibble(
    respondent_id = built$respondent_ids,
    cluster = factor(best_pam$clustering)
  ) %>%
    dplyr::rename(!!cluster_col := cluster)

  cluster_sizes <- clusters %>%
    dplyr::count(.data[[cluster_col]], name = "cluster_size")

  overall_rate <- built$long_df %>%
    dplyr::distinct(.data[[id_col]], .data[[value_col]]) %>%
    dplyr::count(.data[[value_col]], name = "n_overall") %>%
    dplyr::mutate(overall_pct = n_overall / dplyr::n_distinct(built$respondent_ids))

  cluster_lift <- built$long_df %>%
    dplyr::distinct(.data[[id_col]], .data[[value_col]]) %>%
    dplyr::left_join(clusters, by = setNames("respondent_id", id_col)) %>%
    dplyr::count(.data[[cluster_col]], .data[[value_col]], name = "n_cluster") %>%
    dplyr::left_join(cluster_sizes, by = cluster_col) %>%
    dplyr::left_join(overall_rate, by = value_col) %>%
    dplyr::mutate(
      cluster_pct = n_cluster / cluster_size,
      lift = cluster_pct / overall_pct,
      pct_point_diff = cluster_pct - overall_pct,
      signature_score = lift * cluster_pct
    ) %>%
    dplyr::arrange(.data[[cluster_col]], dplyr::desc(signature_score), dplyr::desc(cluster_pct))

  labels_primary <- cluster_lift %>%
    dplyr::filter(
      .data[[value_col]] != exclude_from_labels,
      cluster_pct >= min_label_cluster_pct
    ) %>%
    dplyr::group_by(.data[[cluster_col]]) %>%
    dplyr::slice_max(order_by = signature_score, n = 3, with_ties = FALSE) %>%
    dplyr::summarise(
      label = paste(.data[[value_col]], collapse = " + "),
      .groups = "drop"
    )

  labels_fallback <- cluster_lift %>%
    dplyr::filter(.data[[value_col]] != exclude_from_labels) %>%
    dplyr::group_by(.data[[cluster_col]]) %>%
    dplyr::slice_max(order_by = signature_score, n = 3, with_ties = FALSE) %>%
    dplyr::summarise(
      fallback_label = paste(.data[[value_col]], collapse = " + "),
      .groups = "drop"
    )

  labels <- clusters %>%
    dplyr::distinct(.data[[cluster_col]]) %>%
    dplyr::left_join(labels_primary, by = cluster_col) %>%
    dplyr::left_join(labels_fallback, by = cluster_col) %>%
    dplyr::mutate(
      !!label_col := dplyr::coalesce(label, fallback_label)
    ) %>%
    dplyr::select(dplyr::all_of(cluster_col), dplyr::all_of(label_col))

  list(
    long_df = built$long_df,
    X = built$X,
    dmat = dmat,
    pam_results = pam_results,
    best_k = best_k,
    best_pam = best_pam,
    clusters = clusters,
    cluster_sizes = cluster_sizes,
    cluster_lift = cluster_lift,
    labels = labels
  )
}


# ------------------------------------------------------------
# MaxDiff task construction
#
# NOTE: The Stan model code lives in stan/hb_maxdiff.stan and
# should be loaded with stan_model(file = ...) in the calling
# script rather than being inlined here.
# ------------------------------------------------------------

derive_bw_mid_options <- function(responses, bw_design = default_bw_design()) {
  task_letters <- names(bw_design)

  for (t in seq_along(task_letters)) {
    opts <- bw_design[[task_letters[t]]]
    bestcol <- paste0("bw", t, "_best")
    worstcol <- paste0("bw", t, "_worst")
    midcol <- paste0("bw", t, "_mid")

    responses[[midcol]] <- purrr::map2_chr(
      responses[[bestcol]], responses[[worstcol]],
      ~ {
        b <- stringr::str_trim(as.character(.x))
        w <- stringr::str_trim(as.character(.y))
        if (is.na(b) || is.na(w) || b == "" || w == "") return(NA_character_)
        if (identical(b, w)) return(NA_character_)
        remain <- setdiff(opts, c(b, w))
        if (length(remain) == 1) remain else NA_character_
      }
    )
  }

  responses
}

build_bw_tasks_table <- function(responses, bw_design = default_bw_design()) {
  task_letters <- names(bw_design)
  attrs <- sort(unique(unlist(bw_design)))
  idx <- stats::setNames(seq_along(attrs), attrs)

  tasks <- purrr::map_dfr(seq_len(nrow(responses)), function(r) {
    row <- responses[r, ]
    purrr::map_dfr(seq_along(task_letters), function(t) {
      tibble::tibble(
        id = r,
        best = row[[paste0("bw", t, "_best")]],
        worst = row[[paste0("bw", t, "_worst")]],
        mid = row[[paste0("bw", t, "_mid")]]
      )
    })
  }) %>%
    dplyr::mutate(
      best = dplyr::na_if(stringr::str_trim(as.character(best)), ""),
      worst = dplyr::na_if(stringr::str_trim(as.character(worst)), ""),
      mid = dplyr::na_if(stringr::str_trim(as.character(mid)), "")
    ) %>%
    dplyr::filter(is.na(best) | is.na(worst) | best != worst) %>%
    dplyr::filter(!is.na(best), !is.na(worst), !is.na(mid)) %>%
    dplyr::filter(best %in% names(idx), worst %in% names(idx), mid %in% names(idx)) %>%
    dplyr::mutate(
      best = unname(idx[best]),
      worst = unname(idx[worst]),
      mid = unname(idx[mid])
    )

  list(
    tasks = tasks,
    task_letters = task_letters,
    attrs = attrs,
    idx = idx,
    K_attr = length(attrs)
  )
}


# ------------------------------------------------------------
# HB cluster assignments
# ------------------------------------------------------------

make_hb_cluster_assignments <- function(responses, u_mean, seed = 812) {
  x_centered <- sweep(u_mean, 1, rowMeans(u_mean), "-")
  x_scaled <- scale(x_centered)

  set.seed(seed)
  km_final <- stats::kmeans(x_scaled, centers = 2, nstart = 100)
  raw_clusters <- km_final$cluster

  smoothness_idx <- which(colnames(u_mean) == "Nib/tip Smoothness")
  if (length(smoothness_idx) != 1) {
    stop("Could not identify Nib/tip Smoothness column in u_mean")
  }

  smooth_by_cluster <- tapply(u_mean[, smoothness_idx], raw_clusters, mean)
  smooth_cluster_id <- as.integer(names(which.max(smooth_by_cluster)))
  fit_cluster_id <- setdiff(sort(unique(raw_clusters)), smooth_cluster_id)

  label_map <- c()
  label_map[as.character(smooth_cluster_id)] <- "Smoothness-Focused Writers"
  label_map[as.character(fit_cluster_id)] <- "Fit-Focused Writers"

  cluster_assignments <- tibble::tibble(
    respondent_id = responses$respondent_id,
    survey_id = responses$survey_id,
    cluster_raw = factor(raw_clusters),
    cluster = factor(
      unname(label_map[as.character(raw_clusters)]),
      levels = c("Smoothness-Focused Writers", "Fit-Focused Writers")
    )
  )

  list(
    cluster_assignments = cluster_assignments,
    km_final = km_final,
    x_centered = x_centered,
    x_scaled = x_scaled,
    smooth_by_cluster = smooth_by_cluster
  )
}
