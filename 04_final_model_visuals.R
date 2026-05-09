# ============================================================
# 04_final_model_visuals.R
#
# Purpose:
#   - Load the HB bundle and GLMM review bundle
#   - Recreate the core visuals from analysis_clean_with_plots_2.R
#   - Ensure the visuals reflect the selected final model
#
# Note: predict_fixed_prob, binary_token_set, binary_flag_value,
#       binary_flag_label, binary_flag_value_vec, binary_flag_label_vec,
#       pmax_na, choose_optimal_balance, choose_good_option_cutoff, and
#       build_fountain_catalog_prediction_base all live in pipeline_utils.R.
# ============================================================

project_dir <- if (exists("project_dir")) {
  normalizePath(project_dir, mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(project_dir, "config.R"))
source(file.path(project_dir, "pipeline_utils.R"))

ensure_packages(c("tidyverse", "ggplot2", "ggtext", "tidytext", "scales"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggtext)
  library(tidytext)
  library(scales)
})

paths <- make_pipeline_paths(project_dir)
make_output_dirs(paths)
review_bundle_path <- file.path(paths$out$review, "03_glmm_review_bundle.rds")
maxdiff_bundle_path <- file.path(paths$out$maxdiff, "02_maxdiff_bundle.rds")
clean_bundle_path <- file.path(paths$out$clean, "01_clean_data_bundle.rds")

if (!file.exists(review_bundle_path)) {
  stop("Could not find 03_glmm_review_bundle.rds. Run 03_glmm_review_and_validation.R first.")
}
if (!file.exists(maxdiff_bundle_path)) {
  stop("Could not find 02_maxdiff_bundle.rds. Run 02_maxdiff_hb_and_clusters.R first.")
}
if (!file.exists(clean_bundle_path)) {
  stop("Could not find 01_clean_data_bundle.rds. Run 01_clean_data_and_features.R first.")
}

review_bundle <- readRDS(review_bundle_path)
hb_bundle <- readRDS(maxdiff_bundle_path)
clean_bundle <- readRDS(clean_bundle_path)

final_model <- review_bundle$selected_models$final_model
plot_data <- review_bundle$selected_models$final_model_data %>%
  mutate(hand_size = factor(hand_size, levels = c("Small", "Medium", "Large")))

balance_col <- "darkgrey"
balance_col_no_title <- "steelblue"
base_col <- "steelblue"

# ------------------------------------------------------------
# Figure 1: HB stated-importance segments
# ------------------------------------------------------------
plot_pct_ci <- hb_bundle$share_sum %>%
  left_join(hb_bundle$cluster_sizes, by = "cluster") %>%
  mutate(
    cluster_label = case_when(
      as.character(cluster) == "Smoothness-Focused Writers" ~
        paste0("bold('Smoothness-Focused Writers')~'(' * italic(n) == ", n, " * ')'"),
      as.character(cluster) == "Fit-Focused Writers" ~
        paste0("bold('Fit-Focused Writers')~'(' * italic(n) == ", n, " * ')'"),
      TRUE ~ as.character(cluster)
    ),
    cluster_label = factor(
      cluster_label,
      levels = c(
        unique(cluster_label[cluster == "Smoothness-Focused Writers"]),
        unique(cluster_label[cluster == "Fit-Focused Writers"])
      )
    ),
    is_balance = Attribute == "Balance",
    Attribute_r = reorder_within(Attribute, Mean, cluster_label),
    label = sprintf("%.1f%%", Mean)
  ) %>%
  ggplot(aes(x = Mean, y = Attribute_r, fill = is_balance)) +
  geom_col() +
  geom_text(aes(label = label), hjust = -0.1, size = 3) +
  facet_wrap(
    ~ cluster_label,
    ncol = 1,
    scales = "free_y",
    labeller = label_parsed
  ) +
  scale_y_reordered() +
  scale_fill_manual(
    values = c(`FALSE` = base_col, `TRUE` = balance_col),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "People split into two groups when asked what matters for writing comfort",
    subtitle = paste0(
      "One group prioritized <b>smoothness</b>, while the other spread importance across <b>fit-related</b> features ",
      "<br><br>",
      "For both groups <span style='color:", balance_col, ";'><b>balance</b></span> ranked low on stated importance."
    ),
    x = "Importance Share (%)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = ggtext::element_markdown(size = 12),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 30, 25, 10)
  ) +
  coord_cartesian(clip = "off")

plot_pct_ci_no_title <- hb_bundle$share_sum %>%
  left_join(hb_bundle$cluster_sizes, by = "cluster") %>%
  mutate(
    cluster_label = case_when(
      as.character(cluster) == "Smoothness-Focused Writers" ~
        paste0("bold('Smoothness-Focused Writers')~'(' * italic(n) == ", n, " * ')'"),
      as.character(cluster) == "Fit-Focused Writers" ~
        paste0("bold('Fit-Focused Writers')~'(' * italic(n) == ", n, " * ')'"),
      TRUE ~ as.character(cluster)
    ),
    cluster_label = factor(
      cluster_label,
      levels = c(
        unique(cluster_label[cluster == "Smoothness-Focused Writers"]),
        unique(cluster_label[cluster == "Fit-Focused Writers"])
      )
    ),
    is_balance = Attribute == "Balance",
    Attribute_r = reorder_within(Attribute, Mean, cluster_label),
    label = sprintf("%.1f%%", Mean)
  ) %>%
  ggplot(aes(x = Mean, y = Attribute_r, fill = is_balance)) +
  geom_col() +
  geom_text(aes(label = label), hjust = -0.1, size = 3) +
  facet_wrap(
    ~ cluster_label,
    ncol = 1,
    scales = "free_y",
    labeller = label_parsed
  ) +
  scale_y_reordered() +
  scale_fill_manual(
    values = c(`FALSE` = base_col, `TRUE` = balance_col_no_title),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Importance Share (%)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = ggtext::element_markdown(size = 12),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 30, 25, 10)
  ) +
  coord_cartesian(clip = "off")

# ------------------------------------------------------------
# Predictions from the selected final model
# ------------------------------------------------------------
hand_levels <- c("Small", "Medium", "Large")

diam_mu <- mean(plot_data$diameter, na.rm = TRUE)
diam_sd <- sd(plot_data$diameter, na.rm = TRUE)
diam_rng <- quantile(plot_data$diameter, probs = c(0.02, 0.98), na.rm = TRUE)

cg_mu <- mean(plot_data$cg_rel_used, na.rm = TRUE)
cg_sd <- sd(plot_data$cg_rel_used, na.rm = TRUE)
cg_rng <- quantile(plot_data$cg_rel_used, probs = c(0.02, 0.98), na.rm = TRUE)

hand_props <- plot_data %>%
  count(hand_size) %>%
  mutate(prop = n / sum(n)) %>%
  { setNames(.$prop, .$hand_size) }

fountain_props <- plot_data %>%
  count(is_fountain) %>%
  mutate(prop = n / sum(n)) %>%
  { setNames(.$prop, .$is_fountain) }

# ------------------------------------------------------------
# Figure 2: Relative balance by hand size
# ------------------------------------------------------------
cg_seq <- seq(cg_rng[1], cg_rng[2], length.out = 250)
cg_min <- min(cg_seq, na.rm = TRUE)
cg_max <- max(cg_seq, na.rm = TRUE)
cg_span <- cg_max - cg_min
cg_pad <- 0.03 * cg_span

plot_cg_hand <- expand.grid(
  cg_rel_used = cg_seq,
  hand_size = hand_levels,
  is_fountain = c(0, 1),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    survey_id = plot_data$survey_id[1],
    hand_size = factor(hand_size, levels = hand_levels),
    diameter_z = 0,
    cg_rel_used_z = (cg_rel_used - cg_mu) / cg_sd
  )

plot_cg_hand$prob_high <- predict_fixed_prob(final_model, plot_cg_hand)

plot_cg_hand_marginal <- plot_cg_hand %>%
  mutate(w_f = as.numeric(fountain_props[as.character(is_fountain)])) %>%
  group_by(cg_rel_used, hand_size) %>%
  summarise(
    prob_high = weighted.mean(prob_high, w = w_f),
    .groups = "drop"
  )

p_cg_hand <- ggplot(plot_cg_hand_marginal, aes(x = cg_rel_used, y = prob_high)) +
  geom_line(aes(color = hand_size, linetype = hand_size), linewidth = 1.1, alpha = 0.85) +
  scale_color_manual(
    name = "Hand Size",
    values = c(
      "Small" = "#d95f02",
      "Medium" = "#2596be",
      "Large" = "#7570b3"
    )
  ) +
  scale_linetype_manual(
    name = "Hand Size",
    values = c(
      "Small" = "solid",
      "Medium" = "solid",
      "Large" = "solid"
    )
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 1)
  ) +
  scale_y_continuous(
    limits = c(-0.10, 1),
    labels = label_percent(accuracy = 1)
  ) +
  labs(
    title = "Balance seems to matter differently by hand size",
    subtitle = "Back-heavier balance looked worst for small hands, flatter for medium hands, and moderately worse for large hands",
    x = NULL,
    y = "Probability of High Comfort Rating",
    caption = "Note: model-estimated chance of high comfort, with grip diameter and other model factors held constant"
  ) +
  annotate(
    "text",
    x = cg_min + cg_pad,
    y = -0.06,
    label = "<<- More front-heavy",
    hjust = 0,
    vjust = 1,
    size = 3.8,
    color = "black"
  ) +
  annotate(
    "text",
    x = cg_max - cg_pad,
    y = -0.06,
    label = "More back-heavy ->>",
    hjust = 1,
    vjust = 1,
    size = 3.8,
    color = "black"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.caption = element_text(hjust = 0, size = 9, color = "gray30"),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = margin(10, 20, 25, 10)
  )

p_cg_hand_no_title <- ggplot(plot_cg_hand_marginal, aes(x = cg_rel_used, y = prob_high)) +
  geom_line(aes(color = hand_size, linetype = hand_size), linewidth = 1.1, alpha = 0.85) +
  scale_color_manual(
    name = "Hand Size",
    values = c(
      "Small" = "#d95f02",
      "Medium" = "#2596be",
      "Large" = "#7570b3"
    )
  ) +
  scale_linetype_manual(
    name = "Hand Size",
    values = c(
      "Small" = "solid",
      "Medium" = "solid",
      "Large" = "solid"
    )
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 1)
  ) +
  scale_y_continuous(
    limits = c(-0.10, 1),
    labels = label_percent(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = "Probability of High Comfort Rating",
  ) +
  annotate(
    "text",
    x = cg_min + cg_pad,
    y = -0.06,
    label = "<<- More front-heavy",
    hjust = 0,
    vjust = 1,
    size = 3.8,
    color = "black"
  ) +
  annotate(
    "text",
    x = cg_max - cg_pad,
    y = -0.06,
    label = "More back-heavy ->>",
    hjust = 1,
    vjust = 1,
    size = 3.8,
    color = "black"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.caption = element_text(hjust = 0, size = 9, color = "gray30"),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = margin(10, 20, 25, 10)
  )

# ------------------------------------------------------------
# Figure 3: Diameter effect
# ------------------------------------------------------------
plot_diam <- expand.grid(
  diameter_mm = seq(diam_rng[1], diam_rng[2], length.out = 250),
  hand_size = hand_levels,
  is_fountain = c(0, 1),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    survey_id = plot_data$survey_id[1],
    hand_size = factor(hand_size, levels = hand_levels),
    diameter_z = (diameter_mm - diam_mu) / diam_sd,
    cg_rel_used_z = 0
  )

plot_diam$prob_high <- predict_fixed_prob(final_model, plot_diam)

plot_diam_marginal <- plot_diam %>%
  mutate(
    w_h = as.numeric(hand_props[as.character(hand_size)]),
    w_f = as.numeric(fountain_props[as.character(is_fountain)]),
    w_all = w_h * w_f
  ) %>%
  group_by(diameter_mm) %>%
  summarise(
    prob_high = weighted.mean(prob_high, w = w_all),
    .groups = "drop"
  )

p_diam_final <- ggplot(plot_diam_marginal, aes(x = diameter_mm, y = prob_high)) +
  geom_line(linewidth = 1.1, alpha = 0.9, color = "#2596be") +
  scale_y_continuous(
    limits = c(0, 1),
    labels = label_percent(accuracy = 1)
  ) +
  labs(
    title = "Wider grip sections were linked to higher comfort",
    subtitle = "Estimated comfort increased with diameter after accounting for hand size, balance, and instrument type",
    x = "Grip Diameter (mm)",
    y = "Probability of High Comfort Rating"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    panel.grid.minor = element_blank()
  )

p_diam_final_no_title <- ggplot(plot_diam_marginal, aes(x = diameter_mm, y = prob_high)) +
  geom_line(linewidth = 1.1, alpha = 0.9, color = "#2596be") +
  scale_y_continuous(
    limits = c(0, 1),
    labels = label_percent(accuracy = 1)
  ) +
  labs(
    x = "Grip Diameter (mm)",
    y = "Probability of High Comfort Rating"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------
# Figure 4: Fountain vs non-fountain effect
# ------------------------------------------------------------
plot_fountain <- plot_data %>%
  mutate(
    diameter_z = 0,
    cg_rel_used_z = 0
  )

plot_fountain_cf <- bind_rows(
  plot_fountain %>% mutate(is_fountain = 0, instrument_type = "Non-fountain"),
  plot_fountain %>% mutate(is_fountain = 1, instrument_type = "Fountain")
)

plot_fountain_cf$prob_high <- predict_fixed_prob(final_model, plot_fountain_cf)

fountain_story <- plot_fountain_cf %>%
  group_by(instrument_type) %>%
  summarise(
    prob_high = mean(prob_high, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    instrument_type = factor(instrument_type, levels = c("Non-fountain", "Fountain")),
    label = label_percent(accuracy = 1)(prob_high)
  )

p_fountain_story <- ggplot(fountain_story, aes(x = instrument_type, y = prob_high, fill = instrument_type)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(
    aes(label = label),
    vjust = -0.6,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Non-fountain" = "grey70",
      "Fountain" = "#2596be"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = label_percent(accuracy = 1)
  ) +
  labs(
    title = "Fountain pens are associated with higher comfort in this sample",
    subtitle = "Estimated comfort is higher for fountain pens when diameter and relative balance are held constant",
    x = NULL,
    y = "Probability of High Comfort Rating",
    caption = "Note: bars show fixed-effect predictions averaged across the observed hand-size mix"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.caption = element_text(hjust = 0, size = 9, color = "gray30"),
    panel.grid.minor = element_blank()
  )

p_fountain_story_no_title <- ggplot(fountain_story, aes(x = instrument_type, y = prob_high, fill = instrument_type)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(
    aes(label = label),
    vjust = -0.6,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Non-fountain" = "grey70",
      "Fountain" = "#2596be"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = label_percent(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = "Probability of High Comfort Rating"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    plot.caption = element_text(hjust = 0, size = 9, color = "gray30"),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------
# Figure 5: Default Explorer view
# ------------------------------------------------------------
catalog_prediction_base <- build_fountain_catalog_prediction_base(
  product_lookup = clean_bundle$product_lookup,
  proxy_models = clean_bundle$proxy_models,
  final_model = final_model,
  final_model_data = plot_data
)

catalog_cutoff_eval <- choose_good_option_cutoff(
  fit = final_model,
  data = plot_data,
  outcome = "high6"
)
catalog_good_option_cutoff <- catalog_cutoff_eval$threshold[1]
catalog_total_n <- nrow(catalog_prediction_base)

if (catalog_total_n == 0) {
  stop("No fountain-pen catalog entries with complete prediction inputs are available for the default Explorer view.", call. = FALSE)
}

catalog_lollipop_data <- dplyr::bind_rows(
  catalog_prediction_base %>% dplyr::transmute(Configuration = "Posted", `Hand Size` = "Small", prob_high = prob_small_posted_effective),
  catalog_prediction_base %>% dplyr::transmute(Configuration = "Posted", `Hand Size` = "Medium", prob_high = prob_medium_posted_effective),
  catalog_prediction_base %>% dplyr::transmute(Configuration = "Not Posted", `Hand Size` = "Small", prob_high = prob_small_not_posted),
  catalog_prediction_base %>% dplyr::transmute(Configuration = "Not Posted", `Hand Size` = "Medium", prob_high = prob_medium_not_posted),
  catalog_prediction_base %>% dplyr::transmute(Configuration = "Best Setup", `Hand Size` = "Small", prob_high = optimal_prob_small),
  catalog_prediction_base %>% dplyr::transmute(Configuration = "Best Setup", `Hand Size` = "Medium", prob_high = optimal_prob_medium)
) %>%
  dplyr::mutate(
    Configuration = factor(Configuration, levels = c("Posted", "Not Posted", "Best Setup")),
    `Hand Size` = factor(`Hand Size`, levels = c("Small", "Medium"))
  ) %>%
  dplyr::group_by(Configuration, `Hand Size`) %>%
  dplyr::summarise(
    n_good = sum(prob_high >= catalog_good_option_cutoff, na.rm = TRUE),
    share_good = n_good / catalog_total_n,
    .groups = "drop"
  ) %>%
  dplyr::mutate(label = scales::percent(share_good, accuracy = 1))

p_catalog_lollipop_default <- ggplot(catalog_lollipop_data, aes(x = `Hand Size`, y = share_good, fill = Configuration)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.68) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.74),
    vjust = -0.4,
    size = 3.2
  ) +
  scale_fill_manual(values = c(
    "Posted" = "#8c8c8c",
    "Not Posted" = "#4c78a8",
    "Best Setup" = "#59a14f"
  )) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    title = "Small hands have a narrower pool of likely-fit fountain pens",
    subtitle = "Posting also reduced the likely-fit pool much more for small hands than for medium hands",
    x = NULL,
    y = "% Pens Likely to Work Well",
    fill = NULL,
    caption = "Note: model applied to 1.1k unique pens. For each pen, 'best setup' takes the higher comfort estimate between posted and not posted"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(hjust = 0, size = 9, color = "gray30"),
    plot.margin = margin(10, 16, 20, 10)
  )

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------
ggsave(file.path(paths$out$visuals, "01_hb_cluster_share.png"), plot_pct_ci, width = 12, height = 7, dpi = 300, bg = "white")
ggsave(file.path(paths$out$visuals, "01_hb_cluster_share_no_title.png"), plot_pct_ci_no_title, width = 8, height = 8, dpi = 300)
ggsave(file.path(paths$out$visuals, "02_balance_by_hand_size.png"), p_cg_hand, width = 12, height = 7, dpi = 300, bg = "white")
ggsave(file.path(paths$out$visuals, "02_balance_by_hand_size_no_title.png"), p_cg_hand_no_title, width = 7, height = 6, dpi = 300)
ggsave(file.path(paths$out$visuals, "03_diameter_effect.png"), p_diam_final, width = 12, height = 7, dpi = 300, bg = "white")
ggsave(file.path(paths$out$visuals, "03_diameter_effect_no_title.png"), p_diam_final_no_title, width = 8, height = 6, dpi = 300)
ggsave(file.path(paths$out$visuals, "04_fountain_effect.png"), p_fountain_story, width = 12, height = 7, dpi = 300, bg = "white")
ggsave(file.path(paths$out$visuals, "04_fountain_effect_no_title.png"), p_fountain_story_no_title, width = 7, height = 6, dpi = 300)
ggsave(file.path(paths$out$visuals, "05_explorer_default_view.png"), p_catalog_lollipop_default, width = 12, height = 7, dpi = 300, bg = "white")

pdf(file.path(paths$out$visuals, "final_model_visuals.pdf"), width = 8.5, height = 7)
print(plot_pct_ci)
print(plot_pct_ci_no_title)
print(p_cg_hand)
print(p_cg_hand_no_title)
print(p_diam_final)
print(p_diam_final_no_title)
print(p_fountain_story)
print(p_fountain_story_no_title)
dev.off()

visual_bundle <- list(
  metadata = list(
    created_at = Sys.time(),
    project_dir = project_dir,
    final_model_name = review_bundle$selected_models$final_model_name,
    model_selection_note = review_bundle$metadata$model_selection_note
  ),
  plot_data = plot_data,
  hb_share_plot_data = hb_bundle$share_sum,
  plots = list(
    plot_pct_ci = plot_pct_ci,
    p_cg_hand = p_cg_hand,
    p_diam_final = p_diam_final,
    p_fountain_story = p_fountain_story,
    plot_pct_ci_no_title = plot_pct_ci_no_title,
    p_cg_hand_no_title = p_cg_hand_no_title,
    p_diam_final_no_title = p_diam_final_no_title,
    p_fountain_story_no_title = p_fountain_story_no_title
  )
)

saveRDS(visual_bundle, file.path(paths$out$visuals, "04_visual_bundle.rds"))

cat("\n============================================================\n")
cat("FINAL MODEL VISUALS COMPLETE\n")
cat("============================================================\n")
cat(review_bundle$metadata$model_selection_note, "\n\n")
cat("Outputs written to:\n", paths$out$visuals, "\n", sep = "")
