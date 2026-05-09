# R/utils_catalog.R — Catalog prediction, binary flag helpers, and cutoff selection

# --- Fixed-effects prediction ---

predict_fixed_prob <- function(fit, newdata) {
  predict(
    fit,
    newdata = newdata,
    re.form = NA,
    type = "response",
    allow.new.levels = TRUE
  )
}

# --- Binary flag parsing ---

binary_token_set <- function(x) {
  x_chr <- as.character(x)

  if (length(x_chr) == 0 || all(is.na(x_chr))) {
    return(integer(0))
  }

  bits <- unlist(strsplit(x_chr[1], ","))
  bits <- trimws(bits)
  vals <- parse_capped_flag(bits)
  vals <- vals[!is.na(vals)]

  unique(as.integer(vals))
}

binary_flag_value <- function(x) {
  toks <- binary_token_set(x)

  if (length(toks) == 1 && toks == 1L) return(1)
  if (length(toks) == 1 && toks == 0L) return(0)

  NA_real_
}

binary_flag_label <- function(x) {
  toks <- sort(binary_token_set(x))

  if (length(toks) == 1 && toks == 1L) return("Yes")
  if (length(toks) == 1 && toks == 0L) return("No")

  "Mixed / Unknown"
}

binary_flag_value_vec <- function(x) {
  vapply(x, binary_flag_value, numeric(1))
}

binary_flag_label_vec <- function(x) {
  vapply(x, binary_flag_label, character(1))
}

# --- Balance and threshold helpers ---

pmax_na <- function(a, b) {
  out <- ifelse(
    is.na(a) & is.na(b),
    NA_real_,
    ifelse(is.na(a), b, ifelse(is.na(b), a, pmax(a, b)))
  )

  as.numeric(out)
}

choose_optimal_balance <- function(cg_not_posted, cg_posted,
                                   prob_small_not_posted, prob_small_posted,
                                   prob_medium_not_posted, prob_medium_posted,
                                   prob_large_not_posted, prob_large_posted) {
  small_pref_posted <- !is.na(prob_small_posted) &
    (is.na(prob_small_not_posted) | prob_small_posted > prob_small_not_posted)

  medium_pref_posted <- !is.na(prob_medium_posted) &
    (is.na(prob_medium_not_posted) | prob_medium_posted > prob_medium_not_posted)

  large_pref_posted <- !is.na(prob_large_posted) &
    (is.na(prob_large_not_posted) | prob_large_posted > prob_large_not_posted)

  out <- ifelse(
    small_pref_posted & medium_pref_posted & large_pref_posted,
    cg_posted,
    ifelse(
      !small_pref_posted & !medium_pref_posted & !large_pref_posted,
      cg_not_posted,
      rowMeans(cbind(cg_not_posted, cg_posted), na.rm = TRUE)
    )
  )

  out[!is.finite(out)] <- NA_real_
  out
}

choose_good_option_cutoff <- function(fit, data, outcome = "high6",
                                      threshold_grid = seq(0.05, 0.95, by = 0.005)) {
  y <- data[[outcome]]

  pred <- rep(NA_real_, nrow(data))
  cc <- stats::complete.cases(data[, c("survey_id", "hand_size", "diameter_z", "cg_rel_used_z", "is_fountain", outcome)])

  if (any(cc)) {
    pred[cc] <- as.numeric(
      predict(
        fit,
        newdata = data[cc, , drop = FALSE],
        re.form = NA,
        type = "response",
        allow.new.levels = TRUE
      )
    )
  }

  eval_df <- tibble::tibble(y = y, pred = pred) %>%
    tidyr::drop_na()

  if (nrow(eval_df) == 0) {
    stop("Could not derive a probability cutoff from the final model data.", call. = FALSE)
  }

  tibble::tibble(threshold = threshold_grid) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      tp = sum(eval_df$pred >= threshold & eval_df$y == 1),
      tn = sum(eval_df$pred < threshold & eval_df$y == 0),
      fp = sum(eval_df$pred >= threshold & eval_df$y == 0),
      fn = sum(eval_df$pred < threshold & eval_df$y == 1),
      sensitivity = ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_),
      specificity = ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_),
      balanced_accuracy = (sensitivity + specificity) / 2,
      youden_j = sensitivity + specificity - 1
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(
      dplyr::desc(youden_j),
      dplyr::desc(balanced_accuracy),
      abs(threshold - 0.5)
    )
}

# --- Fountain pen catalog prediction ---

build_fountain_catalog_prediction_base <- function(product_lookup, proxy_models,
                                                   final_model, final_model_data) {
  diameter_mu <- mean(final_model_data$diameter, na.rm = TRUE)
  diameter_sd <- sd(final_model_data$diameter, na.rm = TRUE)
  cg_rel_mu <- mean(final_model_data$cg_rel_used, na.rm = TRUE)
  cg_rel_sd <- sd(final_model_data$cg_rel_used, na.rm = TRUE)

  if (!is.finite(diameter_sd) || diameter_sd <= 0 ||
      !is.finite(cg_rel_sd) || cg_rel_sd <= 0) {
    stop("Could not derive scaling parameters from final_model_data.", call. = FALSE)
  }

  survey_placeholder <- final_model_data$survey_id[1]
  hand_levels_model <- if (is.factor(final_model_data$hand_size)) {
    levels(final_model_data$hand_size)
  } else {
    c("Small", "Medium", "Large")
  }

  predict_catalog_prob <- function(diameter_vec, cg_rel_vec, hand_size_value) {
    newdata <- tibble(
      survey_id = rep(survey_placeholder, length(diameter_vec)),
      hand_size = factor(rep(hand_size_value, length(diameter_vec)), levels = hand_levels_model),
      diameter_z = (diameter_vec - diameter_mu) / diameter_sd,
      cg_rel_used_z = (cg_rel_vec - cg_rel_mu) / cg_rel_sd,
      is_fountain = 1L
    )

    preds <- rep(NA_real_, nrow(newdata))
    cc <- complete.cases(newdata)

    if (any(cc)) {
      preds[cc] <- as.numeric(predict_fixed_prob(final_model, newdata[cc, , drop = FALSE]))
    }

    preds
  }

  catalog <- product_lookup %>%
    filter(product_type == "Fountain Pens") %>%
    mutate(
      catalog_row_id = row_number(),
      manufacturer = dplyr::coalesce(dplyr::na_if(trimws(manufacturer), ""), "Unknown / missing"),
      material_family_raw = dplyr::coalesce(dplyr::na_if(trimws(material_family), ""), "Other / Unknown"),
      material_family = dplyr::case_when(
        material_family_raw %in% c("Celluloid", "Cellulose Acetate", "Ebonite") ~ "Vintage Materials",
        material_family_raw %in% c("Metal", "Metal - Aluminum", "Metal - Steel", "Metal - Brass", "Metal - Magnesium") ~ "General Metals",
        material_family_raw %in% c("Metal - Gold", "Metal - Silver", "Metal - Bronze", "Metal - Copper", "Metal - Titanium", "Metal - Zirconium", "Metal - Platinum") ~ "Premium Metals",
        material_family_raw %in% c("Plastic", "Plastic - Bioplastic") ~ "General Plastics",
        material_family_raw %in% c("Plastic - Resin/Acrylic", "Plastic - Polycarbonate", "Celluloid Derivitave", "Celluloid Derivative") ~ "Premium Plastics",
        material_family_raw %in% c("Stone", "Wood", "Bone", "Carbon Fiber", "Glass", "Lucite", "Micarta", "Enamel", "Ultem", "Delrin", "Torlon", "Other / Unknown", "Unknown / missing") ~ "Other / Unknown",
        TRUE ~ material_family_raw
      ),
      model_name = dplyr::coalesce(dplyr::na_if(trimws(model_name), ""), "Unknown model"),
      item_label = stringr::str_squish(paste(manufacturer, model_name)),
      retractable_flag = binary_flag_value_vec(retractable),
      postable_flag = binary_flag_value_vec(postable),
      retractable_filter = binary_flag_label_vec(retractable),
      postable_filter = binary_flag_label_vec(postable),
      capped_flag_raw = parse_capped_flag(capped),
      price = safe_num(med_price),
      diameter = safe_num(med_diameter_grip),
      body_diameter = safe_num(med_diameter),
      med_length_uncapped = safe_num(med_length_uncapped),
      med_length_posted = safe_num(med_length_posted),
      med_length = safe_num(med_length),
      med_weight = safe_num(med_weight),
      med_weight_cap = safe_num(med_weight_cap),
      med_weight_uncapped = dplyr::coalesce(safe_num(med_weight) - safe_num(med_weight_cap), safe_num(med_weight)),
      capped_flag = case_when(
        retractable_flag == 1 ~ 0,
        capped_flag_raw == 1 ~ 1,
        !is.na(med_weight_cap) | !is.na(med_length_posted) ~ 1,
        TRUE ~ 0
      ),
      length_not_posted = dplyr::coalesce(med_length_uncapped, med_length),
      length_closed_model = case_when(
        capped_flag == 0 ~ dplyr::coalesce(med_length_uncapped, med_length),
        TRUE ~ med_length
      ),
      length_posted_model = case_when(
        capped_flag == 0 ~ dplyr::coalesce(med_length_uncapped, med_length),
        TRUE ~ med_length_posted
      ),
      weight_not_posted = case_when(
        capped_flag == 0 ~ med_weight,
        TRUE ~ med_weight_uncapped
      ),
      cap_weight_model = case_when(
        capped_flag == 0 ~ 0,
        TRUE ~ med_weight_cap
      ),
      price_band = case_when(
        is.na(price) ~ "Unknown / missing",
        price < 25 ~ "Under $25",
        price < 50 ~ "$25-$49",
        price < 100 ~ "$50-$99",
        price < 200 ~ "$100-$199",
        price < 500 ~ "$200-$499",
        TRUE ~ "$500+"
      )
    )

  pred_input_no_wcp <- catalog %>%
    transmute(
      LPn = length_not_posted,
      LCp = length_closed_model,
      LPs = length_posted_model,
      DBa = body_diameter,
      D25 = body_diameter,
      WPn = weight_not_posted
    )

  pred_input_with_wcp <- catalog %>%
    transmute(
      LPn = length_not_posted,
      LCp = length_closed_model,
      LPs = length_posted_model,
      DBa = body_diameter,
      D25 = body_diameter,
      WPn = weight_not_posted,
      WCp = cap_weight_model
    )

  catalog <- catalog %>%
    mutate(
      cg_raw_not_posted = predict_lm_model(proxy_models$fit_cgu, pred_input_with_wcp),
      cg_raw_posted = predict_ridge_model(proxy_models$fit_cgp, pred_input_with_wcp),
      cg_rel_not_posted = predict_lm_model(proxy_models$fit_cgpn, pred_input_no_wcp),
      cg_rel_posted = predict_lm_model(proxy_models$fit_cgps, pred_input_with_wcp)
    ) %>%
    mutate(
      posted_state_available = case_when(
        retractable_flag == 1 ~ TRUE,
        postable_flag == 0 ~ FALSE,
        !is.na(med_length_posted) & !is.na(med_weight) & !is.na(cg_raw_posted) & !is.na(cg_rel_posted) ~ TRUE,
        TRUE ~ FALSE
      ),
      length_posted_effective = case_when(
        retractable_flag == 1 ~ length_not_posted,
        posted_state_available ~ med_length_posted,
        TRUE ~ NA_real_
      ),
      weight_posted_effective = case_when(
        retractable_flag == 1 ~ weight_not_posted,
        posted_state_available ~ med_weight,
        TRUE ~ NA_real_
      ),
      cg_raw_posted_effective = case_when(
        retractable_flag == 1 ~ cg_raw_not_posted,
        posted_state_available ~ cg_raw_posted,
        TRUE ~ NA_real_
      ),
      cg_rel_posted_effective = case_when(
        retractable_flag == 1 ~ cg_rel_not_posted,
        posted_state_available ~ cg_rel_posted,
        TRUE ~ NA_real_
      )
    )

  torque_not_posted <- calc_torque(
    mass_g = catalog$weight_not_posted,
    cg_mm = catalog$cg_raw_not_posted,
    grip_mm = 25
  )

  torque_posted_effective <- calc_torque(
    mass_g = catalog$weight_posted_effective,
    cg_mm = catalog$cg_raw_posted_effective,
    grip_mm = 25
  )

  catalog <- catalog %>%
    mutate(
      torque_not_posted = torque_not_posted$torque_abs_n_m,
      torque_posted_effective = torque_posted_effective$torque_abs_n_m,
      complete_not_posted = !is.na(diameter) & !is.na(cg_rel_not_posted),
      complete_posted_effective = !is.na(diameter) & !is.na(cg_rel_posted_effective),
      prob_small_not_posted = predict_catalog_prob(diameter, cg_rel_not_posted, "Small"),
      prob_medium_not_posted = predict_catalog_prob(diameter, cg_rel_not_posted, "Medium"),
      prob_large_not_posted = predict_catalog_prob(diameter, cg_rel_not_posted, "Large"),
      prob_small_posted_effective = predict_catalog_prob(diameter, cg_rel_posted_effective, "Small"),
      prob_medium_posted_effective = predict_catalog_prob(diameter, cg_rel_posted_effective, "Medium"),
      prob_large_posted_effective = predict_catalog_prob(diameter, cg_rel_posted_effective, "Large")
    ) %>%
    mutate(
      prob_small_not_posted = if_else(complete_not_posted, prob_small_not_posted, NA_real_),
      prob_medium_not_posted = if_else(complete_not_posted, prob_medium_not_posted, NA_real_),
      prob_large_not_posted = if_else(complete_not_posted, prob_large_not_posted, NA_real_),
      prob_small_posted_effective = if_else(complete_posted_effective, prob_small_posted_effective, NA_real_),
      prob_medium_posted_effective = if_else(complete_posted_effective, prob_medium_posted_effective, NA_real_),
      prob_large_posted_effective = if_else(complete_posted_effective, prob_large_posted_effective, NA_real_),
      optimal_prob_small = pmax_na(prob_small_not_posted, prob_small_posted_effective),
      optimal_prob_medium = pmax_na(prob_medium_not_posted, prob_medium_posted_effective),
      optimal_prob_large = pmax_na(prob_large_not_posted, prob_large_posted_effective),
      optimal_cg_rel = choose_optimal_balance(
        cg_not_posted = cg_rel_not_posted,
        cg_posted = cg_rel_posted_effective,
        prob_small_not_posted = prob_small_not_posted,
        prob_small_posted = prob_small_posted_effective,
        prob_medium_not_posted = prob_medium_not_posted,
        prob_medium_posted = prob_medium_posted_effective,
        prob_large_not_posted = prob_large_not_posted,
        prob_large_posted = prob_large_posted_effective
      )
    ) %>%
    filter(complete_not_posted | complete_posted_effective)

  catalog
}
