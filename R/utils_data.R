# R/utils_data.R — Data cleaning, feature engineering, and regression utilities

# --- Numeric helpers ---

safe_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.infinite(x)] <- NA_real_
  x
}

safe_scalar <- function(x) {
  out <- tryCatch(as.numeric(x)[1], error = function(e) NA_real_)
  if (length(out) == 0) NA_real_ else out
}

parse_capped_flag <- function(x) {
  x_chr <- trimws(as.character(x))
  dplyr::case_when(
    x_chr %in% c("1", "1.0", "TRUE", "T", "Yes", "yes") ~ 1,
    x_chr %in% c("0", "0.0", "FALSE", "F", "No", "no", "") ~ 0,
    TRUE ~ NA_real_
  )
}

clamp_to_range <- function(x, rng) {
  pmin(pmax(x, rng[1]), rng[2])
}

impute_with <- function(df, medians) {
  out <- df
  for (nm in names(medians)) {
    out[[nm]] <- ifelse(is.na(out[[nm]]), medians[[nm]], out[[nm]])
  }
  tibble::as_tibble(out)
}

# --- Regression models ---

fit_lm_model <- function(data, target, predictors) {
  train <- data %>%
    dplyr::select(dplyr::all_of(c(target, predictors))) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    tidyr::drop_na(dplyr::all_of(target))

  train_x <- train %>% dplyr::select(dplyr::all_of(predictors))
  medians <- purrr::map_dbl(train_x, ~ stats::median(.x, na.rm = TRUE))
  train_x <- impute_with(train_x, medians)
  train_fit <- dplyr::bind_cols(train %>% dplyr::select(dplyr::all_of(target)), train_x)

  form <- stats::as.formula(paste(target, "~", paste(predictors, collapse = " + ")))
  model <- stats::lm(form, data = train_fit)

  list(
    model = model,
    predictors = predictors,
    medians = medians,
    target_range = range(train[[target]], na.rm = TRUE)
  )
}

predict_lm_model <- function(fit_obj, newdata) {
  x <- newdata %>%
    dplyr::select(dplyr::all_of(fit_obj$predictors)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))

  complete_idx <- stats::complete.cases(x)
  preds <- rep(NA_real_, nrow(x))

  if (any(complete_idx)) {
    preds[complete_idx] <- stats::predict(
      fit_obj$model,
      newdata = x[complete_idx, , drop = FALSE]
    )
    preds[complete_idx] <- clamp_to_range(preds[complete_idx], fit_obj$target_range)
  }

  preds
}

fit_ridge_model <- function(data, target, predictors) {
  train <- data %>%
    dplyr::select(dplyr::all_of(c(target, predictors))) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    tidyr::drop_na(dplyr::all_of(target))

  train_x <- train %>% dplyr::select(dplyr::all_of(predictors))
  medians <- purrr::map_dbl(train_x, ~ stats::median(.x, na.rm = TRUE))
  train_x <- impute_with(train_x, medians)

  x_mat <- as.matrix(train_x)
  y_vec <- train[[target]]

  cv_fit <- glmnet::cv.glmnet(
    x = x_mat,
    y = y_vec,
    alpha = 0,
    standardize = TRUE,
    nfolds = 5
  )

  list(
    model = cv_fit,
    predictors = predictors,
    medians = medians,
    target_range = range(y_vec, na.rm = TRUE)
  )
}

predict_ridge_model <- function(fit_obj, newdata) {
  x <- newdata %>%
    dplyr::select(dplyr::all_of(fit_obj$predictors)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))

  complete_idx <- stats::complete.cases(x)
  preds <- rep(NA_real_, nrow(x))

  if (any(complete_idx)) {
    preds[complete_idx] <- as.numeric(
      stats::predict(
        fit_obj$model,
        newx = as.matrix(x[complete_idx, , drop = FALSE]),
        s = "lambda.min"
      )
    )
    preds[complete_idx] <- clamp_to_range(preds[complete_idx], fit_obj$target_range)
  }

  preds
}

# --- CV helpers ---

rmse_vec <- function(obs, pred) {
  cc <- stats::complete.cases(obs, pred)
  if (!any(cc)) return(NA_real_)
  sqrt(mean((obs[cc] - pred[cc])^2))
}

mae_vec <- function(obs, pred) {
  cc <- stats::complete.cases(obs, pred)
  if (!any(cc)) return(NA_real_)
  mean(abs(obs[cc] - pred[cc]))
}

rsq_vec <- function(obs, pred) {
  cc <- stats::complete.cases(obs, pred)
  if (sum(cc) < 2) return(NA_real_)
  stats::cor(obs[cc], pred[cc])^2
}

cross_validate_regression <- function(data, target, predictors, model_name,
                                      model_type = c("lm", "ridge"),
                                      k = 5, repeats = 3, seed = 812) {
  model_type <- match.arg(model_type)

  d <- data %>%
    dplyr::select(dplyr::all_of(c(target, predictors))) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    tidyr::drop_na()

  k_use <- min(k, nrow(d))
  if (k_use < 2) {
    return(tibble::tibble(
      model = model_name,
      target = target,
      model_type = model_type,
      repeat_id = integer(),
      fold_id = integer(),
      n_test = integer(),
      RMSE = numeric(),
      MAE = numeric(),
      R2 = numeric()
    ))
  }

  out <- vector("list", length = repeats * k_use)
  idx <- 1L

  for (rep_i in seq_len(repeats)) {
    set.seed(seed + rep_i + nchar(target))
    fold_id <- sample(rep(seq_len(k_use), length.out = nrow(d)))

    for (fold_i in seq_len(k_use)) {
      train <- d[fold_id != fold_i, , drop = FALSE]
      test <- d[fold_id == fold_i, , drop = FALSE]

      fit_obj <- if (model_type == "lm") {
        fit_lm_model(train, target, predictors)
      } else {
        fit_ridge_model(train, target, predictors)
      }

      pred <- if (model_type == "lm") {
        predict_lm_model(fit_obj, test)
      } else {
        predict_ridge_model(fit_obj, test)
      }

      out[[idx]] <- tibble::tibble(
        model = model_name,
        target = target,
        model_type = model_type,
        repeat_id = rep_i,
        fold_id = fold_i,
        n_test = sum(stats::complete.cases(test[[target]], pred)),
        RMSE = rmse_vec(test[[target]], pred),
        MAE = mae_vec(test[[target]], pred),
        R2 = rsq_vec(test[[target]], pred)
      )
      idx <- idx + 1L
    }
  }

  dplyr::bind_rows(out)
}

# --- Physics helpers ---

calc_torque <- function(mass_g, cg_mm, grip_mm) {
  mass_kg <- mass_g / 1000
  lever_m <- (cg_mm - grip_mm) / 1000
  force_n <- mass_kg * 9.80665

  tibble::tibble(
    torque_signed_n_m = force_n * lever_m,
    torque_abs_n_m = force_n * abs(lever_m),
    torque_signed_n_mm = force_n * lever_m * 1000,
    torque_abs_n_mm = force_n * abs(lever_m) * 1000
  )
}

# --- Survey response cleaning ---

clean_press_force <- function(x) {
  dplyr::case_when(
    x %in% c("Light", "light", "LIGHT") ~ "Light",
    x %in% c("Moderate", "moderate", "MODERATE", "Medium", "medium", "Heavy", "heavy", "HEAVY") ~ "Moderate/Heavy",
    TRUE ~ NA_character_
  )
}

clean_balance_pref <- function(x) {
  dplyr::case_when(
    x %in% c("Back-heavy") ~ "Back-heavy",
    x %in% c("Balanced") ~ "Balanced",
    x %in% c("Front-heavy") ~ "Front-heavy",
    TRUE ~ NA_character_
  )
}

clean_hand_size <- function(x) {
  dplyr::case_when(
    x %in% c("Extra Small", "Small") ~ "Small",
    x %in% c("Medium") ~ "Medium",
    x %in% c("Large", "Extra Large") ~ "Large",
    TRUE ~ NA_character_
  )
}

clean_handwriting_size <- function(x) {
  dplyr::case_when(
    x %in% c("Small") ~ "Small",
    x %in% c("Medium") ~ "Medium",
    x %in% c("Large") ~ "Large",
    TRUE ~ NA_character_
  )
}

clean_grip_posture <- function(x) {
  dplyr::case_when(
    x %in% c(12, 13) ~ "Tripod",
    x %in% c(19, 20) ~ "Quadrupod",
    x %in% c(1:11, 14:18) ~ "Other",
    TRUE ~ NA_character_
  )
}

clean_longest_session <- function(x) {
  dplyr::case_when(
    x %in% c("<15 min", "_15 min") ~ "<15 min",
    x %in% c("15-60 min", "15–60 min") ~ "15-60 min",
    x %in% c("1-2 h", "1–2 h", ">2 h", "_2 h") ~ ">60 min",
    TRUE ~ NA_character_
  )
}

# --- Product feature lookup ---

build_product_feature_lookup <- function(products_raw) {
  if (!"pocket_pen_bool" %in% names(products_raw)) {
    products_raw$pocket_pen_bool <- NA_real_
  }

  products_raw %>%
    dplyr::arrange(version, product_line_subgroup, product_line, product_family,
                   product_group, product_type, manufacturer) %>%
    dplyr::mutate(
      model_name = dplyr::coalesce(
        dplyr::na_if(trimws(product_line_subgroup), ""),
        dplyr::na_if(trimws(product_line), ""),
        dplyr::na_if(trimws(product_family), "")
      ),
      diameter_body_mm = readr::parse_number(as.character(diameter_body_mm)),
      diameter_grip_mm = readr::parse_number(as.character(diameter_grip_mm)),
      diameter_section_min_mm = readr::parse_number(as.character(diameter_section_min_mm)),
      diameter_section_max_mm = readr::parse_number(as.character(diameter_section_max_mm)),
      diameter_body_max_mm = readr::parse_number(as.character(diameter_body_max_mm)),
      diameter_max_mm = readr::parse_number(as.character(diameter_max_mm)),
      price = readr::parse_number(as.character(price))
    ) %>%
    dplyr::filter(multi_bool == 0 | is.na(multi_bool), product_type != "Wooden Pencils") %>%
    dplyr::mutate(
      price = as.numeric(price),
      price = dplyr::case_when(
        currency == "bps" ~ round(price * 1.30, 2),
        currency == "euro" ~ round(price * 1.15, 2),
        currency == "cad" ~ round(price * 0.75, 2),
        TRUE ~ price
      )
    ) %>%
    dplyr::group_by(
      manufacturer,
      product_category,
      material_family,
      product_type,
      product_group,
      product_family,
      product_line,
      product_line_subgroup,
      model_name
    ) %>%
    dplyr::summarise(
      total_records = dplyr::n(),
      total_sources = length(unique(source)),
      d_body_median = median(unique(as.numeric(diameter_body_mm)), na.rm = TRUE),
      d_grip_median = median(unique(as.numeric(diameter_grip_mm)), na.rm = TRUE),
      d_grip_min = median(unique(as.numeric(diameter_section_min_mm)), na.rm = TRUE),
      d_grip_max = median(unique(as.numeric(diameter_section_max_mm)), na.rm = TRUE),
      d_body_max_median = median(unique(as.numeric(diameter_body_max_mm)), na.rm = TRUE),
      d_max_median = median(unique(as.numeric(diameter_max_mm)), na.rm = TRUE),
      l_capped_median = median(unique(length_capped_mm), na.rm = TRUE),
      l_uncapped_median = median(unique(length_uncapped_mm), na.rm = TRUE),
      l_retracted_median = median(unique(length_retracted_mm), na.rm = TRUE),
      l_body_median = median(unique(length_body_mm), na.rm = TRUE),
      l_posted_median = median(unique(length_posted_mm), na.rm = TRUE),
      w_pen_empty_median = median(unique(weight_whole_pen_empty_g), na.rm = TRUE),
      w_cap_median = median(unique(weight_cap_g), na.rm = TRUE),
      w_pen_ink_median = median(unique(weight_whole_pen_with_ink_g), na.rm = TRUE),
      w_pencil_median = median(unique(weight_whole_pencil_w_lead_g), na.rm = TRUE),
      clippable = paste0(sort(unique(clippable)), collapse = ", "),
      retractable = paste0(sort(unique(retractable_bool)), collapse = ", "),
      retract_mechanism = paste0(sort(unique(retract_mechanism)), collapse = ", "),
      capped = paste0(sort(unique(capped_bool)), collapse = ", "),
      cap_type = paste0(sort(unique(cap_type)), collapse = ", "),
      postable = paste0(sort(unique(postable)), collapse = ", "),
      pocket_pen_bool = {
        pocket_vals <- unique(parse_capped_flag(pocket_pen_bool))
        pocket_vals <- pocket_vals[!is.na(pocket_vals)]
        if (length(pocket_vals) == 0) {
          NA_real_
        } else if (any(pocket_vals == 1)) {
          1
        } else {
          0
        }
      },
      model_key2 = paste0(sort(unique(model_key)), collapse = ", "),
      quantity_in_pack2 = paste0(sort(unique(quantity_in_pack)), collapse = ", "),
      filling_mechanism = paste0(sort(unique(filling_mechanism)), collapse = ", "),
      filler_permament = paste0(sort(unique(filler_permament_bool)), collapse = ", "),
      body_color = paste0(sort(unique(body_color)), collapse = ", "),
      item_number2 = paste0(sort(unique(item_number)), collapse = ", "),
      nib_size = paste0(sort(unique(nib_size)), collapse = ", "),
      tip_size = paste0(sort(unique(tip_size)), collapse = ", "),
      med_price = median(unique(price), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      med_diameter_1 = dplyr::coalesce(d_body_median, d_body_max_median, d_max_median),
      med_d_grip_1 = d_grip_median,
      med_d_grip_min_1 = d_grip_min,
      med_d_grip_max_1 = d_grip_max,
      med_length_1 = dplyr::coalesce(l_retracted_median, l_body_median, l_capped_median),
      med_length_posted_1 = l_posted_median,
      med_length_uncapped_1 = l_uncapped_median,
      med_weight_1 = dplyr::coalesce(w_pen_ink_median, w_pencil_median, w_pen_empty_median),
      med_weight_cap_1 = w_cap_median
    ) %>%
    dplyr::group_by(
      manufacturer, product_category, material_family, product_type,
      product_group, product_family, product_line, product_line_subgroup
    ) %>%
    dplyr::mutate(
      med_diameter_2 = median(unique(as.numeric(med_diameter_1)), na.rm = TRUE),
      med_d_grip_2 = median(unique(as.numeric(med_d_grip_1)), na.rm = TRUE),
      med_d_grip_min_2 = median(unique(as.numeric(med_d_grip_min_1)), na.rm = TRUE),
      med_d_grip_max_2 = median(unique(as.numeric(med_d_grip_max_1)), na.rm = TRUE),
      med_length_2 = median(unique(as.numeric(med_length_1)), na.rm = TRUE),
      med_length_posted_2 = median(unique(as.numeric(med_length_posted_1)), na.rm = TRUE),
      med_length_uncapped_2 = median(unique(as.numeric(med_length_uncapped_1)), na.rm = TRUE),
      med_weight_2 = median(unique(as.numeric(med_weight_1)), na.rm = TRUE),
      med_weight_cap_2 = median(unique(as.numeric(med_weight_cap_1)), na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(
      manufacturer, product_category, material_family, product_type,
      product_group, product_family, product_line
    ) %>%
    dplyr::mutate(
      med_diameter_3 = median(unique(as.numeric(med_diameter_1)), na.rm = TRUE),
      med_d_grip_3 = median(unique(as.numeric(med_d_grip_1)), na.rm = TRUE),
      med_d_grip_min_3 = median(unique(as.numeric(med_d_grip_min_1)), na.rm = TRUE),
      med_d_grip_max_3 = median(unique(as.numeric(med_d_grip_max_1)), na.rm = TRUE),
      med_length_3 = median(unique(as.numeric(med_length_1)), na.rm = TRUE),
      med_length_posted_3 = median(unique(as.numeric(med_length_posted_1)), na.rm = TRUE),
      med_length_uncapped_3 = median(unique(as.numeric(med_length_uncapped_1)), na.rm = TRUE),
      med_weight_3 = median(unique(as.numeric(med_weight_1)), na.rm = TRUE),
      med_weight_cap_3 = median(unique(as.numeric(med_weight_cap_1)), na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(
      manufacturer, product_category, material_family, product_type,
      product_group, product_family
    ) %>%
    dplyr::mutate(
      med_diameter_4 = median(unique(as.numeric(med_diameter_1)), na.rm = TRUE),
      med_d_grip_4 = median(unique(as.numeric(med_d_grip_1)), na.rm = TRUE),
      med_d_grip_min_4 = median(unique(as.numeric(med_d_grip_min_1)), na.rm = TRUE),
      med_d_grip_max_4 = median(unique(as.numeric(med_d_grip_max_1)), na.rm = TRUE),
      med_length_4 = median(unique(as.numeric(med_length_1)), na.rm = TRUE),
      med_length_posted_4 = median(unique(as.numeric(med_length_posted_1)), na.rm = TRUE),
      med_length_uncapped_4 = median(unique(as.numeric(med_length_uncapped_1)), na.rm = TRUE),
      med_weight_4 = median(unique(as.numeric(med_weight_1)), na.rm = TRUE),
      med_weight_cap_4 = median(unique(as.numeric(med_weight_cap_1)), na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      med_diameter = dplyr::coalesce(med_diameter_1, med_diameter_2, med_diameter_3, med_diameter_4),
      med_diameter_grip = dplyr::coalesce(med_d_grip_1, med_d_grip_2, med_d_grip_3, med_d_grip_4),
      med_diameter_grip_min = dplyr::coalesce(med_d_grip_min_1, med_d_grip_min_2, med_d_grip_min_3, med_d_grip_min_4),
      med_diameter_grip_max = dplyr::coalesce(med_d_grip_max_1, med_d_grip_max_2, med_d_grip_max_3, med_d_grip_max_4),
      med_length = dplyr::coalesce(med_length_1, med_length_2, med_length_3, med_length_4),
      med_length_posted = dplyr::coalesce(med_length_posted_1, med_length_posted_2, med_length_posted_3, med_length_posted_4),
      med_length_uncapped = dplyr::coalesce(med_length_uncapped_1, med_length_uncapped_2, med_length_uncapped_3, med_length_uncapped_4),
      med_weight = dplyr::coalesce(med_weight_1, med_weight_2, med_weight_3, med_weight_4),
      med_weight_cap = dplyr::coalesce(med_weight_cap_1, med_weight_cap_2, med_weight_cap_3, med_weight_cap_4)
    ) %>%
    dplyr::ungroup()
}

# --- CG proxy models ---

fit_jrg_proxy_models <- function(jrg) {
  list(
    fit_cgu = fit_lm_model(jrg, "CGu", c("LPn", "LCp", "LPs", "DBa", "D25", "WPn", "WCp")),
    fit_cgp = fit_ridge_model(jrg, "CGp", c("LPn", "LCp", "LPs", "DBa", "D25", "WPn", "WCp")),
    fit_cgpn = fit_lm_model(jrg, "CgPn", c("LPn", "LCp", "LPs", "DBa", "D25", "WPn")),
    fit_cgps = fit_lm_model(jrg, "CgPs", c("LPn", "LCp", "LPs", "DBa", "D25", "WPn", "WCp"))
  )
}

summarize_jrg_proxy_cv <- function(jrg) {
  bind_rows(
    cross_validate_regression(
      data = jrg,
      target = "CGu",
      predictors = c("LPn", "LCp", "LPs", "DBa", "D25", "WPn", "WCp"),
      model_name = "CGu_linear",
      model_type = "lm"
    ),
    cross_validate_regression(
      data = jrg,
      target = "CGp",
      predictors = c("LPn", "LCp", "LPs", "DBa", "D25", "WPn", "WCp"),
      model_name = "CGp_ridge",
      model_type = "ridge"
    ),
    cross_validate_regression(
      data = jrg,
      target = "CgPn",
      predictors = c("LPn", "LCp", "LPs", "DBa", "D25", "WPn"),
      model_name = "CgPn_linear",
      model_type = "lm"
    ),
    cross_validate_regression(
      data = jrg,
      target = "CgPs",
      predictors = c("LPn", "LCp", "LPs", "DBa", "D25", "WPn", "WCp"),
      model_name = "CgPs_linear",
      model_type = "lm"
    )
  ) %>%
    dplyr::group_by(model, target, model_type) %>%
    dplyr::summarise(
      folds_total = dplyr::n(),
      mean_RMSE = mean(RMSE, na.rm = TRUE),
      sd_RMSE = stats::sd(RMSE, na.rm = TRUE),
      mean_MAE = mean(MAE, na.rm = TRUE),
      sd_MAE = stats::sd(MAE, na.rm = TRUE),
      mean_R2 = mean(R2, na.rm = TRUE),
      sd_R2 = stats::sd(R2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(mean_RMSE)
}

# --- Feature building ---

build_features_2_final <- function(comfort_items_raw, product_lookup, proxy_models) {
  ratings <- comfort_items_raw %>%
    dplyr::mutate(
      material = dplyr::na_if(trimws(material), ""),
      instrument = dplyr::na_if(trimws(instrument), ""),
      join_model_name = paste(manufacturer, model_name, product_type)
    )

  c1 <- ratings %>%
    dplyr::left_join(
      product_lookup,
      by = c(
        "instrument" = "model_name",
        "material" = "material_family"
      )
    )

  c2 <- ratings %>%
    dplyr::left_join(
      product_lookup,
      by = c(
        "join_model_name" = "model_name",
        "material_other" = "material_family"
      )
    )

  joined_ratings <- c1 %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(setdiff(names(c2), names(ratings))),
        ~ dplyr::coalesce(.x, c2[[dplyr::cur_column()]])
      )
    ) %>%
    dplyr::select(-join_model_name)

  if (!"pocket_pen_bool" %in% names(joined_ratings)) {
    joined_ratings$pocket_pen_bool <- NA_real_
  }

  features_2 <- joined_ratings %>%
    dplyr::mutate(
      capped_flag = parse_capped_flag(capped),
      med_weight_uncapped = safe_num(med_weight) - safe_num(med_weight_cap)
    ) %>%
    dplyr::transmute(
      survey_id = survey_id,
      created_at_pt = created_at_pt,
      instrument = instrument,
      material = material,
      cap_posting = cap_posting,
      manufacturer = manufacturer.x,
      model_name = model_name,
      product_type_resp = product_type.x,
      material_other = material_other,
      comfort_rating = comfort_rating,
      product_type = product_type.y,
      capped = capped_flag,
      pocket_pen_bool = safe_num(pocket_pen_bool),
      filling_mechanism = filling_mechanism,
      total_records = total_records,
      l_writing = dplyr::case_when(
        capped_flag == 1 & cap_posting == "Yes" ~ safe_num(med_length_posted),
        capped_flag == 1 & cap_posting == "No" ~ safe_num(med_length_uncapped),
        capped_flag == 1 & cap_posting %in% c("", "I don't know", "Not applicable") ~ safe_num(med_length),
        capped_flag == 0 ~ dplyr::coalesce(safe_num(med_length_uncapped), safe_num(med_length)),
        TRUE ~ dplyr::coalesce(safe_num(med_length_uncapped), safe_num(med_length))
      ),
      w_writing = dplyr::case_when(
        capped_flag == 1 & cap_posting == "Yes" ~ safe_num(med_weight),
        capped_flag == 1 & cap_posting == "No" ~ safe_num(med_weight_uncapped),
        capped_flag == 1 & cap_posting %in% c("", "I don't know", "Not applicable") ~ safe_num(med_weight),
        capped_flag == 0 ~ safe_num(med_weight),
        TRUE ~ dplyr::coalesce(safe_num(med_weight_uncapped), safe_num(med_weight))
      ),
      med_diameter_grip = safe_num(med_diameter_grip),
      med_diameter_grip_min = safe_num(med_diameter_grip_min),
      med_diameter_grip_max = safe_num(med_diameter_grip_max),
      med_diameter = safe_num(med_diameter),
      med_length_uncapped = safe_num(med_length_uncapped),
      med_length_posted = safe_num(med_length_posted),
      med_weight_uncapped = safe_num(med_weight_uncapped),
      med_weight_cap = safe_num(med_weight_cap),
      med_length = safe_num(med_length),
      med_weight = safe_num(med_weight)
    )

  features_2_pred <- features_2 %>%
    dplyr::mutate(
      LPn_pred = dplyr::case_when(
        capped == 0 ~ safe_num(l_writing),
        TRUE ~ safe_num(med_length_uncapped)
      ),
      LCp_pred = dplyr::case_when(
        capped == 0 ~ safe_num(l_writing),
        TRUE ~ safe_num(med_length)
      ),
      LPs_pred = dplyr::case_when(
        capped == 0 ~ safe_num(l_writing),
        TRUE ~ safe_num(med_length_posted)
      ),
      DBa_pred = safe_num(med_diameter),
      D25_pred = safe_num(med_diameter),
      WPn_pred = dplyr::case_when(
        capped == 0 ~ safe_num(med_weight),
        TRUE ~ safe_num(med_weight_uncapped)
      ),
      WCp_pred = dplyr::case_when(
        capped == 0 ~ 0,
        TRUE ~ safe_num(med_weight_cap)
      )
    )

  pred_input_no_wcp <- features_2_pred %>%
    dplyr::transmute(
      LPn = LPn_pred,
      LCp = LCp_pred,
      LPs = LPs_pred,
      DBa = DBa_pred,
      D25 = D25_pred,
      WPn = WPn_pred
    )

  pred_input_with_wcp <- features_2_pred %>%
    dplyr::transmute(
      LPn = LPn_pred,
      LCp = LCp_pred,
      LPs = LPs_pred,
      DBa = DBa_pred,
      D25 = D25_pred,
      WPn = WPn_pred,
      WCp = WCp_pred
    )

  features_2_pred <- features_2_pred %>%
    dplyr::mutate(
      cgu_est = predict_lm_model(proxy_models$fit_cgu, pred_input_with_wcp),
      cgp_est = predict_ridge_model(proxy_models$fit_cgp, pred_input_with_wcp),
      cgpn_est = predict_lm_model(proxy_models$fit_cgpn, pred_input_no_wcp),
      cgps_est = predict_lm_model(proxy_models$fit_cgps, pred_input_with_wcp)
    ) %>%
    dplyr::mutate(
      cg_est_used_raw = dplyr::case_when(
        cap_posting == "Yes" ~ cgp_est,
        TRUE ~ cgu_est
      ),
      cg_est_used_rel = dplyr::case_when(
        cap_posting == "Yes" ~ cgps_est,
        TRUE ~ cgpn_est
      )
    )

  torque_uncapped_25 <- calc_torque(
    mass_g = features_2_pred$med_weight_uncapped,
    cg_mm = features_2_pred$cgu_est,
    grip_mm = 25
  )

  torque_uncapped_38 <- calc_torque(
    mass_g = features_2_pred$med_weight_uncapped,
    cg_mm = features_2_pred$cgu_est,
    grip_mm = 38
  )

  torque_posted_25 <- calc_torque(
    mass_g = features_2_pred$med_weight,
    cg_mm = features_2_pred$cgp_est,
    grip_mm = 25
  )

  torque_posted_38 <- calc_torque(
    mass_g = features_2_pred$med_weight,
    cg_mm = features_2_pred$cgp_est,
    grip_mm = 38
  )

  torque_used_25 <- calc_torque(
    mass_g = features_2_pred$w_writing,
    cg_mm = features_2_pred$cg_est_used_raw,
    grip_mm = 25
  )

  torque_used_38 <- calc_torque(
    mass_g = features_2_pred$w_writing,
    cg_mm = features_2_pred$cg_est_used_raw,
    grip_mm = 38
  )

  features_2_final <- features_2_pred %>%
    dplyr::bind_cols(
      torque_uncapped_25 %>% dplyr::rename_with(~ paste0("uncapped_25mm_", .x)),
      torque_uncapped_38 %>% dplyr::rename_with(~ paste0("uncapped_38mm_", .x)),
      torque_posted_25 %>% dplyr::rename_with(~ paste0("posted_25mm_", .x)),
      torque_posted_38 %>% dplyr::rename_with(~ paste0("posted_38mm_", .x)),
      torque_used_25 %>% dplyr::rename_with(~ paste0("used_25mm_", .x)),
      torque_used_38 %>% dplyr::rename_with(~ paste0("used_38mm_", .x))
    ) %>%
    dplyr::mutate(
      dplyr::across(
        c(
          cgp_est,
          cgps_est,
          med_length_posted,
          dplyr::starts_with("posted_25mm_")
        ),
        ~ dplyr::if_else(capped == 0, NA_real_, .x)
      )
    ) %>%
    dplyr::select(
      survey_id,
      created_at_pt,
      instrument,
      material,
      cap_posting,
      manufacturer,
      model_name,
      product_type_resp,
      material_other,
      comfort_rating,
      product_type,
      capped,
      pocket_pen_bool,
      filling_mechanism,
      total_records,
      l_writing,
      w_writing,
      med_diameter_grip,
      med_diameter_grip_min,
      med_diameter_grip_max,
      med_diameter,
      med_length_uncapped,
      med_length_posted,
      med_weight_uncapped,
      med_weight_cap,
      med_length,
      med_weight,
      cgu_est,
      cgp_est,
      cgpn_est,
      cgps_est,
      cg_est_used_raw,
      cg_est_used_rel,
      dplyr::starts_with("uncapped_25mm_"),
      dplyr::starts_with("posted_25mm_"),
      dplyr::starts_with("used_25mm_")
    )

  list(
    joined_ratings = joined_ratings,
    features_2 = features_2,
    features_2_pred = features_2_pred,
    features_2_final = features_2_final
  )
}

# --- Respondent clustering ---

name_activity_clusters <- function(responses) {
  activity_cluster_name_map <- responses %>%
    dplyr::distinct(activity_cluster, activity_cluster_label) %>%
    dplyr::mutate(
      activity_cluster_name = dplyr::case_when(
        stringr::str_detect(activity_cluster_label, "Academic work") ~ "Academic-work oriented",
        stringr::str_detect(activity_cluster_label, "Letter writing") ~ "Letter-writing oriented",
        stringr::str_detect(activity_cluster_label, "Professional writing|Creative / long-form writing") ~ "Expressive / professional multi-activity",
        stringr::str_detect(activity_cluster_label, "Journaling") ~ "Journaling / note-taking core",
        stringr::str_detect(activity_cluster_label, "Note taking") &
          stringr::str_detect(activity_cluster_label, "Planning / task organization") ~ "Functional note-taking / planning",
        TRUE ~ activity_cluster_label
      )
    )

  responses %>%
    dplyr::left_join(
      activity_cluster_name_map %>% dplyr::select(activity_cluster, activity_cluster_name),
      by = "activity_cluster"
    )
}

name_types_clusters <- function(responses) {
  types_cluster_name_map <- responses %>%
    dplyr::distinct(types_cluster, types_cluster_label) %>%
    dplyr::mutate(
      types_cluster_name = dplyr::case_when(
        stringr::str_detect(types_cluster_label, "Mechanical Pencil") &
          stringr::str_detect(types_cluster_label, "Fountain Pen") ~ "Mechanical pencil / precision users",
        stringr::str_detect(types_cluster_label, "Ballpoint Pen") &
          stringr::str_detect(types_cluster_label, "Fountain Pen") ~ "Ballpoint-forward traditional users",
        stringr::str_detect(types_cluster_label, "^Fountain Pen") ~ "Fountain pen focused users",
        stringr::str_detect(types_cluster_label, "Gel Pen") &
          stringr::str_detect(types_cluster_label, "Mechanical Pencil") &
          stringr::str_detect(types_cluster_label, "Ballpoint Pen") ~ "General-purpose multi-instrument users",
        TRUE ~ types_cluster_label
      )
    )

  responses %>%
    dplyr::left_join(
      types_cluster_name_map %>% dplyr::select(types_cluster, types_cluster_name),
      by = "types_cluster"
    )
}

build_responses_enriched <- function(responses_raw) {
  responses <- responses_raw %>%
    dplyr::mutate(respondent_id = dplyr::row_number())

  activities_fit <- fit_multiselect_clusters(
    data = responses,
    id_col = "respondent_id",
    value_col = "activities",
    cluster_col = "activity_cluster",
    label_col = "activity_cluster_label",
    min_prevalence = 2,
    k_values = 2:8,
    min_cluster_size = 8,
    silhouette_tolerance = 0.02,
    forced_k = NULL
  )

  responses <- responses %>%
    dplyr::select(-dplyr::any_of(c("activity_cluster", "activity_cluster_label", "activity_cluster_name"))) %>%
    dplyr::left_join(activities_fit$clusters, by = "respondent_id") %>%
    dplyr::left_join(activities_fit$labels, by = "activity_cluster") %>%
    name_activity_clusters()

  types_fit <- fit_multiselect_clusters(
    data = responses,
    id_col = "respondent_id",
    value_col = "types_used",
    cluster_col = "types_cluster",
    label_col = "types_cluster_label",
    min_prevalence = 2,
    k_values = 2:8,
    min_cluster_size = 8,
    silhouette_tolerance = 0.02,
    forced_k = 4
  )

  responses <- responses %>%
    dplyr::select(-dplyr::any_of(c("types_cluster", "types_cluster_label", "types_cluster_name"))) %>%
    dplyr::left_join(types_fit$clusters, by = "respondent_id") %>%
    dplyr::left_join(types_fit$labels, by = "types_cluster") %>%
    name_types_clusters()

  list(
    responses_enriched = responses,
    activities_fit = activities_fit,
    types_fit = types_fit
  )
}
