# R/utils_models.R — GLMM fitting, summarization, and coefficient helpers


# --- Respondent covariate preparation ---

prepare_respondent_covariates <- function(responses_enriched, cluster_assignments = NULL) {
  out <- responses_enriched %>%
    dplyr::transmute(
      survey_id,
      press_force = clean_press_force(press_force),
      hand_size = clean_hand_size(hand_size),
      handwriting_size = clean_handwriting_size(handwriting_size),
      balance_preference = clean_balance_pref(weight_pref),
      grip_posture = clean_grip_posture(grip_posture),
      longest_session = clean_longest_session(longest_session),
      activity_cluster = activity_cluster_name,
      types_cluster = types_cluster_name
    )

  if (!is.null(cluster_assignments)) {
    out <- out %>%
      dplyr::left_join(cluster_assignments %>% dplyr::select(survey_id, cluster), by = "survey_id")
  } else {
    out$cluster <- NA_character_
  }

  out %>%
    dplyr::mutate(
      press_force = factor(press_force, levels = c("Light", "Moderate/Heavy")),
      hand_size = factor(hand_size, levels = c("Small", "Medium", "Large")),
      handwriting_size = factor(handwriting_size, levels = c("Small", "Medium", "Large")),
      balance_preference = factor(balance_preference, levels = c("Front-heavy", "Balanced", "Back-heavy")),
      grip_posture = factor(grip_posture, levels = c("Tripod", "Quadrupod", "Other")),
      longest_session = factor(longest_session, levels = c("<15 min", "15-60 min", ">60 min")),
      activity_cluster = factor(
        activity_cluster,
        levels = c(
          "Journaling / note-taking core",
          "Expressive / professional multi-activity",
          "Letter-writing oriented",
          "Academic-work oriented",
          "Functional note-taking / planning"
        )
      ),
      types_cluster = factor(
        types_cluster,
        levels = c(
          "General-purpose multi-instrument users",
          "Fountain pen focused users",
          "Ballpoint-forward traditional users",
          "Mechanical pencil / precision users"
        )
      ),
      cluster = factor(
        cluster,
        levels = c("Smoothness-Focused Writers", "Fit-Focused Writers")
      )
    )
}

prep_model_data <- function(data) {
  data %>%
    dplyr::mutate(
      press_force = factor(press_force, levels = c("Light", "Moderate/Heavy")),
      hand_size = factor(hand_size, levels = c("Small", "Medium", "Large")),
      handwriting_size = factor(handwriting_size, levels = c("Small", "Medium", "Large")),
      balance_preference = factor(balance_preference, levels = c("Front-heavy", "Balanced", "Back-heavy")),
      grip_posture = factor(grip_posture, levels = c("Tripod", "Quadrupod", "Other")),
      longest_session = factor(longest_session, levels = c("<15 min", "15-60 min", ">60 min")),
      diameter_z = as.numeric(scale(diameter)),
      weight_z = as.numeric(scale(weight)),
      length_z = as.numeric(scale(length)),
      cg_abs_used_z = as.numeric(scale(cg_abs_used)),
      cg_rel_used_z = as.numeric(scale(cg_rel_used)),
      torque_used_z = as.numeric(scale(torque_used))
    )
}


# --- GLMM fitting ---

fit_glmer_safe <- function(formula, data, control) {
  warnings_seen <- character(0)

  fit <- withCallingHandlers(
    tryCatch(
      lme4::glmer(
        formula = formula,
        data = data,
        family = stats::binomial,
        control = control
      ),
      error = function(e) e
    ),
    warning = function(w) {
      warnings_seen <<- unique(c(warnings_seen, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )

  list(
    fit = fit,
    warnings = warnings_seen,
    formula = paste(deparse(formula), collapse = " ")
  )
}


# --- GLMM metrics and summarization ---

get_random_intercept_var <- function(fit) {
  vc <- tryCatch(as.data.frame(lme4::VarCorr(fit)), error = function(e) NULL)
  if (is.null(vc)) return(NA_real_)
  out <- vc %>%
    dplyr::filter(grp == "survey_id", var1 == "(Intercept)") %>%
    dplyr::pull(vcov)
  if (length(out) == 0) NA_real_ else out[1]
}

binary_glmm_metrics <- function(fit, data, outcome_name, threshold = 0.5, re_form = NA) {
  truth <- data[[outcome_name]]

  prob <- tryCatch(
    stats::predict(
      fit,
      newdata = data,
      type = "response",
      re.form = re_form,
      allow.new.levels = TRUE
    ),
    error = function(e) rep(NA_real_, nrow(data))
  )

  keep <- stats::complete.cases(truth, prob)
  if (!any(keep)) {
    return(tibble::tibble(
      auc = NA_real_,
      accuracy = NA_real_,
      misclassification_rate = NA_real_
    ))
  }

  truth_vec <- as.integer(truth[keep])
  prob_vec <- as.numeric(prob[keep])
  pred_class <- ifelse(prob_vec >= threshold, 1L, 0L)

  auc_val <- tryCatch(
    as.numeric(pROC::roc(response = truth_vec, predictor = prob_vec, quiet = TRUE)$auc),
    error = function(e) NA_real_
  )

  acc_val <- mean(pred_class == truth_vec)
  mis_val <- mean(pred_class != truth_vec)

  tibble::tibble(
    auc = auc_val,
    accuracy = acc_val,
    misclassification_rate = mis_val
  )
}

summarize_glmer_fit <- function(fit_obj, model_name, dataset_name, outcome_name, data) {
  if (!inherits(fit_obj$fit, "merMod")) {
    return(tibble::tibble(
      dataset = dataset_name,
      outcome = outcome_name,
      model = model_name,
      formula = fit_obj$formula,
      n = nrow(data),
      respondents = dplyr::n_distinct(data$survey_id),
      AIC = NA_real_,
      BIC = NA_real_,
      logLik = NA_real_,
      auc = NA_real_,
      accuracy = NA_real_,
      misclassification_rate = NA_real_,
      R2_marginal = NA_real_,
      R2_conditional = NA_real_,
      ICC = NA_real_,
      overdispersion_ratio = NA_real_,
      overdispersion_p = NA_real_,
      random_intercept_var = NA_real_,
      singular = NA,
      has_conv_msg = NA,
      conv_msg = NA_character_,
      warning_msg = if (length(fit_obj$warnings)) paste(fit_obj$warnings, collapse = " | ") else NA_character_,
      error = conditionMessage(fit_obj$fit)
    ))
  }

  fit <- fit_obj$fit
  r2_obj <- tryCatch(performance::r2_nakagawa(fit), error = function(e) NULL)
  icc_obj <- tryCatch(performance::icc(fit), error = function(e) NULL)
  od_obj <- tryCatch(performance::check_overdispersion(fit), error = function(e) NULL)
  conv_msg_raw <- unlist(fit@optinfo$conv$lme4$messages)

  class_metrics <- binary_glmm_metrics(
    fit = fit,
    data = data,
    outcome_name = outcome_name,
    threshold = 0.5,
    re_form = NA
  )

  tibble::tibble(
    dataset = dataset_name,
    outcome = outcome_name,
    model = model_name,
    formula = fit_obj$formula,
    n = stats::nobs(fit),
    respondents = nlevels(lme4::getME(fit, "flist")[[1]]),
    AIC = stats::AIC(fit),
    BIC = stats::BIC(fit),
    logLik = as.numeric(stats::logLik(fit)),
    auc = class_metrics$auc,
    accuracy = class_metrics$accuracy,
    misclassification_rate = class_metrics$misclassification_rate,
    R2_marginal = if (!is.null(r2_obj)) safe_scalar(r2_obj$R2_marginal) else NA_real_,
    R2_conditional = if (!is.null(r2_obj)) safe_scalar(r2_obj$R2_conditional) else NA_real_,
    ICC = if (!is.null(icc_obj)) safe_scalar(icc_obj$ICC_adjusted) else NA_real_,
    overdispersion_ratio = if (!is.null(od_obj)) safe_scalar(od_obj$dispersion_ratio) else NA_real_,
    overdispersion_p = if (!is.null(od_obj)) safe_scalar(od_obj$p_value) else NA_real_,
    random_intercept_var = get_random_intercept_var(fit),
    singular = lme4::isSingular(fit, tol = 1e-4),
    has_conv_msg = length(conv_msg_raw) > 0,
    conv_msg = if (length(conv_msg_raw) > 0) paste(conv_msg_raw, collapse = " | ") else NA_character_,
    warning_msg = if (length(fit_obj$warnings)) paste(fit_obj$warnings, collapse = " | ") else NA_character_,
    error = NA_character_
  )
}


# --- Coefficient extraction ---

extract_glmer_coefs <- function(fit_obj, model_name, dataset_name, outcome_name) {
  if (!inherits(fit_obj$fit, "merMod")) return(tibble::tibble())

  sm <- as.data.frame(stats::coef(summary(fit_obj$fit)))
  sm$term <- rownames(sm)
  rownames(sm) <- NULL

  tibble::tibble(
    dataset = dataset_name,
    outcome = outcome_name,
    model = model_name,
    term = sm$term,
    estimate = as.numeric(sm$Estimate),
    std_error = as.numeric(sm$`Std. Error`),
    statistic = as.numeric(sm$`z value`),
    p_value = as.numeric(sm$`Pr(>|z|)`),
    odds_ratio = exp(as.numeric(sm$Estimate))
  )
}

extract_collinearity_tbl <- function(fit, model_name, dataset_name, outcome_name) {
  out <- tryCatch(as_tibble(performance::check_collinearity(fit)), error = function(e) NULL)
  if (is.null(out)) return(tibble::tibble())

  out %>%
    dplyr::mutate(
      dataset = dataset_name,
      outcome = outcome_name,
      model = model_name,
      .before = 1
    )
}


# --- Model comparison ---

fit_model_set <- function(data, outcome, rhs_list, dataset_name, control) {
  fits <- purrr::imap(rhs_list, function(rhs, nm) {
    fit_glmer_safe(stats::as.formula(paste(outcome, "~", rhs)), data, control)
  })

  comparison <- purrr::imap_dfr(
    fits,
    ~ summarize_glmer_fit(.x, .y, dataset_name, outcome, data)
  ) %>%
    dplyr::arrange(AIC)

  coefficients <- purrr::imap_dfr(
    fits,
    ~ extract_glmer_coefs(.x, .y, dataset_name, outcome)
  )

  list(
    fits = fits,
    comparison = comparison,
    coefficients = coefficients
  )
}

make_lrt_table <- function(fits, model_names, dataset_name, outcome_name, comparison_name) {
  keep_models <- model_names[vapply(model_names, function(m) inherits(fits[[m]]$fit, "merMod"), logical(1))]
  if (length(keep_models) < 2) return(tibble::tibble())

  an <- tryCatch(
    do.call(
      stats::anova,
      c(
        lapply(keep_models, function(m) fits[[m]]$fit),
        list(test = "Chisq")
      )
    ),
    error = function(e) NULL
  )

  if (is.null(an)) return(tibble::tibble())

  as.data.frame(an) %>%
    tibble::rownames_to_column("model") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      dataset = dataset_name,
      outcome = outcome_name,
      comparison = comparison_name,
      .before = 1
    )
}

tidy_lrt <- function(model_names, ...) {
  an <- as.data.frame(do.call(stats::anova, c(list(...), list(test = "Chisq"))))
  tibble::tibble(
    model = model_names[seq_len(nrow(an))],
    npar = an$npar,
    AIC = an$AIC,
    BIC = an$BIC,
    logLik = an$logLik,
    deviance = an$deviance,
    Chisq = if ("Chisq" %in% names(an)) an$Chisq else NA_real_,
    Df = if ("Df" %in% names(an)) an$Df else NA_real_,
    p_value = if ("Pr(>Chisq)" %in% names(an)) an$`Pr(>Chisq)` else NA_real_
  )
}


# --- Simple GLMM helpers (bridge analysis) ---

coef_glmm_simple <- function(fit, model_name) {
  sm <- as.data.frame(coef(summary(fit)))
  sm$term <- rownames(sm)
  rownames(sm) <- NULL
  tibble(
    model = model_name,
    term = sm$term,
    estimate = sm$Estimate,
    std_error = sm$`Std. Error`,
    statistic = sm$`z value`,
    p_value = sm$`Pr(>|z|)`,
    odds_ratio = exp(sm$Estimate)
  )
}

summarize_glmm_simple <- function(fit, data, outcome_name, model_name, panel_name) {
  class_metrics <- binary_glmm_metrics(
    fit = fit,
    data = data,
    outcome_name = outcome_name,
    threshold = 0.5,
    re_form = NA
  )

  tibble(
    panel = panel_name,
    model = model_name,
    n = nobs(fit),
    respondents = nlevels(getME(fit, "flist")[[1]]),
    AIC = AIC(fit),
    BIC = BIC(fit),
    logLik = as.numeric(logLik(fit)),
    auc = class_metrics$auc,
    accuracy = class_metrics$accuracy,
    misclassification_rate = class_metrics$misclassification_rate,
    R2_marginal = safe_scalar(performance::r2_nakagawa(fit)$R2_marginal),
    R2_conditional = safe_scalar(performance::r2_nakagawa(fit)$R2_conditional)
  )
}

extract_glmm_simple_coefs <- function(fit, model_name, panel_name) {
  sm <- as.data.frame(coef(summary(fit)))
  sm$term <- rownames(sm)
  rownames(sm) <- NULL
  tibble(
    panel = panel_name,
    model = model_name,
    term = sm$term,
    estimate = sm$Estimate,
    std_error = sm$`Std. Error`,
    statistic = sm$`z value`,
    p_value = sm$`Pr(>|z|)`,
    odds_ratio = exp(sm$Estimate)
  )
}
