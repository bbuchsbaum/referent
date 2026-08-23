# Run prediction-level concordance against the checked-in PCNtoolkit fixture.

pcn_scenario_spec <- function(scenario) {
  switch(
    scenario,
    null_linear = referent::ref_spec(referent::ref_gaussian(), location = ~ x1, scale = ~ 1),
    linear_gaussian = referent::ref_spec(referent::ref_gaussian(), location = ~ x1, scale = ~ 1),
    multiple_encoded = referent::ref_spec(referent::ref_gaussian(), location = ~ x1 + x2 + sex_M, scale = ~ 1),
    nonlinear_heteroskedastic = referent::ref_spec(
      referent::ref_gaussian(), location = ~ s(x1, bs = "cr", k = 6),
      scale = ~ s(x1, bs = "cr", k = 5)
    ),
    log_linear = referent::ref_spec(
      referent::ref_gaussian(), location = ~ x1, scale = ~ 1, transform = "log"
    ),
    balanced_site = referent::ref_spec(
      referent::ref_gaussian(), location = ~ x1 + s(site, bs = "re"),
      scale = ~ s(site, bs = "re")
    ),
    skew_heavy = referent::ref_spec(
      referent::ref_shash(), location = ~ s(x1, bs = "cr", k = 6),
      scale = ~ s(x1, bs = "cr", k = 5), skew = ~ 1, tail = ~ 1
    ),
    unequal_site = referent::ref_spec(
      referent::ref_gaussian(), location = ~ x1 + s(site, bs = "re"),
      scale = ~ s(site, bs = "re")
    ),
    covariate_shift = referent::ref_spec(
      referent::ref_gaussian(), location = ~ x1, scale = ~ 1
    ),
    stop("unknown PCNtoolkit scenario: ", scenario, call. = FALSE)
  )
}

pcn_referent_predictions <- function(inputs, scenario, n_draw = 2000L,
                                     return_fit = FALSE) {
  dat <- inputs[inputs$scenario == scenario, , drop = FALSE]
  dat$site <- factor(dat$site)
  train <- dat[dat$split == "train", , drop = FALSE]
  test <- dat[dat$split == "test", , drop = FALSE]
  fit <- referent::ref_fit(pcn_scenario_spec(scenario), train, outcomes = "y")
  score <- stats::predict(
    fit, newdata = test, uncertainty = "total", n_draw = n_draw,
    allow_extrapolation = TRUE
  )
  dist <- stats::predict(
    fit, newdata = test, type = "distribution", uncertainty = "total",
    n_draw = n_draw, allow_extrapolation = TRUE
  )$y
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  quantiles <- vapply(
    probs, function(p) referent:::dist_quantile(dist, p), numeric(nrow(test))
  )
  colnames(quantiles) <- paste0("q", c("05", "25", "50", "75", "95"))
  predictions <- data.frame(
    scenario = scenario,
    row_id = test$row_id,
    observed = test$y,
    median = score$median,
    predictive_sd = (quantiles[, "q95"] - quantiles[, "q05"]) /
      (2 * stats::qnorm(0.95)),
    centile = score$centile,
    z = score$z,
    log_density = score$log_density,
    quantiles,
    check.names = FALSE
  )
  if (isTRUE(return_fit)) {
    return(list(predictions = predictions, fit = fit))
  }
  predictions
}

pcn_concordance_registry <- function() {
  data.frame(
    scenario = c(
      "null_linear", "linear_gaussian", "multiple_encoded", "log_linear",
      "covariate_shift", "nonlinear_heteroskedastic", "balanced_site",
      "skew_heavy", "unequal_site"
    ),
    evidence_class = c(
      rep("matched_estimator", 5), rep("same_estimand", 4)
    ),
    reason = c(
      rep("Gaussian BLR and linear mgcv model on the same design", 3),
      "Gaussian linear model on log(y), scored back on the response scale",
      "Gaussian linear model under a locked covariate-shifted test distribution",
      "B-spline and penalised smooth plus variance estimators differ",
      "PCN fixed batch effects and Referent penalised site effects differ",
      "Gaussian BLR is compared with a selected SHASH distributional model",
      "Fixed batch effects and penalised site effects differ under unequal site sizes"
    ),
    stringsAsFactors = FALSE
  )
}

pcn_run_concordance <- function(root, n_draw = 2000L) {
  inputs <- utils::read.csv(file.path(root, "fitted_inputs.csv"), stringsAsFactors = FALSE)
  comparator <- utils::read.csv(
    file.path(root, "fitted_predictions.csv"), stringsAsFactors = FALSE
  )
  registry <- pcn_concordance_registry()
  rows <- lapply(registry$scenario, function(scenario) {
    referent <- pcn_referent_predictions(inputs, scenario, n_draw = n_draw)
    pcn <- comparator[comparator$scenario == scenario, , drop = FALSE]
    pcn <- pcn[match(referent$row_id, pcn$row_id), , drop = FALSE]
    if (anyNA(pcn$row_id) || !identical(referent$row_id, pcn$row_id)) {
      stop("row identity mismatch in scenario ", scenario, call. = FALSE)
    }
    core <- pcn_concordance(referent, pcn, stats::sd(referent$observed))
    q_delta <- as.matrix(referent[paste0("q", c("05", "25", "50", "75", "95"))]) -
      as.matrix(pcn[paste0("q", c("05", "25", "50", "75", "95"))])
    data.frame(
      scenario = scenario,
      evidence_class = registry$evidence_class[registry$scenario == scenario],
      n_test = nrow(referent),
      as.list(core),
      quantile_nrmse = sqrt(mean(q_delta^2)) / stats::sd(referent$observed),
      mean_log_score_difference = mean(referent$log_density - pcn$log_density),
      referent_coverage90 = mean(
        referent$observed >= referent$q05 & referent$observed <= referent$q95
      ),
      pcntoolkit_coverage90 = mean(pcn$observed >= pcn$q05 & pcn$observed <= pcn$q95),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

pcn_concordance_pass <- function(results) {
  matched <- results$evidence_class == "matched_estimator"
  checks <- cbind(
    median_nrmse = results$median_nrmse <= 0.01,
    predictive_sd_nrmse = results$predictive_sd_nrmse <= 0.02,
    z_rmse = results$z_rmse <= 0.03,
    z_abs_q95 = results$z_abs_q95 <= 0.08,
    coverage90_difference = results$coverage90_difference <= 0.01
  )
  rowSums(!checks[matched, , drop = FALSE]) == 0L
}

if (sys.nframe() == 0L) {
  source("tools/pcntoolkit/compare_results.R")
  root <- pcn_fixture_root()
  pcn_validate_fixture(root)
  result <- pcn_run_concordance(root)
  print(result, row.names = FALSE)
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 2L && identical(args[[1L]], "--output")) {
    dir.create(dirname(args[[2L]]), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(result, args[[2L]], row.names = FALSE)
  } else if (length(args)) {
    stop("usage: run_concordance.R [--output FILE]")
  }
  if (!all(pcn_concordance_pass(result))) quit(status = 1L)
}
