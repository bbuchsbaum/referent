#' Assess a fitted reference model
#'
#' Reports overall probabilistic fit, marginal calibration, conditional
#' calibration, and tail calibration on held-out data.
#'
#' @details
#' `overall` maps onto the PCNtoolkit evaluation metrics as follows:
#' `standardized_log_score` is the negative MSLL: the mean log score
#' minus the log score of a Gaussian with the *reference* sample's mean
#' and population standard deviation (`fit$reference_baseline`, recorded
#' by [norm_fit()]). The baseline never uses the held-out sample's own
#' moments.
#'
#' Note that this is *not* the same baseline PCNtoolkit uses. PCNtoolkit
#' fits its baseline Gaussian to the sample being scored, so its MSLL is
#' computed against a held-out oracle and the two numbers differ by the
#' log-score gap between the two baselines (beyond the sign convention).
#' Compare against PCNtoolkit only after rescoring both sides on one
#' explicitly chosen baseline.
#'
#' `cor` is the Pearson correlation of observed and predicted median;
#' PCNtoolkit's `Rho` column is a Spearman correlation, so the two are
#' related but not interchangeable. `smse` is the
#' standardised mean squared error \eqn{\mathrm{MSE}/\mathrm{var}(y)};
#' `ev` is explained variance \eqn{1 - \mathrm{var}(y - \hat y)/\mathrm{var}(y)}.
#' `rmse`, `mae`, `mean_log_score`, and `crps` are the usual proper and
#' point scores.
#'
#' `conditional` fits light diagnostic GAMs of \eqn{z} and \eqn{z^2 - 1}
#' against every numeric covariate the model uses and reports the
#' largest fitted drift with the GAM standard error at that point.
#'
#' @param fit A [norm_fit].
#' @param newdata Held-out validation data (required; training data are
#'   not stored on the fit). For out-of-fold evaluation use
#'   [norm_crossfit()].
#' @param by Optional grouping column; the marginal calibration table is
#'   then reported per group (column `.group`).
#' @param uncertainty Passed to [predict.norm_fit()]. `"conditional"`
#'   (the default, and the historical behaviour) scores the point-estimate
#'   predictive; `"total"` integrates over the coefficient draws, matching
#'   what [predict.norm_fit()] itself returns by default. Without this
#'   argument the scores in `overall` cannot be reproduced from
#'   `predict()`.
#' @return An object of class `norm_assessment` with `overall`,
#'   `marginal`, `conditional`, and `tail` tibbles plus the scores.
#' @export
norm_assess <- function(fit, newdata, by = NULL,
                        uncertainty = c("conditional", "total")) {
  if (missing(newdata) || is.null(newdata)) {
    cli::cli_abort("Supply {.arg newdata}: held-out data or a cross-fitted frame.")
  }
  uncertainty <- match.arg(uncertainty)
  newdata <- tibble::as_tibble(newdata)
  by_vec <- pull_column(newdata, rlang::enquo(by), default = NULL)
  dists <- predict_dists(fit, newdata, uncertainty = uncertainty)
  scores <- scores_from_dists(fit, dists, newdata, allow_extrapolation = FALSE)
  assess_from_scores(fit, scores, dists, newdata, by_vec)
}

assess_from_scores <- function(fit, scores, dists, newdata, by_vec = NULL) {
  overall <- assess_overall(scores, dists, newdata, fit$reference_baseline)
  marginal <- if (is.null(by_vec)) {
    assess_marginal(scores)
  } else {
    grp <- as.character(by_vec)[scores$.row]
    parts <- lapply(split(seq_len(nrow(scores)), grp), function(idx) {
      m <- assess_marginal(scores[idx, , drop = FALSE])
      tibble::tibble(.group = grp[[idx[[1L]]]], m)
    })
    dplyr_bind(parts)
  }
  structure(
    list(
      overall = overall,
      marginal = marginal,
      conditional = assess_conditional(scores, newdata, fit),
      tail = assess_tail(scores),
      scores = scores,
      n = nrow(newdata),
      in_sample = nrow(scores) > 0 && all(scores$.in_sample)
    ),
    class = "norm_assessment"
  )
}

# Apply `fun(outcome_name, outcome_scores)` per outcome and bind the rows.
by_outcome <- function(scores, fun) {
  by_out <- split(scores, scores$.outcome)
  dplyr_bind(lapply(names(by_out), function(nm) fun(nm, by_out[[nm]])))
}

assess_overall <- function(scores, dists, newdata, baseline = NULL) {
  by_outcome(scores, function(nm, sc) {
    ok <- is.finite(sc$log_density) & is.finite(sc$observed)
    log_score <- mean(sc$log_density[ok])
    naive <- naive_log_score(sc$observed[ok], baseline[[nm]])
    crps <- if (is.null(dists[[nm]])) {
      NA_real_
    } else {
      mean(crps_from_dist(dists[[nm]], sc$observed), na.rm = TRUE)
    }
    var_y <- if (sum(ok) > 1) stats::var(sc$observed[ok]) else NA_real_
    tibble::tibble(
      .outcome = nm,
      mean_log_score = log_score,
      standardized_log_score = log_score - naive,
      crps = crps,
      mae = mean(abs(sc$residual[ok])),
      rmse = sqrt(mean(sc$residual[ok]^2)),
      smse = mean(sc$residual[ok]^2) / var_y,
      ev = 1 - stats::var(sc$residual[ok]) / var_y,
      cor = if (sum(ok) > 2) stats::cor(sc$observed[ok], sc$median[ok]) else NA_real_
    )
  })
}

# Log score of the unconditional Gaussian baseline: the reference
# sample's mean and population sd, or NA when the fit has none (an
# outcome without spread in the reference).
naive_log_score <- function(y, base = NULL) {
  if (is.null(base)) {
    return(NA_real_)
  }
  mean(stats::dnorm(y, base$mean, base$sd, log = TRUE), na.rm = TRUE)
}

assess_marginal <- function(scores) {
  by_outcome(scores, function(nm, sc) {
    z <- sc$z[is.finite(sc$z)]
    u <- sc$centile[is.finite(sc$centile)]
    tibble::tibble(
      .outcome = nm,
      n = length(z),
      mean_z = mean(z),
      var_z = stats::var(z),
      skew_z = if (length(z) > 3) std_moment(z, 3) else NA_real_,
      excess_kurtosis = if (length(z) > 4) std_moment(z, 4) - 3 else NA_real_,
      cover_50 = mean(u > 0.25 & u < 0.75, na.rm = TRUE),
      cover_80 = mean(u > 0.10 & u < 0.90, na.rm = TRUE),
      cover_90 = mean(u > 0.05 & u < 0.95, na.rm = TRUE),
      cover_95 = mean(u > 0.025 & u < 0.975, na.rm = TRUE),
      cover_99 = mean(u > 0.005 & u < 0.995, na.rm = TRUE)
    )
  })
}

# k-th standardised moment; 0 when x has no spread.
std_moment <- function(x, k) {
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) {
    return(0)
  }
  mean(((x - mean(x)) / s)^k)
}

assess_conditional <- function(scores, newdata, fit) {
  covs <- intersect(fit$covariates, names(newdata))
  numeric_covs <- covs[vapply(newdata[covs], is.numeric, logical(1))]
  rows <- list()
  for (nm in unique(scores$.outcome)) {
    sc <- scores[scores$.outcome == nm, , drop = FALSE]
    for (cv in numeric_covs) {
      x <- newdata[[cv]][sc$.row]
      ok <- is.finite(sc$z) & is.finite(x)
      dat <- data.frame(z = sc$z[ok], x = x[ok])
      loc <- drift_gam(z ~ s(x, k = 5), dat)
      sc2 <- drift_gam(I(z^2 - 1) ~ s(x, k = 5), dat)
      rows[[length(rows) + 1L]] <- tibble::tibble(
        .outcome = nm, covariate = cv,
        location_drift = loc$drift, location_se = loc$se,
        scale_drift = sc2$drift, scale_se = sc2$se
      )
    }
  }
  if (!length(rows)) {
    return(tibble::tibble(
      .outcome = character(), covariate = character(),
      location_drift = numeric(), location_se = numeric(),
      scale_drift = numeric(), scale_se = numeric()
    ))
  }
  dplyr_bind(rows)
}

# Largest absolute fitted value of a diagnostic GAM and its SE there; NA
# with fewer than 20 usable rows or when the GAM fails.
drift_gam <- function(formula, dat) {
  m <- if (nrow(dat) < 20L) {
    NULL
  } else {
    tryCatch(mgcv::gam(formula, data = dat), error = function(e) NULL)
  }
  if (is.null(m)) {
    return(list(drift = NA_real_, se = NA_real_))
  }
  pr <- stats::predict(m, se.fit = TRUE)
  i <- which.max(abs(pr$fit))
  list(drift = abs(pr$fit[[i]]), se = as.numeric(pr$se.fit[[i]]))
}

assess_tail <- function(scores) {
  levels <- c(0.005, 0.01, 0.025, 0.05)
  by_outcome(scores, function(nm, sc) {
    tp <- sc$tail_prob[is.finite(sc$tail_prob)]
    n <- length(tp)
    tibble::tibble(
      .outcome = nm,
      tail_level = levels,
      expected = levels,
      observed = vapply(levels, function(a) mean(tp <= a), numeric(1)),
      n = n,
      se = sqrt(levels * (1 - levels) / pmax(n, 1))
    )
  })
}

#' @export
print.norm_assessment <- function(x, ...) {
  cli::cli_text("{.cls norm_assessment} n = {x$n}")
  print(x$overall)
  invisible(x)
}
