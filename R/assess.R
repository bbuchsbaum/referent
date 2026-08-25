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
#' by [ref_fit()]). The baseline never uses the held-out sample's own
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
#' `marginal` reports the moments of \eqn{z} (`mean_z`, `var_z`,
#' `skew_z`, `excess_kurtosis_z`), the Shapiro-Wilk statistic `shapiro_w`
#' of \eqn{z} (on a fixed-seed subsample of 5000 when \eqn{n} is larger),
#' the coverage of central intervals, and `mace`, PCNtoolkit's mean
#' absolute centile error: the mean over the centiles 0.05, 0.25, 0.5,
#' 0.75, and 0.95 of the absolute difference between the centile and the
#' fraction of observations at or below its curve.
#'
#' `conditional` checks the scores against every covariate the model
#' uses. For a numeric covariate it fits light diagnostic GAMs of
#' \eqn{z} and \eqn{z^2 - 1} against it and reports the largest absolute
#' fitted value with its pointwise standard error (`location_drift`,
#' `location_se`, `scale_drift`, `scale_se`; descriptive amplitudes, since
#' the maximum is selected) and the approximate p-value of the smooth
#' term (`location_p`, `scale_p`), which is the test of no drift. For a
#' factor covariate it reports one row per `level` with the mean of
#' \eqn{z} and the excess of \eqn{\mathrm{var}(z)} over 1, their standard
#' errors, and two-sided p-values.
#'
#' @param fit A [ref_fit].
#' @param newdata Held-out validation data (required; training data are
#'   not stored on the fit). For out-of-fold evaluation use
#'   [ref_crossfit()].
#' @param by Optional grouping column; the marginal calibration table is
#'   then reported per group (column `.group`).
#' @param uncertainty Passed to [predict.ref_fit()]. `"conditional"`
#'   scores the point-estimate predictive; `"total"` (the default) integrates
#'   over coefficient uncertainty and matches [predict.ref_fit()].
#' @param allow_in_sample Set to `TRUE` only for explicitly labelled training
#'   diagnostics. The default aborts if `newdata` exactly matches the fitting
#'   data; use [ref_crossfit()] for honest reference-sample assessment.
#' @param allow_calibration_reuse Set to `TRUE` only for an explicitly labelled
#'   diagnostic on the same rows used to estimate a calibration map. The
#'   default requires evaluation data independent of fitting and calibration.
#' @return An object of class `ref_assessment` with `overall`,
#'   `marginal`, `conditional`, and `tail` tibbles plus the scores.
#' @export
ref_assess <- function(fit, newdata, by = NULL,
                        uncertainty = c("total", "conditional"),
                        allow_in_sample = FALSE,
                        allow_calibration_reuse = FALSE) {
  if (missing(newdata) || is.null(newdata)) {
    cli::cli_abort("Supply {.arg newdata}: held-out data or a cross-fitted frame.")
  }
  uncertainty <- match.arg(uncertainty)
  newdata <- tibble::as_tibble(newdata)
  in_sample <- is_in_sample_data(fit, newdata)
  if (in_sample && !isTRUE(allow_in_sample)) {
    cli::cli_abort(
      "{.arg newdata} are the model's training data; use {.fn ref_crossfit} or set {.arg allow_in_sample = TRUE} for a labelled diagnostic."
    )
  }
  if (is_calibration_data(fit, newdata) && !isTRUE(allow_calibration_reuse)) {
    cli::cli_abort(
      "{.arg newdata} are the calibration data; use an independent evaluation set or set {.arg allow_calibration_reuse = TRUE} for a labelled diagnostic."
    )
  }
  by_vec <- pull_column(newdata, rlang::enquo(by), default = NULL)
  dists <- predict_dists(fit, newdata, uncertainty = uncertainty)
  scores <- scores_from_dists(
    fit, dists, newdata, allow_extrapolation = FALSE,
    in_sample = in_sample
  )
  assess_from_scores(fit, scores, dists, newdata, by_vec)
}

is_calibration_data <- function(fit, newdata) {
  cal <- fit$calibration
  cols <- unique(c(fit$covariates, fit$outcomes))
  if (is.null(cal) || is.null(cal$data_hash) || !all(cols %in% names(newdata)) ||
      nrow(newdata) != cal$n) {
    return(FALSE)
  }
  hash_version <- as.integer(cal$data_hash_version %||% 3L)
  identical(
    digest_data(newdata[, cols, drop = FALSE], hash_version),
    cal$data_hash
  )
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
    class = "ref_assessment"
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
      excess_kurtosis_z = if (length(z) > 4) std_moment(z, 4) - 3 else NA_real_,
      shapiro_w = shapiro_w(z),
      mace = mace(u),
      cover_50 = mean(u > 0.25 & u < 0.75, na.rm = TRUE),
      cover_80 = mean(u > 0.10 & u < 0.90, na.rm = TRUE),
      cover_90 = mean(u > 0.05 & u < 0.95, na.rm = TRUE),
      cover_95 = mean(u > 0.025 & u < 0.975, na.rm = TRUE),
      cover_99 = mean(u > 0.005 & u < 0.995, na.rm = TRUE)
    )
  })
}

# Mean absolute centile error (PCNtoolkit): the mean over the centile
# grid of |q - (fraction of observations at or below the q-th centile
# curve)|, which for PIT values u is |q - mean(u <= q)|.
mace <- function(u, grid = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  if (!length(u)) {
    return(NA_real_)
  }
  mean(vapply(grid, function(q) abs(q - mean(u <= q)), numeric(1)))
}

# Shapiro-Wilk W of z. stats::shapiro.test() accepts at most 5000
# values, so larger samples are subsampled once with a fixed seed.
shapiro_w <- function(z) {
  if (length(z) < 3L || stats::sd(z) == 0) {
    return(NA_real_)
  }
  if (length(z) > 5000L) {
    z <- withr::with_seed(1L, sample(z, 5000L))
  }
  unname(stats::shapiro.test(z)$statistic)
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
  rows <- list()
  for (nm in unique(scores$.outcome)) {
    sc <- scores[scores$.outcome == nm, , drop = FALSE]
    for (cv in covs) {
      x <- newdata[[cv]][sc$.row]
      ok <- is.finite(sc$z) & !is.na(x)
      z <- sc$z[ok]
      x <- x[ok]
      rows[[length(rows) + 1L]] <- if (is.numeric(x)) {
        tibble::tibble(.outcome = nm, covariate = cv, level = NA_character_,
                       n = length(z), drift_numeric(z, x))
      } else {
        drift_levels(nm, cv, z, as.character(x))
      }
    }
  }
  if (!length(rows)) {
    return(tibble::tibble(
      .outcome = character(), covariate = character(), level = character(),
      n = integer(),
      location_drift = numeric(), location_se = numeric(), location_p = numeric(),
      scale_drift = numeric(), scale_se = numeric(), scale_p = numeric()
    ))
  }
  dplyr_bind(rows)
}

# Numeric covariate: diagnostic GAMs of z and z^2 - 1 on x. The drift is
# the largest absolute fitted value with its pointwise SE (a descriptive
# amplitude; the maximum of a fitted curve is selected, so its SE does
# not give a test), and `p` is the approximate p-value of the smooth term
# (Wood 2013), which is the test of no drift.
drift_numeric <- function(z, x) {
  dat <- data.frame(z = z, x = x)
  loc <- drift_gam(z ~ s(x, k = 5), dat)
  sc2 <- drift_gam(I(z^2 - 1) ~ s(x, k = 5), dat)
  tibble::tibble(
    location_drift = loc$drift, location_se = loc$se, location_p = loc$p,
    scale_drift = sc2$drift, scale_se = sc2$se, scale_p = sc2$p
  )
}

# Factor covariate: one row per level with the mean of z and the excess
# of var(z) over 1, their standard errors, and two-sided normal p-values.
drift_levels <- function(nm, cv, z, lev) {
  parts <- split(z, lev)
  dplyr_bind(lapply(names(parts), function(l) {
    zl <- parts[[l]]
    n <- length(zl)
    if (n < 2L) {
      return(tibble::tibble(
        .outcome = nm, covariate = cv, level = l, n = n,
        location_drift = if (n) mean(zl) else NA_real_, location_se = NA_real_,
        location_p = NA_real_, scale_drift = NA_real_, scale_se = NA_real_,
        scale_p = NA_real_
      ))
    }
    loc <- mean(zl)
    loc_se <- stats::sd(zl) / sqrt(n)
    sq <- zl^2 - 1
    sc <- mean(sq)
    sc_se <- stats::sd(sq) / sqrt(n)
    tibble::tibble(
      .outcome = nm, covariate = cv, level = l, n = n,
      location_drift = loc, location_se = loc_se,
      location_p = 2 * stats::pnorm(-abs(loc / loc_se)),
      scale_drift = sc, scale_se = sc_se,
      scale_p = if (sc_se > 0) 2 * stats::pnorm(-abs(sc / sc_se)) else NA_real_
    )
  }))
}

# Largest absolute fitted value of a diagnostic GAM, its SE there, and
# the p-value of the smooth term; NA with fewer than 20 usable rows or
# when the GAM fails.
drift_gam <- function(formula, dat) {
  m <- if (nrow(dat) < 20L) {
    NULL
  } else {
    tryCatch(mgcv::gam(formula, data = dat), error = function(e) NULL)
  }
  if (is.null(m)) {
    return(list(drift = NA_real_, se = NA_real_, p = NA_real_))
  }
  pr <- stats::predict(m, se.fit = TRUE)
  i <- which.max(abs(pr$fit))
  st <- tryCatch(summary(m)$s.table, error = function(e) NULL)
  p <- if (is.null(st) || !nrow(st)) NA_real_ else unname(st[1L, "p-value"])
  list(drift = abs(pr$fit[[i]]), se = as.numeric(pr$se.fit[[i]]), p = p)
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
print.ref_assessment <- function(x, ...) {
  cli::cli_text("{.cls ref_assessment} n = {x$n}")
  if (isTRUE(x$in_sample)) {
    cli::cli_alert_warning("IN-SAMPLE diagnostic: do not report as held-out performance.")
  }
  print(x$overall)
  invisible(x)
}
