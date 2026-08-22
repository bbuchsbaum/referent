#' Assess a fitted reference model
#'
#' Reports overall probabilistic fit, marginal calibration, conditional
#' calibration, and tail calibration. Phase 1 returns numbers, not plots.
#'
#' @param fit A [norm_fit].
#' @param newdata Held-out validation data. If omitted, training data
#'   scores are used and marked in-sample.
#' @param by Optional grouping column for conditional summaries.
#' @return An object of class `norm_assessment`.
#' @export
norm_assess <- function(fit, newdata = NULL, by = NULL) {
  if (is.null(newdata)) {
    cli::cli_warn("Assessing training data; prefer held-out or cross-fitted scores.")
    newdata <- recover_newdata(fit)
  }
  scores <- predict(fit, newdata = newdata, type = "scores", uncertainty = "conditional")
  dists <- predict(fit, newdata = newdata, type = "distribution", uncertainty = "conditional")
  overall <- assess_overall(scores, dists, newdata)
  marginal <- assess_marginal(scores)
  conditional <- assess_conditional(scores, newdata, fit)
  tail <- assess_tail(scores)
  structure(
    list(
      overall = overall,
      marginal = marginal,
      conditional = conditional,
      tail = tail,
      scores = scores,
      n = nrow(newdata),
      in_sample = isTRUE(attr(scores, "in_sample"))
    ),
    class = "norm_assessment"
  )
}

recover_newdata <- function(fit) {
  cli::cli_abort(
    "Supply {.arg newdata} for assessment. Training frames are not stored on the fit by default."
  )
}

assess_overall <- function(scores, dists, newdata) {
  by_out <- split(scores, scores$.outcome)
  rows <- lapply(names(by_out), function(nm) {
    sc <- by_out[[nm]]
    ok <- is.finite(sc$log_density) & is.finite(sc$observed)
    log_score <- mean(sc$log_density[ok])
    naive <- naive_log_score(sc$observed[ok])
    crps <- mean_crps(dists[[nm]], sc$observed)
    tibble::tibble(
      .outcome = nm,
      mean_log_score = log_score,
      standardized_log_score = log_score - naive,
      crps = crps,
      mae = mean(abs(sc$residual[ok])),
      rmse = sqrt(mean(sc$residual[ok]^2)),
      cor = if (sum(ok) > 2) stats::cor(sc$observed[ok], sc$median[ok]) else NA_real_
    )
  })
  dplyr_bind(rows)
}

naive_log_score <- function(y) {
  mu <- mean(y, na.rm = TRUE)
  s <- stats::sd(y, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) {
    s <- 1
  }
  mean(stats::dnorm(y, mu, s, log = TRUE), na.rm = TRUE)
}

mean_crps <- function(dist, y) {
  if (is.null(dist)) {
    return(NA_real_)
  }
  y <- recycle_to(y, length(dist))
  if (has_pkg("scoringRules") && identical(attr(dist, "family"), "gaussian")) {
    p <- norm_params(dist)
    return(mean(scoringRules::crps_norm(y, mean = p$location, sd = p$scale), na.rm = TRUE))
  }
  mean(crps_from_dist(dist, y), na.rm = TRUE)
}

crps_from_dist <- function(dist, y) {
  if (identical(attr(dist, "family"), "gaussian")) {
    p <- norm_params(dist)
    return(crps_norm(y, p$location, p$scale))
  }
  u <- (seq_len(99L) - 0.5) / 99
  vapply(seq_along(y), function(i) {
    di <- vctrs::vec_slice(dist, i)
    qs <- dist_quantile(di, u)
    mean(abs(qs - y[[i]])) - 0.5 * mean(abs(outer(qs, qs, `-`)))
  }, numeric(1))
}

crps_norm <- function(y, mu, sigma) {
  z <- (y - mu) / sigma
  sigma * (z * (2 * stats::pnorm(z) - 1) + 2 * stats::dnorm(z) - 1 / sqrt(pi))
}

assess_marginal <- function(scores) {
  by_out <- split(scores, scores$.outcome)
  rows <- lapply(names(by_out), function(nm) {
    z <- by_out[[nm]]$z
    z <- z[is.finite(z)]
    u <- by_out[[nm]]$centile
    u <- u[is.finite(u)]
    tibble::tibble(
      .outcome = nm,
      n = length(z),
      mean_z = mean(z),
      var_z = stats::var(z),
      skew_z = if (length(z) > 3) skewness(z) else NA_real_,
      excess_kurtosis = if (length(z) > 4) excess_kurtosis(z) else NA_real_,
      cover_50 = mean(u > 0.25 & u < 0.75, na.rm = TRUE),
      cover_80 = mean(u > 0.10 & u < 0.90, na.rm = TRUE),
      cover_90 = mean(u > 0.05 & u < 0.95, na.rm = TRUE),
      cover_95 = mean(u > 0.025 & u < 0.975, na.rm = TRUE),
      cover_99 = mean(u > 0.005 & u < 0.995, na.rm = TRUE)
    )
  })
  dplyr_bind(rows)
}

skewness <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) {
    return(0)
  }
  mean(((x - m) / s)^3)
}

excess_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) {
    return(0)
  }
  mean(((x - m) / s)^4) - 3
}

assess_conditional <- function(scores, newdata, fit) {
  covs <- intersect(fit$covariate_names, names(newdata))
  numeric_covs <- covs[vapply(newdata[covs], is.numeric, logical(1))]
  if (!length(numeric_covs) || !requireNamespace("mgcv", quietly = TRUE)) {
    return(tibble::tibble(
      .outcome = unique(scores$.outcome),
      location_drift = 0,
      scale_drift = 0
    ))
  }
  rows <- lapply(unique(scores$.outcome), function(nm) {
    sc <- scores[scores$.outcome == nm, , drop = FALSE]
    z <- sc$z
    ok <- is.finite(z)
    if (sum(ok) < 20L) {
      return(tibble::tibble(.outcome = nm, location_drift = NA_real_, scale_drift = NA_real_))
    }
    dat <- cbind(z = z, newdata[sc$.row, numeric_covs[1], drop = FALSE])
    names(dat)[2] <- "x"
    dat <- dat[ok, , drop = FALSE]
    loc <- tryCatch(
      mgcv::gam(z ~ s(x, k = 5), data = dat),
      error = function(e) NULL
    )
    sc2 <- tryCatch(
      mgcv::gam(I(z^2 - 1) ~ s(x, k = 5), data = dat),
      error = function(e) NULL
    )
    tibble::tibble(
      .outcome = nm,
      location_drift = if (is.null(loc)) NA_real_ else max(abs(stats::fitted(loc))),
      scale_drift = if (is.null(sc2)) NA_real_ else max(abs(stats::fitted(sc2)))
    )
  })
  dplyr_bind(rows)
}

assess_tail <- function(scores) {
  levels <- c(0.005, 0.01, 0.025, 0.05)
  by_out <- split(scores, scores$.outcome)
  rows <- lapply(names(by_out), function(nm) {
    tp <- by_out[[nm]]$tail_prob
    tp <- tp[is.finite(tp)]
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
  dplyr_bind(rows)
}

#' @export
print.norm_assessment <- function(x, ...) {
  cli::cli_text("{.cls norm_assessment} n = {x$n}")
  print(x$overall)
  invisible(x)
}

acceptable_calibration <- function(assessment, var_lo = 0.5, var_hi = 1.8,
                                   mean_z_max = 0.35) {
  m <- assessment$marginal
  ok <- is.finite(m$mean_z) & abs(m$mean_z) <= mean_z_max &
    is.finite(m$var_z) & m$var_z >= var_lo & m$var_z <= var_hi
  all(ok)
}
