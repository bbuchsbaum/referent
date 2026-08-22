#' Assess a fitted reference model
#'
#' Reports overall probabilistic fit, marginal calibration, conditional
#' calibration, and tail calibration on held-out data.
#'
#' @details
#' `overall` maps onto the PCNtoolkit evaluation metrics as follows:
#' `standardized_log_score` is the negative MSLL (mean log score minus
#' the log score of an unconditional Gaussian); `cor` is Rho (Pearson
#' correlation of observed and predicted median); `smse` is the
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
#' @return An object of class `norm_assessment` with `overall`,
#'   `marginal`, `conditional`, and `tail` tibbles plus the scores.
#' @export
norm_assess <- function(fit, newdata, by = NULL) {
  if (missing(newdata) || is.null(newdata)) {
    cli::cli_abort("Supply {.arg newdata}: held-out data or a cross-fitted frame.")
  }
  newdata <- tibble::as_tibble(newdata)
  by_vec <- pull_column(newdata, rlang::enquo(by), default = NULL)
  dists <- predict_dists(fit, newdata, uncertainty = "conditional")
  scores <- scores_from_dists(fit, dists, newdata, allow_extrapolation = FALSE)
  assess_from_scores(fit, scores, dists, newdata, by_vec)
}

assess_from_scores <- function(fit, scores, dists, newdata, by_vec = NULL) {
  overall <- assess_overall(scores, dists, newdata)
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
      in_sample = isTRUE(attr(scores, "in_sample"))
    ),
    class = "norm_assessment"
  )
}

assess_overall <- function(scores, dists, newdata) {
  by_out <- split(scores, scores$.outcome)
  rows <- lapply(names(by_out), function(nm) {
    sc <- by_out[[nm]]
    ok <- is.finite(sc$log_density) & is.finite(sc$observed)
    log_score <- mean(sc$log_density[ok])
    naive <- naive_log_score(sc$observed[ok])
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

# Per-observation CRPS: closed form for a Gaussian predictive, otherwise
# the sample CRPS of scoringRules on seeded draws, or a quantile-grid
# approximation when scoringRules is not installed.
crps_from_dist <- function(dist, y, times = 500L, seed = 1L) {
  y <- rep_len(as.numeric(y), length(dist))
  u <- dist_unpack(dist)
  if (inherits(u, "dist_normal")) {
    return(crps_norm(y, u$mu, u$sigma))
  }
  if (has_pkg("scoringRules")) {
    draws <- withr::with_seed(seed, dist_generate(dist, times))
    ok <- is.finite(y) & rowSums(!is.finite(draws)) == 0L
    out <- rep(NA_real_, length(y))
    out[ok] <- scoringRules::crps_sample(y[ok], draws[ok, , drop = FALSE])
    return(out)
  }
  p <- (seq_len(99L) - 0.5) / 99
  vapply(seq_along(y), function(i) {
    qs <- dist_quantile(dist[i], p)
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
  covs <- intersect(fit$covariates, names(newdata))
  numeric_covs <- covs[vapply(newdata[covs], is.numeric, logical(1))]
  empty <- tibble::tibble(
    .outcome = character(), covariate = character(),
    location_drift = numeric(), location_se = numeric(),
    scale_drift = numeric(), scale_se = numeric()
  )
  if (!length(numeric_covs)) {
    return(empty)
  }
  rows <- list()
  for (nm in unique(scores$.outcome)) {
    sc <- scores[scores$.outcome == nm, , drop = FALSE]
    for (cv in numeric_covs) {
      x <- newdata[[cv]][sc$.row]
      ok <- is.finite(sc$z) & is.finite(x)
      if (sum(ok) < 20L) {
        rows[[length(rows) + 1L]] <- tibble::tibble(
          .outcome = nm, covariate = cv,
          location_drift = NA_real_, location_se = NA_real_,
          scale_drift = NA_real_, scale_se = NA_real_
        )
        next
      }
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
  dplyr_bind(rows)
}

# Largest absolute fitted value of a diagnostic GAM and its SE there.
drift_gam <- function(formula, dat) {
  m <- tryCatch(mgcv::gam(formula, data = dat), error = function(e) NULL)
  if (is.null(m)) {
    return(list(drift = NA_real_, se = NA_real_))
  }
  pr <- stats::predict(m, se.fit = TRUE)
  i <- which.max(abs(pr$fit))
  list(drift = abs(pr$fit[[i]]), se = as.numeric(pr$se.fit[[i]]))
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
