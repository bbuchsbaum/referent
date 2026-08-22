#' Held-out probability recalibration
#'
#' Maps PIT values through a monotone map estimated on a calibration
#' sample: \eqn{F^\star(y\mid x)=G(F_0(y\mid x))}. This targets global
#' (or group-specific) calibration. It does not create full conditional
#' calibration, so pre-calibration diagnostics (the assessment tables,
#' without the per-row scores) are kept in `fit$calibration$pre`.
#'
#' @details
#' The map is a sinh-arcsinh reference distribution for the normal score
#' \eqn{z_0=\Phi^{-1}(F_0(y\mid x))}, that is
#' \eqn{G(u)=P_{\mathrm{SHASH}}(\Phi^{-1}(u))}. Its four parameters move
#' the location, spread, skew, and tail weight of \eqn{z}, and it nests
#' the identity (\eqn{\mu=0,\sigma=1,\epsilon=0,\delta=1}), so a model
#' that is already calibrated is left alone. It is smooth and strictly
#' monotone on the whole line, so extreme observations keep extreme,
#' finite, correctly ordered scores.
#'
#' Because \eqn{F^\star=G(F_0)}, the whole predictive moves, not just the
#' centile: the density is \eqn{f^\star=g(F_0(y))f_0(y)}, so `log_density`
#' and every density-based score change, and the calibrated median is
#' \eqn{Q_0(G^{-1}(1/2))}, so `median` and `residual` change as well.
#'
#' The map is kept only where it earns its place: the fitted parameters are
#' used when their log-likelihood beats the identity's by more than the four
#' parameters they cost, which is the AIC comparison of the two, and the map
#' falls back to the exact identity when they do not. On calibration scores
#' that are already uniform the fitted map survives that comparison for about
#' 5 percent of samples, and is close to the identity when it does. A map is
#' absent altogether, leaving `calibrated = FALSE`, only when there are fewer
#' than twenty usable scores or the fit does not converge.
#'
#' @param fit A [ref_fit].
#' @param data Calibration reference data, not reused for evaluation, or
#'   an out-of-fold `ref_scores` table from [ref_crossfit()]. Passing the
#'   cross-fitted scores calibrates on the reference sample itself
#'   without holding any of it out.
#' @param by Optional grouping column (e.g. site). Groups absent from the
#'   calibration data fall back to the pooled map. Not available when
#'   `data` is a score table.
#' @param uncertainty Which predictive the calibration data are scored
#'   under, as in [ref_assess()]. A map corrects the predictive it was
#'   fitted on, so this should match the `uncertainty` that
#'   [predict.ref_fit()] will later be asked for. Ignored when `data` is a
#'   score table, which already carries its own.
#' @return The fit with a `calibration` slot; [predict.ref_fit()] applies
#'   the map to `centile`, `z`, `log_density`, `median`, `residual`, and
#'   the tail columns and sets `calibrated = TRUE`.
#' @export
ref_calibrate <- function(fit, data, by = NULL,
                          uncertainty = c("conditional", "total")) {
  uncertainty <- match.arg(uncertainty)
  by_quo <- rlang::enquo(by)
  base <- fit
  base$calibration <- NULL
  if (inherits(data, "ref_scores")) {
    if (!rlang::quo_is_null(by_quo)) {
      cli::cli_abort("{.arg by} needs calibration {.arg data}, not a score table.")
    }
    scores <- data
    by_vec <- NULL
    pre <- NULL
    n <- length(unique(scores$.row))
  } else {
    data <- tibble::as_tibble(data)
    by_vec <- pull_column(data, by_quo, default = NULL)
    dists <- predict_dists(base, data, uncertainty = uncertainty)
    scores <- scores_from_dists(base, dists, data)
    pre <- assess_from_scores(base, scores, dists, data, by_vec = NULL)
    pre$scores <- NULL
    n <- nrow(data)
  }
  maps <- lapply(split(scores, scores$.outcome), function(sc) {
    out <- list(.global = pit_map(sc$z))
    if (!is.null(by_vec)) {
      groups <- split(sc$z, as.character(by_vec)[sc$.row])
      out <- c(out, lapply(groups, pit_map))
    }
    out
  })
  fit$calibration <- structure(
    list(
      by = if (is.null(by_vec)) NULL else as_col_name(by_quo),
      maps = maps,
      n = n,
      pre = pre
    ),
    class = "ref_calibration"
  )
  fit
}

# Sinh-arcsinh reference distribution for the calibration normal scores,
# fitted by maximum likelihood and kept only when its log-likelihood beats
# the identity's by more than its four parameters (an AIC comparison).
pit_map <- function(z, penalty = 4) {
  z <- z[is.finite(z)]
  n <- length(z)
  if (n < 5L * penalty) {
    return(NULL)
  }
  par <- shash_map_ml(z)
  if (is.null(par)) {
    return(NULL)
  }
  gain <- sum(shash_log_density(z, par[[1L]], par[[2L]], par[[3L]], par[[4L]]) -
    stats::dnorm(z, log = TRUE))
  if (!is.finite(gain) || gain <= penalty) {
    par <- c(0, 1, 0, 1) # the identity: the sample is already calibrated
    gain <- 0
  }
  structure(
    list(mu = par[[1L]], sigma = par[[2L]], eps = par[[3L]], delta = par[[4L]],
         n = n, log_lik_gain = gain),
    class = "ref_pit_map"
  )
}

shash_map_ml <- function(z) {
  nll <- function(p) {
    v <- shash_log_density(z, p[[1L]], exp(p[[2L]]), p[[3L]], exp(p[[4L]]))
    if (!all(is.finite(v))) {
      return(.Machine$double.xmax)
    }
    -sum(v)
  }
  start <- c(mean(z), log(stats::sd(z)), 0, 0)
  o <- stats::optim(start, nll, control = list(maxit = 2000L, reltol = 1e-10))
  o <- stats::optim(o$par, nll, control = list(maxit = 2000L, reltol = 1e-10))
  if (!is.finite(o$value) || o$value >= .Machine$double.xmax) {
    return(NULL)
  }
  c(o$par[[1L]], exp(o$par[[2L]]), o$par[[3L]], exp(o$par[[4L]]))
}

# Calibrated log tails and the log derivative of the map at the observed
# PIT. Everything is evaluated in log space, so very large |z| stay
# distinguishable and ordered.
apply_pit_map <- function(map, z) {
  if (is.null(map)) {
    return(list(
      lower = stats::pnorm(z, log.p = TRUE),
      upper = stats::pnorm(z, log.p = TRUE, lower.tail = FALSE),
      log_jacobian = rep(0, length(z))
    ))
  }
  list(
    lower = shash_log_cdf(z, map$mu, map$sigma, map$eps, map$delta, TRUE),
    upper = shash_log_cdf(z, map$mu, map$sigma, map$eps, map$delta, FALSE),
    log_jacobian = shash_log_density(z, map$mu, map$sigma, map$eps, map$delta) -
      stats::dnorm(z, log = TRUE)
  )
}

# The map's own median on the PIT scale: G^-1(1/2), so that the calibrated
# median is the raw predictive quantile there. The identity gives 1/2.
pit_map_median <- function(map) {
  if (is.null(map)) {
    return(0.5)
  }
  stats::pnorm(shash_quantile(0.5, map$mu, map$sigma, map$eps, map$delta))
}

# Calibration group of every row of `newdata` (".global" when the map is
# pooled or the group is missing), or NULL when the fit is not calibrated.
calibration_groups <- function(cal, newdata) {
  if (is.null(cal)) {
    return(NULL)
  }
  by_nm <- cal$by
  grp <- if (!is.null(by_nm) && by_nm %in% names(newdata)) {
    as.character(newdata[[by_nm]])
  } else {
    rep(".global", nrow(newdata))
  }
  grp[is.na(grp)] <- ".global"
  grp
}

# Map one outcome's score table through its PIT maps, one vectorised pass
# per group. Groups without a map use the pooled map, and rows that no map
# reached keep `calibrated = FALSE`. F* = G(F0) implies f* = g(F0(y)) f0(y)
# and median* = Q0(G^-1(1/2)), so the density and the median move with the
# centile; `d` is the raw predictive the score table came from.
calibrate_scores <- function(cal, outcome, grp, sc, d) {
  maps <- cal$maps[[outcome]]
  sc$calibrated <- rep(FALSE, nrow(sc))
  if (is.null(maps)) {
    return(sc)
  }
  lower <- stats::pnorm(sc$z, log.p = TRUE)
  upper <- stats::pnorm(sc$z, log.p = TRUE, lower.tail = FALSE)
  jac <- rep(0, nrow(sc))
  for (g in unique(grp)) {
    idx <- which(grp == g)
    map <- maps[[g]] %||% maps$.global
    if (is.null(map)) {
      next
    }
    mapped <- apply_pit_map(map, sc$z[idx])
    lower[idx] <- mapped$lower
    upper[idx] <- mapped$upper
    jac[idx] <- mapped$log_jacobian
    p_med <- pit_map_median(map)
    if (p_med != 0.5) {
      sc$median[idx] <- dist_eval(d[idx], quantile, rep(p_med, length(idx)))
      sc$residual[idx] <- sc$observed[idx] - sc$median[idx]
    }
    sc$calibrated[idx] <- TRUE
  }
  out <- scores_from_log_tails(lower, upper)
  sc$centile <- out$centile
  sc$z <- out$z
  sc$tail_prob <- out$tail_prob
  sc$tail_surprisal <- out$tail_surprisal
  ok <- is.finite(jac)
  sc$log_density[ok] <- sc$log_density[ok] + jac[ok]
  sc
}

#' @export
print.ref_calibration <- function(x, ...) {
  cli::cli_text("{.cls ref_calibration} n = {x$n}{if (is.null(x$by)) '' else paste0(', by ', x$by)}")
  invisible(x)
}
