#' History-conditioned forecast
#'
#' Returns a [dist_conditioned] predictive distribution at each requested time. The
#' dynamic process only changes the predictive distribution; scores and
#' graphics then follow from the existing contract.
#'
#' The forecast refuses to run when the process for `outcome` is not
#' identified (see [ref_dynamics()]): a history-conditioned distribution
#' from an unestimated kernel would be silently wrong. History rows with a
#' non-finite time or outcome are dropped and the number actually used is
#' returned as `history_n`. Epistemic uncertainty of the reference location
#' is carried in the marginal distributions but the uncertainty of the
#' kernel parameters is not propagated into the forecast.
#'
#' @param dynamic A [ref_dynamics] object.
#' @param history Subject history data frame. Must contain the time column
#'   used in [ref_dynamics()] and the outcome.
#' @param times Future times.
#' @param outcome Outcome name.
#' @return A `ref_forecast` with `$dist` (a `distribution` vector, one element per time),
#'   `$summary` (median and 90% interval), `$history_n`, and
#'   `$lag_support` (`"in"` or `"extrapolated_lag"` per forecast time,
#'   judged from the last history visit against the reference lag range).
#' @export
ref_forecast <- function(dynamic, history, times, outcome = NULL) {
  outcome <- outcome %||% dynamic$outcomes[[1]]
  if (!outcome %in% dynamic$outcomes) {
    cli::cli_abort("{.val {outcome}} is not an outcome of the dynamics object.")
  }
  pr <- dynamic$processes[[outcome]]
  if (!isTRUE(pr$identified)) {
    cli::cli_abort(c(
      "The longitudinal process for {.field {outcome}} is not identified; no history-conditioned forecast is possible.",
      i = pr$reason %||% "See print(dynamic)."
    ))
  }
  time_nm <- dynamic$time_name
  history <- tibble::as_tibble(history)
  if (!time_nm %in% names(history)) {
    cli::cli_abort("History must contain the time column {.field {time_nm}} used by {.fn ref_dynamics}.")
  }
  if (!outcome %in% names(history)) {
    cli::cli_abort("History must contain the outcome column {.field {outcome}}.")
  }
  times <- as.numeric(times)
  keep <- is.finite(as.numeric(history[[time_nm]])) & is.finite(as.numeric(history[[outcome]]))
  n_dropped <- sum(!keep)
  history <- history[keep, , drop = FALSE]
  if (!nrow(history)) {
    cli::cli_abort("No history rows with finite {.field {time_nm}} and {.field {outcome}}.")
  }
  if (n_dropped > 0L) {
    cli::cli_inform("Dropped {n_dropped} history row{?s} with non-finite time or outcome.")
  }
  history <- history[order(history[[time_nm]]), , drop = FALSE]
  t_hist <- as.numeric(history[[time_nm]])
  y_hist <- as.numeric(history[[outcome]])
  new_grid <- history[rep(nrow(history), length(times)), , drop = FALSE]
  new_grid[[time_nm]] <- times
  new_grid[[outcome]] <- NA_real_
  marg <- predict_dists(dynamic$reference, new_grid, uncertainty = "conditional")[[outcome]]
  hist_marg <- predict_dists(dynamic$reference, history, uncertainty = "conditional")[[outcome]]
  z_hist <- dist_z(hist_marg, y_hist)
  cond <- condition_history(z_hist, t_hist, times, pr)
  dists <- dist_conditioned(marg, m = cond$m, s = cond$s)
  lag_support <- classify_temporal_support(dynamic, rep(max(t_hist), length(times)),
                                           times, rep(length(t_hist), length(times)))
  summary <- tibble::tibble(
    time = times,
    median = dist_quantile(dists, 0.5),
    lower = dist_quantile(dists, 0.05),
    upper = dist_quantile(dists, 0.95),
    support = lag_support
  )
  structure(
    list(
      dynamic = dynamic,
      outcome = outcome,
      times = times,
      dist = dists,
      summary = summary,
      history = history,
      time_name = time_nm,
      history_n = length(z_hist),
      lag_support = lag_support
    ),
    class = "ref_forecast"
  )
}

#' Population-chart derivatives
#'
#' Estimates \eqn{\partial Q_t(p)/\partial t} of the modeled quantile
#' curves by central finite differences, with a delta-method standard
#' error from the fitted coefficient covariance (`Vp`) applied to the
#' finite-difference of the linear-predictor matrix. This is a
#' population-chart derivative, not an individual longitudinal estimate.
#'
#' @param reference A [ref_fit].
#' @param newdata Grid of covariate values.
#' @param with_respect_to The time covariate, as a bare column name or a
#'   string.
#' @param centiles Probability levels.
#' @param h Finite-difference step.
#' @param outcomes Optional subset of outcomes (default: all).
#' @return A `ref_derivative` tibble with columns `.outcome`, `time`,
#'   `centile`, `chart_velocity`, and `chart_velocity_se` (`NA` when the
#'   engine does not expose a coefficient covariance).
#' @examples
#' ref <- ref_simulate(150, seed = 1)
#' fit <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 5) + sex), ref, "y")
#' grid <- data.frame(age = c(30, 50, 70), sex = factor("F", levels = c("F", "M")))
#' ref_derivative(fit, grid, with_respect_to = age, centiles = c(0.1, 0.5, 0.9))
#' @export
ref_derivative <- function(reference,
                            newdata,
                            with_respect_to,
                            centiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
                            h = 1e-2,
                            outcomes = NULL) {
  t_nm <- as_col_name(rlang::enquo(with_respect_to))
  if (is.null(t_nm)) {
    cli::cli_abort("{.arg with_respect_to} must name the time covariate.")
  }
  newdata <- tibble::as_tibble(newdata)
  if (!t_nm %in% names(newdata)) {
    cli::cli_abort("{.arg newdata} must contain {.field {t_nm}}.")
  }
  nms <- outcomes %||% reference$outcomes
  nms <- intersect(as.character(nms), reference$outcomes)
  plus <- newdata
  minus <- newdata
  plus[[t_nm]] <- plus[[t_nm]] + h
  minus[[t_nm]] <- minus[[t_nm]] - h
  rows <- lapply(nms, function(nm) {
    fit_one <- reference$models[[nm]]
    if (!fit_ok(fit_one)) {
      return(NULL)
    }
    d_plus <- predict_mgcv_dist(fit_one, plus)
    d_minus <- predict_mgcv_dist(fit_one, minus)
    vel <- vapply(centiles, function(p) {
      (dist_quantile(d_plus, p) - dist_quantile(d_minus, p)) / (2 * h)
    }, numeric(nrow(newdata)))
    vel <- matrix(vel, nrow(newdata), length(centiles))
    se <- derivative_se(fit_one, plus, minus, centiles, h)
    dplyr_bind(lapply(seq_along(centiles), function(k) {
      tibble::tibble(
        .outcome = nm,
        time = newdata[[t_nm]],
        centile = centiles[[k]],
        chart_velocity = vel[, k],
        chart_velocity_se = se[, k]
      )
    }))
  })
  rows <- Filter(Negate(is.null), rows)
  out <- if (length(rows)) {
    dplyr_bind(rows)
  } else {
    tibble::tibble(.outcome = character(), time = numeric(), centile = numeric(),
                   chart_velocity = numeric(), chart_velocity_se = numeric())
  }
  structure(out, class = c("ref_derivative", class(out)))
}

# Delta-method SE of the chart velocity: the velocity is a function of the
# coefficient vector through the lp matrices at t +/- h; its gradient is
# taken by central differences in beta and combined with Vp.
derivative_se <- function(fit_one, plus, minus, centiles, h) {
  n <- nrow(plus)
  na <- matrix(NA_real_, n, length(centiles))
  model <- fit_one$model
  fam <- fit_one$family$name
  if (is.null(model$Vp)) {
    return(na)
  }
  lp_p <- tryCatch(predict_gam_quiet(model, plus, type = "lpmatrix"), error = function(e) NULL)
  lp_m <- tryCatch(predict_gam_quiet(model, minus, type = "lpmatrix"), error = function(e) NULL)
  beta <- stats::coef(model)
  vp <- as.matrix(model$Vp)
  if (is.null(lp_p) || is.null(lp_m) || !all(is.finite(beta)) || !all(is.finite(vp)) ||
      ncol(vp) != length(beta)) {
    return(na)
  }
  cond_p <- mgcv_parameters(model, plus, fam, fit_one)
  cond_m <- mgcv_parameters(model, minus, fam, fit_one)
  q_at <- function(lp, b, cond) {
    par <- params_from_eta(eta_from_lp(lp, b), model, fam, fit_one, cond)
    d <- make_dist(fam, par)
    vapply(centiles, function(p) dist_quantile(d, p), numeric(n))
  }
  vel_at <- function(b) (q_at(lp_p, b, cond_p) - q_at(lp_m, b, cond_m)) / (2 * h)
  step <- pmax(1e-4 * abs(beta), 1e-6)
  grads <- lapply(seq_along(beta), function(j) {
    e <- replace(numeric(length(beta)), j, step[[j]])
    (vel_at(beta + e) - vel_at(beta - e)) / (2 * step[[j]])
  })
  out <- na
  for (k in seq_along(centiles)) {
    g <- do.call(cbind, lapply(grads, function(gm) matrix(gm, n, length(centiles))[, k]))
    out[, k] <- sqrt(pmax(rowSums((g %*% vp) * g), 0))
  }
  out
}
