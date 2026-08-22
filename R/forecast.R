#' History-conditioned forecast
#'
#' Returns predictive [norm_dist] objects at each requested time. The
#' dynamic process only changes the predictive distribution; scores and
#' graphics then follow from the existing contract.
#'
#' @param dynamic A [norm_dynamics] object.
#' @param history Subject history data frame.
#' @param times Future times.
#' @param outcome Outcome name.
#' @return A `norm_forecast`.
#' @export
norm_forecast <- function(dynamic, history, times, outcome = NULL) {
  outcome <- outcome %||% dynamic$outcomes[[1]]
  time_nm <- dynamic$time_name
  if (!time_nm %in% names(history)) {
    if ("age" %in% names(history)) {
      time_nm <- "age"
    } else {
      cli::cli_abort("History must contain the time column {.field {dynamic$time_name}}.")
    }
  }
  t_hist <- as.numeric(history[[time_nm]])
  y_hist <- history[[outcome]]
  new_grid <- history[rep(nrow(history), length(times)), , drop = FALSE]
  new_grid[[time_nm]] <- times
  if (outcome %in% names(new_grid)) {
    new_grid[[outcome]] <- NA_real_
  }
  marg <- predict(dynamic$reference, newdata = new_grid, type = "distribution",
                  uncertainty = "conditional")[[outcome]]
  hist_marg <- predict(dynamic$reference, newdata = history, type = "distribution",
                       uncertainty = "conditional")[[outcome]]
  z_hist <- stats::qnorm(clamp_prob(cdf(hist_marg, y_hist)))
  pr <- dynamic$processes[[outcome]]
  dists <- lapply(seq_along(times), function(i) {
    cond <- condition_z(z_hist, t_hist, times[[i]], pr$psi, pr$process)
    condition_norm_dist(vctrs::vec_slice(marg, i), cond$m, cond$s)
  })
  summary <- tibble::tibble(
    time = times,
    median = vapply(dists, function(d) center(d), numeric(1)),
    lower = vapply(dists, function(d) quantile(d, 0.05), numeric(1)),
    upper = vapply(dists, function(d) quantile(d, 0.95), numeric(1))
  )
  structure(
    list(
      dynamic = dynamic,
      outcome = outcome,
      times = times,
      dist = dists,
      summary = summary,
      history_n = length(z_hist)
    ),
    class = "norm_forecast"
  )
}

#' Population-chart derivatives
#'
#' Estimates \eqn{\partial Q_t(p)/\partial t} of the modeled quantile
#' curves. This is a population-chart derivative, not an individual
#' longitudinal estimate.
#'
#' @param reference A [norm_fit].
#' @param newdata Grid of covariate values.
#' @param with_respect_to Name of the time covariate.
#' @param centiles Probability levels.
#' @param type Must be `"chart"`.
#' @param h Finite-difference step.
#' @export
norm_derivative <- function(reference,
                            newdata,
                            with_respect_to,
                            centiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
                            type = c("chart"),
                            h = 1e-2) {
  type <- match.arg(type)
  t_nm <- as.character(with_respect_to)
  outcome <- reference$outcomes[[1]]
  plus <- newdata
  minus <- newdata
  plus[[t_nm]] <- plus[[t_nm]] + h
  minus[[t_nm]] <- minus[[t_nm]] - h
  d_plus <- predict(reference, newdata = plus, type = "distribution",
                    uncertainty = "conditional")[[outcome]]
  d_minus <- predict(reference, newdata = minus, type = "distribution",
                     uncertainty = "conditional")[[outcome]]
  rows <- lapply(centiles, function(p) {
    tibble::tibble(
      time = newdata[[t_nm]],
      centile = p,
      chart_velocity = (dist_quantile(d_plus, p) - dist_quantile(d_minus, p)) / (2 * h)
    )
  })
  structure(dplyr_bind(rows), class = c("norm_derivative", "tbl_df", "tbl", "data.frame"))
}
