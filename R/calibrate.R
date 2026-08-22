#' Held-out probability recalibration
#'
#' Maps PIT values through a monotone map estimated on a calibration
#' sample: \eqn{F^\star(y\mid x)=G(F_0(y\mid x))}. This targets global
#' (or group-specific) calibration. It does not create full conditional
#' calibration, so pre-calibration diagnostics (the assessment tables,
#' without the per-row scores) are kept in `fit$calibration$pre`.
#'
#' @details
#' With calibration PITs \eqn{u_{(1)} \le \dots \le u_{(n)}}, the map
#' interpolates linearly between the knots
#' \eqn{(u_{(r)}, r/(n+1))} and continues linearly to \eqn{(0,0)} and
#' \eqn{(1,1)} outside the observed range. Beyond the range the map is
#' therefore a rescaled identity: extreme tails are compressed or
#' stretched but never collapsed to 0 or 1, and the calibrated tail is
#' evaluated in log space so very large `|z|` remain distinguishable.
#'
#' @param fit A [norm_fit].
#' @param data Calibration reference data, not reused for evaluation.
#' @param by Optional grouping column (e.g. site). Groups absent from the
#'   calibration data fall back to the pooled map.
#' @return The fit with a `calibration` slot; [predict.norm_fit()] applies
#'   the map to `centile`, `z`, and the tail columns and sets
#'   `calibrated = TRUE`.
#' @export
norm_calibrate <- function(fit, data, by = NULL) {
  data <- tibble::as_tibble(data)
  by_quo <- rlang::enquo(by)
  by_vec <- pull_column(data, by_quo, default = NULL)
  base <- fit
  base$calibration <- NULL
  dists <- predict_dists(base, data, uncertainty = "conditional")
  scores <- scores_from_dists(base, dists, data)
  pre <- assess_from_scores(base, scores, dists, data, by_vec = NULL)
  pre$scores <- NULL
  maps <- lapply(split(scores, scores$.outcome), function(sc) {
    out <- list(.global = pit_map(sc$centile))
    if (!is.null(by_vec)) {
      groups <- split(sc$centile, as.character(by_vec)[sc$.row])
      out <- c(out, lapply(groups, pit_map))
    }
    out
  })
  fit$calibration <- structure(
    list(
      method = "rank",
      by = if (is.null(by_vec)) NULL else as_col_name(by_quo),
      maps = maps,
      n = nrow(data),
      pre = pre
    ),
    class = "norm_calibration"
  )
  fit
}

# Rank-interpolated PIT map with linear (rescaled identity) tails.
pit_map <- function(u) {
  u <- sort(u[is.finite(u) & u > 0 & u < 1])
  n <- length(u)
  if (n < 2L) {
    return(NULL)
  }
  structure(
    list(x = c(0, u, 1), y = c(0, seq_len(n) / (n + 1), 1), n = n),
    class = "norm_pit_map"
  )
}

# Log lower and upper calibrated tails from raw log tails. Inside the
# knot range the map is interpolated on the probability scale; beyond it
# the tail is scaled by the slope of the boundary segment in log space.
apply_pit_map <- function(map, log_lower, log_upper) {
  if (is.null(map)) {
    return(list(lower = log_lower, upper = log_upper))
  }
  x <- map$x
  y <- map$y
  n <- map$n
  u <- exp(log_lower)
  g <- stats::approx(x, y, xout = u, ties = "ordered", rule = 2)$y
  lower <- log(g)
  upper <- log1p(-g)
  lo_slope <- log(y[2]) - log(x[2])
  hi_slope <- log1p(-y[n + 1]) - log1p(-x[n + 1])
  below <- is.finite(u) & u < x[2]
  above <- is.finite(u) & u > x[n + 1]
  lower[below] <- log_lower[below] + lo_slope
  upper[below] <- log1p(-exp(lower[below]))
  upper[above] <- log_upper[above] + hi_slope
  lower[above] <- log1p(-exp(upper[above]))
  list(lower = lower, upper = upper)
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

# Map the tails of one outcome's score table through its PIT maps, one
# vectorised pass per group. Groups without a map use the pooled map.
calibrate_scores <- function(cal, outcome, grp, sc) {
  maps <- cal$maps[[outcome]]
  if (is.null(maps)) {
    return(sc)
  }
  tails <- log_tails_from_z(sc$z)
  lower <- tails$lower
  upper <- tails$upper
  for (g in unique(grp)) {
    idx <- which(grp == g)
    mapped <- apply_pit_map(maps[[g]] %||% maps$.global, lower[idx], upper[idx])
    lower[idx] <- mapped$lower
    upper[idx] <- mapped$upper
  }
  out <- scores_from_log_tails(lower, upper)
  sc$centile <- out$centile
  sc$z <- out$z
  sc$tail_prob <- out$tail_prob
  sc$tail_surprisal <- out$tail_surprisal
  sc
}

#' @export
print.norm_calibration <- function(x, ...) {
  cli::cli_text("{.cls norm_calibration} method = {x$method}, n = {x$n}{if (is.null(x$by)) '' else paste0(', by ', x$by)}")
  invisible(x)
}
