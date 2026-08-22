#' Held-out probability recalibration
#'
#' Maps PIT values through a monotone map estimated on a calibration
#' sample: \eqn{F^\star(y\mid x)=G(F_0(y\mid x))}. This targets global
#' (or group-specific) calibration. It does not create full conditional
#' calibration, so pre-calibration diagnostics are kept in
#' `fit$calibration$pre`.
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
#' @param method `"rank"` (rank-interpolated empirical CDF).
#' @param by Optional grouping column (e.g. site). Groups absent from the
#'   calibration data fall back to the pooled map.
#' @return The fit with a `calibration` slot; [predict.norm_fit()] applies
#'   the map to `centile`, `z`, and the tail columns and sets
#'   `calibrated = TRUE`.
#' @export
norm_calibrate <- function(fit, data, method = c("rank"), by = NULL) {
  method <- match.arg(method)
  data <- tibble::as_tibble(data)
  by_quo <- rlang::enquo(by)
  by_vec <- pull_column(data, by_quo, default = NULL)
  base <- fit
  base$calibration <- NULL
  scores <- predict(base, newdata = data, type = "scores",
                    uncertainty = "conditional", allow_extrapolation = TRUE)
  pre <- norm_assess(base, newdata = data)
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
      method = method,
      by = if (is.null(by_vec)) NULL else rlang::as_name(by_quo),
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

apply_calibration <- function(fit, scores, newdata = NULL) {
  cal <- fit$calibration
  if (is.null(cal) || !nrow(scores)) {
    return(scores)
  }
  by_nm <- cal$by
  grp <- if (!is.null(by_nm) && !is.null(newdata) && by_nm %in% names(newdata)) {
    as.character(newdata[[by_nm]])[scores$.row]
  } else {
    rep(".global", nrow(scores))
  }
  grp[is.na(grp)] <- ".global"
  tails <- log_tails_from_z(scores$z)
  lower <- tails$lower
  upper <- tails$upper
  key <- paste(scores$.outcome, grp)
  for (k in unique(key)) {
    idx <- which(key == k)
    nm <- scores$.outcome[[idx[[1L]]]]
    maps <- cal$maps[[nm]]
    if (is.null(maps)) {
      next
    }
    map <- maps[[grp[[idx[[1L]]]]]] %||% maps$.global
    mapped <- apply_pit_map(map, lower[idx], upper[idx])
    lower[idx] <- mapped$lower
    upper[idx] <- mapped$upper
  }
  out <- scores_from_log_tails(lower, upper)
  scores$centile <- out$centile
  scores$z <- out$z
  scores$tail_prob <- out$tail_prob
  scores$tail_surprisal <- out$tail_surprisal
  scores$calibrated <- TRUE
  scores
}

#' @export
print.norm_calibration <- function(x, ...) {
  cli::cli_text("{.cls norm_calibration} method = {x$method}, n = {x$n}{if (is.null(x$by)) '' else paste0(', by ', x$by)}")
  invisible(x)
}
