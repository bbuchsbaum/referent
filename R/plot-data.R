#' Tidy plot data for reference charts
#'
#' These helpers return tibbles. `autoplot()` consumes them; users can
#' also rebuild a chart by hand.
#'
#' @name fortify_centiles
#' @param fit A [ref_fit].
#' @param outcome Outcome name.
#' @param x Numeric covariate for the x-axis.
#' @param by Optional grouping factor (column name). One chart is built
#'   per level; other covariates stay at their reference value.
#' @param centiles Probability levels to evaluate.
#' @param n Grid length along `x`.
#' @return A list of class `ref_centile_data` with `lines` and `ribbons` tibbles.
NULL

#' @rdname fortify_centiles
#' @export
fortify_centiles <- function(fit,
                             outcome = NULL,
                             x = NULL,
                             by = NULL,
                             centiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
                             n = 120) {
  outcome <- outcome %||% fit$outcomes[[1]]
  x_nm <- x %||% default_x(fit)
  by_nm <- by
  centiles <- sort(unique(as.numeric(centiles)))
  if (any(!is.finite(centiles) | centiles <= 0 | centiles >= 1)) {
    cli::cli_abort("{.arg centiles} must lie in (0, 1).")
  }
  grid <- expand_centile_grid(fit, x_nm, by_nm, n)
  dists <- predict(fit, newdata = grid, type = "distribution",
                   uncertainty = "conditional")[[outcome]]
  lines <- dplyr_bind(lapply(centiles, function(p) {
    tibble::tibble(
      x = grid[[x_nm]],
      y = as.numeric(dist_quantile(dists, rep(p, length(dists)))),
      centile = p,
      .label = centile_label(p),
      .group = grid$.group,
      .outcome = outcome
    )
  }))
  ribbons <- ribbon_from_lines(lines)
  structure(
    list(lines = lines, ribbons = ribbons, x_name = x_nm, by = by_nm, outcome = outcome),
    class = "ref_centile_data"
  )
}

expand_centile_grid <- function(fit, x_nm, by_nm, n) {
  if (is.null(by_nm)) {
    grid <- centile_grid(fit, x_nm, n = n)
    grid$.group <- NA_character_
    return(grid)
  }
  levs <- fit$support_ref$factor_levels[[by_nm]]
  if (is.null(levs)) {
    cli::cli_abort("{.arg by} must be a factor covariate on the fit, not {.val {by_nm}}.")
  }
  dplyr_bind(lapply(levs, function(lv) {
    extras <- list()
    extras[[by_nm]] <- factor(lv, levels = levs)
    g <- centile_grid(fit, x_nm, n = n, extras = extras)
    g$.group <- as.character(lv)
    g
  }))
}

centile_grid <- function(fit, x_nm, n = 120, extras = NULL) {
  r <- fit$support_ref$numeric[[x_nm]]
  if (is.null(r)) {
    cli::cli_abort("No numeric covariate {.field {x_nm}} for a centile plot.")
  }
  grid <- tibble::tibble(!!x_nm := seq(r$min, r$max, length.out = n))
  for (nm in setdiff(fit$covariates, x_nm)) {
    if (!is.null(extras) && nm %in% names(extras)) {
      grid[[nm]] <- extras[[nm]]
    } else if (nm %in% names(fit$support_ref$numeric)) {
      grid[[nm]] <- fit$support_ref$numeric[[nm]]$mean
    } else if (nm %in% names(fit$support_ref$factor_levels)) {
      grid[[nm]] <- factor(
        fit$support_ref$factor_levels[[nm]][[1]],
        levels = fit$support_ref$factor_levels[[nm]]
      )
    }
  }
  grid
}

ribbon_from_lines <- function(lines) {
  p <- sort(unique(lines$centile))
  if (length(p) < 2L) {
    return(tibble::tibble())
  }
  groups <- unique(lines$.group)
  rows <- list()
  i <- 1L
  j <- length(p)
  band_i <- 1L
  while (i < j) {
    for (g in groups) {
      lo <- lines[lines$centile == p[[i]] & eq_group(lines$.group, g), ]
      hi <- lines[lines$centile == p[[j]] & eq_group(lines$.group, g), ]
      rows[[length(rows) + 1L]] <- tibble::tibble(
        x = lo$x,
        ymin = lo$y,
        ymax = hi$y,
        band = paste(centile_label(p[[i]]), centile_label(p[[j]]), sep = "-"),
        .group = g,
        .rank = band_i
      )
    }
    i <- i + 1L
    j <- j - 1L
    band_i <- band_i + 1L
  }
  dplyr_bind(rows)
}

eq_group <- function(x, g) {
  if (length(g) == 1L && is.na(g)) {
    return(is.na(x))
  }
  x == g
}

fortify_kernel <- function(fit, lags = NULL) {
  if (!inherits(fit, "ref_dynamics")) {
    cli::cli_abort("{.fn fortify_kernel} expects a {.cls ref_dynamics} object.")
  }
  max_lag <- fit$lag_range[[2]]
  hi <- if (is.finite(max_lag) && max_lag > 0) {
    1.5 * max_lag
  } else {
    max(diff(fit$time_range), 1, na.rm = TRUE)
  }
  lags <- lags %||% seq(0, hi, length.out = 80)
  dplyr_bind(lapply(fit$outcomes, function(nm) {
    pr <- fit$processes[[nm]]
    if (!isTRUE(pr$identified)) {
      return(tibble::tibble(lag = lags, correlation = NA_real_, .outcome = nm))
    }
    tibble::tibble(lag = lags, correlation = process_kernel(lags, pr$psi, pr$process),
                   .outcome = nm)
  }))
}

# Conditional-quantile ("thrive") paths.
#
# For each anchor centile `p0` occupied at anchor time `t0`, the path gives
# the `thrive` quantile of the history-conditioned forecast distribution at
# `t0 + h`. On normal scores this is
#
#   z(h) = r(h) z0 + sqrt(1 - r(h)^2) qnorm(thrive),   z0 = qnorm(p0),
#
# so a thrive line is [ref_forecast()] evaluated on a grid of one-visit
# histories, not a separate estimator. Because r(h) comes from the fitted
# continuous-time kernel of [ref_dynamics()] rather than a binned age-by-age
# table, it is defined at every real lag (the path is a curve, not a
# one-step segment) and it carries a standard error from the kernel's
# observed information, which the delta method turns into an interval on the
# path.
#
# `centile` is the primary quantity and `value` is that centile read off the
# marginal reference distribution at `t0 + h`: a monotone reparametrisation
# of the outcome moves `value` and leaves `centile` alone.
fortify_thrive <- function(dynamic,
                           outcome = NULL,
                           anchors = c(0.05, 0.25, 0.5, 0.75, 0.95),
                           from = NULL,
                           horizon = 1,
                           thrive = 0.025,
                           level = 0.9,
                           x = NULL) {
  if (!inherits(dynamic, "ref_dynamics")) {
    cli::cli_abort("{.fn fortify_thrive} expects a {.cls ref_dynamics} object.")
  }
  outcome <- outcome %||% dynamic$outcomes[[1]]
  if (!outcome %in% dynamic$outcomes) {
    cli::cli_abort("{.val {outcome}} is not an outcome of the dynamics object.")
  }
  pr <- dynamic$processes[[outcome]]
  if (!isTRUE(pr$identified)) {
    cli::cli_abort(c(
      "The longitudinal process for {.field {outcome}} is not identified; no thrive line is possible.",
      i = pr$reason %||% "See print(dynamic)."
    ))
  }
  anchors <- sort(unique(as.numeric(anchors)))
  thrive <- sort(unique(as.numeric(thrive)))
  probs <- c(anchors, thrive)
  if (any(!is.finite(probs) | probs <= 0 | probs >= 1)) {
    cli::cli_abort("{.arg anchors} and {.arg thrive} must lie in (0, 1).")
  }
  fit <- dynamic$reference
  x_nm <- x %||% dynamic$time_name
  rng <- fit$support_ref$numeric[[x_nm]]
  if (is.null(rng)) {
    cli::cli_abort("No numeric covariate {.field {x_nm}} on the reference fit.")
  }
  h <- sort(unique(c(0, as.numeric(horizon))))
  h <- h[is.finite(h) & h >= 0]
  if (length(h) < 2L) {
    cli::cli_abort("{.arg horizon} must contain at least one positive step.")
  }
  lo <- max(rng$min, dynamic$time_range[[1]])
  hi <- min(rng$max, dynamic$time_range[[2]]) - max(h)
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    cli::cli_abort("The horizon leaves no anchor times inside the reference range.")
  }
  if (is.null(from)) {
    step <- max(signif((hi - lo) / 12, 1), .Machine$double.eps)
    from <- seq(ceiling(lo / step) * step, hi, by = step)
  }
  from <- sort(unique(as.numeric(from[is.finite(from)])))

  # r(h) and its standard error from the kernel's observed information.
  r_of <- function(theta, lag) process_kernel(lag, process_par(theta, pr$process), pr$process)
  r <- vapply(h, function(u) r_of(pr$theta, u), numeric(1))
  r_se <- rep(NA_real_, length(h))
  if (!is.null(pr$vcov)) {
    r_se <- vapply(h, function(u) {
      g <- numeric_gradient(function(th) r_of(th, u), pr$theta)
      sqrt(max(drop(t(g) %*% pr$vcov %*% g), 0))
    }, numeric(1))
  }

  ix <- expand.grid(h = seq_along(h), a = seq_along(anchors),
                    t = seq_along(from), k = seq_along(thrive))
  z0 <- stats::qnorm(anchors)[ix$a]
  kq <- stats::qnorm(thrive)[ix$k]
  rr <- r[ix$h]
  rse <- r_se[ix$h]
  s <- sqrt(pmax(1 - rr^2, 0))
  z <- rr * z0 + s * kq
  # delta method: dz/dr = z0 - r kq / sqrt(1 - r^2); exactly zero at h = 0.
  dz <- ifelse(s > 1e-8, z0 - rr * kq / s, 0)
  sd_z <- abs(dz) * rse
  q <- if (is.null(level) || !isTRUE(is.finite(level))) {
    NA_real_
  } else {
    stats::qnorm(1 - (1 - level) / 2)
  }
  t0 <- from[ix$t]
  tt <- t0 + h[ix$h]

  grid <- centile_grid(fit, x_nm, n = length(tt))
  grid[[x_nm]] <- tt
  dists <- predict_dists(fit, grid, uncertainty = "conditional")[[outcome]]
  qat <- function(zz) as.numeric(dist_quantile(dists, stats::pnorm(zz)))

  out <- tibble::tibble(
    .outcome = outcome,
    .segment = paste(ix$t, ix$a, ix$k, sep = "-"),
    anchor_time = t0,
    anchor_centile = anchors[ix$a],
    thrive = thrive[ix$k],
    horizon = h[ix$h],
    time = tt,
    r = rr,
    r_se = rse,
    z = z,
    centile = stats::pnorm(z),
    value = qat(z),
    lower = if (is.na(q)) NA_real_ else qat(z - q * sd_z),
    upper = if (is.na(q)) NA_real_ else qat(z + q * sd_z),
    support = classify_temporal_support(dynamic, t0, tt, rep(1L, length(tt)))
  )
  out <- out[order(out$.segment, out$horizon), , drop = FALSE]
  attr(out, "x_name") <- x_nm
  attr(out, "level") <- level
  out
}

