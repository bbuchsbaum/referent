#' Tidy plot data for reference charts
#'
#' These helpers return tibbles. `autoplot()` consumes them; users can
#' also rebuild a chart by hand.
#'
#' @name fortify_centiles
#' @param fit A [norm_fit].
#' @param outcome Outcome name.
#' @param x Numeric covariate for the x-axis.
#' @param by Optional grouping factor (column name). One chart is built
#'   per level; other covariates stay at their reference value.
#' @param centiles Probability levels to evaluate.
#' @param n Grid length along `x`.
#' @return A list of class `norm_centile_data` with `lines` and `ribbons` tibbles.
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
    class = "norm_centile_data"
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
  if (!inherits(fit, "norm_dynamics")) {
    cli::cli_abort("{.fn fortify_kernel} expects a {.cls norm_dynamics} object.")
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
