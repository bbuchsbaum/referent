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

fortify_visits <- function(data, id) {
  data <- tibble::as_tibble(data)
  visit_table(pull_column(data, rlang::enquo(id)))
}

visit_table <- function(id_vec) {
  n_per <- as.integer(table(as.character(id_vec)))
  tab <- as.data.frame(table(n_per), stringsAsFactors = FALSE)
  names(tab) <- c("n_visits", "n_subjects")
  tab$n_visits <- as.integer(as.character(tab$n_visits))
  tibble::as_tibble(tab)
}

fortify_coverage <- function(data, time, by = NULL) {
  coverage_frame(tibble::as_tibble(data), rlang::enquo(time), rlang::enquo(by))
}

coverage_frame <- function(data, time_quo, by_quo) {
  time_vec <- as.numeric(pull_column(data, time_quo))
  if (rlang::quo_is_null(by_quo) || rlang::quo_is_missing(by_quo)) {
    grp <- rep("all", length(time_vec))
    by_name <- NULL
  } else {
    grp <- as.character(pull_column(data, by_quo))
    by_name <- tryCatch(rlang::as_name(by_quo), error = function(e) "group")
  }
  tibble::tibble(
    time = time_vec,
    .group = grp,
    .by = by_name %||% "group"
  )
}

fortify_kernel <- function(fit, lags = NULL) {
  if (!inherits(fit, "norm_dynamics")) {
    cli::cli_abort("{.fn fortify_kernel} expects a {.cls norm_dynamics} object.")
  }
  hi <- max(c(diff(fit$time_range), fit$lag_range[[2]], 1), na.rm = TRUE)
  lags <- lags %||% seq(0, hi, length.out = 80)
  dplyr_bind(lapply(fit$outcomes, function(nm) {
    pr <- fit$processes[[nm]]
    if (!isTRUE(pr$identified)) {
      return(tibble::tibble(lag = lags, correlation = NA_real_, .outcome = nm))
    }
    r <- vapply(lags, function(h) {
      process_correlation(c(0, h), pr$psi, pr$process)[1, 2]
    }, numeric(1))
    tibble::tibble(lag = lags, correlation = r, .outcome = nm)
  }))
}

default_x <- function(fit) {
  if ("age" %in% fit$covariates) {
    return("age")
  }
  fit$support_ref$numeric_names[[1]] %||% fit$covariates[[1]]
}

as_col_name <- function(quo, default = NULL) {
  if (rlang::quo_is_null(quo) || rlang::quo_is_missing(quo)) {
    return(default)
  }
  tryCatch(rlang::as_name(quo), error = function(e) default)
}
