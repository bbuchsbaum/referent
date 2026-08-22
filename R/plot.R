#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' Graphics for reference models
#'
#' Every `autoplot()` method returns a `ggplot` object styled with
#' [theme_referent()]. Individual observations are overlays, not part of
#' the fitted geometry.
#'
#' @param object A fit, assessment, score table, forecast, transition,
#'   or dynamics object.
#' @param type Plot kind.
#' @param outcome Outcome name for multi-outcome fits.
#' @param x Covariate mapped to the x-axis (default `age` if present).
#' @param newdata Grid or observations to overlay.
#' @param by Optional grouping factor for faceted charts.
#' @param data Visit-level data for composition plots.
#' @param id Subject identifier for visits and trajectories.
#' @param time Time variable for coverage and trajectories.
#' @param centiles Probability levels for centile and fan charts.
#' @param ... Unused.
#' @return A ggplot, or a patchwork object when `type = "composition"`
#'   and patchwork is installed.
#' @exportS3Method ggplot2::autoplot
autoplot.norm_fit <- function(object,
                              type = c(
                                "centiles", "visits", "coverage",
                                "composition", "support", "trajectories",
                                "adaptation"
                              ),
                              outcome = NULL,
                              x = NULL,
                              newdata = NULL,
                              by = NULL,
                              data = NULL,
                              id = NULL,
                              time = NULL,
                              centiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
                              ...) {
  type <- match.arg(type)
  outcome <- outcome %||% object$outcomes[[1]]
  by_nm <- as_col_name(rlang::enquo(by))
  switch(
    type,
    centiles = plot_centiles(object, outcome, x, newdata, by_nm, centiles),
    visits = plot_visits(data %||% newdata, rlang::enquo(id)),
    coverage = plot_coverage(data %||% newdata, rlang::enquo(time), rlang::enquo(by)),
    composition = plot_composition(data %||% newdata, rlang::enquo(id), rlang::enquo(time), rlang::enquo(by)),
    support = plot_support(object, newdata %||% data, x),
    trajectories = plot_trajectories(
      object, data %||% newdata, outcome, x, rlang::enquo(id),
      rlang::enquo(time), by_nm, centiles
    ),
    adaptation = plot_adaptation(object)
  )
}

#' @exportS3Method ggplot2::autoplot
autoplot.norm_assessment <- function(object, type = c("calibration", "worm", "conditional"), ...) {
  type <- match.arg(type)
  switch(
    type,
    calibration = plot_calibration(object),
    worm = plot_worm(object),
    conditional = plot_conditional(object)
  )
}

#' @exportS3Method ggplot2::autoplot
autoplot.norm_scores <- function(object,
                                 type = c("profile", "heatmap"),
                                 id = NULL,
                                 ...) {
  type <- match.arg(type)
  if (identical(type, "profile")) {
    return(plot_profile(object, id))
  }
  plot_heatmap(object)
}

#' @exportS3Method ggplot2::autoplot
autoplot.norm_forecast <- function(object,
                                   type = c("fan", "thrive"),
                                   centiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
                                   ...) {
  type <- match.arg(type)
  plot_fan(object, centiles = centiles, thrive = identical(type, "thrive"))
}

#' @exportS3Method ggplot2::autoplot
autoplot.norm_transition <- function(object, type = c("velocity", "innovation", "change"), ...) {
  type <- match.arg(type)
  plot_transition(object, type)
}

#' @exportS3Method ggplot2::autoplot
autoplot.norm_dynamics <- function(object, type = c("kernel", "calibration"), ...) {
  type <- match.arg(type)
  plot_kernel(object)
}

plot_centiles <- function(fit, outcome, x, newdata, by_nm, centiles) {
  built <- fortify_centiles(fit, outcome, x, by_nm, centiles)
  p <- ggplot2::ggplot()
  if (nrow(built$ribbons)) {
    p <- p + ggplot2::geom_ribbon(
      data = built$ribbons,
      ggplot2::aes(x = .data$x, ymin = .data$ymin, ymax = .data$ymax, fill = .data$band),
      colour = NA
    ) +
      ggplot2::scale_fill_manual(
        values = band_fills(length(unique(built$ribbons$band))),
        name = NULL
      )
  }
  other <- built$lines[built$lines$.label != "Median", , drop = FALSE]
  mid <- built$lines[built$lines$.label == "Median", , drop = FALSE]
  if (nrow(other)) {
    p <- p + ggplot2::geom_line(
      data = other,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$.label),
      colour = referent_cols()$ink,
      linewidth = 0.4,
      linetype = "dashed"
    )
  }
  if (nrow(mid)) {
    p <- p + ggplot2::geom_line(
      data = mid,
      ggplot2::aes(x = .data$x, y = .data$y),
      colour = referent_cols()$ink,
      linewidth = 0.95
    )
  }
  if (has_groups(built$lines$.group)) {
    p <- p + ggplot2::facet_wrap(~.group, nrow = 1)
  }
  if (!is.null(newdata) && all(c(built$x_name, outcome) %in% names(newdata)) &&
      isTRUE(nrow(newdata) > 0)) {
    pts <- tibble::as_tibble(newdata)
    if (!is.null(by_nm) && by_nm %in% names(pts) && has_groups(built$lines$.group)) {
      pts$.group <- as.character(pts[[by_nm]])
    }
    p <- p + ggplot2::geom_point(
      data = pts,
      ggplot2::aes(x = .data[[built$x_name]], y = .data[[outcome]]),
      inherit.aes = FALSE,
      colour = referent_cols()$ink,
      alpha = 0.28,
      size = 1.1
    )
  }
  finish_plot(
    p + pretty_x() + pretty_y(),
    title = NULL,
    xlab = built$x_name,
    ylab = outcome
  )
}

band_fills <- function(n) {
  base <- grDevices::col2rgb(referent_cols()$band) / 255
  alphas <- seq(0.10, 0.28, length.out = max(n, 1L))
  vapply(seq_len(n), function(i) {
    grDevices::rgb(base[1], base[2], base[3], alphas[[i]])
  }, character(1))
}

has_groups <- function(x) {
  length(unique(x[!is.na(x)])) > 1L
}

plot_visits <- function(data, id_quo) {
  if (is.null(data)) {
    cli::cli_abort("{.arg data} is required for a visit plot.")
  }
  if (rlang::quo_is_null(id_quo) || rlang::quo_is_missing(id_quo)) {
    cli::cli_abort("{.arg id} is required for a visit plot.")
  }
  tab <- visit_table(pull_column(tibble::as_tibble(data), id_quo))
  p <- ggplot2::ggplot(tab, ggplot2::aes(.data$n_visits, .data$n_subjects)) +
    ggplot2::geom_col(fill = referent_cols()$accent, width = 0.72) +
    ggplot2::scale_x_continuous(breaks = tab$n_visits) +
    pretty_y()
  finish_plot(p, title = "Visits per subject", xlab = "Number of visits", ylab = "Number of subjects")
}

plot_coverage <- function(data, time_quo, by_quo) {
  if (is.null(data)) {
    cli::cli_abort("{.arg data} is required for a coverage plot.")
  }
  if (rlang::quo_is_null(time_quo) || rlang::quo_is_missing(time_quo)) {
    cli::cli_abort("{.arg time} is required for a coverage plot.")
  }
  df <- coverage_frame(tibble::as_tibble(data), time_quo, by_quo)
  n_g <- length(unique(df$.group))
  if (n_g > 8L && has_pkg("ggridges")) {
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$time, .data$.group, fill = .data$.group)) +
      ggridges::geom_density_ridges(
        colour = "white", linewidth = 0.15, scale = 1.05,
        rel_min_height = 0.01, show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = rep(referent_cols()$band, length.out = n_g)) +
      pretty_x()
    return(finish_plot(p, title = NULL, xlab = "years", ylab = NULL))
  }
  grp_lab <- df$.by[[1]] %||% "group"
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$time, fill = .data$.group, colour = .data$.group)) +
    ggplot2::geom_density(alpha = 0.25, linewidth = 0.4) +
    scale_fill_referent(name = grp_lab) +
    scale_colour_referent(name = grp_lab) +
    pretty_x()
  if (n_g > 4L) {
    p <- p + ggplot2::facet_wrap(~.group, ncol = 1, scales = "free_y") +
      ggplot2::theme(legend.position = "none")
  }
  finish_plot(p, title = NULL, xlab = "years", ylab = "Density")
}

plot_composition <- function(data, id_quo, time_quo, by_quo) {
  p_vis <- plot_visits(data, id_quo)
  p_cov <- plot_coverage(data, time_quo, by_quo)
  if (has_pkg("patchwork")) {
    return(patchwork::wrap_plots(p_vis, p_cov, ncol = 2) +
             patchwork::plot_annotation(tag_levels = "A"))
  }
  list(visits = p_vis, coverage = p_cov)
}

plot_support <- function(fit, newdata, x) {
  if (is.null(newdata)) {
    cli::cli_abort("{.arg newdata} is required for a support plot.")
  }
  x_nm <- x %||% default_x(fit)
  df <- tibble::as_tibble(newdata)
  df$support <- classify_support(fit$support_ref, df)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data[[x_nm]], fill = .data$support)) +
    ggplot2::geom_histogram(bins = 30, position = "identity", alpha = 0.7, colour = NA) +
    support_fill_scale() +
    pretty_x() + pretty_y()
  finish_plot(p, title = NULL, xlab = x_nm, ylab = "Count")
}

plot_trajectories <- function(fit, data, outcome, x, id_quo, time_quo, by_nm, centiles) {
  if (is.null(data)) {
    cli::cli_abort("{.arg data} is required for a trajectory plot.")
  }
  x_nm <- x %||% default_x(fit)
  p <- plot_centiles(fit, outcome, x_nm, newdata = NULL, by_nm, centiles)
  df <- tibble::as_tibble(data)
  df$.id <- as.character(pull_column(df, id_quo, default = seq_len(nrow(df))))
  if (!rlang::quo_is_null(time_quo) && !rlang::quo_is_missing(time_quo)) {
    t_nm <- as_col_name(time_quo, x_nm)
    if (t_nm %in% names(df)) {
      x_nm <- t_nm
    }
  }
  if (!is.null(by_nm) && by_nm %in% names(df)) {
    df$.group <- as.character(df[[by_nm]])
  }
  p +
    ggplot2::geom_path(
      data = df,
      ggplot2::aes(x = .data[[x_nm]], y = .data[[outcome]], group = .data$.id),
      inherit.aes = FALSE,
      colour = referent_cols()$ink,
      alpha = 0.45,
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = .data[[x_nm]], y = .data[[outcome]]),
      inherit.aes = FALSE,
      colour = referent_cols()$ink,
      alpha = 0.7,
      size = 1.2
    )
}

plot_calibration <- function(assessment) {
  df <- assessment$marginal
  long <- tibble::tibble(
    outcome = rep(df$.outcome, 5),
    nominal = rep(c(0.5, 0.8, 0.9, 0.95, 0.99), each = nrow(df)),
    observed = c(df$cover_50, df$cover_80, df$cover_90, df$cover_95, df$cover_99)
  )
  p <- ggplot2::ggplot(long, ggplot2::aes(.data$nominal, .data$observed, colour = .data$outcome)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = referent_cols()$muted) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    scale_colour_referent() +
    ggplot2::coord_equal(xlim = c(0.45, 1), ylim = c(0.45, 1)) +
    pretty_x() + pretty_y()
  finish_plot(p, title = NULL, xlab = "Nominal coverage", ylab = "Observed coverage")
}

plot_worm <- function(assessment) {
  sc <- assessment$scores
  sc <- sc[is.finite(sc$z), , drop = FALSE]
  pieces <- lapply(split(sc, sc$.outcome), function(one) {
    one <- one[order(one$z), , drop = FALSE]
    one$expected <- stats::qnorm(stats::ppoints(nrow(one)))
    one$worm <- one$z - one$expected
    one
  })
  sc <- dplyr_bind(pieces)
  p <- ggplot2::ggplot(sc, ggplot2::aes(.data$expected, .data$worm, colour = .data$.outcome)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = referent_cols()$muted) +
    ggplot2::geom_point(alpha = 0.45, size = 1.2) +
    ggplot2::geom_smooth(se = FALSE, linewidth = 0.5, method = "loess", formula = y ~ x, span = 0.7) +
    scale_colour_referent() +
    pretty_x() + pretty_y()
  finish_plot(p, title = NULL, xlab = "Expected Z", ylab = "Observed minus expected")
}

plot_conditional <- function(assessment) {
  p <- ggplot2::ggplot(assessment$conditional, ggplot2::aes(.data$.outcome, .data$location_drift)) +
    ggplot2::geom_col(fill = referent_cols()$band, width = 0.65) +
    pretty_y()
  finish_plot(p, title = NULL, xlab = NULL, ylab = "max |E[z | x]|")
}

plot_profile <- function(scores, id = NULL) {
  if (!is.null(id)) {
    one <- scores[as.character(scores$.id) == as.character(id), , drop = FALSE]
  } else {
    key <- scores$.id[[1]]
    one <- scores[scores$.id == key, , drop = FALSE]
  }
  p <- ggplot2::ggplot(one, ggplot2::aes(.data$.outcome, .data$z)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = referent_cols()$muted) +
    ggplot2::geom_col(fill = referent_cols()$band, width = 0.65) +
    pretty_y()
  finish_plot(p, title = NULL, xlab = NULL, ylab = "z")
}

plot_heatmap <- function(scores) {
  p <- ggplot2::ggplot(scores, ggplot2::aes(.data$.id, .data$.outcome, fill = .data$z)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_referent_diverging()
  finish_plot(p, title = NULL, xlab = "Subject", ylab = NULL)
}

plot_adaptation <- function(object) {
  adapt <- if (inherits(object, "norm_adaptation")) object else object$adaptation
  if (is.null(adapt) || !inherits(adapt, "norm_adaptation")) {
    return(finish_plot(ggplot2::ggplot(), title = "No adaptation attached"))
  }
  rows <- list()
  for (nm in names(adapt$offsets)) {
    for (g in names(adapt$offsets[[nm]])) {
      off <- adapt$offsets[[nm]][[g]]
      rows[[length(rows) + 1L]] <- tibble::tibble(
        outcome = nm, group = g,
        location = off$location, scale = off$scale
      )
    }
  }
  df <- dplyr_bind(rows)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$group, .data$location, fill = .data$outcome)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    scale_fill_referent() +
    pretty_y()
  finish_plot(p, title = NULL, xlab = NULL, ylab = "Location offset")
}

plot_fan <- function(forecast, centiles, thrive = FALSE) {
  times <- forecast$times
  dists <- forecast$dist
  lines <- dplyr_bind(lapply(centiles, function(p) {
    tibble::tibble(
      time = times,
      x = times,
      y = vapply(dists, function(d) as.numeric(quantile(d, p)), numeric(1)),
      centile = p,
      .label = centile_label(p),
      .group = NA_character_
    )
  }))
  ribbons <- ribbon_from_lines(lines)
  p <- ggplot2::ggplot()
  if (!thrive && nrow(ribbons)) {
    p <- p + ggplot2::geom_ribbon(
      data = ribbons,
      ggplot2::aes(x = .data$x, ymin = .data$ymin, ymax = .data$ymax, fill = .data$band),
      colour = NA
    ) +
      ggplot2::scale_fill_manual(values = band_fills(length(unique(ribbons$band))), name = NULL)
  }
  other <- lines[lines$.label != "Median", , drop = FALSE]
  mid <- lines[lines$.label == "Median", , drop = FALSE]
  if (nrow(other)) {
    p <- p + ggplot2::geom_line(
      data = other,
      ggplot2::aes(x = .data$time, y = .data$y, group = .data$.label),
      colour = referent_cols()$ink,
      linewidth = 0.4,
      linetype = "dashed"
    )
  }
  if (nrow(mid)) {
    p <- p + ggplot2::geom_line(
      data = mid,
      ggplot2::aes(x = .data$time, y = .data$y),
      colour = referent_cols()$ink,
      linewidth = 0.95
    )
  }
  hist <- forecast$history
  t_nm <- forecast$time_name %||% "age"
  if (!is.null(hist) && t_nm %in% names(hist) && forecast$outcome %in% names(hist)) {
    p <- p +
      ggplot2::geom_path(
        data = hist,
        ggplot2::aes(x = .data[[t_nm]], y = .data[[forecast$outcome]]),
        inherit.aes = FALSE,
        colour = referent_cols()$accent,
        linewidth = 0.7
      ) +
      ggplot2::geom_point(
        data = hist,
        ggplot2::aes(x = .data[[t_nm]], y = .data[[forecast$outcome]]),
        inherit.aes = FALSE,
        colour = referent_cols()$accent,
        size = 1.8
      )
  }
  finish_plot(
    p + pretty_x() + pretty_y(),
    title = NULL,
    xlab = t_nm,
    ylab = forecast$outcome
  )
}

plot_transition <- function(object, type) {
  y_nm <- switch(
    type,
    velocity = "observed_velocity",
    innovation = "innovation_z",
    change = "change_z"
  )
  ylab <- switch(
    type,
    velocity = "Observed velocity",
    innovation = "Innovation Z",
    change = "Change Z"
  )
  df <- tibble::as_tibble(object)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$.dt, .data[[y_nm]], colour = .data$support)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = referent_cols()$muted) +
    ggplot2::geom_point(alpha = 0.7, size = 1.6) +
    ggplot2::scale_colour_manual(values = referent_cols()$support, name = "Support") +
    pretty_x() + pretty_y()
  finish_plot(p, title = NULL, xlab = "Elapsed time", ylab = ylab)
}

plot_kernel <- function(dyn) {
  df <- fortify_kernel(dyn)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$lag, .data$correlation, colour = .data$.outcome)) +
    ggplot2::geom_hline(yintercept = 0, colour = referent_cols()$muted, linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.7) +
    scale_colour_referent() +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    pretty_x() + pretty_y()
  finish_plot(p, title = NULL, xlab = "Lag", ylab = "Process correlation")
}
