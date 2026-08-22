#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' Graphics for reference models
#'
#' Every `autoplot()` method returns a `ggplot` object styled with
#' [theme_referent()]. Individual observations are overlays, not part of
#' the fitted geometry.
#'
#' @details
#' For a [norm_fit]:
#' * `"centiles"`: the fitted conditional distribution along `x`
#'   (median, dashed outer centiles, and shaded bands), optionally faceted
#'   by a factor covariate and with `newdata` overlaid as points.
#' * `"trajectories"`: the centile chart with subject paths from `data`
#'   drawn over it.
#' * `"support"`: where `newdata` falls relative to the reference support
#'   (histogram of `x` coloured by [norm_support()] status).
#' * `"visits"`: visits per subject in the reference sample (from the `id`
#'   declared at fit time) or in `data`.
#' * `"coverage"`: the distribution of `x` in `data` against the
#'   reference range and the 2nd/98th-percentile edges of the fit.
#' * `"composition"`: visits and coverage side by side (a patchwork when
#'   that package is installed, otherwise a list of two plots).
#' * `"adaptation"`: location offsets (with standard errors) estimated by
#'   [norm_adapt()].
#'
#' For a `norm_dynamics`, `"kernel"` draws the fitted process correlation
#' against lag and `"calibration"` draws a normal Q-Q plot of held-out
#' innovation Z from [norm_transition()] on `data`; the latter needs
#' `data`, `id`, and `time`.
#'
#' @param object A fit, assessment, score table, forecast, transition,
#'   or dynamics object.
#' @param type Plot kind.
#' @param outcome Outcome name for multi-outcome fits.
#' @param x Covariate mapped to the x-axis (default `age` if present).
#' @param newdata Observations to overlay or to classify.
#' @param by Optional grouping factor for faceted charts.
#' @param data Visit-level data for composition, trajectory, and
#'   dynamics-calibration plots.
#' @param id Subject identifier column.
#' @param time Time column.
#' @param centiles Probability levels for centile and fan charts.
#' @param level Coverage of the simulation envelope on the QQ plot.
#' @param ... Unused.
#' @return A ggplot, or a patchwork object when `type = "composition"`
#'   and patchwork is installed.
#' @examples
#' ref <- norm_simulate(150, seed = 1)
#' fit <- norm_fit(norm_spec(norm_gaussian(), ~ s(age, k = 5) + sex), ref, "y")
#' autoplot(fit, type = "centiles", by = sex, newdata = ref[1:20, ])
#' sc <- predict(fit, newdata = ref[1:6, ], uncertainty = "conditional")
#' autoplot(sc, type = "heatmap")
#' @name autoplot.norm_fit
#' @exportS3Method ggplot2::autoplot
autoplot.norm_fit <- function(object,
                              type = c(
                                "centiles", "trajectories", "support",
                                "visits", "coverage", "composition",
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
  by_quo <- rlang::enquo(by)
  by_nm <- as_col_name(by_quo)
  id_quo <- rlang::enquo(id)
  time_quo <- rlang::enquo(time)
  data <- data %||% newdata
  switch(
    type,
    centiles = plot_centiles(object, outcome, x, newdata, by_nm, centiles),
    trajectories = plot_trajectories(
      object, data, outcome, x, id_quo, time_quo, by_nm, centiles
    ),
    support = plot_support(object, data, x),
    visits = plot_visits(object, data, id_quo),
    coverage = plot_coverage(object, data, x, time_quo, by_quo),
    composition = plot_composition(object, data, x, id_quo, time_quo, by_quo),
    adaptation = plot_adaptation(object)
  )
}

#' @rdname autoplot.norm_fit
#' @exportS3Method ggplot2::autoplot
autoplot.norm_assessment <- function(object,
                                     type = c("calibration", "worm", "qq", "conditional"),
                                     level = 0.95,
                                     ...) {
  type <- match.arg(type)
  switch(
    type,
    calibration = plot_calibration(object),
    worm = plot_worm(object),
    qq = plot_qq(object, level = level),
    conditional = plot_conditional(object)
  )
}

#' @rdname autoplot.norm_fit
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

#' @rdname autoplot.norm_fit
#' @exportS3Method ggplot2::autoplot
autoplot.norm_forecast <- function(object,
                                   type = c("fan"),
                                   centiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
                                   ...) {
  type <- match.arg(type)
  plot_fan(object, centiles = centiles)
}

#' @rdname autoplot.norm_fit
#' @exportS3Method ggplot2::autoplot
autoplot.norm_transition <- function(object, type = c("velocity", "innovation", "change"), ...) {
  type <- match.arg(type)
  plot_transition(object, type)
}

#' @rdname autoplot.norm_fit
#' @exportS3Method ggplot2::autoplot
autoplot.norm_dynamics <- function(object,
                                   type = c("kernel", "calibration"),
                                   data = NULL,
                                   id = NULL,
                                   time = NULL,
                                   ...) {
  type <- match.arg(type)
  if (identical(type, "kernel")) {
    return(plot_kernel(object))
  }
  plot_dynamics_calibration(object, data, rlang::enquo(id), rlang::enquo(time))
}

# --- shared centile geometry ----------------------------------------------

# Ribbon and line layers shared by the centile chart and the forecast fan.
# `lines` needs x, y, .label (with "Median" for the centre line); `ribbons`
# needs x, ymin, ymax, band, .rank (outer band first).
centile_layers <- function(lines, ribbons) {
  layers <- list()
  if (nrow(ribbons)) {
    ord <- unique(ribbons[order(ribbons$.rank), "band", drop = TRUE])
    ribbons$band <- factor(ribbons$band, levels = ord)
    fills <- band_fills(length(ord))
    names(fills) <- ord
    layers <- c(layers, list(
      ggplot2::geom_ribbon(
        data = ribbons,
        ggplot2::aes(x = .data$x, ymin = .data$ymin, ymax = .data$ymax, fill = .data$band),
        colour = NA
      ),
      ggplot2::scale_fill_manual(values = fills, breaks = ord, name = NULL)
    ))
  }
  other <- lines[lines$.label != "Median", , drop = FALSE]
  mid <- lines[lines$.label == "Median", , drop = FALSE]
  if (nrow(other)) {
    layers <- c(layers, list(ggplot2::geom_line(
      data = other,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$.label),
      colour = referent_cols()$ink,
      linewidth = 0.4,
      linetype = "dashed"
    )))
  }
  if (nrow(mid)) {
    layers <- c(layers, list(ggplot2::geom_line(
      data = mid,
      ggplot2::aes(x = .data$x, y = .data$y),
      colour = referent_cols()$ink,
      linewidth = 0.95
    )))
  }
  layers
}

# Band fills from the outermost (lightest) to the innermost band; the
# central 25-75 band sits at alpha 0.25.
band_fills <- function(n) {
  base <- grDevices::col2rgb(referent_cols()$band) / 255
  alphas <- if (n <= 1L) 0.25 else seq(0.12, 0.25, length.out = n)
  vapply(seq_len(n), function(i) {
    grDevices::rgb(base[1], base[2], base[3], alphas[[i]])
  }, character(1))
}

has_groups <- function(x) {
  length(unique(x[!is.na(x)])) > 1L
}

# --- norm_fit -------------------------------------------------------------

plot_centiles <- function(fit, outcome, x, newdata, by_nm, centiles) {
  built <- fortify_centiles(fit, outcome, x, by_nm, centiles)
  p <- ggplot2::ggplot() + centile_layers(built$lines, built$ribbons)
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

plot_trajectories <- function(fit, data, outcome, x, id_quo, time_quo, by_nm, centiles) {
  if (is.null(data)) {
    cli::cli_abort("{.arg data} is required for a trajectory plot.")
  }
  x_nm <- x %||% default_x(fit)
  p <- plot_centiles(fit, outcome, x_nm, newdata = NULL, by_nm, centiles)
  df <- tibble::as_tibble(data)
  df$.id <- as.character(pull_column(df, id_quo, default = seq_len(nrow(df))))
  t_nm <- as_col_name(time_quo, x_nm)
  if (t_nm %in% names(df)) {
    x_nm <- t_nm
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

plot_support <- function(fit, newdata, x) {
  if (is.null(newdata)) {
    cli::cli_abort("{.arg newdata} is required for a support plot.")
  }
  x_nm <- x %||% default_x(fit)
  df <- tibble::as_tibble(newdata)
  df$support <- classify_support(fit$support_ref, df)$support
  r <- fit$support_ref$numeric[[x_nm]]
  p <- ggplot2::ggplot(df, ggplot2::aes(.data[[x_nm]], fill = .data$support))
  if (!is.null(r)) {
    p <- p + ggplot2::geom_vline(
      xintercept = c(r$min, r$max), colour = referent_cols()$muted, linewidth = 0.4
    ) +
      ggplot2::geom_vline(
        xintercept = c(r$edge_lo, r$edge_hi), colour = referent_cols()$muted,
        linetype = "dashed", linewidth = 0.3
      )
  }
  p <- p +
    ggplot2::geom_histogram(bins = 30, position = "stack", colour = NA) +
    support_scale(df$support, "fill") +
    pretty_x() + pretty_y()
  single <- single_level_legend(df$support)
  p <- finish_plot(
    p,
    title = NULL,
    subtitle = if (single) paste0("all rows: ", unique(df$support)) else NULL,
    xlab = x_nm,
    ylab = "Count"
  )
  if (single) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p
}

plot_visits <- function(fit, data, id_quo) {
  id_given <- !(rlang::quo_is_null(id_quo) || rlang::quo_is_missing(id_quo))
  if (is.null(data) || !id_given) {
    if (is.null(fit$id_name)) {
      cli::cli_abort(c(
        "No subject identifier was declared at fit time.",
        i = "Pass {.arg id} to {.fn norm_fit}, or supply {.arg data} and {.arg id}."
      ))
    }
    ids <- fit$id
    source <- "reference sample"
  } else {
    ids <- pull_column(tibble::as_tibble(data), id_quo)
    source <- "supplied data"
  }
  tab <- visit_table(ids)
  p <- ggplot2::ggplot(tab, ggplot2::aes(.data$n_visits, .data$n_subjects)) +
    ggplot2::geom_col(fill = referent_cols()$band, width = 0.72) +
    ggplot2::scale_x_continuous(breaks = tab$n_visits) +
    pretty_y()
  finish_plot(
    p,
    title = NULL,
    subtitle = paste0(sum(tab$n_subjects), " subjects, ", source),
    xlab = "Visits per subject",
    ylab = "Subjects"
  )
}

plot_coverage <- function(fit, data, x, time_quo, by_quo) {
  if (is.null(data)) {
    cli::cli_abort("{.arg data} is required for a coverage plot.")
  }
  data <- tibble::as_tibble(data)
  x_nm <- as_col_name(time_quo, x %||% default_x(fit))
  if (!x_nm %in% names(data)) {
    cli::cli_abort("{.arg data} has no column {.field {x_nm}}.")
  }
  df <- coverage_frame(data, rlang::quo(.data[[x_nm]]), by_quo)
  df$time <- as.numeric(data[[x_nm]])
  r <- fit$support_ref$numeric[[x_nm]]
  grp_lab <- df$.by[[1]] %||% "group"
  n_g <- length(unique(df$.group))
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$time))
  if (!is.null(r)) {
    p <- p +
      ggplot2::annotate(
        "rect", xmin = r$edge_lo, xmax = r$edge_hi, ymin = -Inf, ymax = Inf,
        fill = referent_cols()$band, alpha = 0.08
      ) +
      ggplot2::geom_vline(
        xintercept = c(r$min, r$max), colour = referent_cols()$muted, linewidth = 0.4
      )
  }
  if (n_g > 1L) {
    p <- p +
      ggplot2::geom_density(
        ggplot2::aes(fill = .data$.group, colour = .data$.group),
        alpha = 0.25, linewidth = 0.4
      ) +
      scale_fill_referent(name = grp_lab) +
      scale_colour_referent(name = grp_lab)
    if (n_g > 4L) {
      p <- p + ggplot2::facet_wrap(~.group, ncol = 1, scales = "free_y")
    }
  } else {
    p <- p + ggplot2::geom_density(
      fill = referent_cols()$band, colour = referent_cols()$band,
      alpha = 0.25, linewidth = 0.4
    )
  }
  p <- finish_plot(
    p + pretty_x(),
    title = NULL,
    subtitle = if (is.null(r)) NULL else "shaded: central 96% of the reference; lines: reference range",
    xlab = x_nm,
    ylab = "Density"
  )
  if (n_g > 4L) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p
}

plot_composition <- function(fit, data, x, id_quo, time_quo, by_quo) {
  p_vis <- plot_visits(fit, data, id_quo)
  p_cov <- plot_coverage(fit, data, x, time_quo, by_quo)
  if (has_pkg("patchwork")) {
    return(patchwork::wrap_plots(p_vis, p_cov, ncol = 2) +
             patchwork::plot_annotation(tag_levels = "A"))
  }
  list(visits = p_vis, coverage = p_cov)
}

plot_adaptation <- function(object) {
  adapt <- if (inherits(object, "norm_adaptation")) object else object$adaptation
  if (is.null(adapt) || !inherits(adapt, "norm_adaptation")) {
    cli::cli_abort("No adaptation is attached; see {.fn norm_adapt}.")
  }
  rows <- list()
  for (nm in names(adapt$offsets)) {
    for (g in names(adapt$offsets[[nm]])) {
      off <- adapt$offsets[[nm]][[g]]
      rows[[length(rows) + 1L]] <- tibble::tibble(
        outcome = nm, group = g,
        location = off$location, location_se = off$location_se %||% NA_real_,
        scale = off$scale, n = off$n
      )
    }
  }
  df <- dplyr_bind(rows)
  if (!is.null(adapt$by) && any(df$group != ".all")) {
    df <- df[df$group != ".all", , drop = FALSE]
  } else {
    df$group <- "pooled"
  }
  df$outcome <- factor(df$outcome, levels = unique(df$outcome))
  df$group <- factor(df$group, levels = unique(df$group))
  multi <- !single_level_legend(df$outcome)
  dodge <- ggplot2::position_dodge(width = 0.6)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$group, .data$location)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = referent_cols()$muted)
  if (multi) {
    p <- p +
      ggplot2::geom_linerange(
        ggplot2::aes(
          ymin = .data$location - 2 * .data$location_se,
          ymax = .data$location + 2 * .data$location_se,
          colour = .data$outcome
        ),
        position = dodge, linewidth = 0.5, na.rm = TRUE
      ) +
      ggplot2::geom_point(ggplot2::aes(colour = .data$outcome), position = dodge, size = 2.4) +
      scale_colour_referent(name = "Outcome")
  } else {
    p <- p +
      ggplot2::geom_linerange(
        ggplot2::aes(
          ymin = .data$location - 2 * .data$location_se,
          ymax = .data$location + 2 * .data$location_se
        ),
        colour = referent_cols()$band, linewidth = 0.5, na.rm = TRUE
      ) +
      ggplot2::geom_point(colour = referent_cols()$band, size = 2.4)
  }
  finish_plot(
    p + pretty_y(),
    title = NULL,
    subtitle = "location offset, +/- 2 SE",
    xlab = adapt$by %||% NULL,
    ylab = "Location offset"
  )
}

# --- norm_assessment --------------------------------------------------------

plot_calibration <- function(assessment) {
  df <- assessment$marginal
  levels <- c(0.5, 0.8, 0.9, 0.95, 0.99)
  long <- tibble::tibble(
    outcome = rep(df$.outcome, length(levels)),
    nominal = rep(levels, each = nrow(df)),
    observed = unlist(lapply(paste0("cover_", c("50", "80", "90", "95", "99")), function(cl) df[[cl]]))
  )
  long$outcome <- factor(long$outcome, levels = unique(df$.outcome))
  multi <- !single_level_legend(long$outcome)
  p <- ggplot2::ggplot(long, ggplot2::aes(.data$nominal, .data$observed)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = referent_cols()$muted)
  if (multi) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(colour = .data$outcome, group = .data$outcome), linewidth = 0.5) +
      ggplot2::geom_point(ggplot2::aes(colour = .data$outcome), size = 2) +
      scale_colour_referent(name = "Outcome")
  } else {
    p <- p +
      ggplot2::geom_line(colour = referent_cols()$band, linewidth = 0.5) +
      ggplot2::geom_point(colour = referent_cols()$band, size = 2)
  }
  p <- p +
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
  sc$.outcome <- factor(sc$.outcome, levels = unique(assessment$scores$.outcome))
  multi <- !single_level_legend(sc$.outcome)
  p <- ggplot2::ggplot(sc, ggplot2::aes(.data$expected, .data$worm)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = referent_cols()$muted)
  if (multi) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(colour = .data$.outcome), alpha = 0.45, size = 1.2) +
      ggplot2::geom_smooth(
        ggplot2::aes(colour = .data$.outcome), se = FALSE, linewidth = 0.5,
        method = "loess", formula = y ~ x, span = 0.7
      ) +
      scale_colour_referent(name = "Outcome")
  } else {
    p <- p +
      ggplot2::geom_point(colour = referent_cols()$ink, alpha = 0.45, size = 1.2) +
      ggplot2::geom_smooth(
        colour = referent_cols()$band, se = FALSE, linewidth = 0.6,
        method = "loess", formula = y ~ x, span = 0.7
      )
  }
  finish_plot(
    p + pretty_x() + pretty_y(),
    title = NULL, xlab = "Expected Z", ylab = "Observed minus expected Z"
  )
}

# Normal QQ plot of the Z scores with the envelope a calibrated model implies.
#
# The reference is exact rather than nominal. Under a correctly specified
# reference model the held-out centiles are iid uniform, so the i-th ordered
# centile is Beta(i, n - i + 1) and the i-th ordered Z is its normal
# transform. That gives the pointwise band in closed form. The simultaneous
# band uses equal local levels (Aldor-Noiman et al. 2013, Am. Stat. 67:249):
# the common two-sided local level gamma is calibrated by Monte Carlo so that
# the whole ordered sample stays inside with probability `level`.
#
# The corner tally reports both bands. Counting only the simultaneous band
# would let a panel headline "0 outside" while a large share of points sit
# outside the narrower pointwise band that the same panel draws -- the two
# bands answer different questions, so the label names both and states what
# the pointwise rate should be.
plot_qq <- function(assessment, level = 0.95) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
    level <= 0 || level >= 1) {
    cli::cli_abort("{.arg level} must be a single probability in (0, 1).")
  }
  sc <- assessment$scores
  sc <- sc[is.finite(sc$z), , drop = FALSE]
  if (!nrow(sc)) {
    return(finish_plot(ggplot2::ggplot(), title = "No finite Z scores"))
  }
  by_out <- split(sc, sc$.outcome)
  # The envelope depends only on n, so outcomes of equal size share one.
  sizes <- unique(vapply(by_out, nrow, integer(1)))
  envs <- stats::setNames(lapply(sizes, qq_envelope, level = level), sizes)
  pieces <- lapply(by_out, function(one) {
    one <- one[order(one$z), , drop = FALSE]
    env <- envs[[as.character(nrow(one))]]
    one$expected <- env$expected
    one$point_lo <- env$point_lo
    one$point_hi <- env$point_hi
    one$sim_lo <- env$sim_lo
    one$sim_hi <- env$sim_hi
    one$outside <- one$z < env$sim_lo | one$z > env$sim_hi
    one
  })
  df <- dplyr_bind(pieces)

  cols <- referent_cols()
  pct <- format(100 * level, trim = TRUE)
  lab_sim <- paste0(pct, "% simultaneous")
  lab_pt <- paste0(pct, "% pointwise")
  fills <- stats::setNames(band_fills(2L), c(lab_sim, lab_pt))

  # The bands are defined only on the expected quantiles, so the x view is set
  # by those; y follows the data, which can run past them in a heavy tail.
  xlim <- range(df$expected, finite = TRUE)
  ylim <- range(c(df$expected, df$z), finite = TRUE)
  ylim <- ylim + c(-1, 1) * 0.06 * diff(ylim)
  brk <- pretty(ylim, n = 7L)

  # One corner label per panel, naming the band each count belongs to.
  tally <- dplyr_bind(lapply(pieces, function(one) {
    n <- nrow(one)
    k_sim <- sum(one$z < one$sim_lo | one$z > one$sim_hi, na.rm = TRUE)
    k_pt <- sum(one$z < one$point_lo | one$z > one$point_hi, na.rm = TRUE)
    tibble::tibble(
      .outcome = one$.outcome[[1L]],
      expected = xlim[[1L]] + 0.04 * diff(xlim),
      z = ylim[[2L]] - 0.03 * diff(ylim),
      label = sprintf(
        "outside simultaneous: %d of %d\noutside pointwise: %d (%.1f%%)",
        k_sim, n, k_pt, 100 * k_pt / n
      )
    )
  }))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$expected, y = .data$z)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$sim_lo, ymax = .data$sim_hi, fill = lab_sim),
      colour = NA
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$point_lo, ymax = .data$point_hi, fill = lab_pt),
      colour = NA
    ) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      colour = cols$muted, linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        colour = .data$outside, size = .data$outside, alpha = .data$outside
      )
    ) +
    ggplot2::geom_text(
      data = tally, ggplot2::aes(label = .data$label),
      hjust = 0, vjust = 1, size = 2.7, colour = cols$muted,
      inherit.aes = TRUE, show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = fills, breaks = c(lab_pt, lab_sim), name = NULL) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = cols$ink, "TRUE" = cols$accent), guide = "none"
    ) +
    ggplot2::scale_size_manual(values = c("FALSE" = 1.1, "TRUE" = 2.0), guide = "none") +
    ggplot2::scale_alpha_manual(values = c("FALSE" = 0.45, "TRUE" = 0.95), guide = "none") +
    # Both axes are Z, so they share one tick vector and one unit of length.
    ggplot2::scale_x_continuous(breaks = brk) +
    ggplot2::scale_y_continuous(breaks = brk) +
    ggplot2::coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE)

  multi <- has_groups(df$.outcome)
  if (multi) {
    p <- p + ggplot2::facet_wrap(~.outcome, nrow = 1)
  }
  finish_plot(
    p,
    title = if (multi) NULL else df$.outcome[[1L]],
    subtitle = paste0(
      pct, "% pointwise and simultaneous envelopes for a calibrated model; ",
      format(100 * (1 - level), trim = TRUE),
      "% of points are expected outside the pointwise band"
    ),
    xlab = "Expected Z",
    ylab = "Observed Z"
  ) +
    ggplot2::labs(
      caption = if (isTRUE(assessment$in_sample)) {
        "In-sample scores: the envelope assumes held-out data and is optimistic here."
      }
    )
}

# Expected order statistics of Z and the pointwise / simultaneous envelopes
# implied by n iid uniform centiles. Monte Carlo draws calibrate the
# equal-local-level gamma; `reps` follows a fixed work budget so the cost of
# the chart does not grow with the cohort, and the seed keeps it reproducible.
# At the default budget the realised simultaneous level is within about half a
# percentage point of `level`.
qq_envelope <- function(n, level = 0.95, reps = NULL, seed = 20240619L) {
  n <- as.integer(n)
  reps <- as.integer(reps %||% min(10000L, max(1000L, as.integer(2e6 / n))))
  expected <- stats::qnorm(stats::ppoints(n))
  if (n < 2L) {
    na <- rep(NA_real_, n)
    return(list(
      expected = expected, point_lo = na, point_hi = na,
      sim_lo = na, sim_hi = na
    ))
  }
  i <- seq_len(n)
  a <- 1 - level
  q <- function(g) stats::qnorm(stats::qbeta(g, i, n - i + 1))
  gamma_sim <- withr::with_seed(seed, {
    worst <- vapply(seq_len(reps), function(k) {
      pu <- stats::pbeta(sort(stats::runif(n)), i, n - i + 1)
      min(pu, 1 - pu)
    }, numeric(1))
    stats::quantile(worst, probs = a, names = FALSE, type = 1L)
  })
  list(
    expected = expected,
    point_lo = q(a / 2),
    point_hi = q(1 - a / 2),
    sim_lo = q(gamma_sim),
    sim_hi = q(1 - gamma_sim)
  )
}

plot_conditional <- function(assessment) {
  df <- assessment$conditional
  df$.outcome <- factor(df$.outcome, levels = unique(df$.outcome))
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$.outcome, .data$location_drift)) +
    ggplot2::geom_col(fill = referent_cols()$band, width = 0.65, na.rm = TRUE) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        ymin = pmax(.data$location_drift - 2 * .data$location_se, 0),
        ymax = .data$location_drift + 2 * .data$location_se
      ),
      colour = referent_cols()$ink, linewidth = 0.4, na.rm = TRUE
    ) +
    pretty_y()
  if (length(unique(df$covariate)) > 1L) {
    p <- p + ggplot2::facet_wrap(~covariate, nrow = 1)
  }
  finish_plot(
    p,
    title = NULL,
    subtitle = paste0("largest fitted |E[z | covariate]|, +/- 2 SE",
                      if (length(unique(df$covariate)) == 1L) paste0(" (", df$covariate[[1]], ")") else ""),
    xlab = NULL, ylab = "Location drift"
  )
}

# --- norm_scores ------------------------------------------------------------

plot_profile <- function(scores, id = NULL) {
  if (!".id" %in% names(scores) || !".outcome" %in% names(scores)) {
    cli::cli_abort("A profile needs a score table with {.field .id} and {.field .outcome} columns.")
  }
  key <- if (is.null(id)) scores$.id[[1]] else id
  one <- scores[as.character(scores$.id) == as.character(key), , drop = FALSE]
  if (!nrow(one)) {
    cli::cli_abort("No rows for subject {.val {key}}.")
  }
  one$.outcome <- factor(one$.outcome, levels = rev(unique(one$.outcome)))
  one$exceeds <- is.finite(one$z) & abs(one$z) > 2
  lim <- max(2.5, abs(one$z), na.rm = TRUE)
  p <- ggplot2::ggplot(one, ggplot2::aes(y = .data$.outcome, x = .data$z)) +
    ggplot2::geom_vline(xintercept = 0, colour = referent_cols()$ink, linewidth = 0.35) +
    ggplot2::geom_vline(
      xintercept = c(-2, 2), linetype = "dashed", colour = referent_cols()$muted, linewidth = 0.4
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data$z, yend = .data$.outcome),
      colour = referent_cols()$muted, linewidth = 0.5, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$exceeds), size = 2.6, na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = referent_cols()$band, `TRUE` = referent_cols()$accent),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(limits = c(-lim, lim), n.breaks = 7)
  n_na <- sum(!is.finite(one$z))
  finish_plot(
    p,
    title = NULL,
    subtitle = paste0(
      "subject ", key, "; dashed guides at |z| = 2",
      if (n_na) paste0("; ", n_na, " outcome", if (n_na > 1) "s" else "", " without a score") else ""
    ),
    xlab = "z",
    ylab = NULL
  ) + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

plot_heatmap <- function(scores) {
  df <- tibble::as_tibble(scores)
  df$.id <- factor(as.character(df$.id), levels = unique(as.character(df$.id)))
  df$.outcome <- factor(df$.outcome, levels = rev(unique(df$.outcome)))
  lim <- max(2, abs(df$z), na.rm = TRUE)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$.id, .data$.outcome, fill = .data$z)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_referent_diverging(limits = c(-lim, lim), name = "z", na.value = "#d9d9d9") +
    ggplot2::coord_cartesian(expand = FALSE)
  if (nlevels(df$.id) > 20L) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                            axis.ticks.x = ggplot2::element_blank())
  }
  finish_plot(p, title = NULL, xlab = "Subject", ylab = NULL) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank())
}

# --- norm_forecast ----------------------------------------------------------

plot_fan <- function(forecast, centiles) {
  centiles <- sort(unique(as.numeric(centiles)))
  times <- forecast$times
  dists <- forecast$dist
  t_nm <- forecast$time_name %||% "age"
  y_nm <- forecast$outcome
  hist <- forecast$history
  anchor <- NULL
  if (!is.null(hist) && all(c(t_nm, y_nm) %in% names(hist)) && nrow(hist)) {
    last <- hist[which.max(hist[[t_nm]]), , drop = FALSE]
    anchor <- list(t = as.numeric(last[[t_nm]]), y = as.numeric(last[[y_nm]]))
  }
  lines <- dplyr_bind(lapply(centiles, function(p) {
    x <- times
    y <- dist_quantile(dists, rep(p, length(dists)))
    if (!is.null(anchor) && anchor$t < min(times)) {
      x <- c(anchor$t, x)
      y <- c(anchor$y, y)
    }
    tibble::tibble(
      x = x, y = y, centile = p, .label = centile_label(p), .group = NA_character_
    )
  }))
  ribbons <- ribbon_from_lines(lines)
  p <- ggplot2::ggplot() + centile_layers(lines, ribbons)
  if (!is.null(hist) && all(c(t_nm, y_nm) %in% names(hist))) {
    p <- p +
      ggplot2::geom_path(
        data = hist,
        ggplot2::aes(x = .data[[t_nm]], y = .data[[y_nm]]),
        inherit.aes = FALSE,
        colour = referent_cols()$accent,
        linewidth = 0.7
      ) +
      ggplot2::geom_point(
        data = hist,
        ggplot2::aes(x = .data[[t_nm]], y = .data[[y_nm]]),
        inherit.aes = FALSE,
        colour = referent_cols()$accent,
        size = 1.8
      )
  }
  n_extra <- sum(forecast$lag_support %in% "extrapolated_lag")
  finish_plot(
    p + pretty_x() + pretty_y(),
    title = NULL,
    subtitle = if (n_extra) {
      paste0(n_extra, " forecast time", if (n_extra > 1) "s" else "",
             " beyond the reference lag range")
    } else {
      NULL
    },
    xlab = t_nm,
    ylab = y_nm
  )
}

# --- norm_transition --------------------------------------------------------

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
  df <- df[is.finite(df[[y_nm]]), , drop = FALSE]
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$.dt, .data[[y_nm]])) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = referent_cols()$muted)
  if (type != "velocity") {
    p <- p + ggplot2::geom_hline(
      yintercept = c(-2, 2), linetype = "dotted", colour = referent_cols()$muted, linewidth = 0.4
    )
  }
  if (single_level_legend(df$support)) {
    p <- p + ggplot2::geom_point(alpha = 0.7, size = 1.6, colour = referent_cols()$band)
  } else {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(colour = .data$support), alpha = 0.7, size = 1.6) +
      support_scale(df$support, "colour")
  }
  if (length(unique(df$.outcome)) > 1L) {
    p <- p + ggplot2::facet_wrap(~.outcome)
  }
  finish_plot(p + pretty_x() + pretty_y(), title = NULL, xlab = "Elapsed time", ylab = ylab)
}

# --- norm_dynamics ----------------------------------------------------------

plot_kernel <- function(dyn) {
  df <- fortify_kernel(dyn)
  missing <- unique(df$.outcome[!is.finite(df$correlation)])
  df <- df[is.finite(df$correlation), , drop = FALSE]
  df$.outcome <- factor(df$.outcome, levels = dyn$outcomes)
  multi <- !single_level_legend(df$.outcome)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$lag, .data$correlation)) +
    ggplot2::geom_hline(yintercept = 0, colour = referent_cols()$muted, linewidth = 0.3)
  if (multi) {
    p <- p + ggplot2::geom_line(ggplot2::aes(colour = .data$.outcome), linewidth = 0.7) +
      scale_colour_referent(name = "Outcome")
  } else if (nrow(df)) {
    p <- p + ggplot2::geom_line(colour = referent_cols()$band, linewidth = 0.8)
  }
  rng <- dyn$lag_range
  if (!is.null(rng) && all(is.finite(rng)) && rng[[2]] > rng[[1]]) {
    p <- p + ggplot2::annotate(
      "rect", xmin = rng[[1]], xmax = rng[[2]], ymin = -Inf, ymax = Inf,
      fill = referent_cols()$band, alpha = 0.06
    )
  }
  p <- p + ggplot2::coord_cartesian(ylim = c(0, 1)) + pretty_x() + pretty_y()
  finish_plot(
    p,
    title = NULL,
    subtitle = paste0(
      "shaded: observed lag range",
      if (length(missing)) paste0("; not identified: ", paste(missing, collapse = ", ")) else ""
    ),
    xlab = paste0("Lag (", dyn$time_name %||% "time", ")"),
    ylab = "Process correlation"
  )
}

plot_dynamics_calibration <- function(dyn, data, id_quo, time_quo) {
  if (is.null(data)) {
    cli::cli_abort(c(
      "{.arg data} is required for a dynamics calibration plot.",
      i = "Supply held-out visits with {.arg id} and {.arg time}; innovation Z is computed with {.fn norm_transition}."
    ))
  }
  if (rlang::quo_is_null(id_quo) || rlang::quo_is_missing(id_quo) ||
      rlang::quo_is_null(time_quo) || rlang::quo_is_missing(time_quo)) {
    cli::cli_abort("{.arg id} and {.arg time} are required for a dynamics calibration plot.")
  }
  tr <- norm_transition(dyn, data = data, id = !!id_quo, time = !!time_quo)
  df <- tibble::as_tibble(tr)
  df <- df[is.finite(df$innovation_z), , drop = FALSE]
  if (!nrow(df)) {
    cli::cli_abort("No finite innovation Z in {.arg data}; is the process identified?")
  }
  pieces <- lapply(split(df, df$.outcome), function(one) {
    one <- one[order(one$innovation_z), , drop = FALSE]
    one$expected <- stats::qnorm(stats::ppoints(nrow(one)))
    one
  })
  df <- dplyr_bind(pieces)
  df$.outcome <- factor(df$.outcome, levels = dyn$outcomes)
  stats_txt <- vapply(split(df$innovation_z, df$.outcome), function(z) {
    sprintf("mean %.2f, var %.2f, n = %d", mean(z), stats::var(z), length(z))
  }, character(1))
  stats_txt <- stats_txt[nzchar(stats_txt)]
  multi <- !single_level_legend(df$.outcome)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$expected, .data$innovation_z)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = referent_cols()$muted)
  if (multi) {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$.outcome), alpha = 0.5, size = 1.2) +
      scale_colour_referent(name = "Outcome")
  } else {
    p <- p + ggplot2::geom_point(colour = referent_cols()$ink, alpha = 0.5, size = 1.2)
  }
  finish_plot(
    p + ggplot2::coord_equal() + pretty_x() + pretty_y(),
    title = NULL,
    subtitle = paste0("held-out innovation Z: ", paste(
      if (multi) paste0(names(stats_txt), " ", stats_txt) else stats_txt, collapse = "; "
    )),
    xlab = "Expected N(0, 1) quantile",
    ylab = "Innovation Z"
  )
}
