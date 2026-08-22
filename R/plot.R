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
#' * `"adaptation"`: location offsets (with standard errors) estimated by
#'   [norm_adapt()].
#'
#' For a `norm_assessment`, `"calibration"` plots observed against nominal
#' coverage, `"conditional"` the largest fitted drift of `z` against each
#' numeric covariate (factor levels are in the table only), and `"qq"` / `"worm"` the ordered Z scores against their
#' expected normal order statistics, raw or detrended (observed minus
#' expected), inside the envelope a calibrated model implies. Under a
#' correctly specified model the held-out centiles are iid uniform, so
#' the i-th ordered centile is Beta(i, n - i + 1): the pointwise band is
#' its normal transform in closed form, and the simultaneous band uses
#' equal local levels (Aldor-Noiman et al. 2013, Am. Stat. 67:249) with
#' the common local level calibrated by Monte Carlo so the whole ordered
#' sample stays inside with probability `level`. The corner tally names
#' both bands, since about `1 - level` of the points are expected outside
#' the pointwise band even when the model is right.
#'
#' For a `norm_dynamics`, `"kernel"` draws the fitted process correlation
#' against lag and `"calibration"` draws the same Q-Q chart of held-out
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
#' @param data Visit-level data for trajectory and dynamics-calibration
#'   plots.
#' @param id Subject identifier column.
#' @param time Time column.
#' @param centiles Probability levels for centile and fan charts.
#' @param level Coverage of the envelope on the Q-Q and worm plots.
#' @param ... Unused.
#' @return A ggplot.
#' @examples
#' ref <- norm_simulate(150, seed = 1)
#' fit <- norm_fit(norm_spec(norm_gaussian(), ~ s(age, k = 5) + sex), ref, "y")
#' autoplot(fit, type = "centiles", by = sex, newdata = ref[1:20, ])
#' sc <- predict(fit, newdata = ref[1:6, ], uncertainty = "conditional")
#' autoplot(sc, type = "heatmap")
#' @name autoplot.norm_fit
#' @exportS3Method ggplot2::autoplot
autoplot.norm_fit <- function(object,
                              type = c("centiles", "trajectories", "support", "adaptation"),
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
  data <- data %||% newdata
  switch(
    type,
    centiles = plot_centiles(object, outcome, x, newdata, by_nm, centiles),
    trajectories = plot_trajectories(
      object, data, outcome, x, rlang::enquo(id), rlang::enquo(time), by_nm, centiles
    ),
    support = plot_support(object, data, x),
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
    worm = plot_qq(object, level = level, detrend = TRUE),
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
                                   level = 0.95,
                                   ...) {
  type <- match.arg(type)
  if (identical(type, "kernel")) {
    return(plot_kernel(object))
  }
  plot_dynamics_calibration(object, data, rlang::enquo(id), rlang::enquo(time), level)
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

# Colour scale mapped to a grouping variable, with the legend dropped when
# the variable has a single level.
outcome_colour <- function(x, name = "Outcome") {
  scale_colour_referent(name = name, guide = if (has_groups(x)) "legend" else "none")
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
  single <- !has_groups(df$support)
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
  dodge <- ggplot2::position_dodge(width = 0.6)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$group, .data$location, colour = .data$outcome)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = referent_cols()$muted) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        ymin = .data$location - 2 * .data$location_se,
        ymax = .data$location + 2 * .data$location_se
      ),
      position = dodge, linewidth = 0.5, na.rm = TRUE
    ) +
    ggplot2::geom_point(position = dodge, size = 2.4) +
    outcome_colour(df$outcome)
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
  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(.data$nominal, .data$observed, colour = .data$outcome, group = .data$outcome)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = referent_cols()$muted) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    outcome_colour(long$outcome) +
    ggplot2::coord_equal(xlim = c(0.45, 1), ylim = c(0.45, 1)) +
    pretty_x() + pretty_y()
  finish_plot(p, title = NULL, xlab = "Nominal coverage", ylab = "Observed coverage")
}

# Ordered Z scores against their expected normal order statistics, raw
# (`detrend = FALSE`) or as observed minus expected (the worm), inside the
# pointwise and simultaneous envelopes of a calibrated model (see
# `qq_envelope()`). `scores` needs `.outcome` and `z`; `ylab` and `note`
# let the dynamics calibration chart reuse the geometry for innovation Z.
plot_qq <- function(assessment, level = 0.95, detrend = FALSE,
                    ylab = "Observed Z", note = NULL) {
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
  centre <- function(v, expected) if (detrend) v - expected else v
  pieces <- lapply(by_out, function(one) {
    one <- one[order(one$z), , drop = FALSE]
    env <- envs[[as.character(nrow(one))]]
    one$expected <- env$expected
    one$y <- centre(one$z, env$expected)
    one$point_lo <- centre(env$point_lo, env$expected)
    one$point_hi <- centre(env$point_hi, env$expected)
    one$sim_lo <- centre(env$sim_lo, env$expected)
    one$sim_hi <- centre(env$sim_hi, env$expected)
    one$outside <- one$y < one$sim_lo | one$y > one$sim_hi
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
  ylim <- range(c(if (!detrend) df$expected, df$y, df$sim_lo, df$sim_hi), finite = TRUE)
  ylim <- ylim + c(-1, 1) * 0.06 * diff(ylim)

  # One corner label per panel, naming the band each count belongs to.
  tally <- dplyr_bind(lapply(pieces, function(one) {
    n <- nrow(one)
    k_sim <- sum(one$outside, na.rm = TRUE)
    k_pt <- sum(one$y < one$point_lo | one$y > one$point_hi, na.rm = TRUE)
    tibble::tibble(
      .outcome = one$.outcome[[1L]],
      expected = xlim[[1L]] + 0.04 * diff(xlim),
      y = ylim[[2L]] - 0.03 * diff(ylim),
      label = sprintf(
        "outside simultaneous: %d of %d\noutside pointwise: %d (%.1f%%)",
        k_sim, n, k_pt, 100 * k_pt / n
      )
    )
  }))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$expected, y = .data$y)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$sim_lo, ymax = .data$sim_hi, fill = lab_sim),
      colour = NA
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$point_lo, ymax = .data$point_hi, fill = lab_pt),
      colour = NA
    ) +
    ggplot2::geom_abline(
      slope = if (detrend) 0 else 1, intercept = 0, linetype = "dashed",
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
    ggplot2::scale_alpha_manual(values = c("FALSE" = 0.45, "TRUE" = 0.95), guide = "none")
  if (detrend) {
    p <- p + pretty_x() + pretty_y() +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE)
  } else {
    # Both axes are Z, so they share one tick vector and one unit of length.
    brk <- pretty(ylim, n = 7L)
    p <- p + ggplot2::scale_x_continuous(breaks = brk) +
      ggplot2::scale_y_continuous(breaks = brk) +
      ggplot2::coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE)
  }

  multi <- has_groups(df$.outcome)
  if (multi) {
    p <- p + ggplot2::facet_wrap(~.outcome, nrow = 1)
  }
  finish_plot(
    p,
    title = if (multi) NULL else as.character(df$.outcome[[1L]]),
    subtitle = paste0(
      if (!is.null(note)) paste0(note, "; "),
      pct, "% pointwise and simultaneous envelopes for a calibrated model; ",
      format(100 * (1 - level), trim = TRUE),
      "% of points are expected outside the pointwise band"
    ),
    xlab = "Expected Z",
    ylab = if (detrend) paste(ylab, "minus expected") else ylab
  ) +
    ggplot2::labs(
      caption = if (isTRUE(assessment$in_sample)) {
        "In-sample scores: the envelope assumes held-out data and is optimistic here."
      }
    )
}

# Expected order statistics of Z and the pointwise / simultaneous envelopes
# implied by n iid uniform centiles. The pointwise band is the exact Beta
# order-statistic band; the simultaneous band's common local level gamma is
# calibrated by seeded Monte Carlo (2000 ordered uniform samples put the
# realised simultaneous level within about a percentage point of `level`).
qq_envelope <- function(n, level = 0.95, reps = 2000L, seed = 20240619L) {
  n <- as.integer(n)
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
  if ("level" %in% names(df)) {
    df <- df[is.na(df$level), , drop = FALSE]
  }
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
  p <- p +
    ggplot2::geom_point(ggplot2::aes(colour = .data$support), alpha = 0.7, size = 1.6) +
    support_scale(df$support, "colour", guide = if (has_groups(df$support)) "legend" else "none")
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
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$lag, .data$correlation, colour = .data$.outcome)) +
    ggplot2::geom_hline(yintercept = 0, colour = referent_cols()$muted, linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.7) +
    outcome_colour(df$.outcome)
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

plot_dynamics_calibration <- function(dyn, data, id_quo, time_quo, level = 0.95) {
  if (is.null(data)) {
    cli::cli_abort(c(
      "{.arg data} is required for a dynamics calibration plot.",
      i = "Supply held-out visits with {.arg id} and {.arg time}; innovation Z is computed with {.fn norm_transition}."
    ))
  }
  if (is.null(as_col_name(id_quo)) || is.null(as_col_name(time_quo))) {
    cli::cli_abort("{.arg id} and {.arg time} are required for a dynamics calibration plot.")
  }
  tr <- norm_transition(dyn, data = data, id = !!id_quo, time = !!time_quo)
  df <- tibble::tibble(.outcome = tr$.outcome, z = tr$innovation_z)
  df <- df[is.finite(df$z), , drop = FALSE]
  if (!nrow(df)) {
    cli::cli_abort("No finite innovation Z in {.arg data}; is the process identified?")
  }
  df$.outcome <- factor(df$.outcome, levels = dyn$outcomes)
  plot_qq(
    list(scores = df, in_sample = FALSE), level = level,
    ylab = "Innovation Z", note = "held-out innovation Z"
  )
}
