#' Graphics for reference models
#'
#' Every plotting function returns a `ggplot` object. Individual
#' observations are overlays, not part of the fitted geometry.
#'
#' @param object A fit, assessment, score table, forecast, or transition.
#' @param type Plot kind.
#' @param outcome Outcome name for multi-outcome fits.
#' @param x Covariate mapped to the x-axis (default `age` if present).
#' @param newdata Grid or observations to overlay.
#' @param ... Unused.
#' @return A ggplot.
#' @exportS3Method ggplot2::autoplot
autoplot.norm_fit <- function(object,
                              type = c("centiles", "support", "adaptation"),
                              outcome = NULL,
                              x = NULL,
                              newdata = NULL,
                              ...) {
  type <- match.arg(type)
  outcome <- outcome %||% object$outcomes[[1]]
  if (type == "centiles") {
    return(plot_centiles(object, outcome, x, newdata))
  }
  if (type == "support") {
    return(plot_support(object, newdata))
  }
  plot_adaptation(object)
}

#' @exportS3Method ggplot2::autoplot
autoplot.norm_assessment <- function(object, type = c("calibration", "worm", "conditional"), ...) {
  type <- match.arg(type)
  if (type == "calibration") {
    return(plot_calibration(object))
  }
  if (type == "worm") {
    return(plot_worm(object))
  }
  plot_conditional(object)
}

#' @exportS3Method ggplot2::autoplot
autoplot.norm_scores <- function(object, type = c("profile", "heatmap"), ...) {
  type <- match.arg(type)
  if (type == "profile") {
    return(plot_profile(object))
  }
  plot_heatmap(object)
}

#' @exportS3Method ggplot2::autoplot
autoplot.norm_forecast <- function(object, type = c("fan"), ...) {
  plot_fan(object)
}

#' @exportS3Method ggplot2::autoplot
autoplot.norm_transition <- function(object, type = c("velocity", "innovation", "change"), ...) {
  type <- match.arg(type)
  plot_transition(object, type)
}

#' @exportS3Method ggplot2::autoplot
autoplot.norm_dynamics <- function(object, type = c("calibration"), ...) {
  ggplot2::ggplot() +
    ggplot2::labs(title = "Dynamic calibration", subtitle = "See norm_assess() on forecast innovations")
}

plot_centiles <- function(fit, outcome, x, newdata) {
  x_nm <- x %||% default_x(fit)
  grid <- centile_grid(fit, x_nm)
  dists <- predict(fit, newdata = grid, type = "distribution", uncertainty = "conditional")[[outcome]]
  bands <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  qs <- lapply(bands, function(p) dist_quantile(dists, rep(p, length(dists))))
  names(qs) <- paste0("p", bands)
  df <- tibble::as_tibble(c(list(x = grid[[x_nm]]), qs))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$p0.05, ymax = .data$p0.95), alpha = 0.12) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$p0.25, ymax = .data$p0.75), alpha = 0.18) +
    ggplot2::geom_line(ggplot2::aes(y = .data$p0.5)) +
    ggplot2::labs(x = x_nm, y = outcome, title = "Conditional reference distribution")
  if (!is.null(newdata) && all(c(x_nm, outcome) %in% names(newdata))) {
    p <- p + ggplot2::geom_point(
      data = newdata,
      ggplot2::aes(x = .data[[x_nm]], y = .data[[outcome]]),
      inherit.aes = FALSE,
      alpha = 0.35
    )
  }
  p
}

centile_grid <- function(fit, x_nm) {
  r <- fit$support_ref$numeric[[x_nm]]
  if (is.null(r)) {
    cli::cli_abort("No numeric covariate {.field {x_nm}} for a centile plot.")
  }
  grid <- tibble::tibble(!!x_nm := seq(r$min, r$max, length.out = 80))
  for (nm in setdiff(fit$covariate_names, x_nm)) {
    if (nm %in% names(fit$support_ref$numeric)) {
      grid[[nm]] <- fit$support_ref$numeric[[nm]]$mean
    } else if (nm %in% names(fit$support_ref$factor_levels)) {
      grid[[nm]] <- factor(fit$support_ref$factor_levels[[nm]][[1]],
                           levels = fit$support_ref$factor_levels[[nm]])
    }
  }
  grid
}

default_x <- function(fit) {
  if ("age" %in% fit$covariate_names) {
    return("age")
  }
  fit$support_ref$numeric_names[[1]] %||% fit$covariate_names[[1]]
}

plot_calibration <- function(assessment) {
  df <- assessment$marginal
  long <- tibble::tibble(
    outcome = rep(df$.outcome, 5),
    nominal = rep(c(0.5, 0.8, 0.9, 0.95, 0.99), each = nrow(df)),
    observed = c(df$cover_50, df$cover_80, df$cover_90, df$cover_95, df$cover_99)
  )
  ggplot2::ggplot(long, ggplot2::aes(.data$nominal, .data$observed, colour = .data$outcome)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Nominal versus observed coverage", x = "Nominal", y = "Observed")
}

plot_worm <- function(assessment) {
  sc <- assessment$scores
  sc <- sc[is.finite(sc$z), , drop = FALSE]
  sc$expected <- stats::qnorm(ppoints(nrow(sc)))
  sc <- sc[order(sc$z), , drop = FALSE]
  sc$expected <- stats::qnorm(ppoints(nrow(sc)))
  sc$worm <- sc$z - sc$expected
  ggplot2::ggplot(sc, ggplot2::aes(.data$expected, .data$worm, colour = .data$.outcome)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::labs(title = "Worm plot", x = "Expected Z", y = "Observed minus expected")
}

plot_conditional <- function(assessment) {
  ggplot2::ggplot(assessment$conditional, ggplot2::aes(.data$.outcome, .data$location_drift)) +
    ggplot2::geom_col() +
    ggplot2::labs(title = "Conditional location drift", x = NULL, y = "max |E[z|x]|")
}

plot_support <- function(fit, newdata) {
  if (is.null(newdata)) {
    cli::cli_abort("{.arg newdata} is required for a support plot.")
  }
  st <- classify_support(fit$support_ref, newdata)
  x_nm <- default_x(fit)
  df <- newdata
  df$support <- st
  ggplot2::ggplot(df, ggplot2::aes(.data[[x_nm]], fill = .data$support)) +
    ggplot2::geom_histogram(bins = 30, position = "identity", alpha = 0.6) +
    ggplot2::labs(title = "Reference support")
}

plot_profile <- function(scores) {
  one <- scores[scores$.row == scores$.row[[1]], , drop = FALSE]
  ggplot2::ggplot(one, ggplot2::aes(.data$.outcome, .data$z)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_col() +
    ggplot2::labs(title = "Individual deviation profile", y = "z")
}

plot_heatmap <- function(scores) {
  ggplot2::ggplot(scores, ggplot2::aes(.data$.id, .data$.outcome, fill = .data$z)) +
    ggplot2::geom_tile() +
    ggplot2::labs(title = "Subject-by-feature deviations")
}

plot_adaptation <- function(object) {
  if (!inherits(object, "norm_adaptation")) {
    return(ggplot2::ggplot() + ggplot2::labs(title = "No adaptation attached"))
  }
  rows <- list()
  for (nm in names(object$offsets)) {
    for (g in names(object$offsets[[nm]])) {
      off <- object$offsets[[nm]][[g]]
      rows[[length(rows) + 1]] <- tibble::tibble(
        outcome = nm, group = g,
        location = off$location, scale = off$scale
      )
    }
  }
  df <- dplyr_bind(rows)
  ggplot2::ggplot(df, ggplot2::aes(.data$group, .data$location, fill = .data$outcome)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(title = "Adapted location offsets")
}

plot_fan <- function(forecast) {
  df <- forecast$summary
  ggplot2::ggplot(df, ggplot2::aes(.data$time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = .data$median)) +
    ggplot2::labs(title = "History-conditioned forecast", y = "response")
}

plot_transition <- function(object, type) {
  y <- switch(type,
    velocity = object$observed_velocity,
    innovation = object$innovation_z,
    change = object$change_z
  )
  ggplot2::ggplot(object, ggplot2::aes(.data$.dt, y)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_point() +
    ggplot2::labs(title = paste("Transition", type), x = "elapsed time")
}

