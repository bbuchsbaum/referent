#' Publication theme for referent plots
#'
#' A restrained `ggplot2` theme: ink axes, light horizontal grid, no
#' panel box, and a plot-aligned title. Every `autoplot()` method ends
#' here so charts share one visual language.
#'
#' @param base_size Base font size.
#' @param base_family Font family.
#' @param grid Which major grid lines to draw.
#' @return A ggplot2 theme.
#' @export
theme_referent <- function(base_size = 11,
                           base_family = "",
                           grid = c("y", "both", "none")) {
  grid <- match.arg(grid)
  major <- ggplot2::element_line(colour = "#e6e6e6", linewidth = 0.3)
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = if (identical(grid, "both")) major else ggplot2::element_blank(),
      panel.grid.major.y = if (grid %in% c("y", "both")) major else ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      axis.line = ggplot2::element_line(colour = "#222222", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(colour = "#222222", linewidth = 0.3),
      axis.title = ggplot2::element_text(colour = "#222222"),
      axis.text = ggplot2::element_text(colour = "#333333"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.05)),
      plot.subtitle = ggplot2::element_text(colour = "#555555"),
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

referent_cols <- function() {
  list(
    ink = "#1c1c1c",
    band = "#3d5a80",
    accent = "#9b2226",
    muted = "#6e7f80",
    okabe = c(
      "#0072B2", "#D55E00", "#009E73", "#CC79A7",
      "#56B4E9", "#E69F00", "#000000", "#F0E442"
    ),
    support = c(
      `in` = "#4d7c6f",
      edge = "#c48b2a",
      out = "#9b2226",
      new_group = "#6c4d8a",
      unknown = "#9a9a9a",
      missing_predictor = "#9a9a9a",
      insufficient_history = "#7a7a7a",
      unsupported_age = "#d08c60",
      extrapolated_lag = "#b07aa1",
      unidentified = "#5c5c5c"
    )
  )
}

scale_colour_referent <- function(...) {
  ggplot2::scale_colour_manual(values = referent_cols()$okabe, ...)
}

scale_fill_referent <- function(...) {
  ggplot2::scale_fill_manual(values = referent_cols()$okabe, ...)
}

scale_fill_referent_diverging <- function(...) {
  ggplot2::scale_fill_gradient2(
    low = "#0072B2",
    mid = "#f4f4f4",
    high = "#9b2226",
    midpoint = 0,
    ...
  )
}

finish_plot <- function(p, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL) {
  p +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme_referent()
}

pretty_x <- function() {
  ggplot2::scale_x_continuous(n.breaks = 6)
}

pretty_y <- function() {
  ggplot2::scale_y_continuous(n.breaks = 5)
}

centile_label <- function(p) {
  ifelse(abs(p - 0.5) < 1e-8, "Median", paste0(format(100 * p, trim = TRUE), "th"))
}

# Manual support scale restricted to the statuses actually present, in the
# canonical order; unknown statuses get a neutral grey.
support_scale <- function(present, aesthetics = "fill", ...) {
  cols <- referent_cols()$support
  present <- unique(as.character(present[!is.na(present)]))
  extra <- setdiff(present, names(cols))
  if (length(extra)) {
    cols <- c(cols, stats::setNames(rep("#7a7a7a", length(extra)), extra))
  }
  keep <- c(intersect(names(cols), present), extra)
  ggplot2::scale_discrete_manual(
    aesthetics = aesthetics,
    values = cols[keep],
    breaks = keep,
    name = "Support",
    ...
  )
}

# More than one non-missing level?
has_groups <- function(x) {
  length(unique(x[!is.na(x)])) > 1L
}
