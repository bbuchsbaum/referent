#' Specification of a distributional reference model
#'
#' @param family A [norm_family] such as [norm_gaussian()] or [norm_shash()].
#' @param location,scale,skew,tail One-sided formulas for the additive
#'   predictors. Unused shape formulas default to `~ 1`.
#' @param engine Fitting backend. Only `"mgcv"` is supported.
#' @param method `mgcv` smoothness-selection method.
#' @param use_bam Use [mgcv::bam()] when the family supports it and
#'   `n >= bam_min_n`.
#' @param bam_min_n Minimum rows before `bam()` is considered.
#' @param control Extra engine controls (list).
#' @return An object of class `norm_spec`.
#' @examples
#' norm_spec(family = norm_gaussian(), location = ~ s(age, k = 8) + sex)
#' @export
norm_spec <- function(family = norm_gaussian(),
                      location = ~1,
                      scale = ~1,
                      skew = ~1,
                      tail = ~1,
                      engine = "mgcv",
                      method = "REML",
                      use_bam = TRUE,
                      bam_min_n = 20000L,
                      control = list()) {
  engine <- match.arg(engine)
  if (!inherits(family, "norm_family")) {
    cli::cli_abort("{.arg family} must be a {.fn norm_family} constructor result.")
  }
  location <- as_rhs_formula(location)
  scale <- as_rhs_formula(scale)
  skew <- as_rhs_formula(skew)
  tail <- as_rhs_formula(tail)
  structure(
    list(
      family = family,
      location = location,
      scale = scale,
      skew = skew,
      tail = tail,
      engine = engine,
      method = method,
      use_bam = isTRUE(use_bam),
      bam_min_n = as.integer(bam_min_n),
      control = control
    ),
    class = "norm_spec"
  )
}

as_rhs_formula <- function(x) {
  if (inherits(x, "formula")) {
    if (length(x) == 3L) {
      x[[2L]] <- NULL
    }
    return(stats::as.formula(x, env = environment(x)))
  }
  cli::cli_abort("Predictors must be formulas.")
}

formula_is_intercept_only <- function(f) {
  identical(deparse(f[[length(f)]]), "1")
}

#' @export
print.norm_spec <- function(x, ...) {
  cli::cli_text("{.cls norm_spec} {x$family$name} via {x$engine}")
  cli::cli_text("  location: {deparse(x$location)}")
  cli::cli_text("  scale:    {deparse(x$scale)}")
  if (x$family$n_parameter >= 3L) {
    cli::cli_text("  skew:     {deparse(x$skew)}")
  }
  if (x$family$n_parameter >= 4L) {
    cli::cli_text("  tail:     {deparse(x$tail)}")
  }
  invisible(x)
}
