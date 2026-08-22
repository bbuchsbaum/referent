#' Specification of a distributional reference model
#'
#' @param family A [ref_family] such as [ref_gaussian()] or [ref_shash()].
#' @param location,scale,skew,tail One-sided formulas for the additive
#'   predictors. Unused shape formulas default to `~ 1`.
#' @param transform Monotone response transform applied before fitting:
#'   `"identity"` (the default), `"sqrt"`, `"log"`, or a Box-Cox power
#'   \eqn{\lambda \ge 0}. The model is fitted to
#'   \eqn{h_\lambda(y)=(y^\lambda-1)/\lambda} (\eqn{\log y} at
#'   \eqn{\lambda=0}) and every predictive distribution is pushed back to
#'   the response scale, so densities carry the log Jacobian and
#'   quantiles, medians, `hilo()`, and CRPS are in the units of `y`. Any
#'   transform other than the identity requires a positive outcome.
#' @param method `mgcv` smoothness-selection method.
#' @param use_bam Use [mgcv::bam()] when the family supports it and
#'   `n >= bam_min_n`.
#' @param bam_min_n Minimum rows before `bam()` is considered.
#' @param control Extra engine controls (list).
#' @return An object of class `ref_spec`.
#' @examples
#' ref_spec(family = ref_gaussian(), location = ~ s(age, k = 8) + sex)
#' ref_spec(family = ref_shash(), location = ~ s(age, k = 8), transform = "log")
#' @export
ref_spec <- function(family = ref_gaussian(),
                      location = ~1,
                      scale = ~1,
                      skew = ~1,
                      tail = ~1,
                      transform = "identity",
                      method = "REML",
                      use_bam = TRUE,
                      bam_min_n = 20000L,
                      control = list()) {
  if (!inherits(family, "ref_family")) {
    cli::cli_abort("{.arg family} must be a {.fn ref_family} constructor result.")
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
      transform = as_transform(transform),
      engine = "mgcv",
      method = method,
      use_bam = isTRUE(use_bam),
      bam_min_n = as.integer(bam_min_n),
      control = control
    ),
    class = "ref_spec"
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

# --- Response transform ---------------------------------------------------
# One monotone Box-Cox map h_lambda of the outcome, named for the three
# powers that have names. The model is fitted to h(y); the predictive
# distribution is pushed back to the response scale by `dist_warped()`,
# which is where the log Jacobian and the back-transformed quantiles live.

as_transform <- function(x) {
  if (inherits(x, "ref_transform")) {
    return(x)
  }
  bad <- function() {
    cli::cli_abort(
      '{.arg transform} must be "identity", "sqrt", "log", or a Box-Cox power >= 0.'
    )
  }
  lambda <- if (is.character(x) && length(x) == 1L) {
    switch(x, identity = 1, sqrt = 0.5, log = 0, bad())
  } else if (is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0) {
    as.numeric(x)
  } else {
    bad()
  }
  structure(list(name = transform_name(lambda), lambda = lambda),
            class = "ref_transform")
}

transform_name <- function(lambda) {
  if (lambda == 1) {
    "identity"
  } else if (lambda == 0) {
    "log"
  } else if (lambda == 0.5) {
    "sqrt"
  } else {
    paste0("power ", format(lambda))
  }
}

# The transform of a spec, tolerating specs built before the argument.
spec_transform <- function(spec) {
  as_transform(spec$transform %||% "identity")
}

transform_is_identity <- function(tr) {
  identical(tr$lambda, 1)
}

# h(y). Values outside the domain (y <= 0 for every non-identity power)
# become NA, so `ref_fit()` drops them with the other incomplete rows.
transform_apply <- function(tr, y) {
  if (transform_is_identity(tr)) {
    return(y)
  }
  y[is.finite(y) & y <= 0] <- NA_real_
  if (tr$lambda == 0) log(y) else (y^tr$lambda - 1) / tr$lambda
}

# h^-1(t), clamped at the lower end of the image.
transform_invert <- function(tr, t) {
  if (transform_is_identity(tr)) {
    return(t)
  }
  if (tr$lambda == 0) {
    return(exp(t))
  }
  out <- rep(NA_real_, length(t))
  low <- is.finite(t) & t <= transform_t_min(tr)
  ok <- is.finite(t) & !low
  out[ok] <- (1 + tr$lambda * t[ok])^(1 / tr$lambda)
  out[low] <- 0
  out[t == Inf] <- Inf
  out
}

# log h'(y).
transform_log_deriv <- function(tr, y) {
  if (transform_is_identity(tr)) {
    return(rep(0, length(y)))
  }
  (tr$lambda - 1) * log(y)
}

# h of the lower end of the domain. Finite only when the image of h is
# bounded below, which is what makes the renormalisation in `dist_warped()`
# necessary: the fitted distribution puts mass below h(0) that no response
# value can produce.
transform_t_min <- function(tr) {
  if (transform_is_identity(tr) || tr$lambda == 0) -Inf else -1 / tr$lambda
}

#' @export
print.ref_spec <- function(x, ...) {
  cli::cli_text("{.cls ref_spec} {x$family$name} via {x$engine}")
  tr <- spec_transform(x)
  if (!transform_is_identity(tr)) {
    cli::cli_text("  transform: {tr$name}")
  }
  cli::cli_text("  location: {deparse(x$location)}")
  cli::cli_text("  scale:    {deparse(x$scale)}")
  if (identical(x$family$name, "shash")) {
    cli::cli_text("  skew:     {deparse(x$skew)}")
    cli::cli_text("  tail:     {deparse(x$tail)}")
  }
  invisible(x)
}
