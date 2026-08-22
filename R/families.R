#' Family constructors for reference models
#'
#' These objects describe a conditional distribution family. They are not
#' fitted models. Engines use them to choose a likelihood and to construct
#' [distributional][dist_shash] vectors. `norm_gaussian()` models location
#' and scale; `norm_shash()` adds skew and tail (four parameters), with the
#' `mgcv` link functions (`identity` location, `logb` / `logeb` scale).
#'
#' @param min_scale Minimum scale used by the `logb` / `logeb` links.
#' @return An object of class `norm_family` with fields `name` and
#'   `min_scale`.
#' @name norm_family
NULL

#' @rdname norm_family
#' @export
norm_gaussian <- function(min_scale = 0.01) {
  structure(list(name = "gaussian", min_scale = min_scale), class = "norm_family")
}

#' @rdname norm_family
#' @export
norm_shash <- function(min_scale = 0.01) {
  structure(list(name = "shash", min_scale = min_scale), class = "norm_family")
}

#' @export
print.norm_family <- function(x, ...) {
  n_par <- if (identical(x$name, "shash")) 4L else 2L
  cli::cli_text("{.cls norm_family} {x$name} ({n_par} parameter{?s})")
  invisible(x)
}
