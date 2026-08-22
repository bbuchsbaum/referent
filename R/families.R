#' Family constructors for reference models
#'
#' These objects describe a conditional distribution family. They are not
#' fitted models. Engines use them to choose a likelihood and to construct
#' [distributional][dist_shash] vectors. `ref_gaussian()` models location
#' and scale; `ref_shash()` adds skew and tail (four parameters), with the
#' `mgcv` link functions (`identity` location, `logb` / `logeb` scale).
#'
#' @param min_scale Minimum scale used by the `logb` / `logeb` links.
#' @return An object of class `ref_family` with fields `name` and
#'   `min_scale`.
#' @name ref_family
NULL

#' @rdname ref_family
#' @export
ref_gaussian <- function(min_scale = 0.01) {
  structure(list(name = "gaussian", min_scale = min_scale), class = "ref_family")
}

#' @rdname ref_family
#' @export
ref_shash <- function(min_scale = 0.01) {
  structure(list(name = "shash", min_scale = min_scale), class = "ref_family")
}

#' @export
print.ref_family <- function(x, ...) {
  n_par <- if (identical(x$name, "shash")) 4L else 2L
  cli::cli_text("{.cls ref_family} {x$name} ({n_par} parameter{?s})")
  invisible(x)
}
