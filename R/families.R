#' Family constructors for reference models
#'
#' These objects describe a conditional distribution family. They are not
#' fitted models. Engines use them to choose a likelihood and to construct
#' [distributional][dist_shash] vectors.
#'
#' @param link_location,link_scale,link_skew,link_tail Character link
#'   names. Defaults follow `mgcv` (`identity` location; `logb` /
#'   `logeb` scale; identity skew and tail).
#' @param min_scale Minimum scale used by the `logb` / `logeb` links.
#' @return An object of class `norm_family`.
#' @name norm_family
NULL

new_norm_family <- function(name, n_parameter, parameter_names, links,
                            extra = list()) {
  structure(
    c(
      list(
        name = name,
        n_parameter = n_parameter,
        parameter_names = parameter_names,
        links = links
      ),
      extra
    ),
    class = c(paste0("norm_family_", name), "norm_family")
  )
}

#' @rdname norm_family
#' @export
norm_gaussian <- function(link_location = "identity",
                          link_scale = "logb",
                          min_scale = 0.01) {
  new_norm_family(
    name = "gaussian",
    n_parameter = 2L,
    parameter_names = c("location", "scale"),
    links = c(location = link_location, scale = link_scale),
    extra = list(min_scale = min_scale)
  )
}

#' @rdname norm_family
#' @export
norm_shash <- function(link_location = "identity",
                       link_scale = "logeb",
                       link_skew = "identity",
                       link_tail = "identity",
                       min_scale = 1e-2) {
  new_norm_family(
    name = "shash",
    n_parameter = 4L,
    parameter_names = c("location", "scale", "skew", "tail"),
    links = c(
      location = link_location,
      scale = link_scale,
      skew = link_skew,
      tail = link_tail
    ),
    extra = list(min_scale = min_scale)
  )
}

#' @export
print.norm_family <- function(x, ...) {
  cli::cli_text("{.cls norm_family} {x$name} ({x$n_parameter} parameter{?s})")
  invisible(x)
}
