#' Parallel outcome fitting
#'
#' Optional `future` support. Sequential and parallel fits must be
#' statistically and numerically reproducible when seeds are set.
#'
#' @param spec A [norm_spec].
#' @param data Reference data.
#' @param outcomes Outcome names (character).
#' @param seed Base seed.
#' @param ... Passed to [norm_fit()] / the engine.
#' @keywords internal
fit_outcomes <- function(spec, data, outcomes, seed = 1L, ...) {
  if (length(outcomes) == 1L || !has_pkg("future.apply")) {
    return(norm_fit(spec, data = data, outcomes = outcomes, ...))
  }
  fits <- future.apply::future_lapply(
    outcomes,
    function(nm) {
      if (!is.null(seed)) {
        set.seed(as.integer(seed) + sum(utf8ToInt(nm)))
      }
      norm_fit(spec, data = data, outcomes = nm, ...)
    },
    future.seed = TRUE
  )
  merge_fits(fits)
}

merge_fits <- function(fits) {
  out <- fits[[1]]
  for (f in fits[-1]) {
    out$outcomes <- c(out$outcomes, f$outcomes)
    out$models <- c(out$models, f$models)
  }
  out
}
