#' Flag threshold exceedances
#'
#' These are threshold exceedances, not abnormalities. With many
#' well-calibrated independent outcomes, about `0.0455 * p` rows are
#' expected to exceed `|z| > 2` by chance.
#'
#' @param scores A `norm_scores` table.
#' @param threshold Absolute Z threshold.
#' @return The score table with `exceedance` and `fdr` columns, of class
#'   `norm_flags`. The attributes `expected_exceedances` and
#'   `observed_exceedances` give the chance expectation under calibration
#'   and the observed count.
#' @examples
#' sc <- as_scores(distributional::dist_normal(0, 1), c(0.2, 2.6, -3.1))
#' sc$.in_sample <- FALSE
#' norm_flag(sc)
#' @export
norm_flag <- function(scores, threshold = 2) {
  warn_in_sample(scores, used_for = "a group comparison")
  out <- scores
  out$exceedance <- is.finite(out$z) & abs(out$z) > threshold
  out$fdr <- NA_real_
  if (any(is.finite(out$tail_prob))) {
    out$fdr <- stats::p.adjust(out$tail_prob, method = "fdr")
  }
  expected <- 2 * stats::pnorm(-threshold) * sum(is.finite(out$z))
  attr(out, "expected_exceedances") <- expected
  attr(out, "observed_exceedances") <- sum(out$exceedance, na.rm = TRUE)
  class(out) <- unique(c("norm_flags", class(out)))
  out
}

#' @rdname norm_flag
#' @details `flag()` is a deprecated alias for `norm_flag()` and will be
#'   removed in the next release.
#' @export
flag <- function(scores, threshold = 2) {
  cli::cli_warn(c(
    "{.fn flag} was deprecated in referent 0.1.0 and will be removed in the next release.",
    i = "Use {.fn norm_flag} instead."
  ), .frequency = "once", .frequency_id = "referent_flag_deprecated")
  norm_flag(scores, threshold = threshold)
}

#' @export
print.norm_flags <- function(x, ...) {
  cli::cli_text(
    "{.cls norm_flags} {attr(x, 'observed_exceedances')} exceedance{?s} (expected {signif(attr(x, 'expected_exceedances'), 3)})"
  )
  NextMethod()
}
