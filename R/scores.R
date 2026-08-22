#' Flag threshold exceedances
#'
#' These are threshold exceedances, not abnormalities. With many
#' well-calibrated independent outcomes, about `0.0455 * p` rows are
#' expected to exceed `|z| > 2` by chance.
#'
#' @param scores A `norm_scores` table.
#' @param threshold Absolute Z threshold.
#' @return The score table with `exceedance` and `fdr` columns, of class
#'   `norm_flags`.
#' @export
flag <- function(scores, threshold = 2) {
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

#' @export
print.norm_flags <- function(x, ...) {
  cli::cli_text(
    "{.cls norm_flags} {attr(x, 'observed_exceedances')} exceedance{?s} (expected {signif(attr(x, 'expected_exceedances'), 3)})"
  )
  NextMethod()
}
