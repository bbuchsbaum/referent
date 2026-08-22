#' Scores derived from a predictive CDF
#'
#' @param distribution A `distribution` vector (see [dist_shash]).
#' @param y Observed values, recycled against `distribution`.
#' @return A tibble of centiles, Z-scores, tails, residuals, and log
#'   densities. Tails are evaluated in log space so that `z = 8, 10, 40`
#'   remain distinct. There is no abnormality column.
#' @examples
#' as_scores(distributional::dist_normal(0, 1), c(-2, 0, 2))
#' @export
as_scores <- function(distribution, y) {
  n <- max(length(distribution), length(y))
  y <- rep_len(as.numeric(y), n)
  distribution <- vctrs::vec_recycle(distribution, n)
  u <- dist_unpack(distribution)
  tails <- scores_from_log_tails(
    dist_eval(distribution, log_tail, y, lower.tail = TRUE, unpacked = u),
    dist_eval(distribution, log_tail, y, lower.tail = FALSE, unpacked = u)
  )
  med <- dist_eval(distribution, quantile, rep(0.5, n), unpacked = u)
  tibble::tibble(
    observed = y,
    median = med,
    centile = tails$centile,
    z = tails$z,
    tail_prob = tails$tail_prob,
    tail_surprisal = tails$tail_surprisal,
    residual = y - med,
    log_density = dist_eval(distribution, log_dens, y, unpacked = u)
  )
}

#' Flag threshold exceedances
#'
#' These are threshold exceedances, not abnormalities. With many
#' well-calibrated independent outcomes, about `0.0455 * p` rows are
#' expected to exceed `|z| > 2` by chance.
#'
#' @param scores A `ref_scores` table.
#' @param threshold Absolute Z threshold.
#' @return The score table with `exceedance` and `fdr` columns, of class
#'   `ref_flags`. The attributes `expected_exceedances` and
#'   `observed_exceedances` give the chance expectation under calibration
#'   and the observed count.
#' @examples
#' sc <- as_scores(distributional::dist_normal(0, 1), c(0.2, 2.6, -3.1))
#' sc$.in_sample <- FALSE
#' ref_flag(sc)
#' @export
ref_flag <- function(scores, threshold = 2) {
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
  class(out) <- unique(c("ref_flags", class(out)))
  out
}

#' @export
print.ref_flags <- function(x, ...) {
  cli::cli_text(
    "{.cls ref_flags} {attr(x, 'observed_exceedances')} exceedance{?s} (expected {signif(attr(x, 'expected_exceedances'), 3)})"
  )
  NextMethod()
}
