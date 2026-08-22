#' Long-form score tables
#'
#' `predict(..., type = "scores")` returns a long table of class
#' `norm_scores` with one row per subject–outcome pair.
#'
#' @name norm_scores
#' @param scores A `norm_scores` table.
#' @param value Column to spread across outcomes.
#' @param .data Data to filter; typically a `norm_scores` table.
#' @param ... Expressions evaluated in the score table, as in
#'   `abs(z) > 2`.
#' @param threshold Absolute Z threshold for [flag()].
#' @return A wide tibble, a filtered score table, or threshold exceedances.
#' @export
as_wide <- function(scores, value = "z") {
  if (!value %in% names(scores)) {
    cli::cli_abort("Unknown score column {.field {value}}.")
  }
  row_col <- if (".row" %in% names(scores)) ".row" else ".id"
  ids <- unique(scores[[row_col]])
  outs <- unique(scores$.outcome)
  mat <- matrix(NA, length(ids), length(outs),
                dimnames = list(NULL, outs))
  key <- match(scores[[row_col]], ids)
  col <- match(scores$.outcome, outs)
  mat[cbind(key, col)] <- scores[[value]]
  out <- tibble::as_tibble(mat)
  out$.row <- ids
  if (".id" %in% names(scores)) {
    key_id <- if (row_col == ".row") ".row" else ".id"
    map <- scores[!duplicated(scores[[key_id]]), c(key_id, ".id"), drop = FALSE]
    out$.id <- map$.id[match(out$.row, map[[key_id]])]
  }
  out
}

#' @rdname norm_scores
#' @export
filter_scores <- function(.data, ...) {
  warn_in_sample(.data, used_for = "a downstream comparison")
  dots <- rlang::enquos(...)
  keep <- Reduce(`&`, lapply(dots, function(q) as.logical(rlang::eval_tidy(q, data = .data))))
  out <- tibble::as_tibble(.data)[keep, , drop = FALSE]
  if (inherits(.data, "norm_scores")) {
    class(out) <- unique(c("norm_scores", class(out)))
    attr(out, "in_sample") <- attr(.data, "in_sample")
  }
  out
}

#' Flag threshold exceedances
#'
#' These are threshold exceedances, not abnormalities. With many
#' well-calibrated independent outcomes, about `0.0455 * p` rows are
#' expected to exceed `|z| > 2` by chance.
#'
#' @rdname norm_scores
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
