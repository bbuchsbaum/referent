#' Gaussian-copula joint deviation layer
#'
#' 1. Take out-of-fold (or transition) normal scores.
#' 2. Estimate a shrinkage correlation matrix.
#' 3. Compute \eqn{D_i^2 = z_i^\top R^{-1} z_i}.
#' 4. Calibrate \eqn{D_i} empirically on the reference subjects.
#'
#' @param scores A `norm_scores` or `norm_transition` table.
#' @param method Must be `"gaussian_copula"`.
#' @param covariance `"shrinkage"` or `"identity"`.
#' @param value Score column (`"z"` or `"innovation_z"`).
#' @return A `norm_joint` object.
#' @export
norm_joint <- function(scores,
                       method = c("gaussian_copula"),
                       covariance = c("shrinkage", "identity"),
                       value = NULL) {
  method <- match.arg(method)
  covariance <- match.arg(covariance)
  value <- value %||% if ("innovation_z" %in% names(scores)) "innovation_z" else "z"
  wide <- scores_matrix(scores, value = value)
  Z <- wide$matrix
  R <- joint_correlation(Z, covariance)
  d2 <- apply(Z, 1, function(z) {
    ok <- is.finite(z)
    if (!any(ok)) {
      return(NA_real_)
    }
    zz <- z[ok]
    rr <- R[ok, ok, drop = FALSE]
    rr <- rr + diag(1e-6, nrow(rr))
    drop(t(zz) %*% solve(rr) %*% zz)
  })
  emp <- stats::ecdf(d2[is.finite(d2)])
  tibble::tibble(
    .row = wide$.row,
    .id = wide$.id,
    d2 = d2,
    joint_centile = emp(d2),
    joint_z = stats::qnorm(clamp_prob(emp(d2))),
    n_observed = rowSums(is.finite(Z))
  ) |>
    structure(
      class = c("norm_joint", "tbl_df", "tbl", "data.frame"),
      correlation = R,
      method = method,
      value = value
    )
}

joint_correlation <- function(Z, covariance) {
  p <- ncol(Z)
  if (identical(covariance, "identity") || p < 2L) {
    return(diag(p))
  }
  s <- stats::cor(Z, use = "pairwise.complete.obs")
  s[!is.finite(s)] <- 0
  diag(s) <- 1
  n <- nrow(Z)
  lam <- p / (p + n)
  (1 - lam) * s + lam * diag(p)
}

# Pivot a long score table to a rows x outcomes matrix of `value`.
scores_matrix <- function(scores, value = "z") {
  if (!value %in% names(scores)) {
    cli::cli_abort("Unknown score column {.field {value}}.")
  }
  row_col <- if (".row" %in% names(scores)) ".row" else ".id"
  rows <- unique(scores[[row_col]])
  outs <- unique(scores$.outcome)
  mat <- matrix(NA_real_, length(rows), length(outs), dimnames = list(NULL, outs))
  mat[cbind(match(scores[[row_col]], rows), match(scores$.outcome, outs))] <-
    as.numeric(scores[[value]])
  ids <- if (".id" %in% names(scores)) {
    scores$.id[match(rows, scores[[row_col]])]
  } else {
    rows
  }
  list(matrix = mat, .row = rows, .id = ids)
}
