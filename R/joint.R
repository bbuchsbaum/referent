#' Gaussian-copula joint deviation layer
#'
#' Fits a second-stage model on reference normal scores and scores new
#' subjects against it:
#'
#' 1. Take out-of-fold (or transition) normal scores of the reference.
#' 2. Estimate a Ledoit-Wolf / Schäfer-Strimmer shrinkage correlation
#'    matrix \eqn{R}.
#' 3. For a subject with observed outcomes \eqn{o}, compute
#'    \eqn{D^2 = z_o^\top R_{oo}^{-1} z_o} on the observed submatrix.
#' 4. Base PIT \eqn{p = P(\chi^2_{|o|} \le D^2)}, which already accounts
#'    for how many outcomes were observed.
#' 5. Correct that PIT through its empirical distribution on the
#'    reference (leave-one-out ranks, \eqn{(\mathrm{rank}-0.5)/n}), so
#'    the joint centile is calibrated even when the copula is imperfect.
#'
#' @param scores A `ref_scores` or `ref_transition` table of reference
#'   scores (out-of-fold for honest calibration), or the wide output of
#'   [augment()] (its `.z_<outcome>` columns).
#' @param covariance `"shrinkage"` or `"identity"`.
#' @param value Score column (`"z"` or `"innovation_z"`).
#' @return A `ref_joint` model. `predict(joint, scores)` returns one row
#'   per subject with `d2`, `n_observed`, `joint_centile`, and `joint_z`;
#'   `joint$reference` holds the reference rows scored leave-one-out.
#' @export
ref_joint <- function(scores,
                       covariance = c("shrinkage", "identity"),
                       value = NULL) {
  covariance <- match.arg(covariance)
  value <- value %||% if ("innovation_z" %in% names(scores)) "innovation_z" else "z"
  wide <- scores_matrix(scores, value = value)
  Z <- drop_degenerate_columns(wide$matrix)
  R <- joint_correlation(Z, covariance)
  base <- joint_base_pit(Z, R)
  p_ref <- base$p[is.finite(base$p)]
  # Leave-one-out mid-rank: predict() maps a new PIT to
  # (#{reference PITs below it} + 0.5) / (n + 1); for reference row i
  # against the other n - 1 rows that count is rank_i - 1, giving
  # (rank_i - 0.5) / n.
  loo <- (rank(base$p, na.last = "keep", ties.method = "average") - 0.5) / length(p_ref)
  obj <- structure(
    list(
      correlation = R,
      outcomes = colnames(Z),
      covariance = covariance,
      value = value,
      reference_pit = sort(p_ref),
      n_reference = length(p_ref)
    ),
    class = "ref_joint"
  )
  obj$reference <- joint_table(wide, base, loo)
  obj
}

#' @export
#' @rdname ref_joint
#' @param object A `ref_joint`.
#' @param newdata A score table for new subjects (same outcomes).
#' @param ... Unused.
predict.ref_joint <- function(object, newdata, ...) {
  wide <- scores_matrix(newdata, value = object$value)
  Z <- wide$matrix
  outs <- object$outcomes
  full <- matrix(NA_real_, nrow(Z), length(outs), dimnames = list(NULL, outs))
  common <- intersect(colnames(Z), outs)
  full[, common] <- Z[, common, drop = FALSE]
  base <- joint_base_pit(full, object$correlation)
  ref <- object$reference_pit
  n <- length(ref)
  centile <- (findInterval(base$p, ref) + 0.5) / (n + 1)
  centile[!is.finite(base$p)] <- NA_real_
  joint_table(wide, base, centile)
}

joint_table <- function(wide, base, centile) {
  tibble::tibble(
    .row = wide$.row,
    .id = wide$.id,
    d2 = base$d2,
    n_observed = base$n_observed,
    base_centile = base$p,
    joint_centile = centile,
    joint_z = stats::qnorm(centile)
  )
}

# Outcomes with fewer than two finite scores or no spread cannot enter a
# correlation matrix; they are dropped with a message.
drop_degenerate_columns <- function(Z) {
  spread <- apply(Z, 2, function(z) {
    z <- z[is.finite(z)]
    length(z) >= 2L && stats::sd(z) > 0
  })
  if (any(!spread)) {
    cli::cli_inform(
      "Dropping outcome{?s} {.field {colnames(Z)[!spread]}} from the joint model: constant or fewer than two finite scores."
    )
  }
  if (!any(spread)) {
    cli::cli_abort("No outcome has two or more finite, varying scores.")
  }
  Z[, spread, drop = FALSE]
}

joint_base_pit <- function(Z, R) {
  d2 <- apply(Z, 1, function(z) {
    ok <- is.finite(z)
    if (!any(ok)) {
      return(NA_real_)
    }
    zz <- z[ok]
    rr <- R[ok, ok, drop = FALSE]
    drop(t(zz) %*% solve(rr) %*% zz)
  })
  n_obs <- rowSums(is.finite(Z))
  p <- stats::pchisq(d2, df = pmax(n_obs, 1))
  p[n_obs == 0] <- NA_real_
  list(d2 = d2, n_observed = n_obs, p = p)
}

# Ledoit-Wolf / Schaefer-Strimmer shrinkage of the correlation matrix
# toward the identity: lambda* = sum Var(r_ij) / sum r_ij^2 over i != j.
joint_correlation <- function(Z, covariance) {
  p <- ncol(Z)
  if (identical(covariance, "identity") || p < 2L) {
    return(matrix(diag(p), p, p, dimnames = list(colnames(Z), colnames(Z))))
  }
  X <- scale(Z)
  X[!is.finite(X)] <- 0
  n <- nrow(X)
  r <- stats::cor(Z, use = "pairwise.complete.obs")
  r[!is.finite(r)] <- 0
  diag(r) <- 1
  # w_ij = x_i x_j per row; its mean and mean square via crossprod
  w_mean <- crossprod(X) / n
  w2 <- crossprod(X^2)
  var_r <- (n / (n - 1)^3) * (w2 - n * w_mean^2)
  off <- row(r) != col(r)
  lambda <- sum(var_r[off]) / sum(r[off]^2)
  if (!is.finite(lambda)) {
    lambda <- 1
  }
  lambda <- min(max(lambda, 0), 1)
  out <- (1 - lambda) * r
  diag(out) <- 1
  attr(out, "lambda") <- lambda
  out
}

#' @export
print.ref_joint <- function(x, ...) {
  cli::cli_text(
    "{.cls ref_joint} gaussian copula, {length(x$outcomes)} outcome{?s}, n = {x$n_reference} reference subjects"
  )
  lam <- attr(x$correlation, "lambda")
  if (!is.null(lam)) {
    cli::cli_text("shrinkage lambda = {signif(lam, 3)}")
  }
  invisible(x)
}
