#' Out-of-fold reference scores
#'
#' Every reference observation receives a score from a model that did not
#' include that observation. A final deployment model is then refit on all
#' reference rows.
#'
#' @param spec A [ref_spec].
#' @param data Reference data.
#' @param outcomes Outcome selection, as in [ref_fit()].
#' @param folds Number of folds (at least 2).
#' @param strata Optional stratification column. Rows with a missing
#'   stratum form a stratum of their own.
#' @param cluster Optional cluster / subject column. Entire clusters stay
#'   in one fold.
#' @param id Optional identifier stored on scores.
#' @param uncertainty Passed to [predict.ref_fit()] when the held-out
#'   fold is scored, as in [ref_assess()]. Use `"total"` to match the
#'   predictive that [predict.ref_fit()] returns by default, which is
#'   what a map fitted on these scores by [ref_calibrate()] will be
#'   applied to.
#' @param n_draw Number of coefficient draws when `uncertainty = "total"`.
#'   The value is retained as calibration provenance.
#' @param ... Passed to [ref_fit()].
#' @return A `ref_scores` object with `.in_sample = FALSE`, `.row`
#'   indexing rows of `data`, a per-row `crps` column, and a `.fold`
#'   column, plus a `deployment` fit attribute.
#' @export
ref_crossfit <- function(spec,
                          data,
                          outcomes,
                          folds = 5,
                          strata = NULL,
                          cluster = NULL,
                          id = NULL,
                          uncertainty = c("conditional", "total"),
                          n_draw = NULL,
                          ...) {
  uncertainty <- match.arg(uncertainty)
  data <- tibble::as_tibble(data)
  outcome_names <- select_outcomes(rlang::enquo(outcomes), data)
  strata_vec <- pull_column(data, rlang::enquo(strata), default = NULL)
  cluster_vec <- pull_column(data, rlang::enquo(cluster), default = seq_len(nrow(data)))
  id_quo <- rlang::enquo(id)
  id_vec <- pull_column(data, id_quo, default = seq_len(nrow(data)))
  fold_id <- make_folds(cluster_vec, folds, strata_vec)
  pieces <- lapply(seq_len(folds), function(k) {
    train <- data[fold_id != k, , drop = FALSE]
    test <- data[fold_id == k, , drop = FALSE]
    if (!nrow(test) || !nrow(train)) {
      return(NULL)
    }
    fit <- ref_fit(spec, data = train, outcomes = outcome_names, ...)
    dists <- predict_dists(fit, test, uncertainty = uncertainty, n_draw = n_draw)
    sc <- scores_from_dists(fit, dists, test, allow_extrapolation = FALSE)
    sc$crps <- NA_real_
    for (nm in names(dists)) {
      if (!is.null(dists[[nm]]) && nm %in% names(test)) {
        sel <- sc$.outcome == nm
        sc$crps[sel] <- crps_from_dist(dists[[nm]], test[[nm]])
      }
    }
    rows <- which(fold_id == k)
    sc$.id <- id_vec[rows][sc$.row]
    sc$.row <- rows[sc$.row]
    sc$.in_sample <- FALSE
    sc$.fold <- k
    sc
  })
  scores <- dplyr_bind(Filter(Negate(is.null), pieces))
  scores <- scores[order(scores$.outcome, scores$.row), , drop = FALSE]
  deployment <- ref_fit(spec, data = data, outcomes = outcome_names, id = !!id_quo, ...)
  structure(
    scores,
    class = c("ref_scores", class(scores)),
    deployment = deployment,
    folds = fold_id,
    uncertainty = uncertainty,
    n_draw = if (identical(uncertainty, "total")) {
      as.integer(n_draw %||% deployment$spec$control$n_draw %||% 200L)
    } else {
      NA_integer_
    },
    data_hash = deployment$data_hash,
    data_hash_version = as.integer(deployment$data_hash_version %||% 3L)
  )
}

# Fold of every row. Clusters are assigned whole; within each stratum
# (an `NA` stratum is a stratum of its own) clusters are spread evenly.
make_folds <- function(cluster, folds, strata = NULL) {
  folds <- as.integer(folds)
  if (!is.finite(folds) || folds < 2L) {
    cli::cli_abort("{.arg folds} must be at least 2 (got {folds}).")
  }
  cl <- as.character(cluster)
  uniq <- unique(cl)
  if (length(uniq) < 2L) {
    cli::cli_abort(
      "Cross-fitting needs at least two clusters; {.arg cluster} has {length(uniq)}."
    )
  }
  if (length(uniq) < folds) {
    cli::cli_warn("Only {length(uniq)} cluster{?s} for {folds} folds; some folds are empty.")
  }
  if (is.null(strata)) {
    fold_of <- sample(rep(seq_len(folds), length.out = length(uniq)))
    names(fold_of) <- uniq
    return(unname(fold_of[cl]))
  }
  st <- as.character(strata)
  st[is.na(st)] <- ".na"
  fold <- integer(length(cl))
  for (s in unique(st)) {
    u <- unique(cl[st == s])
    fo <- sample(rep(seq_len(folds), length.out = length(u)))
    names(fo) <- u
    fold[cl %in% u] <- unname(fo[cl[cl %in% u]])
  }
  fold
}
