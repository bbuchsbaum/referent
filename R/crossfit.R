#' Out-of-fold reference scores
#'
#' Every reference observation receives a score from a model that did not
#' include that observation. A final deployment model is then refit on all
#' reference rows.
#'
#' @param spec A [norm_spec].
#' @param data Reference data.
#' @param outcomes Outcome selection, as in [norm_fit()].
#' @param folds Number of folds.
#' @param strata Optional stratification column.
#' @param cluster Optional cluster / subject column. Entire clusters stay
#'   in one fold.
#' @param id Optional identifier stored on scores.
#' @param ... Passed to [norm_fit()].
#' @return A `norm_scores` object with `.in_sample = FALSE`, `.row`
#'   indexing rows of `data`, a per-row `crps` column, and a `.fold`
#'   column, plus a `deployment` fit attribute.
#' @export
norm_crossfit <- function(spec,
                          data,
                          outcomes,
                          folds = 5,
                          strata = NULL,
                          cluster = NULL,
                          id = NULL,
                          ...) {
  data <- tibble::as_tibble(data)
  outcome_names <- select_outcomes(rlang::enquo(outcomes), data)
  strata_vec <- pull_column(data, rlang::enquo(strata), default = NULL)
  cluster_vec <- pull_column(data, rlang::enquo(cluster), default = seq_len(nrow(data)))
  id_vec <- pull_column(data, rlang::enquo(id), default = seq_len(nrow(data)))
  fold_id <- make_folds(cluster_vec, folds, strata_vec)
  pieces <- lapply(seq_len(folds), function(k) {
    train <- data[fold_id != k, , drop = FALSE]
    test <- data[fold_id == k, , drop = FALSE]
    if (!nrow(test) || !nrow(train)) {
      return(NULL)
    }
    fit <- norm_fit(spec, data = train, outcomes = outcome_names, ...)
    dists <- predict_dists(fit, test, uncertainty = "conditional")
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
  deployment <- norm_fit(spec, data = data, outcomes = outcome_names, id = id_vec, ...)
  structure(
    scores,
    class = c("norm_scores", class(scores)),
    in_sample = FALSE,
    deployment = deployment,
    folds = fold_id
  )
}

make_folds <- function(cluster, folds, strata = NULL) {
  folds <- as.integer(folds)
  cl <- as.character(cluster)
  uniq <- unique(cl)
  if (is.null(strata)) {
    fold_of <- sample(rep(seq_len(folds), length.out = length(uniq)))
    names(fold_of) <- uniq
    return(unname(fold_of[cl]))
  }
  st <- as.character(strata)
  fold <- integer(length(cl))
  for (s in unique(st)) {
    u <- unique(cl[st == s])
    fo <- sample(rep(seq_len(folds), length.out = length(u)))
    names(fo) <- u
    fold[cl %in% u] <- unname(fo[cl[cl %in% u]])
  }
  fold
}
