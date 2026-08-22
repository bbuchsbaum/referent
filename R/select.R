#' Lexicographic out-of-sample model ladder
#'
#' Every criterion is computed from out-of-fold cross-fitted scores:
#'
#' 1. Drop failed or unstable fits.
#' 2. Drop models whose out-of-fold calibration is unacceptable (see
#'    [acceptable_calibration()]).
#' 3. Rank survivors by out-of-fold mean log score; ties are broken by
#'    out-of-fold CRPS.
#' 4. Keep the candidates within one standard error of the best, where
#'    the standard error is that of the *paired* per-observation
#'    log-density difference against the best model. A candidate whose
#'    CRPS is worse than the best by more than one paired standard error
#'    is dropped as well.
#' 5. Require `shape_threshold` extra mean log score before accepting any
#'    spec with covariate-dependent skew or tail.
#' 6. Choose the simplest remaining model (lowest ladder level).
#'
#' @param specs A named list of [norm_spec] objects, simplest first.
#' @param data Reference data.
#' @param outcomes Outcome selection.
#' @param folds Cross-fit folds.
#' @param shape_threshold Extra log-score gain required to accept a
#'   covariate-dependent skew/tail model.
#' @param ... Passed to [norm_crossfit()].
#' @return A `norm_selection`: the selected spec and name, the selected
#'   model's cross-fit scores, and a `comparison` table with the
#'   out-of-fold log score, CRPS, and calibration gate values per model.
#' @export
norm_select <- function(specs,
                        data,
                        outcomes,
                        folds = 5,
                        shape_threshold = 0.02,
                        ...) {
  if (inherits(specs, "norm_spec")) {
    specs <- list(model = specs)
  }
  if (is.null(names(specs))) {
    names(specs) <- paste0("m", seq_along(specs))
  }
  outcome_names <- select_outcomes(rlang::enquo(outcomes), data)
  fits <- lapply(names(specs), function(nm) {
    tryCatch(
      norm_crossfit(specs[[nm]], data = data, outcomes = outcome_names, folds = folds, ...),
      error = function(e) e
    )
  })
  names(fits) <- names(specs)
  rows <- lapply(names(specs), function(nm) {
    cf <- fits[[nm]]
    level <- match(nm, names(specs))
    shape <- spec_has_covariate_shape(specs[[nm]])
    if (inherits(cf, "error") || !any(cf$status == "ok")) {
      return(tibble::tibble(
        model = nm, level = level, status = "nonconverged", shape = shape,
        mean_log_score = -Inf, se_log_score = NA_real_, crps = Inf,
        n = 0L, mean_z = NA_real_, var_z = NA_real_, cover_95 = NA_real_,
        tail_05 = NA_real_, calibrated = FALSE, selected = FALSE
      ))
    }
    marg <- assess_marginal(cf)
    tail <- assess_tail(cf)
    ok <- is.finite(cf$log_density)
    tibble::tibble(
      model = nm,
      level = level,
      status = if (all(cf$status %in% c("ok", "missing_predictor"))) "ok" else "partial",
      shape = shape,
      mean_log_score = mean(cf$log_density[ok]),
      se_log_score = se_mean(cf$log_density),
      crps = mean(cf$crps, na.rm = TRUE),
      n = sum(ok),
      mean_z = stats::weighted.mean(marg$mean_z, marg$n),
      var_z = stats::weighted.mean(marg$var_z, marg$n),
      cover_95 = stats::weighted.mean(marg$cover_95, marg$n),
      tail_05 = stats::weighted.mean(tail$observed[tail$tail_level == 0.05],
                                     tail$n[tail$tail_level == 0.05]),
      calibrated = acceptable_calibration(list(marginal = marg, tail = tail)),
      selected = FALSE
    )
  })
  tab <- dplyr_bind(rows)
  tab$shape <- as.logical(tab$shape)
  tab$calibrated <- as.logical(tab$calibrated)
  usable <- tab$status %in% c("ok", "partial")
  survivors <- tab[usable & tab$calibrated, , drop = FALSE]
  if (!nrow(survivors)) {
    cli::cli_warn("No ladder candidate passed the calibration gate; ranking all usable fits.")
    survivors <- tab[usable, , drop = FALSE]
  }
  if (!nrow(survivors)) {
    cli::cli_abort("No ladder candidate produced a usable fit.")
  }
  ord <- order(-survivors$mean_log_score, survivors$crps)
  best <- survivors$model[[ord[[1L]]]]
  paired <- paired_se(fits, best, survivors$model)
  within <- survivors$mean_log_score >= tab$mean_log_score[tab$model == best] -
    paired$log_density &
    survivors$crps <= tab$crps[tab$model == best] + paired$crps
  cand <- survivors[within, , drop = FALSE]
  if (any(cand$shape) && any(!cand$shape)) {
    base_best <- max(cand$mean_log_score[!cand$shape])
    keep <- !cand$shape | cand$mean_log_score >= base_best + shape_threshold
    cand <- cand[keep, , drop = FALSE]
  }
  pick <- cand$model[[which.min(cand$level)]]
  tab$selected <- tab$model == pick
  structure(
    list(
      selected = specs[[pick]],
      selected_name = pick,
      comparison = tab,
      crossfit = fits[[pick]]
    ),
    class = "norm_selection"
  )
}

spec_has_covariate_shape <- function(spec) {
  length(all.vars(spec$skew)) > 0L || length(all.vars(spec$tail)) > 0L
}

# Standard errors of the paired per-row differences (best minus model)
# in log density and CRPS, keyed by (.row, .outcome).
paired_se <- function(fits, best, models) {
  ref <- fits[[best]]
  key_ref <- paste(ref$.row, ref$.outcome)
  se_ld <- 0
  se_crps <- 0
  for (nm in setdiff(models, best)) {
    cf <- fits[[nm]]
    m <- match(paste(cf$.row, cf$.outcome), key_ref)
    d_ld <- ref$log_density[m] - cf$log_density
    d_crps <- ref$crps[m] - cf$crps
    se_ld <- max(se_ld, se_mean(d_ld), na.rm = TRUE)
    se_crps <- max(se_crps, se_mean(d_crps), na.rm = TRUE)
  }
  list(log_density = se_ld, crps = se_crps)
}

#' Calibration gate used by the model ladder
#'
#' A model passes when, for every outcome with \eqn{n} out-of-fold scores,
#' \eqn{|\bar z| \le 0.10 + 2/\sqrt{n}},
#' \eqn{0.80 - 2\sqrt{2/n} \le \mathrm{var}(z) \le 1.25 + 2\sqrt{2/n}},
#' the 95% coverage lies in \eqn{[0.93, 0.97] \pm 2\sqrt{0.05 \cdot 0.95/n}},
#' and the observed 5% two-sided tail rate lies in
#' \eqn{0.05 \pm (0.015 + 2\sqrt{0.05 \cdot 0.95/n})}. The fixed margins
#' are practical tolerances; the \eqn{2/\sqrt{n}} terms are two standard
#' errors under perfect calibration, so small samples are not rejected
#' for noise alone.
#'
#' @param assessment A [norm_assess()] result or a list with `marginal`
#'   and `tail` tables.
#' @return `TRUE` or `FALSE`.
#' @export
acceptable_calibration <- function(assessment) {
  m <- assessment$marginal
  if (!nrow(m)) {
    return(FALSE)
  }
  n <- pmax(m$n, 1)
  se_mean <- 1 / sqrt(n)
  se_var <- sqrt(2 / n)
  se_cov <- sqrt(0.05 * 0.95 / n)
  ok <- is.finite(m$mean_z) & abs(m$mean_z) <= 0.10 + 2 * se_mean &
    is.finite(m$var_z) & m$var_z >= 0.80 - 2 * se_var & m$var_z <= 1.25 + 2 * se_var &
    is.finite(m$cover_95) & m$cover_95 >= 0.93 - 2 * se_cov & m$cover_95 <= 0.97 + 2 * se_cov
  t5 <- assessment$tail
  if (!is.null(t5)) {
    t5 <- t5[t5$tail_level == 0.05, , drop = FALSE]
    t5 <- t5[match(m$.outcome, t5$.outcome), , drop = FALSE]
    tol <- 0.015 + 2 * sqrt(0.05 * 0.95 / pmax(t5$n, 1))
    ok <- ok & (is.na(t5$observed) | abs(t5$observed - 0.05) <= tol)
  }
  all(ok)
}

#' @export
print.norm_selection <- function(x, ...) {
  cli::cli_text("{.cls norm_selection} selected {x$selected_name}")
  print(x$comparison[, c("model", "level", "status", "mean_log_score", "crps",
                         "mean_z", "var_z", "cover_95", "calibrated", "selected")])
  invisible(x)
}
