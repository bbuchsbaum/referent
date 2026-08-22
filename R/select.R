#' Lexicographic out-of-sample model ladder
#'
#' Every criterion is computed from out-of-fold cross-fitted scores:
#'
#' 1. Drop failed or unstable fits.
#' 2. Drop models whose out-of-fold calibration is unacceptable (the
#'    calibration gate below).
#' 3. Rank survivors by out-of-fold mean log score; ties are broken by
#'    out-of-fold CRPS.
#' 4. Keep the candidates within one standard error of the best, where
#'    the standard error is that of the *paired* per-observation
#'    log-density difference between the candidate and the best model
#'    (column `se_log_score_paired`; 0 for the best model itself). A
#'    candidate whose CRPS is worse than the best by more than one paired
#'    standard error (`se_crps_paired`) is dropped as well.
#' 5. Require `shape_threshold` extra mean log score before accepting any
#'    spec with covariate-dependent skew or tail.
#' 6. Choose the simplest remaining model (lowest ladder level).
#'
#' @details
#' The calibration gate passes a model when, for every outcome with
#' \eqn{n} out-of-fold scores, \eqn{|\bar z| \le 0.10 + 2/\sqrt{n}},
#' \eqn{0.80 - 2\sqrt{2/n} \le \mathrm{var}(z) \le 1.25 + 2\sqrt{2/n}},
#' the 95% coverage lies in \eqn{[0.93, 0.97] \pm 2\sqrt{0.05 \cdot 0.95/n}},
#' and the observed 5% two-sided tail rate lies in
#' \eqn{0.05 \pm (0.015 + 2\sqrt{0.05 \cdot 0.95/n})}. The fixed margins
#' are practical tolerances; the \eqn{2/\sqrt{n}} terms are two standard
#' errors under perfect calibration, so small samples are not rejected
#' for noise alone. The `calibrated` column of `comparison` records the
#' gate's verdict per model.
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
#'   out-of-fold log score, CRPS, the paired standard errors used by the
#'   one-SE rule, and calibration gate values per model.
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
        mean_log_score = -Inf, crps = Inf,
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
  paired <- paired_se(fits, best, tab$model[usable])
  tab$se_log_score_paired <- unname(paired$log_density[match(tab$model, names(paired$log_density))])
  tab$se_crps_paired <- unname(paired$crps[match(tab$model, names(paired$crps))])
  survivors <- tab[usable & tab$model %in% survivors$model, , drop = FALSE]
  within <- survivors$mean_log_score >= tab$mean_log_score[tab$model == best] -
    survivors$se_log_score_paired &
    survivors$crps <= tab$crps[tab$model == best] + survivors$se_crps_paired
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

# Standard errors, per model, of the paired per-row differences (best
# minus model) in log density and CRPS, keyed by (.row, .outcome). The
# best model's own entries are 0.
paired_se <- function(fits, best, models) {
  ref <- fits[[best]]
  key_ref <- paste(ref$.row, ref$.outcome)
  se_ld <- stats::setNames(rep(0, length(models)), models)
  se_crps <- se_ld
  for (nm in setdiff(models, best)) {
    cf <- fits[[nm]]
    m <- match(paste(cf$.row, cf$.outcome), key_ref)
    se_ld[[nm]] <- se_mean(ref$log_density[m] - cf$log_density)
    se_crps[[nm]] <- se_mean(ref$crps[m] - cf$crps)
  }
  list(log_density = se_ld, crps = se_crps)
}

# Calibration gate of the ladder; see the details of norm_select().
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
