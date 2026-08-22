#' Lexicographic out-of-sample model ladder
#'
#' 1. Drop failed or unstable fits.
#' 2. Drop models with unacceptable calibration.
#' 3. Rank survivors by out-of-sample log score.
#' 4. Use CRPS as a secondary score.
#' 5. Choose the simplest model within one standard error of the best.
#' 6. Require extra evidence before covariate-dependent shape.
#'
#' @param specs A named list of [norm_spec] objects, simplest first.
#' @param data Reference data.
#' @param outcomes Outcome selection.
#' @param folds Cross-fit folds.
#' @param shape_threshold Extra log-score gain required to accept a
#'   covariate-dependent skew/tail model.
#' @param ... Passed to [norm_crossfit()].
#' @return A list with the selected spec, scores, and a comparison table.
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
  rows <- lapply(names(specs), function(nm) {
    spec <- specs[[nm]]
    cf <- tryCatch(
      norm_crossfit(spec, data = data, outcomes = outcome_names, folds = folds, ...),
      error = function(e) e
    )
    if (inherits(cf, "error")) {
      return(tibble::tibble(
        model = nm,
        level = match(nm, names(specs)),
        status = "nonconverged",
        mean_log_score = -Inf,
        se_log_score = NA_real_,
        crps = Inf,
        calibrated = FALSE,
        selected = FALSE
      ))
    }
    fit <- attr(cf, "deployment")
    val <- split_for_assess(data, cf)
    assess <- tryCatch(
      norm_assess(fit, newdata = val),
      error = function(e) NULL
    )
    log_s <- mean(cf$log_density, na.rm = TRUE)
    se_s <- se_mean(cf$log_density)
    crps <- if (is.null(assess)) Inf else mean(assess$overall$crps, na.rm = TRUE)
    cal <- if (is.null(assess)) FALSE else acceptable_calibration(assess)
    tibble::tibble(
      model = nm,
      level = match(nm, names(specs)),
      status = "ok",
      mean_log_score = log_s,
      se_log_score = se_s,
      crps = crps,
      calibrated = cal,
      selected = FALSE,
      .crossfit = list(cf)
    )
  })
  tab <- dplyr_bind(rows)
  survivors <- tab[tab$status == "ok" & tab$calibrated, , drop = FALSE]
  if (!nrow(survivors)) {
    survivors <- tab[tab$status == "ok", , drop = FALSE]
  }
  if (!nrow(survivors)) {
    cli::cli_abort("No ladder candidate produced a usable fit.")
  }
  best <- max(survivors$mean_log_score, na.rm = TRUE)
  best_se <- survivors$se_log_score[which.max(survivors$mean_log_score)] %||% 0
  within <- survivors$mean_log_score >= (best - (best_se %||% 0))
  cand <- survivors[within, , drop = FALSE]
  # extra evidence for the most flexible (last) shape model
  last_name <- names(specs)[length(specs)]
  if (last_name %in% cand$model && nrow(cand) > 1) {
    base_best <- max(cand$mean_log_score[cand$model != last_name], na.rm = TRUE)
    if (is.finite(base_best) &&
        cand$mean_log_score[cand$model == last_name] < base_best + shape_threshold) {
      cand <- cand[cand$model != last_name, , drop = FALSE]
    }
  }
  pick <- cand$model[[which.min(cand$level)]]
  tab$selected <- tab$model == pick
  structure(
    list(
      selected = specs[[pick]],
      selected_name = pick,
      comparison = tab[, c("model", "level", "status", "mean_log_score",
                           "se_log_score", "crps", "calibrated", "selected")],
      crossfit = tab$.crossfit[tab$model == pick][[1]]
    ),
    class = "norm_selection"
  )
}

split_for_assess <- function(data, scores) {
  data
}

#' @export
print.norm_selection <- function(x, ...) {
  cli::cli_text("{.cls norm_selection} selected {x$selected_name}")
  print(x$comparison)
  invisible(x)
}
