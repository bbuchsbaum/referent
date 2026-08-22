#' Frozen reference bundle and model card
#'
#' The bundle omits raw training observations. Each fitted model is
#' reduced to what prediction needs (coefficients, their covariance, the
#' smooth constructions, terms, and factor levels): the model frame,
#' fitted values, residuals, weights, and the mgcv family object (rebuilt
#' when predicting) are dropped, so the bundle is a small fraction of the
#' fit's size. Predictions from the bundle reproduce those of the fit
#' exactly. The bundle keeps what is required for support checks,
#' adaptation, calibration, and provenance. Printing a `norm_reference`
#' shows the model card.
#'
#' @param fit A [norm_fit].
#' @param criteria,units,missing_policy Optional provenance fields.
#' @return An object of class `norm_reference` (also a `norm_fit`).
#' @export
norm_reference <- function(fit,
                           criteria = NULL,
                           units = NULL,
                           missing_policy = "complete-case per outcome; predictors are never imputed") {
  structure(
    list(
      spec = fit$spec,
      outcomes = fit$outcomes,
      models = lapply(fit$models, strip_fit_one),
      covariates = fit$covariates,
      support_ref = fit$support_ref,
      calibration = fit$calibration,
      adaptation = fit$adaptation,
      id_name = fit$id_name,
      data_hash = fit$data_hash,
      reference_baseline = fit$reference_baseline,
      n = fit$n,
      criteria = criteria,
      units = units,
      missing_policy = missing_policy,
      versions = list(
        referent = as.character(utils::packageVersion("referent")),
        mgcv = as.character(utils::packageVersion("mgcv")),
        r = paste(R.version$major, R.version$minor, sep = ".")
      )
    ),
    class = c("norm_reference", "norm_fit")
  )
}

strip_fit_one <- function(fit_one) {
  if (!is.null(fit_one$model)) {
    fit_one$model <- strip_gam(fit_one$model)
  }
  fit_one
}

#' @export
print.norm_reference <- function(x, ...) {
  cli::cli_h1("referent model card")
  cli::cli_text("family: {x$spec$family$name}")
  cli::cli_text("engine: {x$spec$engine}")
  cli::cli_text("outcomes: {.field {x$outcomes}}")
  cli::cli_text("covariates: {.field {x$covariates}}")
  cli::cli_text("n: {x$n}")
  cli::cli_text("missing-data policy: {x$missing_policy}")
  cli::cli_text(
    "package {x$versions$referent}, mgcv {x$versions$mgcv}, R {x$versions$r}"
  )
  st <- fit_statuses(x)
  cli::cli_text("statuses: {paste(paste0(names(st), '=', st), collapse = ', ')}")
  if (!is.null(x$adaptation)) {
    cli::cli_text("adapted: {paste(x$adaptation$parameters, collapse = ', ')} (local n = {x$adaptation$n_local})")
  }
  if (!is.null(x$calibration)) {
    cli::cli_text("calibrated: {x$calibration$method} on n = {x$calibration$n}")
  }
  ranges <- x$support_ref$numeric
  if (length(ranges)) {
    cli::cli_text("Covariate ranges")
    for (nm in names(ranges)) {
      r <- ranges[[nm]]
      cli::cli_text("  {nm}: [{signif(r$min, 4)}, {signif(r$max, 4)}]")
    }
  }
  invisible(x)
}
