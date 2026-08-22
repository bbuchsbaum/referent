#' Frozen reference bundle and model card
#'
#' The bundle omits raw training observations. Each fitted model is
#' reduced to what prediction needs (coefficients, their covariance, the
#' smooth constructions, terms, and factor levels): the model frame,
#' fitted values, residuals, weights, and the mgcv family object (rebuilt
#' when predicting) are dropped, so the bundle is a small fraction of the
#' fit's size. Predictions from the bundle reproduce those of the fit
#' exactly. The bundle keeps what is required for support checks,
#' adaptation, calibration, and provenance. Printing a `ref_freeze`
#' shows the model card.
#'
#' @param fit A [ref_fit].
#' @param criteria,units,missing_policy Optional provenance fields.
#' @return An object of class `ref_freeze` (also a `ref_fit`).
#' @export
ref_freeze <- function(fit,
                           criteria = NULL,
                           units = NULL,
                           missing_policy = "complete-case per outcome; predictors are never imputed") {
  structure(
    list(
      spec = strip_spec_env(fit$spec),
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
    class = c("ref_freeze", "ref_fit")
  )
}

strip_fit_one <- function(fit_one) {
  if (!is.null(fit_one$model)) {
    fit_one$model <- strip_gam(fit_one$model)
  }
  if (!is.null(fit_one$spec)) {
    fit_one$spec <- strip_spec_env(fit_one$spec)
  }
  fit_one
}

# The spec's formulas keep the environment they were created in. When that
# is a function frame or a knitr chunk rather than the global environment,
# `saveRDS()` serialises everything in it along with the bundle, so the
# formulas are given an empty child of the global environment instead,
# as `strip_gam()` does for the model's own formulas.
strip_spec_env <- function(spec) {
  env <- new.env(parent = globalenv())
  for (nm in c("location", "scale", "skew", "tail")) {
    if (inherits(spec[[nm]], "formula")) {
      environment(spec[[nm]]) <- env
    }
  }
  spec
}

#' @export
print.ref_freeze <- function(x, ...) {
  cli::cli_h1("referent model card")
  cli::cli_text("family: {x$spec$family$name}")
  cli::cli_text("engine: {x$spec$engine}")
  tr <- spec_transform(x$spec)
  if (!transform_is_identity(tr)) {
    cli::cli_text("response transform: {tr$name}")
  }
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
    cli::cli_text("calibrated: on n = {x$calibration$n}")
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
