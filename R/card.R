#' Frozen reference bundle and model card
#'
#' The bundle omits raw training observations. It keeps what is required
#' for prediction, support checks, adaptation, calibration, and
#' provenance. Printing a `norm_reference` shows the model card.
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
      models = fit$models,
      covariates = fit$covariates,
      support_ref = fit$support_ref,
      calibration = fit$calibration,
      adaptation = fit$adaptation,
      id_name = fit$id_name,
      data_hash = fit$data_hash,
      statuses = fit_statuses(fit),
      n = fit$n,
      criteria = criteria,
      units = units,
      missing_policy = missing_policy,
      factor_levels = fit$support_ref$factor_levels,
      ranges = fit$support_ref$numeric,
      family = fit$spec$family$name,
      formulas = list(
        location = fit$spec$location,
        scale = fit$spec$scale,
        skew = fit$spec$skew,
        tail = fit$spec$tail
      ),
      versions = list(
        referent = as.character(utils::packageVersion("referent")),
        mgcv = as.character(utils::packageVersion("mgcv")),
        r = paste(R.version$major, R.version$minor, sep = ".")
      )
    ),
    class = c("norm_reference", "norm_fit")
  )
}

#' @export
print.norm_reference <- function(x, ...) {
  cli::cli_h1("referent model card")
  cli::cli_text("family: {x$family}")
  cli::cli_text("engine: {x$spec$engine}")
  cli::cli_text("outcomes: {.field {x$outcomes}}")
  cli::cli_text("covariates: {.field {x$covariates}}")
  cli::cli_text("n: {x$n}")
  cli::cli_text("missing-data policy: {x$missing_policy}")
  cli::cli_text(
    "package {x$versions$referent}, mgcv {x$versions$mgcv}, R {x$versions$r}"
  )
  st <- x$statuses
  cli::cli_text("statuses: {paste(paste0(names(st), '=', st), collapse = ', ')}")
  if (!is.null(x$adaptation)) {
    cli::cli_text("adapted: {paste(x$adaptation$parameters, collapse = ', ')} (local n = {x$adaptation$n_local})")
  }
  if (!is.null(x$calibration)) {
    cli::cli_text("calibrated: {x$calibration$method} on n = {x$calibration$n}")
  }
  if (length(x$ranges)) {
    cli::cli_text("Covariate ranges")
    for (nm in names(x$ranges)) {
      r <- x$ranges[[nm]]
      cli::cli_text("  {nm}: [{signif(r$min, 4)}, {signif(r$max, 4)}]")
    }
  }
  invisible(x)
}

#' @export
predict.norm_reference <- function(object, newdata, ...) {
  fit <- object
  class(fit) <- setdiff(class(fit), "norm_reference")
  if (!inherits(fit, "norm_fit")) {
    class(fit) <- c("norm_fit", class(fit))
  }
  predict.norm_fit(fit, newdata = newdata, ...)
}
