#' Frozen reference bundle and model card
#'
#' The bundle omits raw training observations by default. It keeps what
#' is required for prediction, support checks, calibration, and
#' provenance. Printing a `norm_reference` shows the model card.
#'
#' @param fit A [norm_fit].
#' @param include_data If `TRUE`, store the training frame (not default).
#' @param criteria,units,missing_policy Optional provenance fields.
#' @return An object of class `norm_reference` (also a `norm_fit`).
#' @export
norm_reference <- function(fit,
                           include_data = FALSE,
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
      ),
      data = if (isTRUE(include_data)) fit$data else NULL
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
