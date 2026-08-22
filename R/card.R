#' Frozen reference bundle and model card
#'
#' The bundle omits raw training observations by default. It keeps what
#' is required for prediction, support checks, calibration, and
#' provenance.
#'
#' @param fit A [norm_fit].
#' @param include_data If `TRUE`, store the training frame (not default).
#' @param criteria,units,missing_policy Optional provenance fields.
#' @return `norm_reference` or a printed model card.
#' @export
norm_reference <- function(fit,
                           include_data = FALSE,
                           criteria = NULL,
                           units = NULL,
                           missing_policy = "complete-case per outcome; predictors are never imputed") {
  models <- lapply(fit$models, function(m) {
    if (is.null(m$model)) {
      return(m)
    }
    m
  })
  structure(
    list(
      spec = fit$spec,
      outcomes = fit$outcomes,
      models = models,
      covariate_names = fit$covariate_names,
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

#' @rdname norm_reference
#' @export
norm_card <- function(fit) {
  ref <- if (inherits(fit, "norm_reference")) fit else norm_reference(fit)
  structure(ref, class = c("norm_card", class(ref)))
}

#' @export
print.norm_card <- function(x, ...) {
  cat("referent model card\n")
  cat(sprintf("family: %s\n", x$family))
  cat(sprintf("engine: %s\n", x$spec$engine))
  cat(sprintf("outcomes: %s\n", paste(x$outcomes, collapse = ", ")))
  cat(sprintf("n: %s\n", x$n))
  cat(sprintf("missing-data policy: %s\n", x$missing_policy))
  cat(sprintf("package %s, mgcv %s, R %s\n",
              x$versions$referent, x$versions$mgcv, x$versions$r))
  st <- x$statuses
  cat(sprintf("statuses: %s\n", paste(paste0(names(st), "=", st), collapse = ", ")))
  if (length(x$ranges)) {
    cat("Covariate ranges\n")
    for (nm in names(x$ranges)) {
      r <- x$ranges[[nm]]
      cat(sprintf("  %s: [%s, %s]\n", nm, signif(r$min, 4), signif(r$max, 4)))
    }
  }
  invisible(x)
}

#' @export
predict.norm_reference <- function(object, newdata, ...) {
  fit <- object
  class(fit) <- setdiff(class(fit), c("norm_reference", "norm_card"))
  if (!inherits(fit, "norm_fit")) {
    class(fit) <- c("norm_fit", class(fit))
  }
  predict.norm_fit(fit, newdata = newdata, ...)
}
