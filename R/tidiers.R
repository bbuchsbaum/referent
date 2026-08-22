#' Tidy summaries of reference fits and assessments
#'
#' `tidy()` gives one row per outcome: family, fit status, sample size,
#' total effective degrees of freedom and one `edf_<term>` column per
#' smooth term, and any engine message. `glance()` gives one row per fit
#' (or per assessment, averaging the `overall` and `marginal` tables over
#' outcomes). `augment()` returns `newdata` with `.z_<outcome>` and
#' `.centile_<outcome>` columns and the covariate `.support`; the result
#' can be passed straight to [ref_joint()].
#'
#' @param x A [ref_fit] or `ref_assessment`.
#' @param newdata Data frame of target observations.
#' @param uncertainty Passed to [predict.ref_fit()].
#' @param ... Passed to [predict.ref_fit()] by `augment()`; unused
#'   otherwise.
#' @return A tibble.
#' @examples
#' ref <- ref_simulate(120, seed = 3)
#' fit <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 5) + sex), ref, "y")
#' tidy(fit)
#' glance(fit)
#' augment(fit, ref[1:3, ], uncertainty = "conditional")
#' @name tidy.ref_fit
NULL

#' @importFrom generics tidy glance augment
#' @export
generics::tidy

#' @export
generics::glance

#' @export
generics::augment

#' @rdname tidy.ref_fit
#' @export
tidy.ref_fit <- function(x, ...) {
  rows <- lapply(x$outcomes, function(nm) {
    m <- x$models[[nm]]
    model <- m$model
    edf <- if (is.null(model)) list() else smooth_edf(model)
    tibble::tibble(
      outcome = nm,
      family = x$spec$family$name,
      status = m$status %||% "ok",
      n = if (is.null(model)) NA_integer_ else as.integer(m$n_obs),
      edf = if (is.null(model)) NA_real_ else sum(model$edf),
      !!!edf,
      message = m$message %||% NA_character_
    )
  })
  out <- vctrs::vec_rbind(!!!rows)
  out[, c(setdiff(names(out), "message"), "message")]
}

# Effective degrees of freedom per smooth term, named edf_<label>.
smooth_edf <- function(model) {
  sm <- model$smooth
  if (!length(sm)) {
    return(list())
  }
  out <- lapply(sm, function(s) sum(model$edf[s$first.para:s$last.para]))
  names(out) <- paste0("edf_", vapply(sm, `[[`, "", "label"))
  out
}

#' @rdname tidy.ref_fit
#' @export
glance.ref_fit <- function(x, ...) {
  st <- fit_statuses(x)
  tibble::tibble(
    family = x$spec$family$name,
    engine = x$spec$engine,
    n = x$n,
    n_outcomes = length(x$outcomes),
    n_ok = sum(st == "ok"),
    covariates = paste(x$covariates, collapse = ", "),
    adapted = !is.null(x$adaptation),
    calibrated = !is.null(x$calibration)
  )
}

#' @rdname tidy.ref_fit
#' @export
glance.ref_assessment <- function(x, ...) {
  col_mean <- function(tab, cols) {
    lapply(stats::setNames(cols, cols), function(cl) mean(tab[[cl]], na.rm = TRUE))
  }
  tibble::tibble(
    n = x$n,
    n_outcomes = nrow(x$overall),
    !!!col_mean(x$overall, c("mean_log_score", "crps", "rmse", "smse", "ev", "cor")),
    !!!col_mean(x$marginal, c("mean_z", "var_z", "cover_95")),
    in_sample = x$in_sample
  )
}

#' @rdname tidy.ref_fit
#' @export
augment.ref_fit <- function(x, newdata, uncertainty = c("total", "conditional"), ...) {
  uncertainty <- match.arg(uncertainty)
  newdata <- tibble::as_tibble(newdata)
  scores <- predict(x, newdata, uncertainty = uncertainty, ...)
  out <- newdata
  for (value in c("z", "centile")) {
    m <- scores_matrix(scores, value = value)$matrix
    for (nm in colnames(m)) {
      out[[paste0(".", value, "_", nm)]] <- m[, nm]
    }
  }
  out$.support <- classify_support(x$support_ref, newdata)$support
  out
}
