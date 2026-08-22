#' Fit distributional reference models
#'
#' Fits one model per outcome. A failed outcome does not abort the panel:
#' a non-numeric outcome (factor, character, logical) gets the status
#' `"unsupported_type"`, one with fewer than eight finite values or no
#' spread `"insufficient_variation"`, and an engine failure
#' `"nonconverged"`.
#'
#' @param spec A [ref_spec].
#' @param data A data frame of reference observations.
#' @param outcomes Tidyselect specification or character vector of outcome
#'   columns. The covariate frame is kept wide.
#' @param id Optional subject identifier column. It is carried on scores
#'   but is never a covariate.
#' @param ... Passed to the engine fitter.
#' @details
#' Covariates are the variables named in the spec's `location`, `scale`,
#' `skew`, and `tail` formulas; other columns are ignored. Each outcome is
#' fitted on the rows that are complete for that outcome and those
#' covariates, and a message reports how many rows were dropped.
#'
#' Outcomes are fitted in parallel when the `future.apply` package is
#' installed and a non-sequential [future::plan()] is active.
#' @return An object of class `ref_fit`. `fit$covariates` holds the
#'   covariate names.
#' @examples
#' ref <- ref_simulate(80, seed = 1)
#' spec <- ref_spec(family = ref_gaussian(), location = ~ age + sex)
#' fit <- ref_fit(spec, data = ref, outcomes = "y")
#' predict(fit, newdata = ref[1:3, ], uncertainty = "conditional")
#' @export
ref_fit <- function(spec, data, outcomes, id = NULL, ...) {
  if (!inherits(spec, "ref_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls ref_spec}.")
  }
  data <- tibble::as_tibble(data)
  outcome_names <- select_outcomes(rlang::enquo(outcomes), data)
  id_quo <- rlang::enquo(id)
  id_name <- as_col_name(id_quo)
  if (!is.null(id_name) && !id_name %in% names(data)) {
    cli::cli_abort("{.arg id} column {.field {id_name}} not found in {.arg data}.")
  }
  covariate_names <- spec_covariates(spec)
  missing_cov <- setdiff(covariate_names, names(data))
  if (length(missing_cov)) {
    cli::cli_abort("Covariate{?s} {.field {missing_cov}} not found in {.arg data}.")
  }
  failed <- function(nm, status, message) {
    list(
      engine = spec$engine, family = spec$family, outcome = nm,
      status = status, message = message, model = NULL, spec = spec
    )
  }
  fit_one <- function(nm) {
    y <- data[[nm]]
    if (!is.numeric(y)) {
      return(failed(nm, "unsupported_type", paste0(
        "Outcome is ", class(y)[[1L]], "; only numeric outcomes are supported."
      )))
    }
    if (sum(is.finite(y)) < 8L || stats::sd(y, na.rm = TRUE) < .Machine$double.eps) {
      return(failed(nm, "insufficient_variation", "Outcome has insufficient variation."))
    }
    cc <- stats::complete.cases(data[, c(nm, covariate_names), drop = FALSE])
    n_drop <- sum(!cc)
    if (n_drop > 0L) {
      cli::cli_inform(
        "{.field {nm}}: dropped {n_drop} row{?s} with missing outcome or covariate values."
      )
    }
    fit_engine_mgcv(spec, data[cc, , drop = FALSE], nm, ...)
  }
  models <- outcome_lapply(outcome_names, fit_one)
  names(models) <- outcome_names
  support_ref <- support_reference(data, covariate_names)
  structure(
    list(
      spec = spec,
      outcomes = outcome_names,
      models = models,
      id_name = id_name,
      data_hash = digest_data(data[, c(covariate_names, outcome_names), drop = FALSE]),
      reference_baseline = reference_baseline(data, outcome_names),
      n = nrow(data),
      covariates = covariate_names,
      support_ref = support_ref
    ),
    class = "ref_fit"
  )
}

# Mean and population standard deviation (ddof = 0) of every outcome in
# the reference sample. `ref_assess()` scores a model against this
# unconditional Gaussian so that the baseline is fixed by the reference
# population rather than by whatever sample is being scored. Scoring a
# model against the held-out sample's own moments makes the baseline an
# oracle: it absorbs part of the signal and flatters or penalises the
# model depending on how the held-out sample happens to be spread.
# (PCNtoolkit does exactly that, which is why its MSLL is not directly
# comparable with this column.) Outcomes with fewer than two finite
# values, or no spread, get NULL and an NA standardised log score.
reference_baseline <- function(data, outcomes) {
  out <- lapply(outcomes, function(nm) {
    y <- data[[nm]]
    if (!is.numeric(y)) {
      return(NULL)
    }
    y <- y[is.finite(y)]
    if (length(y) < 2L) {
      return(NULL)
    }
    s <- sqrt(mean((y - mean(y))^2))
    if (!is.finite(s) || s <= 0) {
      return(NULL)
    }
    list(mean = mean(y), sd = s)
  })
  names(out) <- outcomes
  out
}

# Covariates are exactly the data variables named in the spec's formulas,
# as mgcv reads them: inside a smooth call (`s()`, `te()`, `ti()`,
# `t2()`) only the unnamed arguments and `by =` name data; named
# arguments such as `k = kk` or `xt = list(...)` are evaluated in the
# formula environment and are not covariates.
spec_covariates <- function(spec) {
  unique(unlist(lapply(
    spec[c("location", "scale", "skew", "tail")],
    function(f) formula_variables(f[[length(f)]])
  ), use.names = FALSE))
}

formula_variables <- function(expr) {
  if (is.name(expr)) {
    return(as.character(expr))
  }
  if (!is.call(expr)) {
    return(character())
  }
  fn <- as.character(expr[[1L]])[[1L]]
  args <- as.list(expr)[-1L]
  if (fn %in% c("s", "te", "ti", "t2")) {
    nms <- names(args) %||% rep("", length(args))
    args <- args[nms == "" | nms == "by"]
  }
  unique(unlist(lapply(args, formula_variables), use.names = FALSE))
}

# Fit outcomes in parallel when a non-sequential future plan is active.
outcome_lapply <- function(outcomes, fun) {
  parallel <- length(outcomes) > 1L &&
    has_pkg("future.apply") &&
    !inherits(future::plan(), "sequential")
  if (parallel) {
    out <- future.apply::future_lapply(outcomes, fun, future.seed = TRUE)
  } else {
    out <- lapply(outcomes, fun)
  }
  names(out) <- outcomes
  out
}

pull_column <- function(data, quo, default = NULL) {
  if (rlang::quo_is_null(quo) || rlang::quo_is_missing(quo)) {
    return(default)
  }
  nm <- tryCatch(rlang::as_name(quo), error = function(e) NULL)
  if (!is.null(nm) && nm %in% names(data)) {
    return(data[[nm]])
  }
  rlang::eval_tidy(quo, data = data)
}

digest_data <- function(data) {
  nms <- sort(names(data))
  data <- data[, nms, drop = FALSE]
  fingerprint <- paste(vapply(data, function(col) {
    if (is.numeric(col)) {
      sprintf("%.8g", sum(as.numeric(col), na.rm = TRUE))
    } else {
      paste(utils::head(as.character(col), 8L), collapse = ",")
    }
  }, character(1)), collapse = ";")
  paste(nrow(data), paste(nms, collapse = ","), fingerprint, sep = "|")
}

fit_statuses <- function(fit) {
  vapply(fit$models, function(m) m$status %||% "ok", character(1))
}

#' @export
print.ref_fit <- function(x, ...) {
  st <- fit_statuses(x)
  tab <- sort(table(st), decreasing = TRUE)
  status_txt <- paste(paste0(names(tab), "=", as.integer(tab)), collapse = ", ")
  cli::cli_text("{.cls ref_fit} {x$spec$family$name} via {x$spec$engine}")
  cli::cli_text("{length(x$outcomes)} outcome{?s}, n = {x$n}")
  cli::cli_text("covariates: {.field {x$covariates}}")
  cli::cli_text("status: {status_txt}")
  if (!is.null(x$adaptation)) {
    cli::cli_text("adapted: {paste(x$adaptation$parameters, collapse = ', ')} (local n = {x$adaptation$n_local})")
  }
  if (!is.null(x$calibration)) {
    cli::cli_text("calibrated: {x$calibration$method}{if (is.null(x$calibration$by)) '' else paste0(' by ', x$calibration$by)}")
  }
  failed <- names(st)[st != "ok"]
  if (length(failed)) {
    cli::cli_alert_warning("Failed or flagged outcome{?s}: {.field {failed}}")
  }
  invisible(x)
}

fit_ok <- function(fit_one) {
  identical(fit_one$status, "ok") && !is.null(fit_one$model)
}

#' Predict distributions or scores from a reference fit
#'
#' One method serves plain, adapted, calibrated, and frozen fits. The
#' steps are: engine prediction, then site adaptation (shifting location
#' and scale, including any coefficient draws), then scoring, then PIT
#' recalibration, and finally extrapolation masking.
#'
#' @param object A [ref_fit].
#' @param newdata Data frame of target observations.
#' @param type `"scores"` or `"distribution"`. Distributions carry
#'   adaptation but not the calibration map, which acts on probabilities.
#' @param uncertainty `"total"` (default) integrates the scores over
#'   coefficient draws (or uses the analytic Gaussian total when the
#'   location is identity-linked with constant scale); `"conditional"`
#'   uses the point estimates.
#' @param outcomes Optional subset of outcomes.
#' @param allow_extrapolation If `FALSE` (default), rows with support
#'   `"out"` or `"new_group"` have `NA` for `z`, `centile`, `tail_prob`,
#'   `tail_surprisal`, and `log_density` (`median` and `residual` are
#'   kept). An unseen level of a grouping covariate is extrapolation
#'   unless the fit has been adapted to it with [ref_adapt()]; a
#'   parametric factor has no prediction at all for such a row.
#' @param n_draw Number of coefficient draws for `"total"`; defaults to
#'   `spec$control$n_draw` or 200. Draws are seeded from
#'   `spec$control$seed` (default 1) so repeated calls agree exactly.
#'   With 200 draws the Monte Carlo error of a score is about 0.02 in
#'   `z` near the centre and larger in the far tails; it falls as
#'   \eqn{1/\sqrt{n_{draw}}}, so pass `n_draw = 2000` (or set
#'   `control = list(n_draw = 2000)` in [ref_spec()]) for reporting
#'   extreme centiles. The constant-scale Gaussian uses the exact
#'   analytic total and is unaffected.
#' @param ... Unused.
#' @return For `"scores"`, a `ref_scores` tibble with one row per
#'   observation and outcome. `status` is `"ok"`, `"missing_predictor"`
#'   (a covariate is `NA`), `"new_group"` (an unseen factor level), or
#'   the fit status of a failed outcome.
#'   For `"distribution"`, a tibble with `.id` and one
#'   [distribution][dist_shash] column per outcome (missing distributions
#'   for failed outcomes).
#' @export
predict.ref_fit <- function(object,
                             newdata,
                             type = c("scores", "distribution"),
                             uncertainty = c("total", "conditional"),
                             outcomes = NULL,
                             allow_extrapolation = FALSE,
                             n_draw = NULL,
                             ...) {
  type <- match.arg(type)
  uncertainty <- match.arg(uncertainty)
  newdata <- tibble::as_tibble(newdata)
  dists <- predict_dists(object, newdata, uncertainty, outcomes, n_draw)
  if (identical(type, "distribution")) {
    n <- nrow(newdata)
    cols <- lapply(dists, function(d) d %||% distributional::dist_missing(n))
    return(tibble::tibble(.id = score_ids(object, newdata), !!!cols))
  }
  scores_from_dists(object, dists, newdata, allow_extrapolation)
}

# Named list of per-outcome distribution vectors (NULL for failed fits),
# after site adaptation.
predict_dists <- function(object, newdata, uncertainty = "conditional",
                          outcomes = NULL, n_draw = NULL) {
  nms <- object$outcomes
  if (!is.null(outcomes)) {
    nms <- intersect(as.character(outcomes), nms)
  }
  ctrl <- object$spec$control %||% list()
  n_draw <- as.integer(n_draw %||% ctrl$n_draw %||% 200L)
  seed <- ctrl$seed %||% 1L
  dists <- lapply(nms, function(nm) {
    m <- object$models[[nm]]
    if (!fit_ok(m)) {
      return(NULL)
    }
    predict_mgcv_dist(m, newdata, uncertainty, n_draw, seed)
  })
  names(dists) <- nms
  if (!is.null(object$adaptation)) {
    dists <- apply_adaptation(object$adaptation, dists, newdata,
                              total = identical(uncertainty, "total"), seed = seed)
  }
  dists
}

# Assemble the long score table from per-outcome distributions, then
# apply the calibration map and the extrapolation mask.
scores_from_dists <- function(object, dists, newdata, allow_extrapolation = TRUE) {
  n <- nrow(newdata)
  support <- classify_support(object$support_ref, newdata)$support
  in_sample <- is_in_sample_data(object, newdata)
  ids <- score_ids(object, newdata)
  covs <- intersect(object$covariates, names(newdata))
  missing_cov <- if (length(covs)) {
    !stats::complete.cases(newdata[, covs, drop = FALSE])
  } else {
    rep(FALSE, n)
  }
  cal_group <- calibration_groups(object$calibration, newdata)
  rows <- lapply(names(dists), function(nm) {
    d <- dists[[nm]]
    y <- if (nm %in% names(newdata) && is.numeric(newdata[[nm]])) {
      newdata[[nm]]
    } else {
      rep(NA_real_, n)
    }
    if (is.null(d)) {
      sc <- as_scores(distributional::dist_missing(n), y)
      status <- rep(object$models[[nm]]$status %||% "nonconverged", n)
    } else {
      sc <- as_scores(d, y)
      status <- rep("ok", n)
      status[support %in% "new_group"] <- "new_group"
      status[missing_cov] <- "missing_predictor"
    }
    calibrated <- rep(FALSE, n)
    if (!is.null(cal_group) && !is.null(d)) {
      sc <- calibrate_scores(object$calibration, nm, cal_group, sc)
      calibrated <- sc$calibrated
      sc$calibrated <- NULL
    }
    if (!isTRUE(allow_extrapolation)) {
      mask <- support %in% c("out", "new_group")
      for (col in c("z", "centile", "tail_prob", "tail_surprisal", "log_density")) {
        sc[[col]][mask] <- NA_real_
      }
    }
    tibble::tibble(
      .row = seq_len(n),
      .id = ids,
      .outcome = nm,
      sc,
      support = support,
      calibrated = calibrated,
      .in_sample = in_sample,
      status = status
    )
  })
  out <- dplyr_bind(rows)
  if (!nrow(out)) {
    out <- tibble::tibble(
      .row = integer(), .id = ids[0], .outcome = character(),
      as_scores(distributional::dist_missing(0), numeric()),
      support = character(), calibrated = logical(),
      .in_sample = logical(), status = character()
    )
  }
  structure(out, class = c("ref_scores", class(out)))
}

score_ids <- function(object, newdata) {
  if (".id" %in% names(newdata)) {
    return(newdata[[".id"]])
  }
  id_name <- object$id_name
  if (!is.null(id_name) && id_name %in% names(newdata)) {
    return(newdata[[id_name]])
  }
  seq_len(nrow(newdata))
}

dplyr_bind <- function(xs) {
  vctrs::vec_rbind(!!!lapply(xs, tibble::as_tibble))
}

is_in_sample_data <- function(fit, newdata) {
  cols <- intersect(c(fit$covariates, fit$outcomes), names(newdata))
  if (!length(cols) || nrow(newdata) != fit$n) {
    return(FALSE)
  }
  identical(digest_data(newdata[, cols, drop = FALSE]), fit$data_hash)
}

#' @export
print.ref_scores <- function(x, ...) {
  n_in <- if (".in_sample" %in% names(x)) {
    sum(x$.in_sample, na.rm = TRUE)
  } else {
    0L
  }
  cli::cli_text("{.cls ref_scores} {nrow(x)} row{?s}")
  if (n_in > 0) {
    cli::cli_alert_warning("{n_in} in-sample score{?s} (.in_sample = TRUE)")
  }
  NextMethod()
}
