#' Fit distributional reference models
#'
#' Fits one model per outcome. A failed outcome does not abort the panel.
#'
#' @param spec A [norm_spec].
#' @param data A data frame of reference observations.
#' @param outcomes Tidyselect specification or character vector of outcome
#'   columns. The covariate frame is kept wide.
#' @param id Optional subject identifier (column name or vector). It is
#'   carried on scores but is never a covariate.
#' @param ... Passed to the engine fitter.
#' @details
#' Covariates are the variables named in the spec's `location`, `scale`,
#' `skew`, and `tail` formulas; other columns are ignored. Each outcome is
#' fitted on the rows that are complete for that outcome and those
#' covariates, and a message reports how many rows were dropped.
#'
#' Outcomes are fitted in parallel when the `future.apply` package is
#' installed and a non-sequential [future::plan()] is active.
#' @return An object of class `norm_fit`. `fit$covariates` holds the
#'   covariate names.
#' @examples
#' ref <- norm_simulate(80, seed = 1)
#' spec <- norm_spec(family = norm_gaussian(), location = ~ age + sex)
#' fit <- norm_fit(spec, data = ref, outcomes = "y")
#' predict(fit, newdata = ref[1:3, ], uncertainty = "conditional")
#' @export
norm_fit <- function(spec, data, outcomes, id = NULL, ...) {
  if (!inherits(spec, "norm_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls norm_spec}.")
  }
  data <- tibble::as_tibble(data)
  outcome_names <- select_outcomes(rlang::enquo(outcomes), data)
  id_quo <- rlang::enquo(id)
  ids <- pull_column(data, id_quo, default = seq_len(nrow(data)))
  id_name <- tryCatch({
    nm <- rlang::as_name(id_quo)
    if (nm %in% names(data)) nm else NULL
  }, error = function(e) NULL)
  covariate_names <- spec_covariates(spec)
  missing_cov <- setdiff(covariate_names, names(data))
  if (length(missing_cov)) {
    cli::cli_abort("Covariate{?s} {.field {missing_cov}} not found in {.arg data}.")
  }
  fit_one <- function(nm) {
    y <- data[[nm]]
    if (stats::sd(y, na.rm = TRUE) < .Machine$double.eps ||
        sum(is.finite(y)) < 8L) {
      return(list(
        engine = spec$engine,
        family = spec$family,
        outcome = nm,
        status = "insufficient_variation",
        message = "Outcome has insufficient variation.",
        model = NULL,
        spec = spec
      ))
    }
    cc <- stats::complete.cases(data[, c(nm, covariate_names), drop = FALSE])
    n_drop <- sum(!cc)
    if (n_drop > 0L) {
      cli::cli_inform(
        "{.field {nm}}: dropped {n_drop} row{?s} with missing outcome or covariate values."
      )
    }
    fit_engine(spec, data[cc, , drop = FALSE], nm, ...)
  }
  models <- outcome_lapply(outcome_names, fit_one)
  names(models) <- outcome_names
  support_ref <- support_reference(data, covariate_names)
  structure(
    list(
      spec = spec,
      outcomes = outcome_names,
      models = models,
      id = ids,
      id_name = id_name,
      data_hash = digest_data(data[, c(covariate_names, outcome_names), drop = FALSE]),
      n = nrow(data),
      covariates = covariate_names,
      support_ref = support_ref,
      in_sample_rows = nrow(data),
      package_version = as.character(utils::packageVersion("referent"))
    ),
    class = "norm_fit"
  )
}

# Covariates are exactly the variables named in the spec's formulas.
spec_covariates <- function(spec) {
  unique(unlist(lapply(
    spec[c("location", "scale", "skew", "tail")],
    all.vars
  ), use.names = FALSE))
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
print.norm_fit <- function(x, ...) {
  st <- fit_statuses(x)
  tab <- sort(table(st), decreasing = TRUE)
  status_txt <- paste(paste0(names(tab), "=", as.integer(tab)), collapse = ", ")
  cli::cli_text("{.cls norm_fit} {x$spec$family$name} via {x$spec$engine}")
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

#' @export
summary.norm_fit <- function(object, ...) {
  st <- fit_statuses(object)
  tibble::tibble(
    outcome = object$outcomes,
    status = unname(st[object$outcomes]),
    family = object$spec$family$name,
    engine = object$spec$engine
  )
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
#' @param object A [norm_fit].
#' @param newdata Data frame of target observations.
#' @param type `"scores"` or `"distribution"`. Distributions carry
#'   adaptation but not the calibration map, which acts on probabilities.
#' @param uncertainty `"total"` (default) integrates the scores over
#'   coefficient draws (or uses the analytic Gaussian total when the
#'   location is identity-linked with constant scale); `"conditional"`
#'   uses the point estimates.
#' @param outcomes Optional subset of outcomes.
#' @param allow_extrapolation If `FALSE` (default), severe out-of-support
#'   rows have `z = NA`.
#' @param n_draw Number of coefficient draws for `"total"`; defaults to
#'   `spec$control$n_draw` or 200. Draws are seeded from
#'   `spec$control$seed` (default 1) so repeated calls agree exactly.
#' @param ... Unused.
#' @return For `"scores"`, a `norm_scores` tibble with one row per
#'   observation and outcome. `status` is `"ok"`, `"missing_predictor"`
#'   (a covariate is `NA`), or the fit status of a failed outcome.
#' @export
predict.norm_fit <- function(object,
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
    predict_engine_dist(m, newdata, uncertainty = uncertainty,
                        n_draw = n_draw, seed = seed)
  })
  names(dists) <- nms
  if (!is.null(object$adaptation)) {
    dists <- apply_adaptation(object$adaptation, dists, newdata)
  }
  if (identical(type, "distribution")) {
    return(dists)
  }
  scores <- scores_from_dists(object, dists, newdata)
  if (!is.null(object$calibration)) {
    scores <- apply_calibration(object, scores, newdata)
  }
  if (!isTRUE(allow_extrapolation)) {
    scores$z[scores$support %in% "out"] <- NA_real_
  }
  scores
}

# Assemble the long score table from per-outcome distributions.
scores_from_dists <- function(object, dists, newdata) {
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
  rows <- lapply(names(dists), function(nm) {
    d <- dists[[nm]]
    y <- if (nm %in% names(newdata)) newdata[[nm]] else rep(NA_real_, n)
    if (is.null(d)) {
      sc <- tibble::tibble(
        observed = as.numeric(y), median = NA_real_, centile = NA_real_,
        z = NA_real_, tail_prob = NA_real_, tail_surprisal = NA_real_,
        residual = NA_real_, log_density = NA_real_,
        aleatoric_sd = NA_real_, epistemic_sd = NA_real_
      )
      status <- rep(object$models[[nm]]$status %||% "nonconverged", n)
    } else {
      sc <- as_scores(d, y)
      status <- ifelse(missing_cov, "missing_predictor", "ok")
    }
    tibble::tibble(
      .row = seq_len(n),
      .id = ids,
      .outcome = nm,
      sc,
      support = support,
      calibrated = FALSE,
      .in_sample = in_sample,
      status = status
    )
  })
  out <- dplyr_bind(rows)
  if (!nrow(out)) {
    out <- tibble::tibble(
      .row = integer(), .id = ids[0], .outcome = character(),
      observed = numeric(), median = numeric(), centile = numeric(),
      z = numeric(), tail_prob = numeric(), tail_surprisal = numeric(),
      residual = numeric(), log_density = numeric(),
      aleatoric_sd = numeric(), epistemic_sd = numeric(),
      support = character(), calibrated = logical(),
      .in_sample = logical(), status = character()
    )
  }
  structure(out, class = c("norm_scores", class(out)), in_sample = in_sample)
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
  tibble::as_tibble(do.call(rbind, lapply(xs, as.data.frame)))
}

is_in_sample_data <- function(fit, newdata) {
  cols <- intersect(c(fit$covariates, fit$outcomes), names(newdata))
  if (!length(cols) || nrow(newdata) != fit$n) {
    return(FALSE)
  }
  identical(digest_data(newdata[, cols, drop = FALSE]), fit$data_hash)
}

#' @export
`[.norm_scores` <- function(x, i, j, drop = FALSE) {
  cl <- class(x)
  out <- NextMethod("[")
  keep <- is.data.frame(out) && all(c("z", ".in_sample") %in% names(out))
  if (keep) {
    class(out) <- cl
    attr(out, "in_sample") <- attr(x, "in_sample")
  }
  out
}

#' @export
print.norm_scores <- function(x, ...) {
  n_in <- if (".in_sample" %in% names(x)) {
    sum(x$.in_sample, na.rm = TRUE)
  } else {
    0L
  }
  cli::cli_text("{.cls norm_scores} {nrow(x)} row{?s}")
  if (n_in > 0) {
    cli::cli_alert_warning("{n_in} in-sample score{?s} (.in_sample = TRUE)")
  }
  NextMethod()
}
