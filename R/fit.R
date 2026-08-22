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
#' @param object A [norm_fit].
#' @param newdata Data frame of target observations.
#' @param type `"scores"` or `"distribution"`.
#' @param uncertainty `"total"` (default for scores) or `"conditional"`.
#' @param outcomes Optional subset of outcomes.
#' @param allow_extrapolation If `FALSE` (default), severe out-of-support
#'   rows have `z = NA`.
#' @param ... Unused.
#' @export
predict.norm_fit <- function(object,
                             newdata,
                             type = c("scores", "distribution"),
                             uncertainty = c("total", "conditional"),
                             outcomes = NULL,
                             allow_extrapolation = FALSE,
                             ...) {
  type <- match.arg(type)
  uncertainty <- match.arg(uncertainty)
  newdata <- tibble::as_tibble(newdata)
  nms <- object$outcomes
  if (!is.null(outcomes)) {
    nms <- intersect(as.character(outcomes), nms)
  }
  dists <- lapply(nms, function(nm) {
    m <- object$models[[nm]]
    if (!fit_ok(m)) {
      return(NULL)
    }
    predict_engine_dist(m, newdata, uncertainty = uncertainty)
  })
  names(dists) <- nms
  if (identical(type, "distribution")) {
    return(dists)
  }
  support <- classify_support(object$support_ref, newdata)$support
  in_sample <- is_in_sample_data(object, newdata)
  rows <- lapply(nms, function(nm) {
    d <- dists[[nm]]
    y <- if (nm %in% names(newdata)) newdata[[nm]] else rep(NA_real_, nrow(newdata))
    if (is.null(d)) {
      return(tibble::tibble(
        .row = seq_len(nrow(newdata)),
        .id = score_ids(object, newdata),
        .outcome = nm,
        observed = y,
        median = NA_real_,
        centile = NA_real_,
        z = NA_real_,
        tail_prob = NA_real_,
        tail_surprisal = NA_real_,
        residual = NA_real_,
        log_density = NA_real_,
        aleatoric_sd = NA_real_,
        epistemic_sd = NA_real_,
        support = support,
        calibrated = FALSE,
        .in_sample = in_sample,
        status = object$models[[nm]]$status %||% "nonconverged"
      ))
    }
    sc <- as_scores(d, y)
    if (identical(uncertainty, "total") && !is.null(attr(d, "location_draws"))) {
      sc$centile <- cdf_total(d, y)
      sc$z <- stats::qnorm(clamp_prob(sc$centile))
      sc$tail_prob <- 2 * pmin(sc$centile, 1 - sc$centile)
      sc$tail_surprisal <- -safe_log(sc$tail_prob)
    }
    sc$.row <- seq_len(nrow(newdata))
    sc$.id <- score_ids(object, newdata)
    sc$.outcome <- nm
    sc$support <- support
    sc$calibrated <- FALSE
    sc$.in_sample <- in_sample
    sc$status <- "ok"
    if (!isTRUE(allow_extrapolation)) {
      severe <- support %in% c("out")
      sc$z[severe] <- NA_real_
    }
    sc
  })
  out <- dplyr_bind(rows)
  out <- out[, c(
    ".row", ".id", ".outcome", "observed", "median", "centile", "z",
    "tail_prob", "tail_surprisal", "residual", "log_density",
    "aleatoric_sd", "epistemic_sd", "support", "calibrated",
    ".in_sample", "status"
  )]
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
