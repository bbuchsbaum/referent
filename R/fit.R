#' Fit distributional reference models
#'
#' Fits one model per outcome. A failed outcome does not abort the panel.
#'
#' @param spec A [norm_spec].
#' @param data A data frame of reference observations.
#' @param outcomes Tidyselect specification or character vector of outcome
#'   columns. The covariate frame is kept wide.
#' @param id Optional subject identifier (column name or vector).
#' @param ... Passed to the engine fitter.
#' @return An object of class `norm_fit`.
#' @examples
#' ref <- norm_simulate(80, kind = "gaussian_location", seed = 1)
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
  ids <- pull_column(data, rlang::enquo(id), default = seq_len(nrow(data)))
  covariate_names <- setdiff(names(data), outcome_names)
  models <- lapply(outcome_names, function(nm) {
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
    fit_engine(spec, data[cc, , drop = FALSE], nm, ...)
  })
  names(models) <- outcome_names
  support_ref <- support_reference(data, covariate_names)
  structure(
    list(
      spec = spec,
      outcomes = outcome_names,
      models = models,
      id = ids,
      data_hash = digest_data(data[, c(covariate_names, outcome_names), drop = FALSE]),
      n = nrow(data),
      covariate_names = covariate_names,
      support_ref = support_ref,
      in_sample_rows = nrow(data),
      package_version = as.character(utils::packageVersion("referent"))
    ),
    class = "norm_fit"
  )
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
  cat(sprintf("norm_fit %s via %s\n", x$spec$family$name, x$spec$engine))
  cat(sprintf("%d outcomes, n = %d\n", length(x$outcomes), x$n))
  cat(sprintf("status: %s\n", status_txt))
  failed <- names(st)[st != "ok"]
  if (length(failed)) {
    cat(sprintf("Failed or flagged outcomes: %s\n", paste(failed, collapse = ", ")))
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
  support <- classify_support(object$support_ref, newdata)
  in_sample <- is_in_sample_data(object, newdata)
  rows <- lapply(nms, function(nm) {
    d <- dists[[nm]]
    y <- if (nm %in% names(newdata)) newdata[[nm]] else rep(NA_real_, nrow(newdata))
    if (is.null(d)) {
      return(tibble::tibble(
        .row = seq_len(nrow(newdata)),
        .id = if (".id" %in% names(newdata)) newdata[[".id"]] else seq_len(nrow(newdata)),
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
    sc$.id <- if (".id" %in% names(newdata)) newdata[[".id"]] else seq_len(nrow(newdata))
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

dplyr_bind <- function(xs) {
  tibble::as_tibble(do.call(rbind, lapply(xs, as.data.frame)))
}

is_in_sample_data <- function(fit, newdata) {
  cols <- intersect(c(fit$covariate_names, fit$outcomes), names(newdata))
  if (!length(cols) || nrow(newdata) != fit$n) {
    return(FALSE)
  }
  identical(digest_data(newdata[, cols, drop = FALSE]), fit$data_hash)
}

#' @export
print.norm_scores <- function(x, ...) {
  n_in <- sum(x$.in_sample, na.rm = TRUE)
  cli::cli_text("{.cls norm_scores} {nrow(x)} row{?s}")
  if (n_in > 0) {
    cli::cli_alert_warning("{n_in} in-sample score{?s} (.in_sample = TRUE)")
  }
  NextMethod()
}
