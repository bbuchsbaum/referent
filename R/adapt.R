#' Adapt a reference model to a new domain
#'
#' Freezes the shared trajectory and estimates shrunk location (and
#' optionally scale) offsets. This is not recalibration and not a refit
#' of the shared trajectory.
#'
#' @param fit A [norm_fit] or [norm_dynamics] object.
#' @param data Local reference observations.
#' @param by Grouping column (typically site).
#' @param parameters Parameters to adapt. For a static fit: `location`
#'   and optionally `scale`. For dynamics: also `measurement`.
#' @param components Alias of `parameters` used by the dynamics API.
#' @param location_prior_n Ridge strength in observation units.
#' @param scale_prior_n Stronger default shrinkage for scale.
#' @return A `norm_adaptation` wrapping the original fit.
#' @export
norm_adapt <- function(fit,
                       data,
                       by = NULL,
                       parameters = c("location"),
                       components = NULL,
                       location_prior_n = 10,
                       scale_prior_n = 25) {
  parameters <- unique(c(parameters, components))
  data <- tibble::as_tibble(data)
  by_quo <- rlang::enquo(by)
  by_vec <- pull_column(data, by_quo, default = rep(".all", nrow(data)))
  by_name <- if (rlang::quo_is_null(by_quo) || rlang::quo_is_missing(by_quo)) {
    NULL
  } else {
    tryCatch(rlang::as_name(by_quo), error = function(e) NULL)
  }
  if (inherits(fit, "norm_dynamics")) {
    return(adapt_dynamics(fit, data, by_vec, by_name, parameters,
                          location_prior_n, scale_prior_n))
  }
  scores <- predict(fit, newdata = data, type = "scores",
                    uncertainty = "conditional", allow_extrapolation = TRUE)
  offsets <- lapply(unique(scores$.outcome), function(nm) {
    sc <- scores[scores$.outcome == nm, , drop = FALSE]
    groups <- split(seq_len(nrow(sc)), by_vec[sc$.row])
    lapply(groups, function(idx) {
      r <- sc$residual[idx]
      s <- sc$aleatoric_sd[idx]
      n_g <- sum(is.finite(r))
      mu <- shrink_mean(r, n_g, location_prior_n)
      log_s <- if ("scale" %in% parameters) {
        shrink_mean(safe_log(pmax(abs(r), 1e-6)) - safe_log(pmax(s, 1e-6)),
                    n_g, scale_prior_n)
      } else {
        0
      }
      list(
        location = mu,
        scale = log_s,
        n = n_g,
        location_se = se_mean(r[is.finite(r)]),
        parameters = parameters
      )
    })
  })
  names(offsets) <- unique(scores$.outcome)
  structure(
    list(
      fit = fit,
      by = by_name,
      parameters = parameters,
      offsets = offsets,
      n_local = nrow(data),
      effective_n = tapply(by_vec, by_vec, length)
    ),
    class = c("norm_adaptation", "norm_fit")
  )
}

shrink_mean <- function(x, n, prior_n) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(0)
  }
  w <- n / (n + prior_n)
  w * mean(x)
}

#' @export
print.norm_adaptation <- function(x, ...) {
  cli::cli_text("{.cls norm_adaptation} parameters: {paste(x$parameters, collapse = ', ')}")
  cli::cli_text("local n = {x$n_local}")
  invisible(x)
}

#' @export
predict.norm_adaptation <- function(object,
                                    newdata,
                                    type = c("scores", "distribution"),
                                    uncertainty = c("total", "conditional"),
                                    ...) {
  type <- match.arg(type)
  uncertainty <- match.arg(uncertainty)
  base <- object$fit
  class(base) <- setdiff(class(base), "norm_adaptation")
  dists <- predict(base, newdata = newdata, type = "distribution",
                   uncertainty = uncertainty, ...)
  by_nm <- object$by
  grp <- if (!is.null(by_nm) && by_nm %in% names(newdata)) {
    as.character(newdata[[by_nm]])
  } else {
    rep(".all", nrow(newdata))
  }
  adapted <- lapply(names(dists), function(nm) {
    d <- dists[[nm]]
    if (is.null(d)) {
      return(NULL)
    }
    p <- norm_params(d)
    loc <- p$location
    sc <- p$scale
    ep <- field_or(d, "epistemic_sd")
    for (i in seq_along(loc)) {
      off <- object$offsets[[nm]][[grp[[i]]]] %||% object$offsets[[nm]][[".all"]]
      if (is.null(off)) {
        next
      }
      loc[[i]] <- loc[[i]] + off$location
      sc[[i]] <- sc[[i]] * exp(off$scale)
      ep[[i]] <- sqrt(ep[[i]]^2 + (off$location_se %||% 0)^2)
    }
    norm_dist(
      family = attr(d, "family"),
      location = loc,
      scale = sc,
      skew = p$skew,
      tail = p$tail,
      aleatoric_sd = sc,
      epistemic_sd = ep
    )
  })
  names(adapted) <- names(dists)
  if (identical(type, "distribution")) {
    return(adapted)
  }
  support <- classify_support(object$fit$support_ref, newdata)$support
  rows <- lapply(names(adapted), function(nm) {
    d <- adapted[[nm]]
    y <- if (nm %in% names(newdata)) newdata[[nm]] else rep(NA_real_, nrow(newdata))
    if (is.null(d)) {
      return(NULL)
    }
    sc <- as_scores(d, y)
    sc$.row <- seq_len(nrow(newdata))
    sc$.id <- seq_len(nrow(newdata))
    sc$.outcome <- nm
    sc$support <- support
    sc$calibrated <- FALSE
    sc$.in_sample <- FALSE
    sc$status <- "ok"
    sc
  })
  dplyr_bind(Filter(Negate(is.null), rows))
}

adapt_dynamics <- function(fit, data, by_vec, by_name, parameters,
                           location_prior_n, scale_prior_n) {
  static <- norm_adapt(
    fit$reference,
    data = data,
    by = NULL,
    parameters = intersect(parameters, c("location", "scale")),
    location_prior_n = location_prior_n,
    scale_prior_n = scale_prior_n
  )
  meas <- 0
  if ("measurement" %in% parameters && !is.null(fit$process)) {
    # local short-interval residuals inflate measurement variance
    y0 <- data[[fit$outcomes[[1]]]]
    meas <- stats::var(as.numeric(scale(y0)[, 1]), na.rm = TRUE)
    meas <- shrink_mean(rep(meas, nrow(data)), nrow(data), scale_prior_n)
  }
  fit$adaptation <- static
  fit$measurement_offset <- meas
  fit$adapted_parameters <- parameters
  class(fit) <- unique(c("norm_adaptation", class(fit)))
  fit
}
