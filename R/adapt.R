#' Adapt a reference model to a new domain
#'
#' Freezes the shared trajectory and estimates shrunk location (and
#' optionally scale) offsets per group. This is not recalibration (see
#' [norm_calibrate()]) and not a refit of the shared trajectory.
#'
#' @details
#' Offsets \eqn{\delta_\mu} (outcome units) and \eqn{\delta_\sigma}
#' (log scale) are estimated per outcome and group by penalised maximum
#' likelihood on the reference family's own density with the shared
#' trajectory frozen, so the estimates are correct for SHASH as well as
#' Gaussian fits. The ridge penalties are worth `location_prior_n` and
#' `scale_prior_n` pseudo-observations of Fisher information, giving a
#' shrinkage factor of about \eqn{n/(n + \text{prior}_n)}. For a Gaussian
#' fit this is the shrunk mean residual and
#' \eqn{\tfrac12\log \mathrm{mean}(z^2)} with \eqn{z} the normal scores
#' under the location-shifted reference; a site drawn from the reference
#' generator has an expected offset of zero. Standard errors come from
#' the penalised Hessian.
#'
#' Adaptation is applied inside [predict.norm_fit()]: the location of
#' every predictive distribution (including coefficient draws under
#' `uncertainty = "total"`) is shifted, the scale multiplied by
#' \eqn{e^{\delta_\sigma}}, and the standard error of the location offset
#' is added to the epistemic SD.
#'
#' @param fit A [norm_fit] or [norm_dynamics] object.
#' @param data Local reference observations.
#' @param by Grouping column (typically site). Rows whose group is not in
#'   the adaptation data (or when `by` is `NULL`) use the pooled offset.
#' @param parameters Parameters to adapt: `"location"` and optionally
#'   `"scale"`.
#' @param location_prior_n Ridge strength in observation units.
#' @param scale_prior_n Stronger default shrinkage for scale.
#' @return The fit with an `adaptation` slot (class `norm_adaptation`)
#'   recording the offsets, their standard errors, and local sample sizes.
#' @export
norm_adapt <- function(fit,
                       data,
                       by = NULL,
                       parameters = c("location"),
                       location_prior_n = 10,
                       scale_prior_n = 25) {
  parameters <- match.arg(parameters, c("location", "scale"), several.ok = TRUE)
  data <- tibble::as_tibble(data)
  by_quo <- rlang::enquo(by)
  by_vec <- pull_column(data, by_quo, default = rep(".all", nrow(data)))
  by_name <- if (rlang::quo_is_null(by_quo) || rlang::quo_is_missing(by_quo)) {
    NULL
  } else {
    tryCatch(rlang::as_name(by_quo), error = function(e) NULL)
  }
  if (inherits(fit, "norm_dynamics")) {
    fit$reference <- norm_adapt(
      fit$reference, data = data, by = !!by_quo, parameters = parameters,
      location_prior_n = location_prior_n, scale_prior_n = scale_prior_n
    )
    return(fit)
  }
  base <- fit
  base$adaptation <- NULL
  dists <- predict(base, newdata = data, type = "distribution",
                   uncertainty = "conditional")
  by_chr <- as.character(by_vec)
  offsets <- lapply(names(dists), function(nm) {
    d <- dists[[nm]]
    if (is.null(d) || !nm %in% names(data)) {
      return(NULL)
    }
    y <- data[[nm]]
    groups <- split(seq_along(y), by_chr)
    groups$.all <- seq_along(y)
    lapply(groups, function(idx) {
      estimate_offsets(dist_slice(d, idx), y[idx], parameters,
                       location_prior_n, scale_prior_n)
    })
  })
  names(offsets) <- names(dists)
  offsets <- Filter(Negate(is.null), offsets)
  fit$adaptation <- structure(
    list(
      by = by_name,
      parameters = parameters,
      offsets = offsets,
      n_local = nrow(data),
      effective_n = table(by_chr),
      location_prior_n = location_prior_n,
      scale_prior_n = scale_prior_n
    ),
    class = "norm_adaptation"
  )
  fit
}

# Penalised maximum-likelihood offsets on the family's own density. The
# location offset is parameterised in units of the mean aleatoric SD and
# the scale offset on the log scale, with ridge penalties equal to
# `prior_n` pseudo-observations of Fisher information, so the shrinkage
# factor is about n / (n + prior_n) for both. For a Gaussian family the
# estimates are the shrunk mean residual and 0.5 * log(mean(z^2)).
estimate_offsets <- function(d, y, parameters, location_prior_n, scale_prior_n) {
  ok <- is.finite(y)
  n <- sum(ok)
  s_bar <- mean(field_or(d, "aleatoric_sd")[ok])
  empty <- list(location = 0, scale = 0, n = n, location_se = NA_real_,
                scale_se = NA_real_, parameters = parameters)
  if (n < 2L || !is.finite(s_bar) || s_bar <= 0) {
    return(empty)
  }
  d <- dist_slice(d, which(ok))
  y <- y[ok]
  fit_scale <- "scale" %in% parameters
  objective <- function(theta) {
    loc <- theta[[1L]] * s_bar
    log_s <- if (fit_scale) theta[[2L]] else 0
    ll <- sum(log_density(shift_dist(d, loc, log_s), y))
    if (!is.finite(ll)) {
      return(1e100)
    }
    -ll + 0.5 * location_prior_n * theta[[1L]]^2 +
      (if (fit_scale) scale_prior_n * log_s^2 else 0)
  }
  start <- if (fit_scale) c(0, 0) else 0
  opt <- if (fit_scale) {
    stats::optim(start, objective, method = "BFGS", hessian = TRUE)
  } else {
    o <- stats::optimize(objective, interval = c(-20, 20))
    h <- stats::optimHess(o$minimum, objective)
    list(par = o$minimum, hessian = h)
  }
  se <- tryCatch(sqrt(diag(solve(opt$hessian))), error = function(e) rep(NA_real_, 2))
  list(
    location = opt$par[[1L]] * s_bar,
    scale = if (fit_scale) opt$par[[2L]] else 0,
    n = n,
    location_se = se[[1L]] * s_bar,
    scale_se = if (fit_scale) se[[2L]] else NA_real_,
    parameters = parameters
  )
}

# Shift location by `loc` and multiply scale by exp(`log_s`), keeping any
# coefficient draws and the total-uncertainty bookkeeping consistent.
shift_dist <- function(d, loc, log_s, loc_se = 0) {
  n <- length(d)
  loc <- recycle_to(loc, n)
  log_s <- recycle_to(log_s, n)
  loc_se <- recycle_to(loc_se, n)
  p <- norm_params(d)
  alea <- field_or(d, "aleatoric_sd") * exp(log_s)
  ep <- sqrt(field_or(d, "epistemic_sd")^2 + loc_se^2)
  total <- identical(attr(d, "uncertainty"), "total")
  draws <- dist_draws(d)
  scale <- if (total && is.null(draws)) {
    sqrt(alea^2 + ep^2)
  } else {
    p$scale * exp(log_s)
  }
  out <- norm_dist(
    family = attr(d, "family"),
    location = p$location + loc,
    scale = scale,
    skew = p$skew,
    tail = p$tail,
    aleatoric_sd = alea,
    epistemic_sd = ep
  )
  if (!is.null(draws)) {
    attr(out, "location_draws") <- draws$location + loc
    attr(out, "scale_draws") <- draws$scale * exp(log_s)
    attr(out, "skew_draws") <- draws$skew
    attr(out, "tail_draws") <- draws$tail
  }
  if (total) {
    attr(out, "uncertainty") <- "total"
  }
  out
}

apply_adaptation <- function(adaptation, dists, newdata) {
  by_nm <- adaptation$by
  grp <- if (!is.null(by_nm) && by_nm %in% names(newdata)) {
    as.character(newdata[[by_nm]])
  } else {
    rep(".all", nrow(newdata))
  }
  grp[is.na(grp)] <- ".all"
  out <- lapply(names(dists), function(nm) {
    d <- dists[[nm]]
    offs <- adaptation$offsets[[nm]]
    if (is.null(d) || is.null(offs) || !length(d)) {
      return(d)
    }
    pick <- function(g, what) {
      off <- offs[[g]] %||% offs[[".all"]]
      off[[what]] %||% 0
    }
    loc <- vapply(grp, pick, numeric(1), what = "location")
    log_s <- vapply(grp, pick, numeric(1), what = "scale")
    loc_se <- vapply(grp, pick, numeric(1), what = "location_se")
    loc_se[!is.finite(loc_se)] <- 0
    shift_dist(d, loc, log_s, loc_se)
  })
  names(out) <- names(dists)
  out
}

#' @export
print.norm_adaptation <- function(x, ...) {
  cli::cli_text("{.cls norm_adaptation} parameters: {paste(x$parameters, collapse = ', ')}")
  cli::cli_text("local n = {x$n_local}")
  for (nm in names(x$offsets)) {
    for (g in setdiff(names(x$offsets[[nm]]), ".all")) {
      off <- x$offsets[[nm]][[g]]
      cli::cli_text(
        "  {nm} / {g}: location {signif(off$location, 3)}, scale x{signif(exp(off$scale), 3)} (n = {off$n})"
      )
    }
  }
  invisible(x)
}
