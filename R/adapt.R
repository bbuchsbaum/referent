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
#' `uncertainty = "total"`) is shifted and the scale multiplied by
#' \eqn{e^{\delta_\sigma}}. Under `uncertainty = "total"` the standard
#' error of the location offset is propagated as well: analytically for a
#' Gaussian predictive and as an extra seeded location perturbation per
#' coefficient draw for a draw mixture.
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
  by_name <- as_col_name(by_quo)
  if (inherits(fit, "norm_dynamics")) {
    fit$reference <- norm_adapt(
      fit$reference, data = data, by = !!by_quo, parameters = parameters,
      location_prior_n = location_prior_n, scale_prior_n = scale_prior_n
    )
    return(fit)
  }
  base <- fit
  base$adaptation <- NULL
  dists <- predict_dists(base, data, uncertainty = "conditional")
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
      estimate_offsets(d[idx], y[idx], parameters,
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
      location_prior_n = location_prior_n,
      scale_prior_n = scale_prior_n
    ),
    class = "norm_adaptation"
  )
  fit
}

# Penalised maximum-likelihood offsets on the family's own density. The
# location offset is parameterised in units of the mean predictive SD and
# the scale offset on the log scale, with ridge penalties equal to
# `prior_n` pseudo-observations of Fisher information, so the shrinkage
# factor is about n / (n + prior_n) for both. For a Gaussian family the
# estimates are the shrunk mean residual and 0.5 * log(mean(z^2)).
estimate_offsets <- function(d, y, parameters, location_prior_n, scale_prior_n) {
  ok <- is.finite(y)
  n <- sum(ok)
  empty <- list(location = 0, scale = 0, n = n, location_se = NA_real_,
                scale_se = NA_real_, parameters = parameters)
  u <- if (n >= 2L) dist_unpack(d[ok]) else NULL
  if (is.null(u)) {
    return(empty)
  }
  s_bar <- sqrt(mean(variance(d[ok])))
  if (!is.finite(s_bar) || s_bar <= 0) {
    return(empty)
  }
  y <- y[ok]
  fit_scale <- "scale" %in% parameters
  objective <- function(theta) {
    loc <- theta[[1L]] * s_bar
    log_s <- if (fit_scale) theta[[2L]] else 0
    ll <- sum(shifted_log_density(u, y, loc, log_s))
    if (!is.finite(ll)) {
      return(1e100)
    }
    -ll + 0.5 * location_prior_n * theta[[1L]]^2 +
      (if (fit_scale) scale_prior_n * log_s^2 else 0)
  }
  opt <- stats::optim(if (fit_scale) c(0, 0) else 0, objective,
                      method = "BFGS", hessian = TRUE)
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

# Log density of `y` under an unpacked distribution vector whose location
# is shifted by `loc` and scale multiplied by exp(`log_s`). The optimiser
# calls this many times, so the draw mixture is evaluated on its fields
# rather than rebuilt.
shifted_log_density <- function(u, y, loc, log_s) {
  check_adaptable(u)
  u$mu <- u$mu + loc
  u$sigma <- u$sigma * exp(log_s)
  log_dens(u, y)
}

check_adaptable <- function(u) {
  if (!inherits(u, c("dist_normal", "dist_shash", "dist_shash_draws"))) {
    cli::cli_abort("Cannot adapt a {.cls {class(u)[[1L]]}} distribution.")
  }
}

# Rebuild an unpacked distribution vector with new location and scale.
dist_rebuild <- function(u, mu, sigma) {
  check_adaptable(u)
  if (inherits(u, "dist_normal")) {
    return(distributional::dist_normal(mu, sigma))
  }
  if (inherits(u, "dist_shash")) {
    return(dist_shash(mu, sigma, u$eps, u$delta))
  }
  dist_shash_draws(mu, sigma, u$eps, u$delta)
}

# Shift the location of every element of a distribution vector by `loc`
# and multiply its scale by exp(`log_s`). A location-offset standard error
# `loc_se` is folded into a Gaussian predictive analytically and into a
# draw mixture as a standardised N(0, loc_se^2) perturbation per draw,
# seeded independently of the coefficient draws.
shift_params <- function(d, loc, log_s, loc_se = 0, seed = 1L) {
  n <- length(d)
  if (!n || all(vapply(vctrs::vec_data(d), is.null, logical(1)))) {
    return(d)
  }
  loc <- rep_len(loc, n)
  log_s <- rep_len(log_s, n)
  loc_se <- rep_len(loc_se, n)
  u <- dist_unpack(d)
  if (is.null(u)) {
    cli::cli_abort("Cannot adapt a distribution vector with missing or mixed elements.")
  }
  mu <- u$mu + loc
  sigma <- u$sigma * exp(log_s)
  if (inherits(u, "dist_normal")) {
    sigma <- sqrt(sigma^2 + loc_se^2)
  } else if (inherits(u, "dist_shash_draws")) {
    k <- ncol(u$mu)
    jitter <- withr::with_seed(seed + 1L, stats::rnorm(k))
    jitter <- if (k > 1L) as.numeric(scale(jitter)) else 0
    mu <- mu + outer(loc_se, jitter)
  }
  dist_rebuild(u, mu, sigma)
}

apply_adaptation <- function(adaptation, dists, newdata, total = FALSE, seed = 1L) {
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
    groups <- unique(grp)
    pick <- function(what) {
      v <- vapply(groups, function(g) (offs[[g]] %||% offs[[".all"]])[[what]] %||% 0,
                  numeric(1))
      unname(v[match(grp, groups)])
    }
    loc <- pick("location")
    log_s <- pick("scale")
    loc_se <- pick("location_se")
    loc_se[!is.finite(loc_se) | !total] <- 0
    shift_params(d, loc, log_s, loc_se, seed = seed)
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
