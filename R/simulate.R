#' Simulate reference-model data
#'
#' Synthetic data with known generative structure for examples and law
#' tests. Every kind returns `age`, `sex`, `site`, an outcome `y`, and
#' three correlated Gaussian markers `marker_01`..`marker_03`. The true
#' generative parameters are stored in `attr(x, "truth")`.
#'
#' * `"gaussian"`: `y` is Gaussian with location
#'   \eqn{\mu(age, sex)} and either a constant or an age-dependent scale.
#' * `"shash"`: `y` is sinh-arcsinh with the same location and scale and
#'   constant skew and tail.
#' * `"longitudinal"`: subjects have 2 to `visits + 1` visits with
#'   irregular lags (or exactly two visits a fixed `lag` apart when `lag`
#'   is given). The normal score of `y` follows a stable subject
#'   component plus a Matern-3/2 process in time plus measurement noise,
#'   \eqn{z_i(t) = \tau_b b_i + \tau_g g_i(t) + \sigma_e e_{it}}, with
#'   \eqn{\tau_b^2 + \tau_g^2 + \sigma_e^2 = 1} so that marginal scores are
#'   standard normal. Adds `participant_id` and `visit` columns.
#'
#' @param n Number of rows (for `"longitudinal"`, approximate).
#' @param kind Generator name.
#' @param sites Number of sites.
#' @param site_shift Location shift per site in outcome units, recycled
#'   to `sites`.
#' @param site_log_scale Log-scale shift per site, recycled to `sites`.
#' @param scale `"constant"` or `"age"`-dependent outcome scale.
#' @param skew,tail SHASH skew and tail (`"shash"` only).
#' @param visits Mean visits per subject (`"longitudinal"` only).
#' @param tau_b,tau_g,sigma_e Standard deviations of the stable, dynamic,
#'   and measurement components (`"longitudinal"` only). They are
#'   rescaled to sum to unit variance.
#' @param ell Matern-3/2 length-scale in time units (`"longitudinal"`
#'   only).
#' @param lag Optional fixed lag (`"longitudinal"` only): every subject
#'   then has exactly two visits this far apart, the design under which
#'   the Matern length-scale is not identified.
#' @param seed Optional seed; the global RNG state is restored afterwards.
#' @return A data frame with a `truth` attribute listing the generative
#'   parameters. For `"longitudinal"`, `truth$correlation(lag)` gives the
#'   implied correlation of normal scores at a given lag.
#' @examples
#' ref <- norm_simulate(200, seed = 1)
#' attr(ref, "truth")$sigma
#' long <- norm_simulate(300, kind = "longitudinal", seed = 2)
#' attr(long, "truth")$correlation(2)
#' @export
norm_simulate <- function(n = 400,
                          kind = c("gaussian", "shash", "longitudinal"),
                          sites = 4,
                          site_shift = 0,
                          site_log_scale = 0,
                          scale = c("constant", "age"),
                          skew = 0.6,
                          tail = 0.85,
                          visits = 3,
                          tau_b = sqrt(0.5),
                          tau_g = sqrt(0.35),
                          sigma_e = sqrt(0.15),
                          ell = 5,
                          lag = NULL,
                          seed = NULL) {
  kind <- match.arg(kind)
  scale <- match.arg(scale)
  run <- function() {
    simulate_impl(n, kind, sites, site_shift, site_log_scale, scale,
                  skew, tail, visits, tau_b, tau_g, sigma_e, ell, lag)
  }
  if (is.null(seed)) run() else withr::with_seed(seed, run())
}

sim_location <- function(age, sex) {
  10 + 0.08 * (age - 50) - 0.001 * (age - 50)^2 + 0.4 * (sex == "M")
}

sim_scale <- function(age, scale) {
  if (identical(scale, "age")) 1.2 + 0.02 * pmax(age - 40, 0) else rep(1.3, length(age))
}

sim_markers <- function(age, z) {
  n <- length(age)
  lambda <- c(0.8, 0.6, 0.3)
  mk <- lapply(seq_along(lambda), function(j) {
    loc <- 10 - 0.2 * (j - 1) * (age - 50) / 10
    loc + 1.1 * (lambda[[j]] * z + sqrt(1 - lambda[[j]]^2) * stats::rnorm(n))
  })
  names(mk) <- sprintf("marker_%02d", seq_along(lambda))
  mk
}

simulate_impl <- function(n, kind, sites, site_shift, site_log_scale, scale,
                          skew, tail, visits, tau_b, tau_g, sigma_e, ell, lag = NULL) {
  site_levels <- LETTERS[seq_len(sites)]
  site_shift <- rep_len(site_shift, sites)
  site_log_scale <- rep_len(site_log_scale, sites)
  names(site_shift) <- names(site_log_scale) <- site_levels
  truth <- list(
    kind = kind,
    location = sim_location,
    scale = scale,
    site_shift = site_shift,
    site_log_scale = site_log_scale
  )
  if (identical(kind, "longitudinal")) {
    fixed <- !is.null(lag)
    n_subject <- max(ceiling(n / if (fixed) 2 else visits), 2L)
    n_visit <- if (fixed) rep(2L, n_subject) else sample(2:(visits + 1L), n_subject, replace = TRUE)
    id <- rep(seq_len(n_subject), n_visit)
    age0 <- stats::runif(n_subject, 20, 72)
    lags <- if (fixed) rep(lag, length(id)) else stats::runif(length(id), 0.5, 3)
    age <- age0[id] + stats::ave(lags, id, FUN = function(l) cumsum(l) - l[[1L]])
    age <- if (fixed) age else pmin(age, 80)
    sex <- factor(sample(c("F", "M"), n_subject, replace = TRUE))[id]
    site <- factor(sample(site_levels, n_subject, replace = TRUE), levels = site_levels)[id]
    tot <- sqrt(tau_b^2 + tau_g^2 + sigma_e^2)
    tau_b <- tau_b / tot
    tau_g <- tau_g / tot
    sigma_e <- sigma_e / tot
    b <- stats::rnorm(n_subject)[id]
    g <- unlist(lapply(split(age, id), function(t) {
      k <- matern32_kernel(abs(outer(t, t, `-`)), ell) + diag(1e-8, length(t))
      drop(t(chol(k)) %*% stats::rnorm(length(t)))
    }), use.names = FALSE)
    g <- g[order(order(id))]
    z <- tau_b * b + tau_g * g + sigma_e * stats::rnorm(length(id))
    truth$tau_b <- tau_b
    truth$tau_g <- tau_g
    truth$sigma_e <- sigma_e
    truth$ell <- ell
    truth$correlation <- function(lag) {
      tau_b^2 + tau_g^2 * matern32_kernel(lag, ell) + sigma_e^2 * (lag == 0)
    }
    visit <- stats::ave(seq_along(id), id, FUN = seq_along)
  } else {
    age <- stats::runif(n, 20, 80)
    sex <- factor(sample(c("F", "M"), n, replace = TRUE))
    site <- factor(sample(site_levels, n, replace = TRUE), levels = site_levels)
    z <- stats::rnorm(n)
    id <- NULL
    visit <- NULL
  }
  mu <- sim_location(age, sex) + site_shift[as.character(site)]
  sigma <- sim_scale(age, scale) * exp(site_log_scale[as.character(site)])
  truth$sigma <- sim_scale
  if (identical(kind, "shash")) {
    truth$skew <- skew
    truth$tail <- tail
    y <- shash_quantile(stats::pnorm(z), mu, sigma, skew, tail)
  } else {
    y <- mu + sigma * z
  }
  out <- data.frame(age = age, sex = sex, site = site, y = as.numeric(y))
  out <- cbind(out, as.data.frame(sim_markers(age, z)))
  if (!is.null(id)) {
    out <- cbind(data.frame(participant_id = id), out, data.frame(visit = visit))
  }
  attr(out, "truth") <- truth
  out
}
