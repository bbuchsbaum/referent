#' Predictive distribution objects
#'
#' A vectorized conditional predictive distribution. Every engine must
#' return an object that supports this contract: `cdf()`, `quantile()`,
#' `log_density()`, `draw()`, `center()`, and `variance()`.
#'
#' Parameters follow the `mgcv` SHASH / Gaussian location-scale
#' conventions:
#'
#' * `location` (\eqn{\mu})
#' * `scale` (\eqn{\sigma > 0})
#' * `skew` (\eqn{\epsilon}; SHASH only)
#' * `tail` (\eqn{\delta > 0}; SHASH only)
#'
#' For SHASH, \eqn{z = (y-\mu)/(\sigma\delta)} and
#' \eqn{F(y) = \Phi(\sinh(\delta\,\mathrm{asinh}(z)-\epsilon))}.
#'
#' @param family Family name or a [norm_family] object.
#' @param location,scale,skew,tail Numeric parameter vectors. Recycled to
#'   a common length. `NA` parameters give `NA` scores (for example when a
#'   predictor is missing).
#' @param aleatoric_sd,epistemic_sd Optional uncertainty summaries,
#'   recycled to the same length.
#' @return An object of class `norm_dist`.
#' @examples
#' d <- norm_dist("gaussian", location = 0, scale = 1)
#' cdf(d, 1.96)
#' as_scores(d, c(-2, 0, 2))
#' @export
norm_dist <- function(family,
                      location,
                      scale,
                      skew = NULL,
                      tail = NULL,
                      aleatoric_sd = NULL,
                      epistemic_sd = NULL) {
  fam <- if (inherits(family, "norm_family")) family$name else as.character(family)
  fam <- fam[[1L]]
  n <- max(length(location), length(scale), length(skew %||% 1), length(tail %||% 1))
  if (!length(location) && !length(scale)) {
    n <- 0L
  }
  location <- if (n == 0L) numeric() else recycle_to(location, n)
  scale <- if (n == 0L) numeric() else recycle_to(scale, n)
  if (n > 0L && any(scale <= 0, na.rm = TRUE)) {
    cli::cli_abort("{.arg scale} must be positive.")
  }
  if (fam == "shash") {
    skew <- recycle_to(skew %||% 0, n)
    tail <- recycle_to(tail %||% 1, n)
    if (any(tail <= 0, na.rm = TRUE)) {
      cli::cli_abort("{.arg tail} must be positive.")
    }
  } else {
    skew <- recycle_to(0, n)
    tail <- recycle_to(1, n)
  }
  aleatoric_sd <- recycle_to(aleatoric_sd %||% scale, n)
  epistemic_sd <- recycle_to(epistemic_sd %||% 0, n)
  vctrs::new_rcrd(
    list(
      location = location,
      scale = scale,
      skew = skew,
      tail = tail,
      aleatoric_sd = aleatoric_sd,
      epistemic_sd = epistemic_sd
    ),
    family = fam,
    class = "norm_dist"
  )
}

#' @exportS3Method vctrs::vec_ptype_abbr
vec_ptype_abbr.norm_dist <- function(x, ...) {
  paste0("ndist<", attr(x, "family"), ">")
}

#' @export
format.norm_dist <- function(x, ...) {
  paste0(
    attr(x, "family"),
    "(mu=", signif(field_or(x, "location"), 4),
    ", sd=", signif(field_or(x, "scale"), 4), ")"
  )
}

field_or <- function(x, name) {
  vctrs::field(x, name)
}

norm_params <- function(x) {
  list(
    location = field_or(x, "location"),
    scale = field_or(x, "scale"),
    skew = field_or(x, "skew"),
    tail = field_or(x, "tail")
  )
}

#' Cumulative distribution function
#'
#' @param distribution A [norm_dist].
#' @param y Numeric observations, recycled to the distribution length.
#' @param x A [norm_dist] (for [stats::quantile()]).
#' @param probs Probabilities in \eqn{(0,1)}.
#' @param ... Unused.
#' @export
cdf <- function(distribution, y) {
  UseMethod("cdf")
}

#' @export
cdf.norm_dist <- function(distribution, y) {
  exp(log_cdf(distribution, y, lower.tail = TRUE))
}

# Log CDF (or log survival function) computed without passing through the
# probability scale, so that tails of 1e-40 remain distinguishable.
# A distribution carrying coefficient draws (uncertainty = "total") is a
# mixture over the draws; its tails are log-mean-exp over the components.
log_cdf <- function(distribution, y, lower.tail = TRUE) {
  UseMethod("log_cdf")
}

#' @export
#' @noRd
log_cdf.norm_dist <- function(distribution, y, lower.tail = TRUE) {
  rec <- recycle_pair(distribution, y)
  distribution <- rec$distribution
  y <- rec$y
  draws <- dist_draws(distribution)
  if (is.null(draws)) {
    p <- norm_params(distribution)
    return(family_log_cdf(attr(distribution, "family"), y, p, lower.tail))
  }
  fam <- attr(distribution, "family")
  k <- ncol(draws$location)
  acc <- matrix(NA_real_, length(y), k)
  for (j in seq_len(k)) {
    acc[, j] <- family_log_cdf(fam, y, draw_params(draws, j), lower.tail)
  }
  log_mean_exp_rows(acc)
}

family_log_cdf <- function(fam, y, p, lower.tail = TRUE) {
  switch(
    fam,
    gaussian = stats::pnorm(y, mean = p$location, sd = p$scale,
                            lower.tail = lower.tail, log.p = TRUE),
    shash = shash_log_cdf(y, p$location, p$scale, p$skew, p$tail, lower.tail),
    cli::cli_abort("Unknown family {.val {fam}}.")
  )
}

log_mean_exp_rows <- function(m) {
  mx <- apply(m, 1, max)
  ok <- is.finite(mx)
  out <- rep(NA_real_, nrow(m))
  out[ok] <- mx[ok] + log(rowMeans(exp(m[ok, , drop = FALSE] - mx[ok])))
  out[!ok & !is.na(mx) & mx == -Inf] <- -Inf
  out
}

# Coefficient-draw parameter matrices attached by the engine for
# uncertainty = "total". Rows are observations, columns are draws.
dist_draws <- function(distribution) {
  loc <- attr(distribution, "location_draws")
  if (is.null(loc)) {
    return(NULL)
  }
  p <- norm_params(distribution)
  k <- ncol(loc)
  fill <- function(nm, default) {
    m <- attr(distribution, nm)
    if (is.null(m)) matrix(default, nrow(loc), k) else m
  }
  list(
    location = loc,
    scale = fill("scale_draws", p$scale),
    skew = fill("skew_draws", p$skew),
    tail = fill("tail_draws", p$tail)
  )
}

draw_params <- function(draws, j) {
  list(
    location = draws$location[, j],
    scale = draws$scale[, j],
    skew = draws$skew[, j],
    tail = draws$tail[, j]
  )
}

# Subset the draw matrices along with the distribution.
dist_slice <- function(distribution, i) {
  out <- vctrs::vec_slice(distribution, i)
  for (nm in c("location_draws", "scale_draws", "skew_draws", "tail_draws")) {
    m <- attr(distribution, nm)
    if (!is.null(m)) {
      attr(out, nm) <- m[i, , drop = FALSE]
    }
  }
  for (nm in setdiff(names(attributes(distribution)), names(attributes(out)))) {
    attr(out, nm) <- attr(distribution, nm)
  }
  out
}

#' @export
#' @rdname cdf
quantile.norm_dist <- function(x, probs = seq(0, 1, 0.25), ...) {
  rec <- recycle_pair(x, probs)
  x <- rec$distribution
  p <- clamp_prob(rec$y)
  par <- norm_params(x)
  draws <- dist_draws(x)
  if (is.null(draws)) {
    return(family_quantile(attr(x, "family"), p, par))
  }
  mixture_quantile(x, p, draws)
}

family_quantile <- function(fam, p, par) {
  switch(
    fam,
    gaussian = stats::qnorm(p, mean = par$location, sd = par$scale),
    shash = shash_quantile(p, par$location, par$scale, par$skew, par$tail),
    cli::cli_abort("Unknown family {.val {fam}}.")
  )
}

# Quantile of the draw mixture by vectorised bisection on the mixture CDF,
# bracketed by the extreme component quantiles.
mixture_quantile <- function(x, p, draws, iter = 60L) {
  fam <- attr(x, "family")
  k <- ncol(draws$location)
  lo <- rep(Inf, length(p))
  hi <- rep(-Inf, length(p))
  for (j in seq_len(k)) {
    qj <- family_quantile(fam, p, draw_params(draws, j))
    lo <- pmin(lo, qj)
    hi <- pmax(hi, qj)
  }
  ok <- is.finite(lo) & is.finite(hi)
  out <- rep(NA_real_, length(p))
  if (!any(ok)) {
    return(out)
  }
  for (it in seq_len(iter)) {
    mid <- (lo + hi) / 2
    u <- exp(log_cdf(x, mid, lower.tail = TRUE))
    below <- u < p
    lo <- ifelse(below, mid, lo)
    hi <- ifelse(below, hi, mid)
    if (all(abs(hi - lo)[ok] < 1e-9 * pmax(1, abs(mid[ok])))) {
      break
    }
  }
  out[ok] <- ((lo + hi) / 2)[ok]
  out
}

# n x length(probs) matrix of quantiles, one column per probability.
quantile_matrix <- function(distribution, probs) {
  n <- length(distribution)
  out <- matrix(NA_real_, n, length(probs))
  for (k in seq_along(probs)) {
    out[, k] <- quantile(distribution, probs = rep(probs[[k]], n))
  }
  out
}

dist_quantile <- function(distribution, p) {
  quantile(distribution, probs = p)
}

#' Log density
#'
#' @inheritParams cdf
#' @export
log_density <- function(distribution, y) {
  UseMethod("log_density")
}

#' @export
log_density.norm_dist <- function(distribution, y) {
  rec <- recycle_pair(distribution, y)
  distribution <- rec$distribution
  y <- rec$y
  draws <- dist_draws(distribution)
  if (is.null(draws)) {
    return(family_log_density(attr(distribution, "family"), y, norm_params(distribution)))
  }
  fam <- attr(distribution, "family")
  k <- ncol(draws$location)
  acc <- matrix(NA_real_, length(y), k)
  for (j in seq_len(k)) {
    acc[, j] <- family_log_density(fam, y, draw_params(draws, j))
  }
  log_mean_exp_rows(acc)
}

family_log_density <- function(fam, y, p) {
  switch(
    fam,
    gaussian = stats::dnorm(y, mean = p$location, sd = p$scale, log = TRUE),
    shash = shash_log_density(y, p$location, p$scale, p$skew, p$tail),
    cli::cli_abort("Unknown family {.val {fam}}.")
  )
}

#' Simulate from a predictive distribution
#'
#' @inheritParams cdf
#' @param n Number of draws *per* observation. Returns a matrix with
#'   `length(distribution)` rows when `n > 1`.
#' @export
draw <- function(distribution, n = 1L) {
  UseMethod("draw")
}

#' @export
draw.norm_dist <- function(distribution, n = 1L) {
  n <- as.integer(n)
  m <- length(distribution)
  u <- matrix(stats::runif(m * n), nrow = m, ncol = n)
  out <- matrix(NA_real_, nrow = m, ncol = n)
  for (j in seq_len(n)) {
    out[, j] <- dist_quantile(distribution, u[, j])
  }
  if (n == 1L) {
    drop(out)
  } else {
    out
  }
}

#' Predictive center (median)
#'
#' @inheritParams cdf
#' @export
center <- function(distribution) {
  UseMethod("center")
}

#' @export
center.norm_dist <- function(distribution) {
  dist_quantile(distribution, rep(0.5, length(distribution)))
}

#' Predictive variance
#'
#' @inheritParams cdf
#' @export
variance <- function(distribution) {
  UseMethod("variance")
}

#' @export
variance.norm_dist <- function(distribution) {
  p <- norm_params(distribution)
  fam <- attr(distribution, "family")
  if (identical(fam, "gaussian") && is.null(dist_draws(distribution))) {
    return(p$scale^2)
  }
  u <- (seq_len(199L) - 0.5) / 199
  q <- quantile_matrix(distribution, u)
  apply(q, 1, stats::var)
}

#' Scores derived from a predictive CDF
#'
#' @inheritParams cdf
#' @param y Observed values.
#' @return A data frame of centiles, Z-scores, tails, residuals, and
#'   log densities. There is no abnormality column.
#' @export
as_scores <- function(distribution, y) {
  rec <- recycle_pair(distribution, y)
  distribution <- rec$distribution
  y <- rec$y
  log_lower <- log_cdf(distribution, y, lower.tail = TRUE)
  log_upper <- log_cdf(distribution, y, lower.tail = FALSE)
  tails <- scores_from_log_tails(log_lower, log_upper)
  med <- center(distribution)
  tibble::tibble(
    observed = y,
    median = med,
    centile = tails$centile,
    z = tails$z,
    tail_prob = tails$tail_prob,
    tail_surprisal = tails$tail_surprisal,
    residual = y - med,
    log_density = log_density(distribution, y),
    aleatoric_sd = field_or(distribution, "aleatoric_sd"),
    epistemic_sd = field_or(distribution, "epistemic_sd")
  )
}

# Centile, z, and two-sided tail from log lower/upper tail probabilities.
# z is taken from whichever tail is smaller, so z = 8, 10, 40 are distinct.
scores_from_log_tails <- function(log_lower, log_upper) {
  use_lower <- is.na(log_upper) | (!is.na(log_lower) & log_lower <= log_upper)
  z <- ifelse(
    use_lower,
    stats::qnorm(log_lower, log.p = TRUE, lower.tail = TRUE),
    stats::qnorm(log_upper, log.p = TRUE, lower.tail = FALSE)
  )
  log_tail <- log(2) + pmin(log_lower, log_upper)
  log_tail <- pmin(log_tail, 0)
  list(
    centile = ifelse(use_lower, exp(log_lower), -expm1(log_upper)),
    z = z,
    tail_prob = exp(log_tail),
    tail_surprisal = -log_tail
  )
}

# Log tails recovered from a z-score (exact inverse of scores_from_log_tails).
log_tails_from_z <- function(z) {
  list(
    lower = stats::pnorm(z, log.p = TRUE, lower.tail = TRUE),
    upper = stats::pnorm(z, log.p = TRUE, lower.tail = FALSE)
  )
}

# --- SHASH (mgcv / Jones-Pewsey parameterization) -------------------------

shash_z <- function(y, mu, sigma, delta) {
  (y - mu) / (sigma * delta)
}

shash_cdf <- function(y, mu, sigma, eps, delta) {
  z <- shash_z(y, mu, sigma, delta)
  s <- sinh(delta * asinh(z) - eps)
  stats::pnorm(s)
}

shash_log_cdf <- function(y, mu, sigma, eps, delta, lower.tail = TRUE) {
  z <- shash_z(y, mu, sigma, delta)
  s <- sinh(delta * asinh(z) - eps)
  stats::pnorm(s, lower.tail = lower.tail, log.p = TRUE)
}

shash_quantile <- function(p, mu, sigma, eps, delta) {
  z <- sinh((asinh(stats::qnorm(p)) + eps) / delta)
  mu + sigma * delta * z
}

shash_log_density <- function(y, mu, sigma, eps, delta) {
  z <- shash_z(y, mu, sigma, delta)
  s <- sinh(delta * asinh(z) - eps)
  c_z <- sqrt(1 + s^2)
  log(c_z) - 0.5 * s^2 - 0.5 * log(2 * pi) - 0.5 * log1p(z^2) - log(sigma)
}

#' Mix a history-conditioned Gaussian score into a marginal CDF
#'
#' Implements \eqn{F_*(y\mid H)=\Phi((z(y)-m_*)/s_*)} from the velocity
#' design memo.
#'
#' @param distribution Marginal [norm_dist] at the forecast time.
#' @param m,s Conditional mean and SD of the latent normal score.
#' @keywords internal
#' @noRd
condition_norm_dist <- function(distribution, m, s) {
  m <- recycle_to(m, length(distribution))
  s <- recycle_to(s, length(distribution))
  if (any(s <= 0, na.rm = TRUE)) {
    cli::cli_abort("{.arg s} must be positive.")
  }
  structure(
    distribution,
    history_m = m,
    history_s = s,
    class = c("norm_dist_conditional", class(distribution))
  )
}

#' @export
cdf.norm_dist_conditional <- function(distribution, y) {
  exp(log_cdf(distribution, y, lower.tail = TRUE))
}

#' @export
#' @noRd
log_cdf.norm_dist_conditional <- function(distribution, y, lower.tail = TRUE) {
  rec <- recycle_pair(distribution, y)
  distribution <- rec$distribution
  y <- rec$y
  bare <- distribution
  class(bare) <- setdiff(class(bare), "norm_dist_conditional")
  z <- scores_from_log_tails(
    log_cdf(bare, y, lower.tail = TRUE),
    log_cdf(bare, y, lower.tail = FALSE)
  )$z
  m <- recycle_to(attr(distribution, "history_m") %||% 0, length(y))
  s <- recycle_to(attr(distribution, "history_s") %||% 1, length(y))
  stats::pnorm((z - m) / s, lower.tail = lower.tail, log.p = TRUE)
}

#' @export
quantile.norm_dist_conditional <- function(x, probs = seq(0, 1, 0.25), ...) {
  rec <- recycle_pair(x, probs)
  x <- rec$distribution
  p <- rec$y
  m <- recycle_to(attr(x, "history_m") %||% 0, length(p))
  s <- recycle_to(attr(x, "history_s") %||% 1, length(p))
  p_marg <- stats::pnorm(m + s * stats::qnorm(clamp_prob(p)))
  class(x) <- setdiff(class(x), "norm_dist_conditional")
  quantile(x, probs = p_marg)
}

#' @export
log_density.norm_dist_conditional <- function(distribution, y) {
  y <- recycle_to(y, length(distribution))
  class_bare <- distribution
  class(class_bare) <- setdiff(class(class_bare), "norm_dist_conditional")
  z <- scores_from_log_tails(
    log_cdf(class_bare, y, lower.tail = TRUE),
    log_cdf(class_bare, y, lower.tail = FALSE)
  )$z
  m <- attr(distribution, "history_m")
  s <- attr(distribution, "history_s")
  log_f <- log_density(class_bare, y)
  log_f + stats::dnorm((z - m) / s, log = TRUE) - log(s) - stats::dnorm(z, log = TRUE)
}

#' @export
center.norm_dist_conditional <- function(distribution) {
  quantile(distribution, probs = 0.5)
}

#' @export
draw.norm_dist_conditional <- function(distribution, n = 1L) {
  n <- as.integer(n)
  m_len <- length(distribution)
  out <- matrix(NA_real_, m_len, n)
  for (j in seq_len(n)) {
    out[, j] <- quantile(distribution, probs = stats::runif(m_len))
  }
  if (n == 1L) drop(out) else out
}
