#' Predictive distributions
#'
#' Reference models return predictive distributions as
#' [distributional](https://pkg.mitchelloharawild.com/distributional/)
#' vectors, so `cdf()`, `quantile()`, `density()`, `generate()`, `mean()`,
#' `variance()`, `hilo()`, tibble columns, and `ggdist` layers all work.
#' Gaussian fits use [distributional::dist_normal()]; SHASH fits use
#' `dist_shash()`; `uncertainty = "total"` predictions from location-scale
#' and SHASH fits are equal-weight mixtures over coefficient draws (an
#' internal `dist_shash_draws` class, one SHASH parameter set per draw and
#' observation); history-conditioned forecasts are `dist_conditioned()`; a
#' spec with a response transform returns the push-forward of the fitted
#' distribution onto the response scale (an internal `dist_warped` class).
#'
#' @details
#' `dist_shash()` follows the `mgcv::shash()` parameterisation. With
#' \eqn{z = (y-\mu)/(\sigma\delta)},
#' \eqn{F(y) = \Phi(\sinh(\delta\,\mathrm{asinh}(z)-\epsilon))}; `eps = 0`,
#' `delta = 1` is the Gaussian \eqn{N(\mu, \sigma^2)}.
#'
#' `dist_conditioned()` maps a marginal predictive distribution through a
#' conditional normal score: \eqn{F_*(y) = \Phi((z(y)-m)/s)} with
#' \eqn{z(y) = \Phi^{-1}(F(y))}, so `m = 0`, `s = 1` recovers the marginal.
#'
#' Tail probabilities of every class are evaluated in log space, so scores
#' for observations 40 standard deviations out remain finite and distinct
#' (see [as_scores()]).
#'
#' @param mu,sigma,eps,delta Location, scale (> 0), skewness, and tail
#'   weight (> 0), recycled to a common length.
#' @param marginal A `distribution` vector of marginal predictive
#'   distributions.
#' @param m,s Conditional mean and standard deviation (> 0) of the latent
#'   normal score, recycled to `length(marginal)`.
#' @return A [distributional][distributional::distributional-package]
#'   `distribution` vector.
#' @examples
#' d <- dist_shash(mu = 0.2, sigma = 1.1, eps = 0.3, delta = 0.9)
#' d
#' cdf(d, 1.5)
#' quantile(d, c(0.05, 0.5, 0.95))
#' mean(d)
#' as_scores(d, c(-2, 0, 2))
#' dist_conditioned(d, m = 0.4, s = 0.7)
#' @name dist_shash
NULL

#' @rdname dist_shash
#' @export
dist_shash <- function(mu, sigma, eps = 0, delta = 1) {
  check_shash_params(sigma, delta)
  distributional::new_dist(
    mu = vctrs::vec_cast(mu, double()),
    sigma = vctrs::vec_cast(sigma, double()),
    eps = vctrs::vec_cast(eps, double()),
    delta = vctrs::vec_cast(delta, double()),
    class = "dist_shash"
  )
}

# Equal-weight mixture over coefficient draws: every parameter is a matrix
# with one row per observation and one column per draw (a vector is
# draw-invariant). Its CDF and density are averages over the draws in log
# space, the quantile is found by safeguarded Newton iteration on the
# mixture CDF, and generate() samples a draw and then the component.
dist_shash_draws <- function(mu, sigma, eps = 0, delta = 1) {
  check_shash_params(sigma, delta)
  n <- max(NROW(mu), NROW(sigma), NROW(eps), NROW(delta))
  k <- max(NCOL(mu), NCOL(sigma), NCOL(eps), NCOL(delta))
  rows <- function(x) {
    x <- as.matrix(x)
    if (nrow(x) == 1L && n > 1L) {
      x <- x[rep(1L, n), , drop = FALSE]
    }
    if (ncol(x) == 1L && k > 1L) {
      x <- x[, rep(1L, k), drop = FALSE]
    }
    lapply(seq_len(nrow(x)), function(i) as.numeric(x[i, ]))
  }
  distributional::new_dist(
    mu = rows(mu), sigma = rows(sigma), eps = rows(eps), delta = rows(delta),
    class = "dist_shash_draws"
  )
}

#' @rdname dist_shash
#' @export
dist_conditioned <- function(marginal, m, s) {
  if (!distributional::is_distribution(marginal)) {
    cli::cli_abort("{.arg marginal} must be a {.cls distribution} vector.")
  }
  if (any(s <= 0, na.rm = TRUE)) {
    cli::cli_abort("{.arg s} must be positive.")
  }
  distributional::new_dist(
    dist = vctrs::vec_data(marginal),
    m = vctrs::vec_cast(m, double()),
    s = vctrs::vec_cast(s, double()),
    class = "dist_conditioned"
  )
}

check_shash_params <- function(sigma, delta) {
  if (any(sigma <= 0, na.rm = TRUE)) {
    cli::cli_abort("{.arg sigma} must be positive.")
  }
  if (any(delta <= 0, na.rm = TRUE)) {
    cli::cli_abort("{.arg delta} must be positive.")
  }
}

#' @importFrom distributional cdf generate hilo variance
#' @export
distributional::cdf

#' @export
distributional::generate

#' @export
distributional::hilo

#' @export
distributional::variance

# --- SHASH kernels (mgcv / Jones-Pewsey parameterisation) -----------------

shash_log_cdf <- function(y, mu, sigma, eps, delta, lower.tail = TRUE) {
  if (all(eps == 0, na.rm = TRUE) && all(delta == 1, na.rm = TRUE)) {
    # sinh(asinh(z)) = z: the Gaussian case skips the hyperbolic functions.
    return(stats::pnorm((y - mu) / sigma, lower.tail = lower.tail, log.p = TRUE))
  }
  z <- (y - mu) / (sigma * delta)
  stats::pnorm(sinh(delta * asinh(z) - eps), lower.tail = lower.tail, log.p = TRUE)
}

shash_quantile <- function(p, mu, sigma, eps, delta) {
  mu + sigma * delta * sinh((asinh(stats::qnorm(p)) + eps) / delta)
}

shash_log_density <- function(y, mu, sigma, eps, delta) {
  z <- (y - mu) / (sigma * delta)
  s <- sinh(delta * asinh(z) - eps)
  0.5 * log1p(s^2) - 0.5 * s^2 - 0.5 * log(2 * pi) - 0.5 * log1p(z^2) - log(sigma)
}

# Jones & Pewsey (2009) moments of S = sinh((asinh(W) + eps) / delta).
shash_pq <- function(q) {
  exp(0.25) / sqrt(8 * pi) * (besselK(0.25, (q + 1) / 2) + besselK(0.25, (q - 1) / 2))
}

shash_mean <- function(mu, sigma, eps, delta) {
  mu + sigma * delta * sinh(eps / delta) * shash_pq(1 / delta)
}

shash_variance <- function(mu, sigma, eps, delta) {
  e_s <- sinh(eps / delta) * shash_pq(1 / delta)
  e_s2 <- 0.5 * (cosh(2 * eps / delta) * shash_pq(2 / delta) - 1)
  (sigma * delta)^2 * (e_s2 - e_s^2)
}

# --- Log tails and log densities -------------------------------------------
# distributional's cdf() has no log.p argument and its own density methods
# ignore `log`, so tails and log densities come from these internal
# generics, which dispatch on an element (or an unpacked pseudo-element,
# see below). Every method is vectorised in its fields.

log_tail <- function(x, q, lower.tail = TRUE) {
  UseMethod("log_tail")
}

#' @export
log_tail.default <- function(x, q, lower.tail = TRUE) {
  p <- cdf(x, q)
  if (lower.tail) log(p) else log1p(-p)
}

#' @export
log_tail.dist_normal <- function(x, q, lower.tail = TRUE) {
  stats::pnorm(q, x[["mu"]], x[["sigma"]], lower.tail = lower.tail, log.p = TRUE)
}

#' @export
log_tail.dist_shash <- function(x, q, lower.tail = TRUE) {
  shash_log_cdf(q, x[["mu"]], x[["sigma"]], x[["eps"]], x[["delta"]], lower.tail)
}

#' @export
log_tail.dist_shash_draws <- function(x, q, lower.tail = TRUE) {
  draws_log_cdf(draws_fields(x, q), lower.tail)
}

#' @export
log_tail.dist_conditioned <- function(x, q, lower.tail = TRUE) {
  z <- elem_z(x[["dist"]], q)
  stats::pnorm((z - x[["m"]]) / x[["s"]], lower.tail = lower.tail, log.p = TRUE)
}

#' @export
log_tail.dist_warped <- function(x, q, lower.tail = TRUE) {
  w <- warped_parts(x)
  t <- transform_apply(w$tr, q)
  outside <- !is.na(q) & is.na(t)
  if (!lower.tail) {
    upper <- log_tail(w$inner, t, FALSE) - w$log_denom
    upper[outside] <- 0
    return(upper)
  }
  lower <- log_tail(w$inner, t, TRUE)
  if (is.finite(w$t_min)) {
    share <- exp(log_tail(w$inner, w$t_min, TRUE) - lower)
    share[!is.finite(share)] <- 0
    lower <- lower + log1p(-pmin(share, 1)) - w$log_denom
  }
  lower[outside] <- -Inf
  lower
}

log_dens <- function(x, at) {
  UseMethod("log_dens")
}

#' @export
log_dens.default <- function(x, at) {
  log(density(x, at))
}

#' @export
log_dens.dist_normal <- function(x, at) {
  stats::dnorm(at, x[["mu"]], x[["sigma"]], log = TRUE)
}

#' @export
log_dens.dist_shash <- function(x, at) {
  shash_log_density(at, x[["mu"]], x[["sigma"]], x[["eps"]], x[["delta"]])
}

#' @export
log_dens.dist_shash_draws <- function(x, at) {
  draws_log_density(draws_fields(x, at))
}

#' @export
log_dens.dist_conditioned <- function(x, at) {
  z <- elem_z(x[["dist"]], at)
  log_dens(x[["dist"]], at) +
    stats::dnorm((z - x[["m"]]) / x[["s"]], log = TRUE) - log(x[["s"]]) -
    stats::dnorm(z, log = TRUE)
}

#' @export
log_dens.dist_warped <- function(x, at) {
  w <- warped_parts(x)
  t <- transform_apply(w$tr, at)
  outside <- !is.na(at) & is.na(t)
  at[outside] <- NA_real_
  out <- log_dens(w$inner, t) + transform_log_deriv(w$tr, at) - w$log_denom
  out[outside] <- -Inf
  out
}

# Normal score of `q` under an element (or unpacked vector) `x`, taken from
# whichever tail is smaller so that z = 8, 10, 40 stay distinct.
elem_z <- function(x, q) {
  scores_from_log_tails(log_tail(x, q, TRUE), log_tail(x, q, FALSE))$z
}

# cdf() and density() of the three classes are the exponentiated log
# tails and log densities; the same bodies serve every class.
cdf_from_log_tail <- function(x, q, ...) {
  exp(log_tail(x, q))
}

density_from_log_dens <- function(x, at, ..., log = FALSE) {
  out <- log_dens(x, at)
  if (log) out else exp(out)
}

# --- dist_shash element methods (vectorised in every argument) ------------

#' @export
format.dist_shash <- function(x, digits = 2, ...) {
  sprintf(
    "SHASH(%s, %s, %s, %s)",
    format(x[["mu"]], digits = digits, ...), format(x[["sigma"]], digits = digits, ...),
    format(x[["eps"]], digits = digits, ...), format(x[["delta"]], digits = digits, ...)
  )
}

#' @method cdf dist_shash
#' @export
cdf.dist_shash <- cdf_from_log_tail

#' @method density dist_shash
#' @export
density.dist_shash <- density_from_log_dens

#' @export
quantile.dist_shash <- function(x, p, ...) {
  shash_quantile(p, x[["mu"]], x[["sigma"]], x[["eps"]], x[["delta"]])
}

#' @export
generate.dist_shash <- function(x, times, ...) {
  quantile(x, stats::runif(times))
}

#' @export
mean.dist_shash <- function(x, ...) {
  shash_mean(x[["mu"]], x[["sigma"]], x[["eps"]], x[["delta"]])
}

#' @export
variance.dist_shash <- function(x, ...) {
  shash_variance(x[["mu"]], x[["sigma"]], x[["eps"]], x[["delta"]])
}

# --- dist_shash_draws: equal-weight mixture over coefficient draws --------
# On an element every field is a draw vector of length K; on an unpacked
# vector every field is an n x K matrix. `draws_fields()` returns the
# fields as matrices with one row per evaluation point so the kernels
# apply to both cases: rows x draws, with `q` recycled down the rows.

draws_fields <- function(x, q) {
  f <- unclass(x)[c("mu", "sigma", "eps", "delta")]
  if (is.matrix(f$mu)) {
    q <- rep_len(q, nrow(f$mu))
  } else {
    k <- length(f$mu)
    f <- lapply(f, function(v) matrix(v, length(q), k, byrow = TRUE))
  }
  c(list(q = q), f)
}

#' @export
format.dist_shash_draws <- function(x, digits = 2, ...) {
  gaussian <- isTRUE(all(x[["eps"]] == 0) && all(x[["delta"]] == 1))
  sprintf(
    "%s[%d](%s, %s)", if (gaussian) "N" else "SHASH", length(x[["mu"]]),
    format(mean(x[["mu"]]), digits = digits, ...),
    format(mean(x[["sigma"]]), digits = digits, ...)
  )
}

# Row-wise log-mean-exp of a matrix.
log_mean_exp <- function(m) {
  mx <- m[cbind(seq_len(nrow(m)), max.col(m, ties.method = "first"))]
  ok <- is.finite(mx)
  out <- mx
  out[ok] <- mx[ok] + log(rowMeans(exp(m[ok, , drop = FALSE] - mx[ok])))
  out
}

draws_log_cdf <- function(a, lower.tail = TRUE) {
  log_mean_exp(shash_log_cdf(a$q, a$mu, a$sigma, a$eps, a$delta, lower.tail))
}

draws_log_density <- function(a) {
  log_mean_exp(shash_log_density(a$q, a$mu, a$sigma, a$eps, a$delta))
}

#' @method cdf dist_shash_draws
#' @export
cdf.dist_shash_draws <- cdf_from_log_tail

#' @method density dist_shash_draws
#' @export
density.dist_shash_draws <- density_from_log_dens

# Quantile by safeguarded Newton iteration on the mixture CDF, started
# from the mean component quantile and bracketed by the extreme component
# quantiles; only unconverged rows are re-evaluated.
#' @export
quantile.dist_shash_draws <- function(x, p, ..., iter = 30L) {
  a <- draws_fields(x, p)
  p <- a$q
  q_comp <- shash_quantile(p, a$mu, a$sigma, a$eps, a$delta)
  i <- seq_len(nrow(q_comp))
  lo <- q_comp[cbind(i, max.col(-q_comp, ties.method = "first"))]
  hi <- q_comp[cbind(i, max.col(q_comp, ties.method = "first"))]
  cur <- rowMeans(q_comp)
  out <- ifelse(p <= 0, lo, ifelse(p >= 1, hi, NA_real_))
  active <- is.finite(lo) & is.finite(hi) & is.na(out)
  for (it in seq_len(iter)) {
    idx <- which(active)
    if (!length(idx)) {
      break
    }
    sub <- lapply(a[c("mu", "sigma", "eps", "delta")], function(f) f[idx, , drop = FALSE])
    sub$q <- cur[idx]
    err <- exp(draws_log_cdf(sub)) - p[idx]
    below <- err < 0
    lo[idx] <- ifelse(below, cur[idx], lo[idx])
    hi[idx] <- ifelse(below, hi[idx], cur[idx])
    nxt <- cur[idx] - err / exp(draws_log_density(sub))
    bad <- !is.finite(nxt) | nxt < lo[idx] | nxt > hi[idx]
    nxt[bad] <- ((lo[idx] + hi[idx]) / 2)[bad]
    tol <- 1e-9 * pmax(1, abs(nxt))
    done <- err == 0 | abs(nxt - cur[idx]) < tol | (hi[idx] - lo[idx]) < tol
    cur[idx] <- nxt
    active[idx[done]] <- FALSE
  }
  out[is.na(out)] <- cur[is.na(out)]
  out
}

#' @export
generate.dist_shash_draws <- function(x, times, ...) {
  a <- draws_fields(x, numeric(1))
  n <- if (is.matrix(x[["mu"]])) nrow(x[["mu"]]) else 1L
  k <- ncol(a$mu)
  j <- sample.int(k, n * times, replace = TRUE)
  pick <- function(f) f[cbind(rep_len(seq_len(n), n * times), j)]
  out <- shash_quantile(stats::runif(n * times), pick(a$mu), pick(a$sigma),
                        pick(a$eps), pick(a$delta))
  if (n == 1L) out else matrix(out, n, times)
}

#' @export
mean.dist_shash_draws <- function(x, ...) {
  a <- draws_fields(x, numeric(1))
  rowMeans(shash_mean(a$mu, a$sigma, a$eps, a$delta))
}

#' @export
variance.dist_shash_draws <- function(x, ...) {
  a <- draws_fields(x, numeric(1))
  m <- shash_mean(a$mu, a$sigma, a$eps, a$delta)
  v <- shash_variance(a$mu, a$sigma, a$eps, a$delta)
  rowMeans(v + m^2) - rowMeans(m)^2
}

# --- dist_conditioned: marginal mapped through a conditional normal score --

#' @export
format.dist_conditioned <- function(x, digits = 2, ...) {
  sprintf(
    "%s | N(%s, %s)", format(x[["dist"]], digits = digits, ...),
    format(x[["m"]], digits = digits, ...), format(x[["s"]], digits = digits, ...)
  )
}

#' @method cdf dist_conditioned
#' @export
cdf.dist_conditioned <- cdf_from_log_tail

#' @method density dist_conditioned
#' @export
density.dist_conditioned <- density_from_log_dens

#' @export
quantile.dist_conditioned <- function(x, p, ...) {
  quantile(x[["dist"]], stats::pnorm(x[["m"]] + x[["s"]] * stats::qnorm(p)))
}

#' @export
generate.dist_conditioned <- function(x, times, ...) {
  quantile(x, stats::runif(times))
}

# --- dist_warped: a fit on h(y), read on the y scale -----------------------
# `ref_spec(transform = ...)` fits the model to t = h(y) for a strictly
# increasing Box-Cox map h. The predictive for Y is the push-forward of the
# predictive for T: tails carry over unchanged, the density picks up
# log h'(y), and quantiles are back-transformed, so a warped fit is scored
# in the units of `y` and its log density is comparable with an untransformed
# fit's. When the image of h is bounded below (any power between 0 and 1),
# the fitted distribution puts mass below h(0) that no response value can
# produce; that mass is renormalised away, which keeps the predictive a
# proper density on (0, Inf) rather than one that integrates to slightly
# less than one.

dist_warped <- function(dist, lambda) {
  distributional::new_dist(
    dist = vctrs::vec_data(dist),
    lambda = rep_len(as.numeric(lambda), length(dist)),
    class = "dist_warped"
  )
}

# The transform, the inner (unwarped) distribution, and the log of the
# probability the inner distribution puts on the image of h.
warped_parts <- function(x) {
  tr <- as_transform(x[["lambda"]][[1L]])
  inner <- x[["dist"]]
  t_min <- transform_t_min(tr)
  list(
    tr = tr, inner = inner, t_min = t_min,
    log_denom = if (is.finite(t_min)) log_tail(inner, t_min, FALSE) else 0
  )
}

#' @export
format.dist_warped <- function(x, digits = 2, ...) {
  sprintf(
    "%s{%s}", as_transform(x[["lambda"]][[1L]])$name,
    format(x[["dist"]], digits = digits, ...)
  )
}

#' @method cdf dist_warped
#' @export
cdf.dist_warped <- cdf_from_log_tail

#' @method density dist_warped
#' @export
density.dist_warped <- density_from_log_dens

#' @export
quantile.dist_warped <- function(x, p, ...) {
  w <- warped_parts(x)
  if (is.finite(w$t_min)) {
    p <- exp(log_tail(w$inner, w$t_min, TRUE)) + p * exp(w$log_denom)
  }
  transform_invert(w$tr, quantile(w$inner, p))
}

#' @export
generate.dist_warped <- function(x, times, ...) {
  quantile(x, stats::runif(times))
}

# The response-scale moments have no closed form; they come from the same
# K-atom quantile discretisation `crps_from_dist()` uses.
warped_atoms <- function(x, K = 199L) {
  n <- max(length(x[["lambda"]]), 1L)
  p <- (seq_len(K) - 0.5) / K
  matrix(vapply(p, function(pk) as.numeric(quantile(x, rep(pk, n))), numeric(n)),
         n, K)
}

#' @export
mean.dist_warped <- function(x, ...) {
  rowMeans(warped_atoms(x))
}

#' @export
variance.dist_warped <- function(x, ...) {
  a <- warped_atoms(x)
  rowMeans(a^2) - rowMeans(a)^2
}

# --- Vectorised evaluation of a whole distribution vector -----------------
# A homogeneous distribution vector is unpacked into one pseudo-element
# whose fields are vectors (or n x K matrices, for draw mixtures). The element
# methods above are vectorised, so calling them on the pseudo-element
# scores every observation at once instead of dispatching per element.
# Missing elements (dist_missing) and heterogeneous vectors fall back to
# per-element dispatch.

dist_unpack <- function(d) {
  el <- vctrs::vec_data(d)
  if (!length(el) || any(vapply(el, is.null, logical(1)))) {
    return(NULL)
  }
  cls <- unique(vapply(el, function(e) class(e)[[1L]], ""))
  if (length(cls) != 1L) {
    return(NULL)
  }
  # Some distributions carry functions rather than parameters: a
  # `dist_transformed` holds its transform and inverse. Those fields cannot be
  # stacked into parallel vectors, so unpacking them corrupts the result. Fall
  # back to the elementwise path, which dispatches on the distribution itself.
  if (any(vapply(el[[1L]], is.function, logical(1)))) {
    return(NULL)
  }
  fields <- list()
  for (f in names(el[[1L]])) {
    v <- lapply(el, `[[`, f)
    fields[[f]] <- if (is.list(v[[1L]])) {
      inner <- dist_unpack(vctrs::new_vctr(v, vars = NULL, class = "distribution"))
      if (is.null(inner)) {
        return(NULL)
      }
      inner
    } else if (length(v[[1L]]) == 1L && !inherits(el[[1L]], "dist_shash_draws")) {
      unlist(v, use.names = FALSE)
    } else {
      do.call(rbind, v)
    }
  }
  structure(fields, class = class(el[[1L]]))
}

# Apply an element function pairwise to a distribution vector and a vector
# `arg` (recycled to a common length). `unpacked` may supply
# `dist_unpack(d)` when the caller evaluates several functions on the
# same vector.
dist_eval <- function(d, fun, arg, ..., unpacked = NULL) {
  n <- max(length(d), length(arg))
  if (!n) {
    return(numeric())
  }
  d <- vctrs::vec_recycle(d, n)
  arg <- rep_len(as.numeric(arg), n)
  u <- unpacked %||% dist_unpack(d)
  if (!is.null(u)) {
    return(as.numeric(fun(u, arg, ...)))
  }
  el <- vctrs::vec_data(d)
  mapply(function(e, a) if (is.null(e)) NA_real_ else fun(e, a, ...), el, arg)
}

dist_quantile <- function(d, p) {
  dist_eval(d, quantile, p)
}

# Normal score of y under d, taken from the smaller tail.
dist_z <- function(d, y) {
  dist_eval(d, elem_z, y)
}

# --- Scores from log tails -------------------------------------------------

# Centile, z, and two-sided tail from log lower/upper tail probabilities.
# z is taken from whichever tail is smaller, so z = 8, 10, 40 are distinct.
scores_from_log_tails <- function(log_lower, log_upper) {
  use_lower <- is.na(log_upper) | (!is.na(log_lower) & log_lower <= log_upper)
  z <- ifelse(
    use_lower,
    stats::qnorm(log_lower, log.p = TRUE, lower.tail = TRUE),
    stats::qnorm(log_upper, log.p = TRUE, lower.tail = FALSE)
  )
  log_tail <- pmin(log(2) + pmin(log_lower, log_upper), 0)
  list(
    centile = as.numeric(ifelse(use_lower, exp(log_lower), -expm1(log_upper))),
    z = as.numeric(z),
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

# Per-observation CRPS: closed form for a Gaussian predictive, otherwise
# the CRPS of the K-atom quantile discretisation of the predictive (atoms
# at the p = (k - 1/2)/K quantiles, equal weights). For a distribution with
# finite first moment the discretisation error is O(1/K); K = 199 keeps it
# well below the sampling noise of any held-out comparison.
crps_from_dist <- function(dist, y, K = 199L) {
  y <- rep_len(as.numeric(y), length(dist))
  u <- dist_unpack(dist)
  if (inherits(u, "dist_normal")) {
    return(crps_norm(y, u$mu, u$sigma))
  }
  n <- length(dist)
  p <- (seq_len(K) - 0.5) / K
  q <- vapply(p, function(pk) dist_eval(dist, quantile, rep(pk, n), unpacked = u), numeric(n))
  q <- matrix(q, n, K)
  # sum_{j,k} |q_j - q_k| over sorted atoms is 2 * sum_k (2k - K - 1) q_(k)
  w <- 2 * seq_len(K) - K - 1
  rowMeans(abs(q - y)) - drop(q %*% w) / K^2
}

crps_norm <- function(y, mu, sigma) {
  z <- (y - mu) / sigma
  sigma * (z * (2 * stats::pnorm(z) - 1) + 2 * stats::dnorm(z) - 1 / sqrt(pi))
}
