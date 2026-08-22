#' Predictive distributions
#'
#' Reference models return predictive distributions as
#' [distributional](https://pkg.mitchelloharawild.com/distributional/)
#' vectors, so `cdf()`, `quantile()`, `density()`, `generate()`, `mean()`,
#' `variance()`, `hilo()`, tibble columns, and `ggdist` layers all work.
#' Gaussian fits use [distributional::dist_normal()]; SHASH fits use
#' `dist_shash()`; `uncertainty = "total"` predictions from location-scale
#' and SHASH fits are equal-weight mixtures over coefficient draws
#' (`dist_shash_mc()`); history-conditioned forecasts are
#' `dist_conditioned()`.
#'
#' @details
#' `dist_shash()` follows the `mgcv::shash()` parameterisation. With
#' \eqn{z = (y-\mu)/(\sigma\delta)},
#' \eqn{F(y) = \Phi(\sinh(\delta\,\mathrm{asinh}(z)-\epsilon))}; `eps = 0`,
#' `delta = 1` is the Gaussian \eqn{N(\mu, \sigma^2)}.
#'
#' `dist_shash_mc()` holds, for every observation, one SHASH parameter set
#' per coefficient draw (a matrix column per draw). Its CDF and density are
#' averages over the draws computed in log space, its quantile is found by
#' safeguarded Newton iteration on the mixture CDF, and `generate()` samples a draw and then
#' the component. A Gaussian draw mixture is the special case
#' `eps = 0`, `delta = 1`.
#'
#' `dist_conditioned()` maps a marginal predictive distribution through a
#' conditional normal score: \eqn{F_*(y) = \Phi((z(y)-m)/s)} with
#' \eqn{z(y) = \Phi^{-1}(F(y))}, so `m = 0`, `s = 1` recovers the marginal.
#'
#' Tail probabilities of all three classes are evaluated in log space, so
#' scores for observations 40 standard deviations out remain finite and
#' distinct (see [as_scores()]).
#'
#' @param mu,sigma,eps,delta Location, scale (> 0), skewness, and tail
#'   weight (> 0). For `dist_shash()` numeric vectors recycled to a common
#'   length; for `dist_shash_mc()` matrices with one row per observation
#'   and one column per draw (a vector is treated as draw-invariant).
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

#' @rdname dist_shash
#' @export
dist_shash_mc <- function(mu, sigma, eps = 0, delta = 1) {
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
    class = "dist_shash_mc"
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
log_tail.dist_shash_mc <- function(x, q, lower.tail = TRUE) {
  mc_log_cdf(mc_fields(x, q), lower.tail)
}

#' @export
log_tail.dist_conditioned <- function(x, q, lower.tail = TRUE) {
  z <- scores_from_log_tails(log_tail(x[["dist"]], q, TRUE), log_tail(x[["dist"]], q, FALSE))$z
  stats::pnorm((z - x[["m"]]) / x[["s"]], lower.tail = lower.tail, log.p = TRUE)
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
log_dens.dist_shash_mc <- function(x, at) {
  mc_log_density(mc_fields(x, at))
}

#' @export
log_dens.dist_conditioned <- function(x, at) {
  z <- scores_from_log_tails(log_tail(x[["dist"]], at, TRUE), log_tail(x[["dist"]], at, FALSE))$z
  log_dens(x[["dist"]], at) +
    stats::dnorm((z - x[["m"]]) / x[["s"]], log = TRUE) - log(x[["s"]]) -
    stats::dnorm(z, log = TRUE)
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

#' @export
cdf.dist_shash <- function(x, q, ...) {
  exp(log_tail(x, q))
}

#' @export
density.dist_shash <- function(x, at, ..., log = FALSE) {
  out <- log_dens(x, at)
  if (log) out else exp(out)
}

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

# --- dist_shash_mc: equal-weight mixture over coefficient draws -----------
# On an element every field is a draw vector of length K; on an unpacked
# vector every field is an n x K matrix. `mc_fields()` returns the fields
# as matrices with one row per evaluation point so the kernels apply to
# both cases: rows x draws, with `q` recycled down the rows.

mc_fields <- function(x, q) {
  f <- unclass(x)[c("mu", "sigma", "eps", "delta")]
  if (is.matrix(f$mu)) {
    q <- rep_len(q, nrow(f$mu))
  } else {
    k <- length(f$mu)
    f <- lapply(f, function(v) matrix(v, length(q), k, byrow = TRUE))
  }
  c(list(q = q), f)
}

mc_draws <- function(x) {
  if (is.matrix(x[["mu"]])) ncol(x[["mu"]]) else length(x[["mu"]])
}

#' @export
format.dist_shash_mc <- function(x, digits = 2, ...) {
  gaussian <- isTRUE(all(x[["eps"]] == 0) && all(x[["delta"]] == 1))
  sprintf(
    "%s[%d](%s, %s)", if (gaussian) "N" else "SHASH", mc_draws(x),
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

mc_log_cdf <- function(a, lower.tail = TRUE) {
  log_mean_exp(shash_log_cdf(a$q, a$mu, a$sigma, a$eps, a$delta, lower.tail))
}

mc_log_density <- function(a) {
  log_mean_exp(shash_log_density(a$q, a$mu, a$sigma, a$eps, a$delta))
}

#' @export
cdf.dist_shash_mc <- function(x, q, ...) {
  exp(log_tail(x, q))
}

#' @export
density.dist_shash_mc <- function(x, at, ..., log = FALSE) {
  out <- log_dens(x, at)
  if (log) out else exp(out)
}

# Quantile by safeguarded Newton iteration on the mixture CDF, started
# from the mean component quantile and bracketed by the extreme component
# quantiles; only unconverged rows are re-evaluated.
#' @export
quantile.dist_shash_mc <- function(x, p, ..., iter = 30L) {
  a <- mc_fields(x, p)
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
    err <- exp(mc_log_cdf(sub)) - p[idx]
    below <- err < 0
    lo[idx] <- ifelse(below, cur[idx], lo[idx])
    hi[idx] <- ifelse(below, hi[idx], cur[idx])
    nxt <- cur[idx] - err / exp(mc_log_density(sub))
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
generate.dist_shash_mc <- function(x, times, ...) {
  a <- mc_fields(x, numeric(1))
  n <- if (is.matrix(x[["mu"]])) nrow(x[["mu"]]) else 1L
  k <- ncol(a$mu)
  j <- sample.int(k, n * times, replace = TRUE)
  pick <- function(f) f[cbind(rep_len(seq_len(n), n * times), j)]
  out <- shash_quantile(stats::runif(n * times), pick(a$mu), pick(a$sigma),
                        pick(a$eps), pick(a$delta))
  if (n == 1L) out else matrix(out, n, times)
}

#' @export
mean.dist_shash_mc <- function(x, ...) {
  a <- mc_fields(x, numeric(1))
  rowMeans(shash_mean(a$mu, a$sigma, a$eps, a$delta))
}

#' @export
variance.dist_shash_mc <- function(x, ...) {
  a <- mc_fields(x, numeric(1))
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

#' @export
cdf.dist_conditioned <- function(x, q, ...) {
  exp(log_tail(x, q))
}

#' @export
quantile.dist_conditioned <- function(x, p, ...) {
  quantile(x[["dist"]], stats::pnorm(x[["m"]] + x[["s"]] * stats::qnorm(p)))
}

#' @export
density.dist_conditioned <- function(x, at, ..., log = FALSE) {
  out <- log_dens(x, at)
  if (log) out else exp(out)
}

#' @export
generate.dist_conditioned <- function(x, times, ...) {
  quantile(x, stats::runif(times))
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
  fields <- list()
  for (f in names(el[[1L]])) {
    v <- lapply(el, `[[`, f)
    fields[[f]] <- if (is.list(v[[1L]])) {
      inner <- dist_unpack(vctrs::new_vctr(v, vars = NULL, class = "distribution"))
      if (is.null(inner)) {
        return(NULL)
      }
      inner
    } else if (length(v[[1L]]) == 1L && !inherits(el[[1L]], "dist_shash_mc")) {
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

dist_log_tail <- function(d, q, lower.tail = TRUE) {
  dist_eval(d, log_tail, q, lower.tail = lower.tail)
}

dist_cdf <- function(d, q) {
  exp(dist_log_tail(d, q))
}

dist_log_density <- function(d, y) {
  dist_eval(d, log_dens, y)
}

dist_quantile <- function(d, p) {
  dist_eval(d, quantile, p)
}

# n x times matrix of simulated values.
dist_generate <- function(d, times) {
  u <- dist_unpack(d)
  if (is.null(u) || length(d) <= 1L) {
    return(do.call(rbind, lapply(generate(d, times), function(g) g %||% rep(NA_real_, times))))
  }
  if (inherits(u, "dist_shash_mc")) {
    return(generate(u, times))
  }
  matrix(quantile(u, stats::runif(length(d) * times)), length(d), times)
}

# --- Scores ---------------------------------------------------------------

#' Scores derived from a predictive CDF
#'
#' @param distribution A `distribution` vector (see [dist_shash]).
#' @param y Observed values, recycled against `distribution`.
#' @return A tibble of centiles, Z-scores, tails, residuals, and log
#'   densities. Tails are evaluated in log space so that `z = 8, 10, 40`
#'   remain distinct. There is no abnormality column.
#' @examples
#' as_scores(distributional::dist_normal(0, 1), c(-2, 0, 2))
#' @export
as_scores <- function(distribution, y) {
  n <- max(length(distribution), length(y))
  y <- rep_len(as.numeric(y), n)
  distribution <- vctrs::vec_recycle(distribution, n)
  u <- dist_unpack(distribution)
  tails <- scores_from_log_tails(
    dist_eval(distribution, log_tail, y, lower.tail = TRUE, unpacked = u),
    dist_eval(distribution, log_tail, y, lower.tail = FALSE, unpacked = u)
  )
  med <- dist_eval(distribution, quantile, rep(0.5, n), unpacked = u)
  tibble::tibble(
    observed = y,
    median = med,
    centile = tails$centile,
    z = tails$z,
    tail_prob = tails$tail_prob,
    tail_surprisal = tails$tail_surprisal,
    residual = y - med,
    log_density = dist_eval(distribution, log_dens, y, unpacked = u)
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
  log_tail <- pmin(log(2) + pmin(log_lower, log_upper), 0)
  list(
    centile = ifelse(use_lower, exp(log_lower), -expm1(log_upper)),
    z = z,
    tail_prob = exp(log_tail),
    tail_surprisal = -log_tail
  )
}

# Normal score of y under d, taken from the smaller tail.
dist_z <- function(d, y) {
  scores_from_log_tails(dist_log_tail(d, y, TRUE), dist_log_tail(d, y, FALSE))$z
}

# Log tails recovered from a z-score (exact inverse of scores_from_log_tails).
log_tails_from_z <- function(z) {
  list(
    lower = stats::pnorm(z, log.p = TRUE, lower.tail = TRUE),
    upper = stats::pnorm(z, log.p = TRUE, lower.tail = FALSE)
  )
}
