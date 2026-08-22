test_that("Gaussian F(Q(p)) recovers p and scores are (y - mu) / sigma", {
  d <- distributional::dist_normal(c(-1, 0, 2), c(0.5, 1, 1.5))
  p <- c(0.05, 0.5, 0.95)
  expect_equal(dist_cdf(d, dist_quantile(d, p)), p, tolerance = 1e-8)
  y <- c(-1, 0, 0.5, 2)
  sc <- as_scores(distributional::dist_normal(0.25, 1.5), y)
  expect_equal(sc$z, (y - 0.25) / 1.5, tolerance = 1e-10)
  expect_equal(sc$centile, stats::pnorm(y, 0.25, 1.5), tolerance = 1e-10)
})

test_that("SHASH F(Q(p)) recovers p and matches mgcv::shash", {
  d <- dist_shash(0, 1.2, 0.4, 0.8)
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  q <- quantile(d, p)[[1]]
  expect_equal(cdf(d, q)[[1]], p, tolerance = 1e-7)
  fam <- mgcv::shash()
  y <- c(-3, -1, 0, 0.7, 2.5)
  # mgcv's cdf takes (mu, tau = log sigma, eps, phi = log delta)
  expect_equal(
    cdf(d, y)[[1]],
    fam$cdf(y, matrix(c(0, log(1.2), 0.4, log(0.8)), 1, 4), 1, 1, FALSE),
    tolerance = 1e-7
  )
  expect_equal(
    quantile(d, p)[[1]],
    fam$qf(p, matrix(c(0, log(1.2), 0.4, log(0.8)), 1, 4), 1, 1),
    tolerance = 1e-7
  )
})

test_that("SHASH reduces to Gaussian when eps = 0 and delta = 1", {
  y <- c(-2, 0, 1.3)
  g <- distributional::dist_normal(0.2, 1.1)
  s <- dist_shash(0.2, 1.1, 0, 1)
  expect_equal(cdf(s, y), cdf(g, y), tolerance = 1e-10)
  expect_equal(density(s, y, log = TRUE), density(g, y, log = TRUE), tolerance = 1e-10)
  expect_equal(median(s), median(g), tolerance = 1e-10)
  expect_equal(mean(s), 0.2, tolerance = 1e-10)
  expect_equal(variance(s), 1.1^2, tolerance = 1e-10)
})

test_that("SHASH mean and variance agree with simulation and integration", {
  d <- dist_shash(1, 0.8, -0.5, 1.2)
  ys <- seq(-12, 12, length.out = 8001)
  dens <- density(d, ys)[[1]]
  integ <- function(f) sum(diff(ys) * (f[-1] + f[-length(f)]) / 2)
  expect_equal(integ(dens), 1, tolerance = 1e-4)
  m <- integ(dens * ys)
  expect_equal(mean(d), m, tolerance = 1e-4)
  expect_equal(variance(d), integ(dens * (ys - m)^2), tolerance = 1e-3)
  set.seed(3)
  g <- generate(d, 40000)[[1]]
  expect_equal(mean(g), mean(d), tolerance = 0.05)
  expect_equal(stats::var(g), variance(d), tolerance = 0.1)
})

test_that("draw mixture cdf is the mean of component cdfs and quantile inverts it", {
  set.seed(4)
  mu <- matrix(stats::rnorm(12, sd = 0.5), 3, 4)
  sig <- matrix(stats::runif(12, 0.8, 1.2), 3, 4)
  m <- dist_shash_draws(mu, sig, eps = 0.3, delta = 0.9)
  y <- c(-1, 0.2, 1.5)
  comp <- sapply(1:4, function(j) cdf(dist_shash(mu[, j], sig[, j], 0.3, 0.9), y[1]))
  expect_equal(cdf(m, y[1]), rowMeans(comp), tolerance = 1e-12)
  p <- c(0.02, 0.5, 0.97)
  expect_equal(dist_cdf(m, dist_quantile(m, p)), p, tolerance = 1e-8)
  # moments: mixture mean is the mean of component means
  mm <- sapply(1:4, function(j) mean(dist_shash(mu[, j], sig[, j], 0.3, 0.9)))
  expect_equal(mean(m), rowMeans(mm), tolerance = 1e-10)
  set.seed(5)
  g <- generate(dist_unpack(m), 20000)
  expect_equal(rowMeans(g), mean(m), tolerance = 0.05)
  expect_equal(apply(g, 1, stats::var), variance(m), tolerance = 0.1)
  # a Gaussian mixture with identical draws is the Gaussian
  g1 <- dist_shash_draws(matrix(0.3, 2, 5), 1.5)
  expect_equal(cdf(g1, 1), rep(stats::pnorm(1, 0.3, 1.5), 2))
  expect_equal(quantile(g1, 0.9), rep(stats::qnorm(0.9, 0.3, 1.5), 2), tolerance = 1e-8)
})

test_that("dist_conditioned at m = 0, s = 1 is the marginal and obeys F(Q(p)) = p", {
  marg <- c(distributional::dist_normal(0, 1), dist_shash(1, 2, 0.2, 1.1))
  same <- dist_conditioned(marg, 0, 1)
  y <- c(-0.5, 2.5)
  expect_equal(dist_cdf(same, y), dist_cdf(marg, y), tolerance = 1e-10)
  expect_equal(dist_log_density(same, y), dist_log_density(marg, y), tolerance = 1e-8)
  d <- dist_conditioned(distributional::dist_normal(0, 1), m = 0.4, s = 0.7)
  p <- c(0.1, 0.5, 0.9)
  expect_equal(dist_cdf(d, dist_quantile(d, p)), p, tolerance = 1e-8)
  expect_equal(as.numeric(quantile(d, 0.5)), stats::qnorm(stats::pnorm(0.4)))
  ys <- seq(-8, 8, length.out = 4001)
  dens <- density(d, ys)[[1]]
  expect_equal(sum(diff(ys) * (dens[-1] + dens[-length(dens)]) / 2), 1, tolerance = 1e-4)
  expect_error(dist_conditioned(marg, 0, -1), "positive")
})

test_that("log-space tails keep z = 40 exact for every class", {
  y <- c(-40, -10, 0, 10, 40)
  dists <- list(
    normal = distributional::dist_normal(0, 1),
    shash = dist_shash(0, 1, 0, 1),
    mc = dist_shash_draws(matrix(0, 1, 3), 1),
    cond = dist_conditioned(distributional::dist_normal(0, 1), 0, 1)
  )
  for (d in dists) {
    expect_equal(as_scores(d, y)$z, y, tolerance = 1e-8)
  }
  sc <- as_scores(dists$normal, y)
  expect_true(all(diff(sc$tail_surprisal[3:5]) > 0))
  expect_equal(sc$tail_surprisal[[5]], -(log(2) + stats::pnorm(-40, log.p = TRUE)))
  # a skewed, heavy-tailed SHASH keeps far-tail scores finite and monotone
  scs <- as_scores(dist_shash(0, 1, 0.6, 0.85), c(-60, -20, 0, 20, 60, 200))
  expect_true(all(is.finite(scs$z)))
  expect_true(all(diff(scs$z) > 0))
  expect_gt(scs$z[[6]], scs$z[[5]] + 0.5)
})

test_that("distribution vectors slice, print, and sit in tibbles", {
  d <- c(dist_shash(0, 1, 0.2, 0.9), distributional::dist_normal(1, 2))
  expect_match(format(d)[[1]], "SHASH")
  tb <- tibble::tibble(d = d, y = c(0, 1))
  expect_equal(nrow(tb), 2L)
  expect_equal(dist_cdf(d[2], 1), 0.5)
  m <- dist_shash_draws(matrix(c(0, 1), 2, 3), 1)
  expect_match(format(m)[[1]], "N\\[3\\]")
  expect_equal(length(m[2]), 1L)
  expect_equal(cdf(m[2], 1), 0.5)
  na <- distributional::dist_missing(2)
  expect_true(all(is.na(as_scores(na, c(0, 1))$z)))
  expect_true(all(is.na(dist_cdf(c(d, na), c(0, 0, 0, 0))[3:4])))
})

test_that("ggdist draws slabs from dist_shash", {
  skip_if_not_installed("ggdist")
  df <- data.frame(g = c("a", "b"), d = dist_shash(c(0, 1), c(1, 1.5), c(0.3, -0.3), 1))
  p <- ggplot2::ggplot(df, ggplot2::aes(xdist = d, y = g)) + ggdist::stat_slab()
  built <- ggplot2::ggplot_build(p)
  expect_gt(nrow(built$data[[1]]), 10L)
  expect_true(all(is.finite(built$data[[1]]$x)))
})
