test_that("Gaussian F(Q(p)) recovers p", {
  d <- norm_dist("gaussian", location = c(-1, 0, 2), scale = c(0.5, 1, 1.5))
  p <- c(0.05, 0.5, 0.95)
  q <- dist_quantile(d, p)
  expect_equal(cdf(d, q), p, tolerance = 1e-8)
})

test_that("SHASH F(Q(p)) recovers p", {
  d <- norm_dist("shash", location = 0, scale = 1.2, skew = 0.4, tail = 0.8)
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  q <- dist_quantile(d, p)
  expect_equal(as.numeric(cdf(d, q)), p, tolerance = 1e-7)
})

test_that("SHASH reduces to Gaussian when skew=0 and tail=1", {
  y <- c(-2, 0, 1.3)
  g <- norm_dist("gaussian", location = 0.2, scale = 1.1)
  s <- norm_dist("shash", location = 0.2, scale = 1.1, skew = 0, tail = 1)
  expect_equal(cdf(s, y), cdf(g, y), tolerance = 1e-10)
  expect_equal(log_density(s, y), log_density(g, y), tolerance = 1e-10)
  expect_equal(center(s), center(g), tolerance = 1e-10)
})

test_that("simulated draws recover declared quantiles", {
  set.seed(12)
  d <- norm_dist("gaussian", location = 3, scale = 2)
  draws <- draw(d, 4000)
  expect_equal(stats::median(draws), 3, tolerance = 0.15)
  expect_equal(stats::sd(draws), 2, tolerance = 0.15)
})

test_that("density integrates to one for SHASH", {
  d <- norm_dist("shash", location = 1, scale = 0.8, skew = -0.5, tail = 1.2)
  ys <- seq(-8, 10, length.out = 2001)
  dens <- exp(log_density(d, ys))
  integ <- sum(diff(ys) * (dens[-1] + dens[-length(dens)]) / 2)
  expect_equal(integ, 1, tolerance = 0.02)
})

test_that("conditional forecast distribution obeys F(Q(p))=p", {
  d <- condition_norm_dist(
    norm_dist("gaussian", location = 0, scale = 1),
    m = 0.4,
    s = 0.7
  )
  p <- c(0.1, 0.5, 0.9)
  expect_equal(as.numeric(cdf(d, quantile(d, p))), p, tolerance = 1e-8)
})
