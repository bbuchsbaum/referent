test_that("scores work on a distributional::dist_transformed", {
  # dist_transformed carries functions in its fields, which cannot be stacked
  # into the vectorised pseudo-element; unpacking it used to abort with
  # "invalid type/length (symbol/0)".
  d <- distributional::dist_transformed(
    distributional::dist_normal(log(1000), 0.4),
    transform = exp, inverse = log
  )
  expect_null(dist_unpack(d))

  crps <- crps_from_dist(d, 1200)
  expect_true(is.finite(crps))
  # Monte Carlo reference on the response scale.
  set.seed(1)
  s <- exp(stats::rnorm(2e5, log(1000), 0.4))
  ref <- mean(abs(s - 1200)) - 0.5 * mean(abs(sample(s) - sample(s)))
  expect_equal(crps, ref, tolerance = 0.02)

  # the Gaussian fast path is untouched
  expect_equal(
    crps_from_dist(distributional::dist_normal(0, 1), 0.5),
    crps_norm(0.5, 0, 1)
  )
})
