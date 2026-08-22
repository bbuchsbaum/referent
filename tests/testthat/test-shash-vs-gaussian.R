test_that("SHASH beats Gaussian on strongly skewed held-out data", {
  skip_on_cran()
  train <- simulate_skewed(280, seed = 61, skew = 1.4, tail = 0.7)
  test <- simulate_skewed(220, seed = 62, skew = 1.4, tail = 0.7)

  fit_g <- norm_fit(simple_spec(scale = TRUE), data = train, outcomes = "y")
  fit_s <- norm_fit(simple_spec("shash", scale = TRUE), data = train, outcomes = "y")
  skip_if(
    identical(fit_s$models$y$status, "nonconverged"),
    "SHASH fit did not converge"
  )

  a_g <- norm_assess(fit_g, newdata = test)
  a_s <- norm_assess(fit_s, newdata = test)

  expect_gt(a_s$overall$mean_log_score, a_g$overall$mean_log_score)
  expect_lt(pit_ks(a_s$scores$centile), pit_ks(a_g$scores$centile))
  expect_lt(worm_rmse(a_s$scores$z), worm_rmse(a_g$scores$z))
  expect_lt(abs(a_s$marginal$skew_z), abs(a_g$marginal$skew_z))
})

test_that("Gaussian remains competitive when the truth is Gaussian", {
  skip_on_cran()
  train <- norm_simulate(240, kind = "gaussian", scale = "age", seed = 63)
  test <- norm_simulate(180, kind = "gaussian", scale = "age", seed = 64)
  fit_g <- norm_fit(simple_spec(scale = TRUE), data = train, outcomes = "y")
  fit_s <- norm_fit(simple_spec("shash", scale = TRUE), data = train, outcomes = "y")
  skip_if(
    identical(fit_s$models$y$status, "nonconverged"),
    "SHASH fit did not converge"
  )
  a_g <- norm_assess(fit_g, newdata = test)
  a_s <- norm_assess(fit_s, newdata = test)
  # Extra SHASH shape parameters should not win by a large margin.
  expect_gt(a_g$overall$mean_log_score, a_s$overall$mean_log_score - 0.04)
  expect_lt(abs(a_g$marginal$var_z - 1), 0.35)
})
