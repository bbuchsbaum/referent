test_that("Gaussian location-scale fit produces calibrated Z on new data", {
  set.seed(21)
  train <- norm_simulate(350, seed = 21)
  test <- norm_simulate(200, seed = 22)
  fit <- norm_fit(simple_spec(), data = train, outcomes = "y")
  expect_equal(fit$models$y$status, "ok")
  sc <- predict(fit, newdata = test, type = "scores", uncertainty = "conditional")
  expect_s3_class(sc, "norm_scores")
  expect_lt(abs(mean(sc$z, na.rm = TRUE)), 0.25)
  expect_lt(abs(stats::var(sc$z, na.rm = TRUE) - 1), 0.35)
  expect_false(any(sc$.in_sample))
})

test_that("failed outcomes do not abort the panel", {
  set.seed(3)
  dat <- norm_simulate(80, seed = 3)
  dat$marker_bad <- 1
  fit <- norm_fit(simple_spec(), data = dat, outcomes = c("marker_01", "marker_bad"))
  expect_equal(fit$models$marker_bad$status, "insufficient_variation")
  expect_equal(fit$models$marker_01$status, "ok")
  out <- paste(cli::cli_fmt(print(fit)), collapse = "\n")
  expect_match(out, "insufficient_variation")
  expect_match(out, "marker_bad")
})

test_that("SHASH engine returns a norm_dist", {
  skip_on_cran()
  set.seed(5)
  dat <- norm_simulate(250, kind = "shash", seed = 5)
  fit <- norm_fit(simple_spec("shash"), data = dat, outcomes = "y")
  expect_equal(fit$models$y$status, "ok")
  d <- predict(fit, newdata = dat[1:5, ], type = "distribution")$y
  expect_s3_class(d, "norm_dist")
  expect_equal(attr(d, "family"), "shash")
  expect_true(all(field_or(d, "scale") > 0))
})

test_that("total uncertainty inflates epistemic sd", {
  set.seed(31)
  dat <- norm_simulate(120, seed = 31)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "y")
  cond <- predict(fit, newdata = dat[1:8, ], type = "distribution",
                  uncertainty = "conditional")$y
  tot <- predict(fit, newdata = dat[1:8, ], type = "distribution",
                 uncertainty = "total")$y
  expect_gt(mean(field_or(tot, "epistemic_sd")), 0)
  # identity-location, constant-scale Gaussian: analytic total N(mu, s^2 + se^2)
  expect_true(is.null(attr(tot, "location_draws")))
  expect_equal(attr(tot, "uncertainty"), "total")
  expect_equal(
    field_or(tot, "scale"),
    sqrt(field_or(cond, "scale")^2 + field_or(tot, "epistemic_sd")^2)
  )
  y <- dat$y[1:8]
  expect_true(all(cdf(tot, y) >= 0 & cdf(tot, y) <= 1))
  expect_true(all(abs(cdf(tot, y) - 0.5) <= abs(cdf(cond, y) - 0.5) + 1e-12))
})

test_that("save/read round trip reproduces scores", {
  set.seed(8)
  dat <- norm_simulate(120, seed = 8)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "y")
  tmp <- tempfile(fileext = ".rds")
  suppressWarnings(saveRDS(fit, tmp))
  fit2 <- readRDS(tmp)
  a <- predict(fit, newdata = dat[1:10, ], uncertainty = "conditional")
  b <- predict(fit2, newdata = dat[1:10, ], uncertainty = "conditional")
  expect_equal(a$z, b$z, tolerance = 1e-10)
})
