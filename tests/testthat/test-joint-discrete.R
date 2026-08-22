test_that("joint scores handle missing outcomes via a submatrix", {
  set.seed(26)
  dat <- norm_simulate(100, kind = "panel", seed = 26)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = starts_with("marker_"))
  sc <- predict(fit, newdata = dat, uncertainty = "conditional", allow_extrapolation = TRUE)
  sc$z[sc$.outcome == "marker_03" & sc$.row <= 10] <- NA_real_
  jt <- norm_joint(sc, covariance = "shrinkage")
  expect_s3_class(jt, "norm_joint")
  expect_true(all(jt$n_observed[1:10] == 2))
  expect_true(all(is.finite(jt$d2[1:10])))
})

test_that("discrete interval PIT is between neighbouring CDFs", {
  pit <- discrete_interval_pit(y = c(0, 1, 2), location = 1, scale = 1)
  expect_true(all(pit$lower <= pit$mid & pit$mid <= pit$upper))
  expect_true(all(pit$mid >= 0 & pit$mid <= 1))
})

test_that("gamlss engine reports unsupported when the package is missing", {
  spec <- norm_spec(family = norm_gaussian(), location = ~ age, engine = "gamlss")
  dat <- norm_simulate(40, kind = "gaussian_location", seed = 1)
  fit <- norm_fit(spec, data = dat, outcomes = "y")
  if (!requireNamespace("gamlss", quietly = TRUE)) {
    expect_equal(fit$models$y$status, "unsupported_family")
  } else {
    expect_true(fit$models$y$status %in% c("ok", "nonconverged"))
  }
})
