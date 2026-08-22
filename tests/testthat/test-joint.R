test_that("joint scores handle missing outcomes via a submatrix", {
  set.seed(26)
  dat <- norm_simulate(100, seed = 26)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = starts_with("marker_"))
  sc <- predict(fit, newdata = dat, uncertainty = "conditional", allow_extrapolation = TRUE)
  sc$z[sc$.outcome == "marker_03" & sc$.row <= 10] <- NA_real_
  jt <- norm_joint(sc, covariance = "shrinkage")
  expect_s3_class(jt, "norm_joint")
  expect_true(all(jt$n_observed[1:10] == 2))
  expect_true(all(is.finite(jt$d2[1:10])))
})
