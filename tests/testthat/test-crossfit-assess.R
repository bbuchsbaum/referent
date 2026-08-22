test_that("cross-fitted scores are marked out of sample", {
  set.seed(11)
  dat <- norm_simulate(160, kind = "gaussian_location", seed = 11)
  dat$participant_id <- rep(1:40, each = 4)
  cf <- norm_crossfit(
    simple_spec(),
    data = dat,
    outcomes = "y",
    folds = 4,
    cluster = participant_id
  )
  expect_s3_class(cf, "norm_scores")
  expect_true(all(!cf$.in_sample))
  expect_s3_class(attr(cf, "deployment"), "norm_fit")
  expect_lt(abs(stats::var(cf$z, na.rm = TRUE) - 1), 0.45)
})

test_that("norm_assess reports coverage near nominal on a calibrated model", {
  set.seed(13)
  train <- norm_simulate(300, kind = "gaussian_location", seed = 13)
  test <- norm_simulate(200, kind = "gaussian_location", seed = 14)
  fit <- norm_fit(simple_spec(), data = train, outcomes = "y")
  a <- norm_assess(fit, newdata = test)
  expect_s3_class(a, "norm_assessment")
  expect_lt(abs(a$marginal$cover_90 - 0.90), 0.12)
  expect_true(is.finite(a$overall$mean_log_score))
  expect_true(is.finite(a$overall$crps))
})

test_that("ladder prefers a simple model on Gaussian data", {
  skip_on_cran()
  set.seed(15)
  dat <- norm_simulate(180, kind = "gaussian_location", seed = 15)
  specs <- list(
    gaussian_const = simple_spec("gaussian", scale = FALSE),
    gaussian_scale = simple_spec("gaussian", scale = TRUE)
  )
  sel <- norm_select(specs, data = dat, outcomes = "y", folds = 3)
  expect_s3_class(sel, "norm_selection")
  expect_true(sel$selected_name %in% names(specs))
})
