test_that("PIT recalibration remains a valid CDF map", {
  set.seed(18)
  train <- norm_simulate(200, kind = "gaussian_location", seed = 18)
  cal <- norm_simulate(120, kind = "gaussian_location", seed = 19)
  test <- norm_simulate(80, kind = "gaussian_location", seed = 20)
  fit <- norm_fit(simple_spec(), data = train, outcomes = "y")
  cal_fit <- norm_calibrate(fit, data = cal)
  expect_s3_class(cal_fit, "norm_calibrated")
  sc <- predict(cal_fit, newdata = test, uncertainty = "conditional")
  expect_true(all(sc$calibrated))
  expect_true(all(sc$centile >= 0 & sc$centile <= 1, na.rm = TRUE))
})

test_that("adaptation recovers a location shift", {
  set.seed(23)
  dat <- norm_simulate(280, kind = "shift", seed = 23)
  train <- dat[dat$site != "B", ]
  local <- dat[dat$site == "B", ]
  local <- local[seq_len(min(40, nrow(local))), ]
  target <- dat[dat$site == "B", ]
  target <- target[setdiff(seq_len(nrow(target)), seq_len(nrow(local))), ]
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ s(age, k = 5) + sex, scale = ~1),
    data = train,
    outcomes = "y"
  )
  ad <- norm_adapt(fit, data = local, by = site, parameters = c("location", "scale"))
  expect_s3_class(ad, "norm_adaptation")
  off <- ad$offsets$y$B$location
  expect_gt(off, 0.3)
  sc <- predict(ad, newdata = target, uncertainty = "conditional")
  expect_lt(abs(mean(sc$z, na.rm = TRUE)), 0.6)
})
