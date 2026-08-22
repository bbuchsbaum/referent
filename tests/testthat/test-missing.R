test_that("missing outcomes do not change the other fit", {
  set.seed(16)
  dat <- norm_simulate(120, seed = 16)
  fit_full <- norm_fit(simple_spec(), data = dat, outcomes = c("marker_01", "marker_02"))
  dat2 <- dat
  dat2$marker_02[seq_len(40)] <- NA
  expect_message(
    fit_miss <- norm_fit(simple_spec(), data = dat2, outcomes = c("marker_01", "marker_02")),
    "marker_02: dropped 40 rows"
  )
  a <- predict(fit_full, newdata = dat[1:8, ], uncertainty = "conditional")
  b <- predict(fit_miss, newdata = dat[1:8, ], uncertainty = "conditional")
  z1a <- a$z[a$.outcome == "marker_01"]
  z1b <- b$z[b$.outcome == "marker_01"]
  expect_equal(z1a, z1b, tolerance = 1e-6)
})

test_that("predictor missingness is not imputed", {
  set.seed(17)
  dat <- norm_simulate(80, seed = 17)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "y")
  extra <- dat[1:2, ]
  extra$age[1] <- NA
  sc <- predict(fit, newdata = extra, uncertainty = "conditional")
  expect_true(is.na(sc$z[[1]]) || is.na(sc$centile[[1]]) || is.na(sc$median[[1]]))
})
