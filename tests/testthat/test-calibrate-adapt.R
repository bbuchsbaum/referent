test_that("PIT recalibration is a valid map that does not collapse tails", {
  train <- norm_simulate(200, seed = 18)
  cal <- norm_simulate(120, seed = 19)
  test <- norm_simulate(80, seed = 20)
  fit <- norm_fit(simple_spec(), data = train, outcomes = "y")
  cal_fit <- norm_calibrate(fit, data = cal)
  expect_s3_class(cal_fit, "norm_fit")
  expect_s3_class(cal_fit$calibration, "norm_calibration")
  sc <- predict(cal_fit, newdata = test, uncertainty = "conditional")
  expect_true(all(sc$calibrated))
  expect_true(all(sc$centile > 0 & sc$centile < 1, na.rm = TRUE))
  # well-behaved data: calibrated z stays in a sane range (no +/-7.03)
  expect_lt(max(abs(sc$z), na.rm = TRUE), 4.5)
  raw <- predict(fit, newdata = test, uncertainty = "conditional")
  expect_gt(stats::cor(raw$z, sc$z), 0.99)
  # an extreme observation keeps an extreme, finite, monotone calibrated z
  far <- test[1:3, ]
  far$y <- far$y + c(10, 20, 60) * 1.3
  scf <- predict(cal_fit, newdata = far, uncertainty = "conditional",
                 allow_extrapolation = TRUE)
  expect_true(all(is.finite(scf$z)))
  expect_true(all(diff(scf$z) > 0))
  expect_gt(scf$z[[3]], 30)
})

test_that("calibration by group uses the group map and pools unknown groups", {
  dat <- norm_simulate(400, site_shift = c(0, 2, 0, 0), seed = 21)
  train <- dat[dat$site %in% c("A", "C", "D"), ]
  cal <- dat[dat$site == "B", ][1:60, ]
  test <- dat[dat$site == "B", ][61:100, ]
  fit <- norm_fit(simple_spec(), data = train, outcomes = "y")
  cal_fit <- norm_calibrate(fit, data = rbind(cal, train[1:60, ]), by = site)
  expect_equal(cal_fit$calibration$by, "site")
  raw <- predict(fit, newdata = test, uncertainty = "conditional")
  sc <- predict(cal_fit, newdata = test, uncertainty = "conditional")
  # site B is shifted by +2 units (~1.5 SD); the B map recentres it
  expect_gt(mean(raw$z), 1)
  expect_lt(abs(mean(sc$z)), 0.35)
})

test_that("adapting to a site from the reference generator changes nothing", {
  dat <- norm_simulate(4000, site_log_scale = 0, seed = 23)
  train <- dat[dat$site != "B", ]
  local <- dat[dat$site == "B", ][1:150, ]
  held <- dat[dat$site == "B", ][151:900, ]
  fit <- norm_fit(simple_spec(), data = train, outcomes = "y")
  ad <- norm_adapt(fit, data = local, by = site, parameters = c("location", "scale"))
  expect_s3_class(ad, "norm_fit")
  expect_s3_class(ad$adaptation, "norm_adaptation")
  off <- ad$adaptation$offsets$y$B
  expect_lt(abs(off$location), 0.25)
  expect_lt(abs(off$scale), 0.1)
  sc <- predict(ad, newdata = held, uncertainty = "conditional", allow_extrapolation = TRUE)
  expect_s3_class(sc, "norm_scores")
  expect_lt(abs(mean(sc$z) - 0), 0.1)
  expect_lt(abs(stats::var(sc$z) - 1), 0.1)
})

test_that("adaptation recovers a scale ratio of 1.5 and a location shift", {
  dat <- norm_simulate(4000, site_shift = c(0, 1.2, 0, 0),
                       site_log_scale = c(0, log(1.5), 0, 0), seed = 24)
  train <- dat[dat$site != "B", ]
  loc_b <- dat[dat$site == "B", ]
  local <- loc_b[1:500, ]
  held <- loc_b[501:nrow(loc_b), ]
  fit <- norm_fit(simple_spec(), data = train, outcomes = "y")
  raw <- predict(fit, newdata = held, uncertainty = "conditional", allow_extrapolation = TRUE)
  expect_gt(stats::var(raw$z), 1.8)
  ad <- norm_adapt(fit, data = local, by = site, parameters = c("location", "scale"))
  off <- ad$adaptation$offsets$y$B
  expect_lt(abs(exp(off$scale) / 1.5 - 1), 0.10)
  expect_lt(abs(off$location - 1.2), 0.25)
  expect_lt(off$location_se, 0.1)
  sc <- predict(ad, newdata = held, uncertainty = "conditional", allow_extrapolation = TRUE)
  expect_lt(abs(mean(sc$z)), 0.15)
  expect_lt(abs(stats::var(sc$z) - 1), 0.15)
})

test_that("adaptation works for SHASH fits via normal scores", {
  dat <- norm_simulate(3200, kind = "shash", site_log_scale = c(0, log(1.5), 0, 0),
                       seed = 25)
  train <- dat[dat$site != "B", ]
  loc_b <- dat[dat$site == "B", ]
  fit <- norm_fit(simple_spec("shash"), data = train, outcomes = "y")
  ad <- norm_adapt(fit, data = loc_b[1:300, ], by = site, parameters = c("location", "scale"))
  expect_lt(abs(exp(ad$adaptation$offsets$y$B$scale) / 1.5 - 1), 0.12)
  sc <- predict(ad, newdata = loc_b[301:nrow(loc_b), ], uncertainty = "conditional",
                allow_extrapolation = TRUE)
  expect_lt(abs(stats::var(sc$z) - 1), 0.15)
})

test_that("adapted predict honours extrapolation, ids, and total uncertainty", {
  dat <- norm_simulate(300, seed = 26)
  dat$subject <- sprintf("s%03d", seq_len(nrow(dat)))
  fit <- norm_fit(simple_spec(), data = dat[1:220, ], outcomes = "y", id = subject)
  ad <- norm_adapt(fit, data = dat[221:260, ], parameters = "location")
  new <- dat[261:270, ]
  new$age[1] <- 140
  sc <- predict(ad, newdata = new, uncertainty = "total")
  expect_s3_class(sc, "norm_scores")
  expect_equal(sc$.id, new$subject)
  expect_true(is.na(sc$z[[1]]))
  expect_equal(sc$support[[1]], "out")
  expect_true(all(is.finite(sc$z[-1])))
  expect_true(all(sc$epistemic_sd > 0))
  sc2 <- predict(ad, newdata = new, uncertainty = "total", allow_extrapolation = TRUE)
  expect_true(is.finite(sc2$z[[1]]))
  # the adaptation offset shifts the total-uncertainty distribution too
  d_base <- predict(fit, newdata = new[-1, ], type = "distribution", uncertainty = "total")$y
  d_ad <- predict(ad, newdata = new[-1, ], type = "distribution", uncertainty = "total")$y
  off <- ad$adaptation$offsets$y$.all$location
  expect_equal(field_or(d_ad, "location"), field_or(d_base, "location") + off)
  expect_true(all(field_or(d_ad, "epistemic_sd") >= field_or(d_base, "epistemic_sd")))
})

test_that("reference bundles keep calibration and adaptation", {
  dat <- norm_simulate(300, seed = 27)
  fit <- norm_fit(simple_spec(), data = dat[1:150, ], outcomes = "y")
  fit <- norm_adapt(fit, data = dat[151:200, ], parameters = "location")
  fit <- norm_calibrate(fit, data = dat[201:250, ])
  ref <- norm_reference(fit)
  expect_s3_class(ref$calibration, "norm_calibration")
  expect_s3_class(ref$adaptation, "norm_adaptation")
  a <- predict(fit, newdata = dat[251:300, ], uncertainty = "conditional")
  b <- predict(ref, newdata = dat[251:300, ], uncertainty = "conditional")
  expect_true(all(b$calibrated))
  expect_equal(a$z, b$z)
})
