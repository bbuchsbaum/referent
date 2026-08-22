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
  sc2 <- predict(ad, newdata = new, uncertainty = "total", allow_extrapolation = TRUE)
  expect_true(is.finite(sc2$z[[1]]))
  # the adaptation offset shifts the total-uncertainty distribution too
  d_base <- predict(fit, newdata = new[-1, ], type = "distribution", uncertainty = "total")$y
  d_ad <- predict(ad, newdata = new[-1, ], type = "distribution", uncertainty = "total")$y
  off <- ad$adaptation$offsets$y$.all$location
  expect_equal(mean(d_ad), mean(d_base) + off)
  expect_gt(ad$adaptation$offsets$y$.all$location_se, 0)
  expect_true(all(variance(d_ad) > variance(d_base)))
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

test_that("total-uncertainty adaptation propagates the offset SE exactly for a draw mixture", {
  dat <- norm_simulate(260, kind = "gaussian", scale = "age", seed = 71)
  spec <- norm_spec(family = norm_gaussian(), location = ~ s(age, k = 5) + sex,
                    scale = ~ s(age, k = 4))
  fit <- norm_fit(spec, data = dat[1:200, ], outcomes = "y")
  ad <- norm_adapt(fit, data = dat[201:240, ], parameters = "location")
  off <- ad$adaptation$offsets$y$.all
  new <- dat[241:260, ]
  d_base <- predict(fit, newdata = new, type = "distribution", uncertainty = "total")$y
  d_ad <- predict(ad, newdata = new, type = "distribution", uncertainty = "total")$y
  expect_s3_class(vctrs::vec_data(d_ad)[[1]], "dist_shash_draws")
  expect_equal(mean(d_ad), mean(d_base) + off$location, tolerance = 1e-10)
  # additive up to the chance covariance of the jitter with each row's draws
  expect_equal(variance(d_ad), variance(d_base) + off$location_se^2, tolerance = 0.02)
  expect_gt(off$location_se, 0)
  # conditional predictions carry the offset but not its uncertainty
  c_base <- predict(fit, newdata = new, type = "distribution", uncertainty = "conditional")$y
  c_ad <- predict(ad, newdata = new, type = "distribution", uncertainty = "conditional")$y
  expect_equal(variance(c_ad), variance(c_base))
})

test_that("joint location+scale adaptation is not biased by feedback between the offsets", {
  # A site shifted by 6 units (~4.6 SD) with log-scale 0.2. The shrunk
  # targets are n / (n + prior_n) * 6 for the location and about
  # 2n / (2n + prior_n) * 0.2 for the log scale; with zero priors the
  # truth itself. Under the old joint fit the location came out at a
  # third of the truth at n = 20 and the log scale at 0.7.
  ref <- norm_simulate(1000, seed = 100)
  fit <- norm_fit(norm_spec(norm_gaussian(), ~ s(age, k = 5) + sex), data = ref, outcomes = "y")
  one <- function(seed, n_local, lp = 10, sp = 25) {
    site <- norm_simulate(n_local + 1000, site_shift = 6, site_log_scale = 0.2, seed = seed)
    loc <- site[seq_len(n_local), ]
    held <- site[-seq_len(n_local), ]
    ad <- norm_adapt(fit, loc, parameters = c("location", "scale"),
                     location_prior_n = lp, scale_prior_n = sp)
    off <- ad$adaptation$offsets$y$.all
    both <- norm_adapt(fit, loc, parameters = "location", location_prior_n = lp)
    sc <- predict(ad, held, uncertainty = "conditional", allow_extrapolation = TRUE)
    c(location = off$location, location_se = off$location_se, scale = off$scale,
      scale_se = off$scale_se, loc_only = both$adaptation$offsets$y$.all$location,
      var_z = stats::var(sc$z))
  }
  small <- t(sapply(1:6, one, n_local = 20))
  # the scale estimate must not drag the location away from the
  # location-only estimate, and the log scale must not explode
  expect_equal(small[, "location"], small[, "loc_only"], tolerance = 1e-3)
  expect_lt(abs(mean(small[, "location"]) - 20 / 30 * 6), 0.4)
  expect_lt(mean(small[, "scale"]), 0.45)

  mid <- t(sapply(1:6, one, n_local = 200))
  expect_lt(abs(mean(mid[, "location"]) - 200 / 210 * 6), 2 * mean(mid[, "location_se"]))
  expect_lt(abs(mean(mid[, "scale"]) - 400 / 425 * 0.2), 2 * mean(mid[, "scale_se"]))
  expect_true(all(abs(mid[, "location"] - 200 / 210 * 6) < 3 * mid[, "location_se"]))
  expect_lt(abs(mean(mid[, "var_z"]) - 1), 0.1)
  expect_true(all(abs(mid[, "var_z"] - 1) < 0.3))

  free <- t(sapply(1:6, one, n_local = 200, lp = 0, sp = 0))
  expect_lt(abs(mean(free[, "location"]) - 6), 0.15)
  expect_lt(abs(mean(free[, "scale"]) - 0.2), 0.05)
})
