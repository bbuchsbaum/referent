test_that("missing outcomes do not change the other fit", {
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
  expect_equal(tidy(fit_miss)$n, c(120L, 80L))
})

test_that("a missing predictor is never imputed: constant-scale and gaulss Gaussian", {
  dat <- norm_simulate(150, kind = "gaussian", scale = "age", seed = 17)
  specs <- list(
    constant = norm_spec(norm_gaussian(), location = ~ s(age, k = 5) + sex),
    gaulss = norm_spec(norm_gaussian(), location = ~ s(age, k = 5) + sex, scale = ~ s(age, k = 4))
  )
  new <- dat[1:3, ]
  new$age[1] <- NA
  for (nm in names(specs)) {
    fit <- norm_fit(specs[[nm]], data = dat, outcomes = "y")
    if (nm == "gaulss") {
      expect_true(inherits(fit$models$y$model$family, "general.family"), info = nm)
    }
    for (unc in c("conditional", "total")) {
      sc <- predict(fit, newdata = new, uncertainty = unc)
      expect_equal(sc$status, c("missing_predictor", "ok", "ok"), info = paste(nm, unc))
      expect_equal(sc$support, c("unknown", "in", "in"), info = paste(nm, unc))
      expect_true(is.na(sc$z[[1]]) && is.na(sc$centile[[1]]) && is.na(sc$median[[1]]),
                  info = paste(nm, unc))
      expect_true(all(is.finite(sc$z[2:3])), info = paste(nm, unc))
      # the other rows are scored exactly as they would be without the NA row
      ref <- predict(fit, newdata = new[2:3, ], uncertainty = unc)
      expect_equal(sc$z[2:3], ref$z, tolerance = 1e-10, info = paste(nm, unc))
    }
    d <- predict(fit, newdata = new, type = "distribution", uncertainty = "conditional")$y
    expect_true(is.na(mean(d)[[1]]))
    expect_true(all(is.finite(mean(d)[2:3])))
  }
})

test_that("SHASH assessment, crossfit, and printing tolerate NA outcomes and predictors", {
  skip_on_cran()
  dat <- norm_simulate(300, kind = "shash", seed = 11)
  sfit <- norm_fit(simple_spec("shash"), data = dat, outcomes = "y")
  test <- norm_simulate(40, kind = "shash", seed = 12)
  test$y[5] <- NA
  test$age[3] <- NA
  a <- norm_assess(sfit, newdata = test)
  expect_true(is.finite(a$overall$crps))
  expect_equal(a$marginal$n, 38L)
  dt <- predict(sfit, newdata = test[1:6, ], type = "distribution", uncertainty = "total")
  expect_no_error(out <- format(dt$y))
  expect_equal(length(out), 6L)
  cf_dat <- dat[1:200, ]
  cf_dat$y[4] <- NA
  cf <- norm_crossfit(simple_spec("shash"), data = cf_dat, outcomes = "y", folds = 2)
  expect_true(is.na(cf$crps[cf$.row == 4]))
  expect_true(mean(is.finite(cf$crps)) > 0.9)
})
