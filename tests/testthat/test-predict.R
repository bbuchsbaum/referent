test_that("empty newdata returns an empty distribution", {
  set.seed(40)
  dat <- norm_simulate(60, seed = 40)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  d <- predict(fit, newdata = dat[0, ], type = "distribution")
  expect_equal(nrow(d), 0L)
  expect_equal(length(d$y), 0L)
  sc <- predict(fit, newdata = dat[0, ], uncertainty = "conditional")
  expect_equal(nrow(sc), 0L)
})

test_that("total uncertainty works for one row and gaulss", {
  set.seed(41)
  dat <- norm_simulate(80, kind = "gaussian", scale = "age", seed = 41)
  fit <- norm_fit(
    norm_spec(
      family = norm_gaussian(),
      location = ~ s(age, k = 5) + sex,
      scale = ~ s(age, k = 4)
    ),
    data = dat,
    outcomes = "y"
  )
  d <- predict(fit, newdata = dat[1, ], type = "distribution", uncertainty = "total")$y
  expect_equal(length(d), 1L)
  u <- dist_unpack(d)
  expect_s3_class(u, "dist_shash_draws")
  expect_equal(dim(u$sigma), dim(u$mu))
  expect_gt(stats::sd(u$mu), 0)
})

test_that("unseen groups do not abort prediction", {
  set.seed(42)
  dat <- norm_simulate(50, seed = 42)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  extra <- dat[1, ]
  extra$sex <- factor("X", levels = c(levels(dat$sex), "X"))
  sc <- predict(fit, newdata = extra, uncertainty = "conditional")
  expect_equal(nrow(sc), 1L)
  expect_equal(sc$support, "new_group")
})

test_that("NA predictors give NA scores with status missing_predictor (gaulss and shash)", {
  dat <- norm_simulate(150, kind = "gaussian", scale = "age", seed = 51)
  spec <- norm_spec(family = norm_gaussian(), location = ~ s(age, k = 5) + sex,
                    scale = ~ s(age, k = 4))
  fit <- norm_fit(spec, data = dat, outcomes = "y")
  new <- dat[1:4, ]
  new$age[2] <- NA
  new$sex[3] <- NA
  for (unc in c("conditional", "total")) {
    sc <- predict(fit, newdata = new, uncertainty = unc)
    expect_equal(sc$status, c("ok", "missing_predictor", "missing_predictor", "ok"))
    expect_true(all(is.na(sc$z[2:3])))
    expect_true(all(is.na(sc$centile[2:3])))
    expect_true(all(is.finite(sc$z[c(1, 4)])))
    expect_equal(sc$support[2:3], c("unknown", "unknown"))
  }
  sfit <- norm_fit(simple_spec("shash"), data = dat, outcomes = "y")
  ssc <- predict(sfit, newdata = new, uncertainty = "conditional")
  expect_equal(ssc$status[2], "missing_predictor")
  expect_true(is.na(ssc$z[[2]]))
})

test_that("total uncertainty is deterministic and the returned distribution is the mixture", {
  dat <- norm_simulate(150, kind = "gaussian", scale = "age", seed = 52)
  spec <- norm_spec(family = norm_gaussian(), location = ~ s(age, k = 5) + sex,
                    scale = ~ s(age, k = 4))
  fit <- norm_fit(spec, data = dat, outcomes = "y")
  new <- dat[1:12, ]
  a <- predict(fit, newdata = new, uncertainty = "total")
  b <- predict(fit, newdata = new, uncertainty = "total")
  expect_identical(a, b)
  dt <- predict(fit, newdata = new, type = "distribution", uncertainty = "total")
  expect_s3_class(dt, "tbl_df")
  expect_named(dt, c(".id", "y"))
  d <- dt$y
  expect_s3_class(d, "distribution")
  expect_equal(ncol(dist_unpack(d)$mu), 200L)
  expect_equal(dist_cdf(d, new$y), a$centile)
  expect_equal(dist_quantile(d, 0.5), a$median)
  expect_equal(dist_log_density(d, new$y), a$log_density)
  # the mixture median is the median of the mixture cdf
  expect_equal(dist_cdf(d, a$median), rep(0.5, 12), tolerance = 1e-6)
  # mixture log density is the log-mean over draws, not the plug-in density
  cond <- predict(fit, newdata = new, uncertainty = "conditional")
  expect_false(isTRUE(all.equal(cond$log_density, a$log_density)))
  # a different seed gives different, a different n_draw gives fewer draws
  spec2 <- spec
  spec2$control <- list(seed = 99, n_draw = 50)
  fit2 <- fit
  fit2$spec <- spec2
  d2 <- predict(fit2, newdata = new, type = "distribution", uncertainty = "total")$y
  expect_equal(ncol(dist_unpack(d2)$mu), 50L)
  expect_false(isTRUE(all.equal(dist_unpack(d2)$mu[, 1], dist_unpack(d)$mu[, 1])))
  d3 <- predict(fit, newdata = new, type = "distribution", uncertainty = "total", n_draw = 20)$y
  expect_equal(ncol(dist_unpack(d3)$mu), 20L)
})

test_that("non-syntactic outcome names fit and predict", {
  dat <- norm_simulate(120, seed = 53)
  dat[["brain volume"]] <- dat$y
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "brain volume")
  expect_equal(fit$models[["brain volume"]]$status, "ok")
  sc <- predict(fit, newdata = dat[1:5, ], uncertainty = "conditional")
  ref <- predict(norm_fit(simple_spec(), data = dat, outcomes = "y"),
                 newdata = dat[1:5, ], uncertainty = "conditional")
  expect_equal(sc$.outcome, rep("brain volume", 5))
  expect_equal(sc$z, ref$z)
})

test_that("norm_assess requires newdata, honours by=, checks every covariate, and reports ev/smse", {
  train <- norm_simulate(300, seed = 54)
  test <- norm_simulate(200, seed = 55)
  fit <- norm_fit(simple_spec(), data = train, outcomes = "y")
  expect_error(norm_assess(fit), "newdata")
  a <- norm_assess(fit, newdata = test, by = site)
  expect_true(".group" %in% names(a$marginal))
  expect_setequal(a$marginal$.group, levels(test$site))
  expect_equal(sum(a$marginal$n), nrow(test))
  expect_equal(sort(unique(a$conditional$covariate)), "age")
  expect_true(all(c("location_se", "scale_se") %in% names(a$conditional)))
  expect_true(all(a$conditional$location_se > 0))
  ov <- a$overall
  expect_true(all(c("ev", "smse") %in% names(ov)))
  expect_equal(ov$ev, 1 - stats::var(a$scores$residual) / stats::var(a$scores$observed))
  expect_equal(ov$smse, mean(a$scores$residual^2) / stats::var(a$scores$observed))
  expect_gt(ov$ev, 0.3)
  expect_lt(ov$smse, 0.7)
  # two numeric covariates -> two conditional rows
  spec2 <- norm_spec(family = norm_gaussian(), location = ~ s(age, k = 5) + marker_01)
  fit2 <- norm_fit(spec2, data = train, outcomes = "y")
  a2 <- norm_assess(fit2, newdata = test)
  expect_setequal(a2$conditional$covariate, c("age", "marker_01"))
})

test_that("norm_support reports unknown for NA covariates", {
  dat <- norm_simulate(80, seed = 56)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "y")
  new <- dat[1:3, ]
  new$age[2] <- NA
  expect_equal(norm_support(fit, new)$support, c("in", "unknown", "in"))
})

test_that("augment matches the long score table", {
  d <- perf_data()
  fit <- norm_calibrate(norm_fit(perf_specs()$gaulss, d$ref, c("y", "marker_01")),
                        d$new, by = site)
  long <- predict(fit, d$new, uncertainty = "conditional")
  wide <- augment(fit, d$new, uncertainty = "conditional")
  for (nm in c("y", "marker_01")) {
    rows <- long[long$.outcome == nm, ]
    expect_equal(wide[[paste0(".z_", nm)]], rows$z)
    expect_equal(wide[[paste0(".centile_", nm)]], rows$centile)
  }
  expect_equal(wide$.support, long$support[long$.outcome == "y"])
})

test_that("draw parameters from the lp matrix match per-draw evaluation", {
  d <- perf_data()
  fit <- norm_fit(perf_specs()$shash, d$ref, "y")
  m <- fit$models$y
  lp <- predict(m$model, d$new[1:5, ], type = "lpmatrix")
  beta <- coef(m$model)
  draws <- cbind(beta, beta * 1.01, beta * 0.99)
  joint <- params_from_eta(eta_from_lp(lp, draws), m$model, "shash", m)
  for (j in 1:3) {
    single <- params_from_eta(eta_from_lp(lp, draws[, j]), m$model, "shash", m)
    for (p in c("location", "scale", "skew", "tail")) {
      expect_equal(joint[[p]][, j], single[[p]])
    }
  }
})

test_that("shifted log density agrees with the distribution methods", {
  d <- perf_data()
  fit <- norm_fit(perf_specs()$shash, d$ref, "y")
  dist <- predict(fit, d$new, type = "distribution", uncertainty = "conditional")$y
  y <- d$new$y
  u <- dist_unpack(dist)
  expect_equal(
    shifted_log_density(u, y, 0.3, 0.1),
    dist_log_density(shift_params(dist, 0.3, 0.1), y)
  )
  total <- predict(fit, d$new[1:10, ], type = "distribution", uncertainty = "total")$y
  expect_equal(
    shifted_log_density(dist_unpack(total), y[1:10], -0.2, 0.05),
    dist_log_density(shift_params(total, -0.2, 0.05), y[1:10])
  )
})
