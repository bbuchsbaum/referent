test_that("tails are computed in log space: z = 8, 10, 40 are distinct and monotone", {
  d <- distributional::dist_normal(0, 1)
  sc <- as_scores(d, c(-40, -10, -8, 0, 8, 10, 40))
  expect_true(all(is.finite(sc$z)))
  expect_equal(sc$z, c(-40, -10, -8, 0, 8, 10, 40), tolerance = 1e-8)
  expect_true(all(diff(sc$z) > 0))
  expect_true(all(diff(sc$tail_surprisal[5:7]) > 0))
  expect_equal(sc$tail_surprisal[[7]], -(log(2) + stats::pnorm(-40, log.p = TRUE)))
  s <- dist_shash(0, 1, 0.6, 0.85)
  scs <- as_scores(s, c(-60, -20, 0, 20, 60, 200))
  expect_true(all(is.finite(scs$z)))
  expect_true(all(diff(scs$z) > 0))
  expect_gt(scs$z[[6]], scs$z[[5]] + 0.5)
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
  expect_equal(mc_draws(vctrs::vec_data(d)[[1]]), 200L)
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
  expect_equal(mc_draws(vctrs::vec_data(d2)[[1]]), 50L)
  expect_false(isTRUE(all.equal(dist_unpack(d2)$mu[, 1], dist_unpack(d)$mu[, 1])))
  d3 <- predict(fit, newdata = new, type = "distribution", uncertainty = "total", n_draw = 20)$y
  expect_equal(mc_draws(vctrs::vec_data(d3)[[1]]), 20L)
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
