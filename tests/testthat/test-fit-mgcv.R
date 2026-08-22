test_that("Gaussian location-scale fit produces calibrated Z on new data", {
  set.seed(21)
  train <- norm_simulate(350, seed = 21)
  test <- norm_simulate(200, seed = 22)
  fit <- norm_fit(simple_spec(), data = train, outcomes = "y")
  expect_equal(fit$models$y$status, "ok")
  sc <- predict(fit, newdata = test, type = "scores", uncertainty = "conditional")
  expect_s3_class(sc, "norm_scores")
  expect_lt(abs(mean(sc$z, na.rm = TRUE)), 0.25)
  expect_lt(abs(stats::var(sc$z, na.rm = TRUE) - 1), 0.35)
  expect_false(any(sc$.in_sample))
})

test_that("failed outcomes do not abort the panel", {
  set.seed(3)
  dat <- norm_simulate(80, seed = 3)
  dat$marker_bad <- 1
  fit <- norm_fit(simple_spec(), data = dat, outcomes = c("marker_01", "marker_bad"))
  expect_equal(fit$models$marker_bad$status, "insufficient_variation")
  expect_equal(fit$models$marker_01$status, "ok")
  out <- paste(cli::cli_fmt(print(fit)), collapse = "\n")
  expect_match(out, "insufficient_variation")
  expect_match(out, "marker_bad")
})

test_that("SHASH engine returns a dist_shash vector", {
  skip_on_cran()
  set.seed(5)
  dat <- norm_simulate(250, kind = "shash", seed = 5)
  fit <- norm_fit(simple_spec("shash"), data = dat, outcomes = "y")
  expect_equal(fit$models$y$status, "ok")
  d <- predict(fit, newdata = dat[1:5, ], type = "distribution",
               uncertainty = "conditional")$y
  expect_s3_class(d, "distribution")
  expect_s3_class(vctrs::vec_data(d)[[1]], "dist_shash")
  expect_true(all(dist_unpack(d)$sigma > 0))
})

test_that("total uncertainty widens the Gaussian predictive analytically", {
  set.seed(31)
  dat <- norm_simulate(120, seed = 31)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "y")
  cond <- predict(fit, newdata = dat[1:8, ], type = "distribution",
                  uncertainty = "conditional")$y
  tot <- predict(fit, newdata = dat[1:8, ], type = "distribution",
                 uncertainty = "total")$y
  # identity-location, constant-scale Gaussian: analytic total N(mu, s^2 + se^2)
  expect_s3_class(vctrs::vec_data(tot)[[1]], "dist_normal")
  se <- stats::predict(fit$models$y$model, newdata = dat[1:8, ], se.fit = TRUE)$se.fit
  expect_true(all(se > 0))
  expect_equal(mean(tot), mean(cond))
  expect_equal(variance(tot), variance(cond) + as.numeric(se)^2)
  y <- dat$y[1:8]
  pt <- dist_cdf(tot, y)
  pc <- dist_cdf(cond, y)
  expect_true(all(pt >= 0 & pt <= 1))
  expect_true(all(abs(pt - 0.5) <= abs(pc - 0.5) + 1e-12))
})

test_that("save/read round trip reproduces scores", {
  set.seed(8)
  dat <- norm_simulate(120, seed = 8)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "y")
  tmp <- tempfile(fileext = ".rds")
  suppressWarnings(saveRDS(fit, tmp))
  fit2 <- readRDS(tmp)
  a <- predict(fit, newdata = dat[1:10, ], uncertainty = "conditional")
  b <- predict(fit2, newdata = dat[1:10, ], uncertainty = "conditional")
  expect_equal(a$z, b$z, tolerance = 1e-10)
})

test_that("predict carries the declared subject id", {
  dat <- norm_simulate(40, seed = 74)
  dat$participant_id <- paste0("S", seq_len(nrow(dat)))
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y",
    id = participant_id
  )
  expect_equal(fit$id_name, "participant_id")
  expect_false("participant_id" %in% fit$covariates)
  sc <- predict(fit, newdata = dat[1:5, ], uncertainty = "conditional")
  expect_equal(as.character(sc$.id), dat$participant_id[1:5])
})
