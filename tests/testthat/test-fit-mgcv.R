test_that("Gaussian location-scale fit produces calibrated Z on new data", {
  set.seed(21)
  train <- ref_simulate(350, seed = 21)
  test <- ref_simulate(200, seed = 22)
  fit <- ref_fit(simple_spec(), data = train, outcomes = "y")
  expect_equal(fit$models$y$status, "ok")
  sc <- predict(fit, newdata = test, type = "scores", uncertainty = "conditional")
  expect_s3_class(sc, "ref_scores")
  expect_lt(abs(mean(sc$z, na.rm = TRUE)), 0.25)
  expect_lt(abs(stats::var(sc$z, na.rm = TRUE) - 1), 0.35)
  expect_false(any(sc$.in_sample))
})

test_that("failed outcomes do not abort the panel", {
  set.seed(3)
  dat <- ref_simulate(80, seed = 3)
  dat$marker_bad <- 1
  fit <- ref_fit(simple_spec(), data = dat, outcomes = c("marker_01", "marker_bad"))
  expect_equal(fit$models$marker_bad$status, "insufficient_variation")
  expect_equal(fit$models$marker_01$status, "ok")
  out <- paste(cli::cli_fmt(print(fit)), collapse = "\n")
  expect_match(out, "insufficient_variation")
  expect_match(out, "marker_bad")
})

test_that("SHASH engine returns a dist_shash vector", {
  skip_on_cran()
  set.seed(5)
  dat <- ref_simulate(250, kind = "shash", seed = 5)
  fit <- ref_fit(simple_spec("shash"), data = dat, outcomes = "y")
  expect_equal(fit$models$y$status, "ok")
  d <- predict(fit, newdata = dat[1:5, ], type = "distribution",
               uncertainty = "conditional")$y
  expect_s3_class(d, "distribution")
  expect_s3_class(vctrs::vec_data(d)[[1]], "dist_shash")
  expect_true(all(dist_unpack(d)$sigma > 0))
})

test_that("total uncertainty widens the Gaussian predictive analytically", {
  set.seed(31)
  dat <- ref_simulate(120, seed = 31)
  fit <- ref_fit(simple_spec(), data = dat, outcomes = "y")
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
  dat <- ref_simulate(120, seed = 8)
  fit <- ref_fit(simple_spec(), data = dat, outcomes = "y")
  tmp <- tempfile(fileext = ".rds")
  suppressWarnings(saveRDS(fit, tmp))
  fit2 <- readRDS(tmp)
  a <- predict(fit, newdata = dat[1:10, ], uncertainty = "conditional")
  b <- predict(fit2, newdata = dat[1:10, ], uncertainty = "conditional")
  expect_equal(a$z, b$z, tolerance = 1e-10)
})

test_that("predict carries the declared subject id", {
  dat <- ref_simulate(40, seed = 74)
  dat$participant_id <- paste0("S", seq_len(nrow(dat)))
  fit <- ref_fit(
    ref_spec(family = ref_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y",
    id = participant_id
  )
  expect_equal(fit$id_name, "participant_id")
  expect_false("participant_id" %in% fit$covariates)
  sc <- predict(fit, newdata = dat[1:5, ], uncertainty = "conditional")
  expect_equal(as.character(sc$.id), dat$participant_id[1:5])
})

test_that("non-numeric outcomes are reported as unsupported_type without aborting the panel", {
  dat <- ref_simulate(100, seed = 3)
  dat$f <- factor(sample(c("a", "b"), 100, TRUE))
  dat$lg <- dat$y > 10
  dat$ch <- as.character(dat$f)
  spec <- ref_spec(ref_gaussian(), ~ s(age, k = 5))
  expect_silent(fit <- suppressMessages(ref_fit(spec, dat, c("y", "f", "lg", "ch"))))
  expect_equal(unname(fit_statuses(fit)),
               c("ok", "unsupported_type", "unsupported_type", "unsupported_type"))
  expect_match(fit$models$f$message, "factor")
  expect_null(fit$reference_baseline$f)
  sc <- predict(fit, dat[1:3, ], uncertainty = "conditional")
  expect_equal(unique(sc$status[sc$.outcome == "lg"]), "unsupported_type")
  expect_true(all(is.finite(sc$z[sc$.outcome == "y"])))
})

test_that("ref_fit aborts when the id column is missing", {
  dat <- ref_simulate(60, seed = 4)
  spec <- ref_spec(ref_gaussian(), ~ age)
  expect_error(ref_fit(spec, dat, "y", id = "nope"), "nope")
  expect_error(ref_fit(spec, dat, "y", id = nope), "nope")
})

test_that("predicting on 0-row newdata with ok and failed outcomes gives a typed 0-row table", {
  dat <- ref_simulate(100, seed = 3)
  dat$const <- 1
  fit <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 5)), dat, c("y", "const"))
  expect_equal(unname(fit_statuses(fit)), c("ok", "insufficient_variation"))
  for (unc in c("conditional", "total")) {
    sc <- predict(fit, dat[0, ], uncertainty = unc)
    expect_s3_class(sc, "ref_scores")
    expect_equal(nrow(sc), 0L)
    full <- predict(fit, dat[1:2, ], uncertainty = unc)
    expect_equal(names(sc), names(full))
    expect_equal(vapply(sc, class, ""), vapply(full, class, ""))
  }
})
