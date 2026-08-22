test_that("covariates come from the spec formulas, not the data frame", {
  dat <- norm_simulate(120, seed = 80)
  spec <- norm_spec(family = norm_gaussian(), location = ~ s(age, k = 5) + sex, scale = ~1)
  fit <- norm_fit(spec, data = dat, outcomes = "y")
  expect_equal(fit$covariates, c("age", "sex"))
  expect_false("site" %in% fit$covariates)
  expect_false("marker_01" %in% fit$covariates)
  spec4 <- norm_spec(
    family = norm_shash(),
    location = ~ s(age, k = 5) + sex,
    scale = ~ s(age, k = 4),
    skew = ~ s(site, bs = "re"),
    tail = ~1
  )
  expect_equal(spec_covariates(spec4), c("age", "sex", "site"))
})

test_that("NA in an unused column drops no rows and emits no message", {
  dat <- norm_simulate(120, seed = 81)
  spec <- norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1)
  fit_clean <- norm_fit(spec, data = dat, outcomes = "y")
  dat_junk <- dat
  dat_junk$junk <- NA_real_
  dat_junk$marker_01[1:60] <- NA
  fit_junk <- NULL
  expect_no_message(fit_junk <- norm_fit(spec, data = dat_junk, outcomes = "y"))
  expect_equal(fit_junk$models$y$model$df.null, fit_clean$models$y$model$df.null)
  expect_equal(nrow(fit_junk$models$y$model$model), nrow(dat))
  a <- predict(fit_clean, newdata = dat[1:10, ], uncertainty = "conditional")
  b <- predict(fit_junk, newdata = dat[1:10, ], uncertainty = "conditional")
  expect_equal(a$z, b$z, tolerance = 1e-10)
})

test_that("NA in a used covariate drops rows with a message", {
  dat <- norm_simulate(120, seed = 82)
  dat$age[1:7] <- NA
  spec <- norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1)
  expect_message(
    fit <- norm_fit(spec, data = dat, outcomes = "y"),
    "dropped 7 rows"
  )
  expect_equal(nrow(fit$models$y$model$model), nrow(dat) - 7L)
})

test_that("unused columns may be absent from newdata", {
  dat <- norm_simulate(100, seed = 83)
  spec <- norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1)
  fit <- norm_fit(spec, data = dat, outcomes = "y")
  slim <- dat[1:5, c("age", "sex", "y")]
  sc <- predict(fit, newdata = slim, uncertainty = "conditional")
  expect_equal(nrow(sc), 5L)
  expect_true(all(is.finite(sc$z)))
  full <- predict(fit, newdata = dat[1:5, ], uncertainty = "conditional")
  expect_equal(sc$z, full$z)
  expect_error(norm_fit(spec, data = dat[, c("age", "y")], outcomes = "y"), "sex")
})

test_that("support uses only formula covariates", {
  dat <- norm_simulate(150, seed = 84)
  spec <- norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1)
  fit <- norm_fit(spec, data = dat, outcomes = "y")
  expect_equal(names(fit$support_ref$numeric), "age")
  expect_equal(names(fit$support_ref$factor_levels), "sex")
  tgt <- dat[1:5, ]
  tgt$marker_01 <- 1e6
  tgt$marker_02 <- -1e6
  tgt$site <- factor("Q")
  st <- norm_support(fit, tgt)
  expect_s3_class(st, "tbl_df")
  expect_equal(names(st), c(".row", "support", "d2"))
  expect_true(all(st$support == "in"))
  tgt$age <- 150
  expect_true(all(norm_support(fit, tgt)$support == "out"))
})

test_that("the model card lists only formula covariates", {
  dat <- norm_simulate(100, seed = 85)
  spec <- norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1)
  fit <- norm_fit(spec, data = dat, outcomes = "y")
  ref <- norm_reference(fit)
  expect_equal(ref$covariates, c("age", "sex"))
  expect_equal(names(ref$support_ref$numeric), "age")
  out <- paste(cli::cli_fmt(print(ref)), collapse = "\n")
  expect_match(out, "covariates: age")
  expect_no_match(out, "marker")
  expect_no_match(out, "site")
})

test_that("fits are reproducible under a parallel future plan", {
  skip_if_not_installed("future.apply")
  skip_on_cran()
  skip_if_not(future::supportsMulticore(), "multicore futures unsupported")
  dat <- norm_simulate(100, seed = 86)
  spec <- norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1)
  seq_fit <- norm_fit(spec, data = dat, outcomes = c("y", "marker_01"))
  old <- suppressWarnings(future::plan(future::multicore, workers = 2))
  on.exit(future::plan(old), add = TRUE)
  par_fit <- norm_fit(spec, data = dat, outcomes = c("y", "marker_01"))
  a <- predict(seq_fit, newdata = dat[1:5, ], uncertainty = "conditional")
  b <- predict(par_fit, newdata = dat[1:5, ], uncertainty = "conditional")
  expect_equal(a$z, b$z, tolerance = 1e-10)
})

test_that("named smooth arguments such as k = kk are not covariates", {
  kk <- 6
  spec <- norm_spec(norm_gaussian(), ~ s(age, k = kk) + sex + s(site, bs = "re"),
                    scale = ~ s(age, k = kk, by = sex))
  expect_setequal(spec_covariates(spec), c("age", "sex", "site"))
  dat <- norm_simulate(150, seed = 14)
  fit <- norm_fit(spec, dat, "y")
  expect_equal(fit$covariates, c("age", "sex", "site"))
  expect_equal(unname(fit_statuses(fit)), "ok")
  sc <- predict(fit, dat[1:3, ], uncertainty = "conditional")
  expect_true(all(is.finite(sc$z)))
})
