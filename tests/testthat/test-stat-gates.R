test_that("held-out Gaussian coverage stays near nominal", {
  set.seed(50)
  train <- norm_simulate(280, seed = 50)
  test <- norm_simulate(220, seed = 51)
  fit <- norm_fit(simple_spec(), data = train, outcomes = "y")
  a <- norm_assess(fit, newdata = test)
  expect_lt(abs(a$marginal$cover_50 - 0.50), 0.12)
  expect_lt(abs(a$marginal$cover_80 - 0.80), 0.10)
  expect_lt(abs(a$marginal$cover_90 - 0.90), 0.10)
  expect_lt(abs(a$marginal$cover_95 - 0.95), 0.08)
  expect_lt(abs(a$marginal$var_z - 1), 0.35)
})

test_that("in-sample Z scores are over-shrunk relative to cross-fit", {
  set.seed(52)
  dat <- norm_simulate(160, seed = 52)
  dat$participant_id <- rep(seq_len(40), each = 4)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "y")
  ins <- predict(fit, newdata = dat, uncertainty = "conditional",
                 allow_extrapolation = TRUE)
  expect_true(all(ins$.in_sample))
  cf <- norm_crossfit(simple_spec(), data = dat, outcomes = "y",
                      folds = 4, cluster = participant_id)
  expect_false(any(cf$.in_sample))
  expect_lte(stats::var(ins$z, na.rm = TRUE), stats::var(cf$z, na.rm = TRUE) + 0.05)
})

test_that("dynamics refuse unusual-change labels without repeats", {
  set.seed(54)
  dat <- norm_simulate(40, seed = 54)
  dat$participant_id <- seq_len(nrow(dat))
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age)
  expect_false(isTRUE(dyn$identifiability$change))
  # two fake visits per person still inherit the unidentified process
  dat2 <- rbind(dat, dat)
  dat2$age[seq_len(nrow(dat)) + nrow(dat)] <- dat$age + 1
  dat2$participant_id <- c(dat$participant_id, dat$participant_id)
  tr <- norm_transition(dyn, data = dat2, id = participant_id, time = age)
  expect_true(all(is.na(tr$innovation_z)))
  expect_true(all(is.na(tr$change_z)))
  expect_true(all(tr$support == "unidentified"))
})

test_that("SHASH beats Gaussian on strongly skewed held-out data", {
  skip_on_cran()
  train <- simulate_skewed(280, seed = 61, skew = 1.4, tail = 0.7)
  test <- simulate_skewed(220, seed = 62, skew = 1.4, tail = 0.7)

  fit_g <- norm_fit(simple_spec(scale = TRUE), data = train, outcomes = "y")
  fit_s <- norm_fit(simple_spec("shash", scale = TRUE), data = train, outcomes = "y")
  skip_if(
    identical(fit_s$models$y$status, "nonconverged"),
    "SHASH fit did not converge"
  )

  a_g <- norm_assess(fit_g, newdata = test)
  a_s <- norm_assess(fit_s, newdata = test)

  expect_gt(a_s$overall$mean_log_score, a_g$overall$mean_log_score)
  expect_lt(pit_ks(a_s$scores$centile), pit_ks(a_g$scores$centile))
  expect_lt(worm_rmse(a_s$scores$z), worm_rmse(a_g$scores$z))
  expect_lt(abs(a_s$marginal$skew_z), abs(a_g$marginal$skew_z))
})

test_that("Gaussian remains competitive when the truth is Gaussian", {
  skip_on_cran()
  train <- norm_simulate(240, kind = "gaussian", scale = "age", seed = 63)
  test <- norm_simulate(180, kind = "gaussian", scale = "age", seed = 64)
  fit_g <- norm_fit(simple_spec(scale = TRUE), data = train, outcomes = "y")
  fit_s <- norm_fit(simple_spec("shash", scale = TRUE), data = train, outcomes = "y")
  skip_if(
    identical(fit_s$models$y$status, "nonconverged"),
    "SHASH fit did not converge"
  )
  a_g <- norm_assess(fit_g, newdata = test)
  a_s <- norm_assess(fit_s, newdata = test)
  # Extra SHASH shape parameters should not win by a large margin.
  expect_gt(a_g$overall$mean_log_score, a_s$overall$mean_log_score - 0.04)
  expect_lt(abs(a_g$marginal$var_z - 1), 0.35)
})
