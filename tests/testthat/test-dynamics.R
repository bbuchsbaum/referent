long_spec <- function() {
  norm_spec(family = norm_gaussian(), location = ~ s(age, k = 5) + sex, scale = ~1)
}

test_that("paper-reduction law: innovation_z equals Z-gain", {
  r <- 0.6
  z1 <- 1.2
  z2 <- -0.4
  expect_gain <- (z2 - r * z1) / sqrt(1 - r^2)
  process <- norm_process("stable")
  psi <- list(tau_b = sqrt(r), tau_g = 0, sigma_e = sqrt(1 - r), ell = Inf)
  R <- process_correlation(c(0, 1), psi, process)
  expect_equal(R[1, 2], r, tolerance = 1e-8)
  cond <- condition_history(z1, 0, 1, list(psi = psi, process = process))
  innov <- (z2 - cond$m) / cond$s
  expect_equal(innov, expect_gain, tolerance = 1e-8)
})

test_that("difference-reduction law: change_z equals standardized difference", {
  r <- 0.5
  z1 <- 0.8
  z2 <- -0.2
  expect_d <- (z2 - z1) / sqrt(2 * (1 - r))
  process <- norm_process(ell = Inf)
  expect_identical(process$name, "stable")
  psi <- list(tau_b = sqrt(r), tau_g = 0, sigma_e = sqrt(1 - r), ell = Inf)
  r12 <- process_kernel(1, psi, process)
  expect_equal((z2 - z1) / sqrt(2 * (1 - r12)), expect_d, tolerance = 1e-8)
})

test_that("history law: added history cannot increase conditional variance", {
  fitted <- list(process = norm_process(ell = 4),
                 psi = list(tau_b = 0.5, tau_g = 0.4, sigma_e = 0.2, ell = 4))
  s1 <- condition_history(0.3, 10, 12, fitted)$s
  s2 <- condition_history(c(0.1, 0.3), c(8, 10), 12, fitted)$s
  expect_lte(s2, s1 + 1e-10)
  # vectorised over forecast times
  both <- condition_history(c(0.1, 0.3), c(8, 10), c(12, 15), fitted)
  expect_equal(both$s[[1]], s2)
  expect_equal(both$m[[1]], condition_history(c(0.1, 0.3), c(8, 10), 12, fitted)$m)
  expect_gt(both$s[[2]], both$s[[1]])
})

test_that("correlation matrices are positive semidefinite", {
  process <- norm_process(ell = 3)
  psi <- list(tau_b = 0.4, tau_g = 0.5, sigma_e = 0.2, ell = 3)
  times <- c(0, 0.4, 1.7, 6)
  k <- process_correlation(times, psi, process)
  expect_equal(diag(k), rep(1, 4))
  ev <- eigen(k, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev >= -1e-10))
})

test_that("batched likelihood equals the per-subject Gaussian likelihood", {
  set.seed(3)
  m <- 6
  n <- 4
  r <- array(0, c(m, n, n))
  z <- matrix(stats::rnorm(m * n), m, n)
  ref <- 0
  for (i in seq_len(m)) {
    t <- sort(stats::runif(n, 0, 10))
    R <- 0.5 + 0.4 * matern32_kernel(abs(outer(t, t, `-`)), 3)
    diag(R) <- 1
    r[i, , ] <- R
    ch <- chol(R)
    w <- backsolve(ch, z[i, ], transpose = TRUE)
    ref <- ref - 0.5 * (n * log(2 * pi) + 2 * sum(log(diag(ch))) + sum(w^2))
  }
  expect_equal(batched_loglik(z, r), ref, tolerance = 1e-10)
})

test_that("simplex parameterisation gives unit total variance", {
  psi <- process_par(c(0.3, -1.2, 0.5), norm_process())
  expect_equal(psi$tau_b^2 + psi$tau_g^2 + psi$sigma_e^2, 1, tolerance = 1e-12)
  psi_s <- process_par(0.7, norm_process("stable"))
  expect_equal(psi_s$tau_b^2 + psi_s$sigma_e^2, 1, tolerance = 1e-12)
  expect_identical(psi_s$ell, Inf)
})

test_that("kernel recovery: r(lag) at the median lag within 0.08 of truth on irregular lags", {
  skip_on_cran()
  dat <- norm_simulate(1800, kind = "longitudinal", seed = 11)
  truth <- attr(dat, "truth")
  fit <- norm_fit(long_spec(), data = dat, outcomes = "y")
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age)
  pr <- dyn$processes$y
  expect_true(pr$identified)
  expect_true(pr$ell_identified)
  expect_identical(pr$process$name, "matern32")
  expect_identical(dyn$z_source, "out_of_fold")
  expect_equal(sum(pr$components), 1, tolerance = 1e-10)
  expect_true(all(is.finite(pr$se)))
  expect_true(is.finite(pr$r_median_se) && pr$r_median_se > 0)
  r_true <- truth$correlation(pr$median_lag)
  expect_lt(abs(pr$r_median - r_true), 0.08)
  comp <- dyn$components
  expect_equal(comp$r_median_lag, pr$r_median)
  expect_match(print_text(dyn), "matern32")
})

test_that("out-of-sample calibration of innovation_z and change_z on irregular 3-visit data", {
  skip_on_cran()
  train <- norm_simulate(1800, kind = "longitudinal", seed = 12)
  test <- norm_simulate(3000, kind = "longitudinal", seed = 13)
  expect_gt(attr(train, "truth")$sigma_e, 0)
  fit <- norm_fit(long_spec(), data = train, outcomes = "y")
  dyn <- norm_dynamics(fit, data = train, id = participant_id, time = age)
  tr <- norm_transition(dyn, data = test, id = participant_id, time = age)
  innov <- tr$innovation_z[is.finite(tr$innovation_z)]
  change <- tr$change_z[is.finite(tr$change_z)]
  expect_gt(length(innov), 1800)
  expect_lt(abs(stats::var(innov) - 1), 0.1)
  expect_lt(abs(stats::var(change) - 1), 0.1)
  expect_lt(abs(mean(innov)), 0.08)
  expect_lt(abs(false_positive_rate(change) - 0.05), 0.015)
  expect_lt(abs(false_positive_rate(innov) - 0.05), 0.015)
  expect_equal(tr$change_centile, stats::pnorm(tr$change_z))
  expect_equal(tr$velocity_centile, stats::pnorm(tr$innovation_z))
  expect_false("z_gain" %in% names(tr))
  expect_equal(tr$measurement_sd[[1]], dyn$processes$y$psi$sigma_e)
  expect_true(all(tr$history_n[tr$history_n > 1] >= 2))
})

test_that("stable model on fixed-lag two-visit data recovers r with ell unidentified", {
  skip_on_cran()
  dat <- simulate_two_visit(600, r = 0.6, lag = 2, seed = 5)
  fit <- norm_fit(long_spec(), data = dat, outcomes = "y")
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age)
  pr <- dyn$processes$y
  expect_true(dyn$identifiability$fixed_lag)
  expect_true(pr$identified)
  expect_false(pr$ell_identified)
  expect_identical(pr$process$name, "stable")
  expect_true(is.na(dyn$components$ell))
  expect_lt(abs(pr$r_median - 0.6), 0.05)
  expect_lt(abs(process_kernel(2, pr$psi, pr$process) - 0.6), 0.05)
  expect_true(is.finite(pr$r_median_se))
  expect_false(dyn$components$ell_identified)
  expect_match(print_text(dyn), "stable")
})

test_that("stable and Matern kernels agree on r(lag) for near-fixed lags", {
  skip_on_cran()
  dat <- simulate_two_visit(600, r = 0.55, lag = 2, seed = 6)
  set.seed(61)
  jitter <- stats::runif(nrow(dat) / 2, -0.4, 0.4)
  dat$age[seq(2, nrow(dat), by = 2)] <- dat$age[seq(2, nrow(dat), by = 2)] + jitter
  fit <- norm_fit(long_spec(), data = dat, outcomes = "y")
  dyn_m <- norm_dynamics(fit, data = dat, id = participant_id, time = age, crossfit = 0)
  dyn_s <- norm_dynamics(fit, data = dat, id = participant_id, time = age,
                         process = norm_process("stable"), crossfit = 0)
  expect_identical(dyn_m$processes$y$process$name, "matern32")
  expect_identical(dyn_s$processes$y$process$name, "stable")
  expect_identical(dyn_m$z_source, "in_sample")
  r_m <- process_kernel(2, dyn_m$processes$y$psi, dyn_m$processes$y$process)
  r_s <- process_kernel(2, dyn_s$processes$y$psi, dyn_s$processes$y$process)
  expect_lt(abs(r_m - r_s), 0.03)
  expect_lt(abs(r_s - 0.55), 0.06)
})

test_that("unidentified dynamics give NA history-conditioned quantities and no forecast", {
  dat <- norm_simulate(200, seed = 77)
  dat$participant_id <- seq_len(nrow(dat))
  # three subjects with a repeat: not enough to identify dependence
  dat$participant_id[1:3] <- dat$participant_id[4:6]
  fit <- norm_fit(norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
                  data = dat, outcomes = "y")
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age, crossfit = 0)
  pr <- dyn$processes$y
  expect_false(pr$identified)
  expect_null(pr$psi)
  expect_match(pr$reason, "too few")
  expect_true(is.na(dyn$components$stable))
  tr <- norm_transition(dyn, data = dat, id = participant_id, time = age)
  expect_equal(nrow(tr), 3L)
  for (col in c("innovation_z", "change_z", "velocity_centile", "change_centile",
                "expected_velocity", "velocity_lower", "velocity_upper", "measurement_sd")) {
    expect_true(all(is.na(tr[[col]])), info = col)
  }
  expect_true(all(is.finite(tr$observed_velocity)))
  expect_true(all(tr$support == "unidentified"))
  expect_false(any(tr$calibrated))
  hist <- dat[dat$participant_id == dat$participant_id[[1]], ]
  expect_error(norm_forecast(dyn, history = hist, times = max(hist$age) + 1), "not identified")
  expect_true(all(is.na(fortify_kernel(dyn)$correlation)))
  expect_match(print_text(dyn), "not identified")
})

test_that("norm_transition handles no transitions and duplicate visit times", {
  dat <- norm_simulate(300, kind = "longitudinal", seed = 31)
  fit <- norm_fit(long_spec(), data = dat, outcomes = "y")
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age, crossfit = 0)
  single <- dat[!duplicated(dat$participant_id), ]
  tr0 <- norm_transition(dyn, data = single, id = participant_id, time = age)
  expect_s3_class(tr0, "norm_transition")
  expect_equal(nrow(tr0), 0L)
  expect_identical(names(tr0), names(transition_template()))
  dup <- dat[dat$participant_id %in% dat$participant_id[1:2], ]
  dup$age[2] <- dup$age[1]
  expect_warning(
    tr_d <- norm_transition(dyn, data = dup, id = participant_id, time = age),
    "duplicate visit time"
  )
  zero <- tr_d$.dt == 0
  expect_true(any(zero))
  expect_true(all(is.na(tr_d$observed_velocity[zero])))
  expect_true(all(is.na(tr_d$innovation_z[zero])))
  expect_true(all(is.finite(tr_d$observed_velocity[!zero])))
})

test_that("temporal support flags extrapolated lags under a fixed-lag reference", {
  dat <- simulate_two_visit(200, r = 0.6, lag = 2, seed = 8)
  fit <- norm_fit(long_spec(), data = dat, outcomes = "y")
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age, crossfit = 0)
  expect_identical(classify_temporal_support(dyn, 40, 42, 1), "in")
  expect_identical(classify_temporal_support(dyn, 40, 50, 1), "extrapolated_lag")
  expect_identical(classify_temporal_support(dyn, 40, 40.5, 1), "extrapolated_lag")
  hist <- dat[dat$participant_id == 1, ]
  fc <- norm_forecast(dyn, history = hist, times = max(hist$age) + c(2, 10))
  expect_equal(fc$history_n, 2L)
  expect_identical(fc$lag_support, c("in", "extrapolated_lag"))
  expect_identical(fc$summary$support, fc$lag_support)
})

test_that("norm_forecast requires the time column and drops non-finite history", {
  dat <- norm_simulate(300, kind = "longitudinal", seed = 32)
  fit <- norm_fit(long_spec(), data = dat, outcomes = "y")
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age, crossfit = 0)
  hist <- dat[dat$participant_id == dat$participant_id[[1]], ]
  bad <- hist
  names(bad)[names(bad) == "age"] <- "years"
  bad$age_at_scan <- bad$years
  expect_error(norm_forecast(dyn, history = bad, times = 60), "time column")
  hist$y[[1]] <- NA_real_
  expect_message(fc <- norm_forecast(dyn, history = hist, times = max(hist$age) + 1),
                 "Dropped 1 history row")
  expect_equal(fc$history_n, nrow(hist) - 1L)
  expect_true(all(is.finite(fc$summary$median)))
  expect_error(norm_forecast(dyn, history = hist, times = 60, outcome = "nope"), "not an outcome")
})

test_that("time-unit law: velocity rescales, innovation does not", {
  skip_on_cran()
  dat <- norm_simulate(300, kind = "longitudinal", seed = 24)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "y")
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age, crossfit = 0)
  tr <- norm_transition(dyn, data = dat, id = participant_id, time = age)
  expect_equal(tr$observed_velocity * tr$.dt, tr$observed_change, tolerance = 1e-10)
  expect_equal(tr$expected_velocity * tr$.dt, tr$expected_change, tolerance = 1e-10)
  dat_m <- dat
  dat_m$age_months <- dat$age * 12
  dyn_m <- dyn
  dyn_m$time_name <- "age_months"
  dyn_m$time_range <- dyn$time_range * 12
  dyn_m$lag_range <- dyn$lag_range * 12
  for (nm in names(dyn_m$processes)) {
    dyn_m$processes[[nm]]$psi$ell <- dyn$processes[[nm]]$psi$ell * 12
  }
  tr_m <- norm_transition(dyn_m, data = dat_m, id = participant_id, time = age_months)
  expect_equal(tr_m$innovation_z, tr$innovation_z, tolerance = 1e-8)
  expect_equal(tr_m$change_z, tr$change_z, tolerance = 1e-8)
  expect_equal(tr_m$observed_velocity * 12, tr$observed_velocity, tolerance = 1e-8)
})

test_that("forecast density law holds after conditioning", {
  dat <- norm_simulate(240, kind = "longitudinal", seed = 25)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age, crossfit = 0)
  hist <- dat[dat$participant_id == dat$participant_id[[1]], ]
  fc <- norm_forecast(dyn, history = hist, times = max(hist$age) + 1)
  d <- fc$dist[1]
  expect_s3_class(vctrs::vec_data(d)[[1]], "dist_conditioned")
  p <- c(0.2, 0.5, 0.8)
  expect_equal(dist_cdf(d, dist_quantile(d, p)), p, tolerance = 1e-6)
})

test_that("norm_derivative reproduces a linear slope and its standard error", {
  dat <- norm_simulate(300, seed = 41)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = c("y", "marker_01")
  )
  grid <- data.frame(age = c(30, 50, 70), sex = factor("F", levels = c("F", "M")))
  dv <- norm_derivative(fit, grid, with_respect_to = "age", centiles = c(0.1, 0.5))
  expect_s3_class(dv, "norm_derivative")
  expect_setequal(unique(dv$.outcome), c("y", "marker_01"))
  expect_equal(nrow(dv), 2 * 3 * 2)
  m <- fit$models$y$model
  slope <- unname(stats::coef(m)[["age"]])
  j <- which(names(stats::coef(m)) == "age"); slope_se <- sqrt(m$Vp[j, j])
  got <- dv[dv$.outcome == "y", ]
  expect_equal(got$chart_velocity, rep(slope, nrow(got)), tolerance = 1e-6)
  expect_equal(got$chart_velocity_se, rep(slope_se, nrow(got)), tolerance = 1e-4)
  dv1 <- norm_derivative(fit, grid, with_respect_to = "age", outcomes = "marker_01")
  expect_true(all(dv1$.outcome == "marker_01"))
  expect_true(all(dv1$chart_velocity_se > 0))
})

test_that("norm_derivative gives finite SEs for smooth SHASH fits", {
  skip_on_cran()
  dat <- norm_simulate(500, kind = "shash", seed = 42)
  fit <- norm_fit(simple_spec("shash"), data = dat, outcomes = "y")
  grid <- data.frame(age = c(30, 50, 70), sex = factor("M", levels = c("F", "M")))
  dv <- norm_derivative(fit, grid, with_respect_to = "age")
  expect_true(all(is.finite(dv$chart_velocity)))
  expect_true(all(is.finite(dv$chart_velocity_se) & dv$chart_velocity_se > 0))
  # the upper centile is steeper where the scale grows with age only if scale varies;
  # here the median derivative must equal the location derivative
  med <- dv[dv$centile == 0.5, ]
  expect_equal(length(unique(round(med$chart_velocity, 8))), 3L)
})
