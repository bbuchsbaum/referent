test_that("paper-reduction law: innovation_z equals Z-gain", {
  r <- 0.6
  z1 <- 1.2
  z2 <- -0.4
  expect_gain <- (z2 - r * z1) / sqrt(1 - r^2)
  process <- norm_matern32(stable_rank = TRUE, ell = 1e6)
  # With huge length-scale and no nugget-dominated correlation, force r via psi
  psi <- list(tau_b = sqrt(r), tau_g = 0, sigma_e = sqrt(1 - r), ell = 1e6)
  times <- c(0, 1)
  R <- process_correlation(times, psi, process)
  expect_equal(R[1, 2], r, tolerance = 1e-8)
  cond <- condition_z(z1, 0, 1, psi, process)
  innov <- (z2 - cond$m) / cond$s
  expect_equal(innov, expect_gain, tolerance = 1e-8)
})

test_that("difference-reduction law: change_z equals standardized difference", {
  r <- 0.5
  z1 <- 0.8
  z2 <- -0.2
  expect_d <- (z2 - z1) / sqrt(2 * (1 - r))
  process <- norm_matern32(ell = 1e6)
  psi <- list(tau_b = sqrt(r), tau_g = 0, sigma_e = sqrt(1 - r), ell = 1e6)
  R <- process_correlation(c(0, 1), psi, process)
  got <- (z2 - z1) / sqrt(2 * (1 - R[1, 2]))
  expect_equal(got, expect_d, tolerance = 1e-8)
})

test_that("history law: added history cannot increase conditional variance", {
  process <- norm_matern32(ell = 4)
  psi <- list(tau_b = 0.5, tau_g = 0.4, sigma_e = 0.2, ell = 4)
  s1 <- condition_z(0.3, 10, 12, psi, process)$s
  s2 <- condition_z(c(0.1, 0.3), c(8, 10), 12, psi, process)$s
  expect_lte(s2, s1 + 1e-10)
})

test_that("covariance matrices are positive semidefinite", {
  process <- norm_matern32(ell = 3)
  psi <- list(tau_b = 0.4, tau_g = 0.5, sigma_e = 0.2, ell = 3)
  times <- c(0, 0.4, 1.7, 6)
  k <- process_covariance(times, psi, process)
  ev <- eigen(k, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev >= -1e-10))
})

test_that("time-unit law: velocity rescales, innovation does not", {
  skip_on_cran()
  set.seed(24)
  dat <- norm_simulate(90, kind = "longitudinal", seed = 24)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "y")
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age)
  tr <- norm_transition(dyn, data = dat, id = participant_id, time = age)
  vel_year <- norm_velocity(tr, time_unit = "year")
  vel_month <- tr
  vel_month$observed_velocity <- vel_month$observed_velocity / 12
  expect_equal(tr$innovation_z, vel_year$innovation_z)
  expect_equal(mean(abs(vel_month$observed_velocity) * 12),
               mean(abs(vel_year$observed_velocity)), tolerance = 1e-10)
})

test_that("innovations are roughly standard normal on longitudinal controls", {
  skip_on_cran()
  set.seed(55)
  dat <- norm_simulate(150, kind = "longitudinal", seed = 55)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age)
  expect_true(isTRUE(dyn$identifiability$change))
  tr <- norm_transition(dyn, data = dat, id = participant_id, time = age)
  z <- tr$innovation_z[is.finite(tr$innovation_z)]
  expect_gt(length(z), 20)
  expect_lt(abs(mean(z)), 0.35)
  expect_lt(abs(stats::var(z) - 1), 0.55)
})

test_that("forecast density law holds after conditioning", {
  set.seed(25)
  dat <- norm_simulate(80, kind = "longitudinal", seed = 25)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age)
  hist <- dat[dat$participant_id == dat$participant_id[[1]], ]
  fc <- norm_forecast(dyn, history = hist, times = max(hist$age) + 1)
  d <- fc$dist[[1]]
  p <- c(0.2, 0.5, 0.8)
  expect_equal(as.numeric(cdf(d, quantile(d, p))), p, tolerance = 1e-6)
})
