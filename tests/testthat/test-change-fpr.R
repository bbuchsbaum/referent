clamp_prob <- function(p, eps = 1e-12) pmin(pmax(as.numeric(p), eps), 1 - eps)

naive_change_z <- function(transition) {
  z1 <- stats::qnorm(clamp_prob(transition$start_centile))
  z2 <- stats::qnorm(clamp_prob(transition$end_centile))
  z2 - z1
}

expected_naive_fpr <- function(r, threshold = stats::qnorm(0.975)) {
  sd_diff <- sqrt(pmax(2 * (1 - r), 1e-8))
  2 * stats::pnorm(-threshold / sd_diff)
}

fit_two_visit_transition <- function(dat) {
  fit <- ref_fit(
    ref_spec(family = ref_gaussian(), location = ~ s(age, k = 5) + sex, scale = ~1),
    data = dat,
    outcomes = "y",
    id = participant_id
  )
  dyn <- ref_dynamics(fit, data = dat, id = participant_id, time = age)
  list(
    dyn = dyn,
    transition = ref_transition(dyn, data = dat, id = participant_id, time = age)
  )
}

test_that("out-of-sample change_z FPR stays near nominal across correlation and lag", {
  skip_on_cran()
  settings <- list(
    list(r = 0.75, lag = 1.5, seed = 71),
    list(r = 0.25, lag = 6, seed = 72),
    list(r = 0.53, lag = 2, seed = 73)
  )
  rows <- lapply(settings, function(cfg) {
    train <- do.call(simulate_two_visit, c(list(n_id = 500), cfg))
    test <- do.call(simulate_two_visit, c(list(n_id = 1000), modifyList(cfg, list(seed = cfg$seed + 100))))
    got <- fit_two_visit_transition(train)
    pr <- got$dyn$processes$y
    expect_true(pr$identified)
    expect_identical(pr$process$name, "stable")
    tr <- ref_transition(got$dyn, data = test, id = participant_id, time = age)
    naive <- naive_change_z(tr)
    tibble::tibble(
      r = cfg$r,
      lag = cfg$lag,
      n = sum(is.finite(tr$change_z)),
      var_change = stats::var(tr$change_z, na.rm = TRUE),
      var_innov = stats::var(tr$innovation_z, na.rm = TRUE),
      fpr_change = false_positive_rate(tr$change_z),
      fpr_naive = false_positive_rate(naive)
    )
  })
  tab <- dplyr_bind(rows)
  expect_true(all(tab$n >= 1000))
  expect_true(all(abs(tab$var_change - 1) < 0.1))
  expect_true(all(abs(tab$var_innov - 1) < 0.1))
  expect_true(all(abs(tab$fpr_change - 0.05) < 0.02))
  expect_gt(max(abs(tab$fpr_naive - 0.05)), abs(tab$fpr_change[which.max(abs(tab$fpr_naive - 0.05))] - 0.05))
  expect_gt(tab$fpr_naive[tab$r == 0.25], 0.07)
  expect_lt(tab$fpr_naive[tab$r == 0.75], 0.04)
})
