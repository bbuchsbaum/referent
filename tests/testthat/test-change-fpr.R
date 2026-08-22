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
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y",
    id = participant_id
  )
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age)
  list(
    dyn = dyn,
    transition = norm_transition(dyn, data = dat, id = participant_id, time = age)
  )
}

test_that("analytic change_z holds 5% FPR while naive Z2-Z1 does not", {
  set.seed(70)
  threshold <- stats::qnorm(0.975)
  n <- 6000
  rows <- lapply(c(0.1, 0.4, 0.7, 0.9), function(r) {
    z1 <- stats::rnorm(n)
    z2 <- r * z1 + sqrt(1 - r^2) * stats::rnorm(n)
    naive <- z2 - z1
    change <- naive / sqrt(2 * (1 - r))
    tibble::tibble(
      r = r,
      fpr_naive = false_positive_rate(naive, threshold),
      fpr_change = false_positive_rate(change, threshold),
      expected_naive = expected_naive_fpr(r, threshold)
    )
  })
  tab <- dplyr_bind(rows)
  expect_true(all(abs(tab$fpr_change - 0.05) < 0.015))
  expect_true(all(abs(tab$fpr_naive - tab$expected_naive) < 0.02))
  expect_gt(max(abs(tab$fpr_naive - 0.05)), 0.04)
})

test_that("out-of-sample change_z FPR stays near nominal across correlation and lag", {
  skip_on_cran()
  settings <- list(
    list(r = 0.75, lag = 1.5, sigma_e = 0, seed = 71),
    list(r = 0.25, lag = 6, sigma_e = 0, seed = 72),
    list(r = 0.6, lag = 2, sigma_e = 0.35, seed = 73)
  )
  rows <- lapply(settings, function(cfg) {
    train <- do.call(simulate_two_visit, c(list(n_id = 500), cfg))
    test <- do.call(simulate_two_visit, c(list(n_id = 1000), modifyList(cfg, list(seed = cfg$seed + 100))))
    got <- fit_two_visit_transition(train)
    pr <- got$dyn$processes$y
    expect_true(pr$identified)
    expect_identical(pr$process$name, "stable")
    tr <- norm_transition(got$dyn, data = test, id = participant_id, time = age)
    naive <- naive_change_z(tr)
    tibble::tibble(
      r = cfg$r,
      lag = cfg$lag,
      sigma_e = cfg$sigma_e,
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

test_that("new subject ids are not new_group support", {
  dat <- norm_simulate(50, seed = 75)
  dat$participant_id <- paste0("S", seq_len(nrow(dat)))
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y",
    id = participant_id
  )
  tgt <- dat[1:8, ]
  tgt$participant_id <- paste0("T", seq_len(nrow(tgt)))
  sc <- predict(fit, newdata = tgt, uncertainty = "conditional")
  expect_false(any(sc$support == "new_group"))
})
