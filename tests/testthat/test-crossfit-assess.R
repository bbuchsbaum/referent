test_that("cross-fitted scores are out of sample, cover every row once, and are calibrated", {
  dat <- norm_simulate(400, seed = 11)
  dat$participant_id <- rep(1:100, each = 4)
  withr::with_seed(11, {
    cf <- norm_crossfit(
      simple_spec(),
      data = dat,
      outcomes = "y",
      folds = 4,
      cluster = participant_id
    )
  })
  expect_s3_class(cf, "norm_scores")
  expect_false(any(cf$.in_sample))
  expect_s3_class(attr(cf, "deployment"), "norm_fit")
  expect_equal(sort(cf$.row), seq_len(nrow(dat)))
  expect_equal(sort(unique(cf$.fold)), 1:4)
  # clusters stay together
  fold_of <- attr(cf, "folds")
  expect_true(all(tapply(fold_of, dat$participant_id, function(f) length(unique(f))) == 1L))
  # the generator is Gaussian, so out-of-fold scores must be calibrated
  m <- assess_marginal(cf)
  expect_lt(abs(m$mean_z), 0.15)
  expect_lt(abs(m$var_z - 1), 0.2)
  expect_lt(abs(m$cover_95 - 0.95), 0.04)
  expect_true(acceptable_calibration(list(marginal = m, tail = assess_tail(cf))))
  # per-row CRPS is the closed-form Gaussian CRPS of the out-of-fold predictive,
  # whose expectation under calibration is sigma / sqrt(pi) with sigma = 1.3
  expect_true(all(is.finite(cf$crps)))
  expect_lt(abs(mean(cf$crps) - 1.3 / sqrt(pi)), 0.1)
})

test_that("norm_assess recovers nominal coverage and a positive log-score gain", {
  train <- norm_simulate(300, seed = 13)
  test <- norm_simulate(400, seed = 14)
  fit <- norm_fit(simple_spec(), data = train, outcomes = "y")
  a <- norm_assess(fit, newdata = test)
  expect_s3_class(a, "norm_assessment")
  expect_false(a$in_sample)
  expect_equal(a$n, 400L)
  m <- a$marginal
  expect_lt(abs(m$cover_50 - 0.50), 0.07)
  expect_lt(abs(m$cover_90 - 0.90), 0.05)
  expect_lt(abs(m$cover_95 - 0.95), 0.04)
  expect_lt(abs(m$var_z - 1), 0.2)
  ov <- a$overall
  # the conditional model must beat an unconditional Gaussian in log score
  expect_gt(ov$standardized_log_score, 0.1)
  # the expected CRPS of a calibrated Gaussian predictive is sigma / sqrt(pi)
  expect_lt(abs(ov$crps - 1.3 / sqrt(pi)), 0.08)
  expect_gt(ov$cor, 0.5)
  tl <- a$tail[a$tail$.outcome == "y", ]
  expect_true(all(abs(tl$observed - tl$expected) <= 3 * tl$se + 0.005))
})

# The unconditional Gaussian that `standardized_log_score` is measured against
# must come from the reference sample, not from the sample being scored.
# Otherwise the column is an oracle score: it moves when the held-out set
# changes, and it is not the MSLL that the normative-modelling literature (and
# PCNtoolkit's `get_statistics_df()`) reports.

msll_longhand <- function(scores, y_ref) {
  sd0 <- sqrt(mean((y_ref - mean(y_ref))^2)) # population sd, denominator n
  -mean(scores$log_density) +
    mean(stats::dnorm(scores$observed, mean(y_ref), sd0, log = TRUE))
}

test_that("standardized_log_score is the negative MSLL against the reference sample", {
  ref <- norm_simulate(300, seed = 21)
  new <- norm_simulate(120, seed = 22)
  fit <- norm_fit(norm_spec(norm_gaussian(), ~ s(age, k = 5) + sex), ref, "y")

  for (rows in list(seq_len(120), seq_len(40))) {
    a <- norm_assess(fit, new[rows, ])
    sc <- predict(fit, new[rows, ], type = "scores", uncertainty = "conditional")
    expect_equal(
      -a$overall$standardized_log_score,
      msll_longhand(sc, ref$y),
      tolerance = 1e-12
    )
  }
})

test_that("the reference baseline survives freezing to a bundle", {
  ref <- norm_simulate(200, seed = 23)
  new <- norm_simulate(80, seed = 24)
  fit <- norm_fit(norm_spec(norm_gaussian(), ~ s(age, k = 5) + sex), ref, "y")
  thawed <- norm_reference(fit)
  expect_equal(norm_assess(thawed, new)$overall, norm_assess(fit, new)$overall,
    tolerance = 1e-12
  )
})

test_that("norm_assess() honours the uncertainty argument", {
  ref <- norm_simulate(200, seed = 25)
  new <- norm_simulate(80, seed = 26)
  fit <- norm_fit(norm_spec(norm_gaussian(), ~ s(age, k = 5) + sex), ref, "y")

  expect_equal(norm_assess(fit, new)$overall, norm_assess(fit, new, uncertainty = "conditional")$overall)
  total <- norm_assess(fit, new, uncertainty = "total")
  expect_equal(
    -total$overall$standardized_log_score,
    msll_longhand(predict(fit, new, type = "scores", uncertainty = "total"), ref$y),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(
    total$overall$mean_log_score,
    norm_assess(fit, new)$overall$mean_log_score
  )))
})
