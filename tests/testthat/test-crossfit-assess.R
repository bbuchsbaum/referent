test_that("cross-fitted scores are out of sample, cover every row once, and are calibrated", {
  dat <- ref_simulate(400, seed = 11)
  dat$participant_id <- rep(1:100, each = 4)
  withr::with_seed(11, {
    cf <- ref_crossfit(
      simple_spec(),
      data = dat,
      outcomes = "y",
      folds = 4,
      cluster = participant_id
    )
  })
  expect_s3_class(cf, "ref_scores")
  expect_false(any(cf$.in_sample))
  expect_s3_class(attr(cf, "deployment"), "ref_fit")
  expect_identical(attr(cf, "uncertainty"), "conditional")
  expect_true(is.character(attr(cf, "data_hash")))
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

test_that("assessment requires held-out provenance and hashes full row content", {
  train <- ref_simulate(140, seed = 81)
  fit <- ref_fit(simple_spec(), train, "y")

  expect_error(ref_assess(fit, train), "training data")
  allowed <- ref_assess(fit, train, allow_in_sample = TRUE)
  expect_true(allowed$in_sample)
  expect_match(paste(cli::cli_fmt(print(allowed)), collapse = "\n"), "IN-SAMPLE")

  changed <- train
  changed$age[c(1, 2)] <- rev(changed$age[c(1, 2)])
  expect_false(identical(digest_data(changed), digest_data(train)))
  sc <- predict(fit, changed, uncertainty = "conditional", allow_extrapolation = TRUE)
  expect_false(any(sc$.in_sample))
})

test_that("ref_assess recovers nominal coverage and a positive log-score gain", {
  train <- ref_simulate(300, seed = 13)
  test <- ref_simulate(400, seed = 14)
  fit <- ref_fit(simple_spec(), data = train, outcomes = "y")
  a <- ref_assess(fit, newdata = test)
  expect_s3_class(a, "ref_assessment")
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
  scores <- scores[is.finite(scores$log_density), ] # out-of-support rows are masked
  -mean(scores$log_density) +
    mean(stats::dnorm(scores$observed, mean(y_ref), sd0, log = TRUE))
}

test_that("standardized_log_score is the negative MSLL against the reference sample", {
  ref <- ref_simulate(300, seed = 21)
  new <- ref_simulate(120, seed = 22)
  fit <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 5) + sex), ref, "y")

  for (rows in list(seq_len(120), seq_len(40))) {
    a <- ref_assess(fit, new[rows, ])
    sc <- predict(fit, new[rows, ], type = "scores", uncertainty = "total")
    expect_equal(
      -a$overall$standardized_log_score,
      msll_longhand(sc, ref$y),
      tolerance = 1e-12
    )
  }
})

test_that("the reference baseline survives freezing to a bundle", {
  ref <- ref_simulate(200, seed = 23)
  new <- ref_simulate(80, seed = 24)
  fit <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 5) + sex), ref, "y")
  thawed <- ref_freeze(fit)
  expect_equal(ref_assess(thawed, new)$overall, ref_assess(fit, new)$overall,
    tolerance = 1e-12
  )
})

test_that("ref_assess() honours the uncertainty argument", {
  ref <- ref_simulate(200, seed = 25)
  new <- ref_simulate(80, seed = 26)
  fit <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 5) + sex), ref, "y")

  total <- ref_assess(fit, new, uncertainty = "total")
  conditional <- ref_assess(fit, new, uncertainty = "conditional")
  expect_equal(ref_assess(fit, new)$overall, total$overall)
  expect_equal(
    -total$overall$standardized_log_score,
    msll_longhand(predict(fit, new, type = "scores", uncertainty = "total"), ref$y),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(
    conditional$overall$mean_log_score,
    ref_assess(fit, new)$overall$mean_log_score
  )))
})

test_that("ref_crossfit rejects one fold or one cluster clearly and scores NA strata", {
  dat <- ref_simulate(90, seed = 3)
  spec <- ref_spec(ref_gaussian(), ~ s(age, k = 5))
  expect_error(ref_crossfit(spec, dat, "y", folds = 1), "at least 2")
  dat$cl <- "one"
  expect_error(ref_crossfit(spec, dat, "y", folds = 3, cluster = cl), "two clusters")
  dat$st <- dat$site
  dat$st[1:10] <- NA
  cf <- ref_crossfit(spec, dat, "y", folds = 3, strata = st)
  expect_equal(nrow(cf), 90L)
  expect_equal(sort(unique(cf$.row)), 1:90)
  expect_true(all(attr(cf, "folds") %in% 1:3))
  # each NA row has a fold, spread across folds like any stratum
  expect_gt(length(unique(attr(cf, "folds")[1:10])), 1L)
})

test_that("ref_crossfit can score the held-out fold with total uncertainty", {
  dat <- ref_simulate(200, seed = 91)
  set.seed(92)
  cond <- ref_crossfit(simple_spec(), data = dat, outcomes = "y", folds = 4)
  set.seed(92)
  tot <- ref_crossfit(simple_spec(), data = dat, outcomes = "y", folds = 4,
                       uncertainty = "total")
  expect_equal(cond$.row, tot$.row)
  # the total predictive is wider, so out-of-fold z shrink toward zero
  expect_lt(stats::var(tot$z, na.rm = TRUE), stats::var(cond$z, na.rm = TRUE))
  expect_gt(mean(tot$log_density, na.rm = TRUE), mean(cond$log_density, na.rm = TRUE))
  expect_error(ref_crossfit(simple_spec(), data = dat, outcomes = "y",
                             uncertainty = "epistemic"), "arg")
})
