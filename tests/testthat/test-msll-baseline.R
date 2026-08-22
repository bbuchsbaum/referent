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
