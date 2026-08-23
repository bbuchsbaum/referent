test_that("benchmark registry locks all prespecified challenge families", {
  registry <- pcn_benchmark_registry()
  expect_setequal(
    registry$scenario,
    c("linear_gaussian", "nonlinear_mean", "heteroskedastic", "skew_heavy",
      "unequal_sites", "site_slopes", "covariate_shift", "unseen_site",
      "small_n", "large_n")
  )
  expect_true(all(registry$n_train > 0 & registry$n_validation > 0 & registry$n_test > 0))
  expect_false(any(grepl("test", registry$selection_rule, fixed = TRUE)))

  unseen <- pcn_simulate_benchmark("unseen_site", seed = 7)
  expect_false("site-5" %in% unseen$site[unseen$split == "train"])
  expect_true("site-5" %in% unseen$site[unseen$split == "test"])
  skew <- pcn_simulate_benchmark("skew_heavy", seed = 8)
  expect_true(all(skew$truth_epsilon == 0.9))
  expect_true(all(skew$truth_delta == 0.7))
  expect_equal(length(unique(skew$row_id)), nrow(skew))
})

test_that("the locked skew-heavy test supports a guarded superiority claim", {
  root <- testthat::test_path("..", "fixtures", "pcntoolkit", "v1.3.0")
  inputs <- utils::read.csv(file.path(root, "fitted_inputs.csv"), stringsAsFactors = FALSE)
  comparator <- utils::read.csv(
    file.path(root, "fitted_predictions.csv"), stringsAsFactors = FALSE
  )
  fitted <- pcn_referent_predictions(
    inputs, "skew_heavy", n_draw = 2000L, return_fit = TRUE
  )
  expect_identical(fitted$fit$models$y$status, "ok")
  referent <- fitted$predictions
  comparator <- comparator[comparator$scenario == "skew_heavy", , drop = FALSE]
  comparator <- comparator[match(referent$row_id, comparator$row_id), , drop = FALSE]
  result <- pcn_compare_predictions(
    referent, comparator, cluster = referent$row_id,
    scenario = "skew_heavy", B = 999L, seed = 91L
  )
  expect_identical(result$classification, "superior")
  expect_gt(result$ci_lower, 0)
  expect_true(result$calibration_pass)
  expect_lte(
    abs(result$referent_coverage90 - 0.9),
    abs(result$pcntoolkit_coverage90 - 0.9) + 0.01
  )
})

test_that("negative controls prevent false superiority and detect tradeoffs", {
  n <- 300L
  rows <- sprintf("row-%03d", seq_len(n))
  observed <- stats::qnorm((seq_len(n) - 0.5) / n)
  base <- data.frame(
    row_id = rows, observed = observed, median = 0,
    centile = stats::pnorm(observed), log_density = stats::dnorm(observed, log = TRUE),
    q05 = stats::qnorm(0.05), q95 = stats::qnorm(0.95)
  )
  identical_result <- pcn_compare_predictions(base, base, rows, B = 199L, seed = 1L)
  expect_identical(identical_result$classification, "equivalent")
  expect_equal(identical_result$log_score_difference, 0)

  expect_identical(
    pcn_classify_benchmark(0.08, 0.03, 0.12, calibration_pass = FALSE),
    "tradeoff"
  )
  expect_identical(
    pcn_classify_benchmark(-0.08, -0.12, -0.03, calibration_pass = TRUE),
    "inferior"
  )
  expect_identical(
    pcn_classify_benchmark(0, 0, 0, calibration_pass = TRUE, conformance_pass = TRUE),
    "conformant"
  )
})

test_that("controlled location perturbations worsen proper score monotonically", {
  y <- stats::qnorm((1:999 - 0.5) / 999)
  shifts <- c(0, 0.2, 0.5, 1)
  score <- vapply(shifts, function(shift) mean(stats::dnorm(y, shift, 1, log = TRUE)), numeric(1))
  expect_true(all(diff(score) < 0))
})

test_that("release superiority requires replicated convergence and calibration", {
  results <- data.frame(
    replicate = seq_len(8),
    log_score_difference = seq(0.12, 0.26, length.out = 8),
    referent_coverage90 = 0.90,
    pcntoolkit_coverage90 = 0.92,
    referent_mace = 0.01,
    pcntoolkit_mace = 0.05,
    referent_tail05 = 0.05,
    pcntoolkit_tail05 = 0.06,
    fit_status = "ok",
    fit_converged = TRUE
  )
  passed <- pcn_release_superiority_summary(results, B = 999L, seed = 14L)
  expect_identical(passed$classification, "superior")
  expect_true(passed$pass)
  expect_gt(passed$ci_lower, 0)
  expect_true(passed$calibration_noninferior)
  expect_true(passed$no_critical_replicate_regression)

  failed_fit <- results
  failed_fit$fit_status[[3L]] <- "nonconverged"
  fit_summary <- pcn_release_superiority_summary(failed_fit, B = 999L, seed = 14L)
  expect_identical(fit_summary$classification, "fit_failure")
  expect_false(fit_summary$pass)

  failed_calibration <- results
  failed_calibration$referent_coverage90[[5L]] <- 0.70
  failed_calibration$pcntoolkit_coverage90[[5L]] <- 0.90
  calibration_summary <- pcn_release_superiority_summary(
    failed_calibration, B = 999L, seed = 14L
  )
  expect_identical(calibration_summary$classification, "tradeoff")
  expect_false(calibration_summary$pass)
  expect_error(
    pcn_release_superiority_summary(results[1:4, ], B = 999L),
    "at least five unique replicates"
  )
})
