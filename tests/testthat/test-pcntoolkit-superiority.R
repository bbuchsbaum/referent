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
  referent <- pcn_referent_predictions(inputs, "skew_heavy", n_draw = 2000L)
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
