test_that("matched PCNtoolkit fitted models satisfy prediction-level margins", {
  root <- testthat::test_path("..", "fixtures", "pcntoolkit", "v1.3.0")
  result <- pcn_run_concordance(root, n_draw = 2000L)
  matched <- result$evidence_class == "matched_estimator"

  expect_true(all(pcn_concordance_pass(result)))
  expect_true(all(is.finite(as.matrix(result[, 4:ncol(result)]))))
  expect_setequal(
    result$scenario[matched],
    c("null_linear", "linear_gaussian", "multiple_encoded", "log_linear",
      "covariate_shift")
  )
  expect_setequal(
    result$scenario[!matched],
    c("nonlinear_heteroskedastic", "balanced_site", "skew_heavy", "unequal_site")
  )
})

test_that("estimator-divergent fitted cases remain diagnostics, not fake parity", {
  registry <- pcn_concordance_registry()
  divergent <- registry$evidence_class == "same_estimand"
  expect_true(all(nzchar(registry$reason[divergent])))
  expect_false(any(registry$evidence_class[divergent] == "matched_estimator"))
})
