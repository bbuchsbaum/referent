fixture_root <- testthat::test_path("..", "fixtures", "pcntoolkit", "v1.3.0")

test_that("PCNtoolkit fixtures are pinned, complete, and unmodified", {
  manifest <- pcn_validate_fixture(fixture_root)
  expect_identical(manifest$schema_version, "1.1.0")
  expect_identical(manifest$dependencies$pcntoolkit, "1.3.0")
  expect_identical(manifest$generator, "tools/pcntoolkit/generate_fixtures.py")
  expect_setequal(
    names(manifest$scenario_classes),
    c("null_linear", "linear_gaussian", "multiple_encoded",
      "nonlinear_heteroskedastic", "log_linear", "balanced_site",
      "skew_heavy", "unequal_site", "covariate_shift")
  )
  expect_identical(manifest$optimizer_controls$skew_heavy$optimizer, "l-bfgs-b")
  expect_equal(manifest$optimizer_controls$skew_heavy$l_bfgs_b_epsilon, 0.01)
})

test_that("fixture validation fails closed on content drift", {
  scratch <- withr::local_tempdir()
  file.copy(list.files(fixture_root, full.names = TRUE), scratch)
  write("scientific drift", file.path(scratch, "semantic_gaussian.csv"), append = TRUE)
  expect_error(pcn_validate_fixture(scratch), "hash mismatch")
})

test_that("PCNtoolkit logp is normalised to a response-scale density", {
  predicted <- utils::read.csv(
    file.path(fixture_root, "fitted_predictions.csv"), stringsAsFactors = FALSE
  )
  gaussian <- predicted$scenario != "log_linear"
  analytic <- stats::dnorm(
    predicted$observed[gaussian], predicted$median[gaussian],
    predicted$predictive_sd[gaussian], log = TRUE
  )
  expect_equal(predicted$log_density[gaussian], analytic, tolerance = 1e-10)

  logged <- predicted[predicted$scenario == "log_linear", ]
  log_sd <- (log(logged$q95) - log(logged$q05)) / (2 * stats::qnorm(0.95))
  expect_equal(
    logged$log_density,
    stats::dlnorm(logged$observed, log(logged$median), log_sd, log = TRUE),
    tolerance = 1e-10
  )
})

test_that("installed validation helpers mirror their source-tree runners", {
  root <- testthat::test_path("..", "..")
  paired <- c(
    "compare_results.R", "run_concordance.R", "benchmark_design.R",
    "benchmark_helpers.R", "site_evidence.R", "benchmark_scenarios.csv",
    "run_hbr_site_evidence.py"
  )
  source_files <- file.path(root, "tools", "pcntoolkit", paired)
  skip_if_not(all(file.exists(source_files)), "source-tree tools are not installed")
  installed_files <- file.path(root, "inst", "pcntoolkit", paired)
  expect_true(all(file.exists(installed_files)))
  for (i in seq_along(paired)) {
    expect_identical(readBin(source_files[[i]], "raw", file.info(source_files[[i]])$size),
                     readBin(installed_files[[i]], "raw", file.info(installed_files[[i]])$size),
                     info = paired[[i]])
  }
})
