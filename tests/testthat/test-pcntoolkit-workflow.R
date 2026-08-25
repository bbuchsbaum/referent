test_that("weekly evidence is rendered from the regenerated fixture", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  workflow_path <- file.path(root, ".github", "workflows", "pcntoolkit-validation.yml")
  renderer_path <- file.path(root, "tools", "pcntoolkit", "render_evidence.R")
  skip_if_not(file.exists(workflow_path) && file.exists(renderer_path))
  workflow <- paste(readLines(workflow_path, warn = FALSE), collapse = "\n")
  expect_match(
    workflow,
    "--fixture-root /tmp/pcntoolkit-v1.3.0",
    fixed = TRUE
  )

  renderer <- new.env(parent = globalenv())
  withr::local_dir(root)
  sys.source(renderer_path, envir = renderer)
  missing_fixture <- file.path(withr::local_tempdir(), "not-the-checked-in-fixture")
  expect_error(
    renderer$render_pcntoolkit_evidence(withr::local_tempdir(), missing_fixture),
    "missing PCNtoolkit manifest"
  )
})

test_that("release workflow runs every declared evidence lane", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  workflow_path <- file.path(root, ".github", "workflows", "pcntoolkit-validation.yml")
  skip_if_not(file.exists(workflow_path))
  workflow <- paste(readLines(workflow_path, warn = FALSE), collapse = "\n")
  expect_match(workflow, "--nuts-sampler nutpie", fixed = TRUE)
  expect_match(workflow, "--tune 1000", fixed = TRUE)
  expect_match(workflow, "--target-accept 0.99", fixed = TRUE)
  expect_match(workflow, "validate_hbr_site_evidence.R", fixed = TRUE)
  expect_match(workflow, "generate_release_superiority.py", fixed = TRUE)
  expect_match(workflow, "--replicates 20", fixed = TRUE)
  expect_match(workflow, "--seed-start 20260924", fixed = TRUE)
  expect_match(workflow, "run_release_superiority.R", fixed = TRUE)
  expect_match(
    workflow,
    "--scenarios linear_gaussian,nonlinear_heteroskedastic,skew_heavy,unequal_site,covariate_shift",
    fixed = TRUE
  )
  expect_match(workflow, "github.event_name == 'release'", fixed = TRUE)
  expect_match(
    workflow,
    "Release replicated scenario-matrix evidence",
    fixed = TRUE
  )
})

test_that("retained release matrix has a valid comparator in every replicate", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  evidence <- file.path(root, "docs", "evidence", "pcntoolkit", "v1.3.0")
  generation_path <- file.path(evidence, "release_generation_receipt.json")
  matrix_path <- file.path(evidence, "release_matrix_receipt.json")
  summary_path <- file.path(evidence, "release_superiority_summary.csv")
  skip_if_not(
    all(file.exists(generation_path, matrix_path, summary_path)),
    "release receipts are not bundled"
  )
  generation <- jsonlite::read_json(generation_path, simplifyVector = TRUE)
  matrix <- jsonlite::read_json(matrix_path, simplifyVector = TRUE)
  summary <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
  skew <- summary[summary$scenario == "skew_heavy", , drop = FALSE]
  expect_identical(generation$schema_version, "1.2.0")
  expect_true(generation$comparator_valid)
  expect_equal(
    generation$optimizer_controls$skew_heavy$l_bfgs_b_epsilon, 0.01
  )
  skew_warnings <- generation$fit_warnings[
    generation$fit_warnings$scenario == "skew_heavy", , drop = FALSE
  ]
  expect_equal(nrow(skew_warnings), 20L)
  expect_true(all(skew_warnings$valid))
  expect_true(all(skew_warnings$critical_count == 0L))
  expect_true(matrix$all_comparator_fits_valid)
  expect_true(matrix$pass)
  expect_identical(skew$classification, "superior")
  expect_true(skew$pass)
})
