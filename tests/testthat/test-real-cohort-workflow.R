nhanes_expect_bound_receipt <- function(evidence) {
  receipt <- jsonlite::read_json(
    file.path(evidence, "receipt.json"), simplifyVector = FALSE
  )
  for (name in names(receipt$output_files)) {
    path <- file.path(evidence, paste0(name, ".csv"))
    table <- utils::read.csv(path, check.names = FALSE)
    record <- receipt$output_files[[name]]
    expect_identical(digest::digest(file = path, algo = "sha256"), record$sha256)
    expect_identical(nrow(table), as.integer(record$rows))
    expect_identical(names(table), unlist(record$columns, use.names = FALSE))
  }
  gate <- utils::read.csv(file.path(evidence, "gate.csv"))
  expect_identical(isTRUE(receipt$calibration_pass), isTRUE(gate$calibration_pass))
  expect_identical(isTRUE(receipt$pass), isTRUE(gate$pass))
  receipt
}

test_that("real-cohort lane separates development from untouched confirmation", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  runner_path <- file.path(root, "tools", "validation",
                           "run_nhanes_external_validation.R")
  workflow_path <- file.path(root, ".github", "workflows",
                             "real-cohort-validation.yml")
  skip_if_not(all(file.exists(runner_path, workflow_path)),
              "repository-only validation assets are excluded from the package tarball")
  runner <- paste(readLines(runner_path, warn = FALSE), collapse = "\n")
  workflow <- paste(readLines(workflow_path, warn = FALSE), collapse = "\n")
  expect_match(runner, "NHANES 2015-2016", fixed = TRUE)
  expect_match(runner, 'cohort = "2017-2018"', fixed = TRUE)
  expect_match(runner, 'cohort = "2013-2014"', fixed = TRUE)
  expect_match(runner, "post_hoc_model_development", fixed = TRUE)
  expect_match(runner, "untouched_confirmation", fixed = TRUE)
  expect_match(runner, "b1301a4e0426eb5b5a83935059053131331670d9", fixed = TRUE)
  expect_match(runner, "acceptable_calibration", fixed = TRUE)
  expect_match(
    runner,
    "scale = ~ s(age, k = 5) + sex + race_ethnicity",
    fixed = TRUE
  )
  expect_match(runner, "uncertainty = \"total\"", fixed = TRUE)
  expect_match(workflow, "if: always()", fixed = TRUE)
  expect_match(workflow, "confirmation", fixed = TRUE)
  expect_match(workflow, "nhanes-external-validation-evidence", fixed = TRUE)
})

test_that("retained 2017-2018 result is labeled post-hoc development evidence", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  evidence <- file.path(root, "docs", "evidence", "nhanes")
  skip_if_not(file.exists(file.path(evidence, "receipt.json")), "receipt is not bundled")
  receipt <- nhanes_expect_bound_receipt(evidence)
  expect_identical(receipt$evidence_role, "post_hoc_model_development")
  expect_identical(receipt$evaluation_cohort, "NHANES 2017-2018")
  expect_identical(receipt$design_registration$status, "post-hoc model development")
})

test_that("untouched NHANES confirmation is bound to the prior Mote contract", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  development <- file.path(root, "docs", "evidence", "nhanes")
  evidence <- file.path(development, "confirmation-2013-2014")
  skip_if_not(file.exists(file.path(evidence, "receipt.json")), "confirmation is not bundled")
  receipt <- nhanes_expect_bound_receipt(evidence)
  development_receipt <- jsonlite::read_json(
    file.path(development, "receipt.json"), simplifyVector = FALSE
  )
  expect_identical(receipt$evidence_role, "untouched_confirmation")
  expect_identical(receipt$evaluation_cohort, "NHANES 2013-2014")
  expect_identical(receipt$design_registration$tracker, "mote")
  expect_identical(receipt$design_registration$issue, "referent-32n")
  expect_identical(
    receipt$design_registration$commit,
    "b1301a4e0426eb5b5a83935059053131331670d9"
  )
  expect_false(identical(
    receipt$input_files$demo_test$sha256,
    development_receipt$input_files$demo_test$sha256
  ))
})
