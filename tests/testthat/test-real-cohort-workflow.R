test_that("real-cohort lane declares independent NHANES cohorts and durable evidence", {
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
  expect_match(runner, "NHANES 2017-2018", fixed = TRUE)
  expect_match(runner, "acceptable_calibration", fixed = TRUE)
  expect_match(
    runner,
    "scale = ~ s(age, k = 5) + sex + race_ethnicity",
    fixed = TRUE
  )
  expect_match(runner, "uncertainty = \"total\"", fixed = TRUE)
  expect_match(workflow, "if: always()", fixed = TRUE)
  expect_match(workflow, "nhanes-external-validation-evidence", fixed = TRUE)
})

test_that("retained NHANES receipt passes the full conditional contract", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  evidence <- file.path(root, "docs", "evidence", "nhanes")
  receipt_path <- file.path(evidence, "receipt.json")
  gate_path <- file.path(evidence, "gate.csv")
  skip_if_not(all(file.exists(receipt_path, gate_path)), "cohort receipts are not bundled")
  receipt <- jsonlite::read_json(receipt_path, simplifyVector = TRUE)
  gate <- utils::read.csv(gate_path, stringsAsFactors = FALSE)
  expect_true(receipt$calibration_pass)
  expect_true(receipt$pass)
  expect_true(gate$calibration_pass)
  expect_true(gate$pass)
})
