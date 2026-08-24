test_that("quality workflow covers supported operating systems and scale budgets", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  workflow_path <- file.path(root, ".github", "workflows", "quality.yml")
  benchmark_path <- file.path(root, "tools", "benchmarks", "run_scale_benchmarks.R")
  skip_if_not(all(file.exists(workflow_path, benchmark_path)),
              "repository-only quality assets are excluded from the package tarball")
  workflow <- paste(readLines(workflow_path,
                              warn = FALSE), collapse = "\n")
  benchmark <- paste(readLines(benchmark_path, warn = FALSE),
                     collapse = "\n")
  expect_match(workflow, "ubuntu-latest", fixed = TRUE)
  expect_match(workflow, "macos-latest", fixed = TRUE)
  expect_match(workflow, "windows-latest", fixed = TRUE)
  expect_match(workflow, "oldrel-1", fixed = TRUE)
  expect_match(workflow, "--as-cran", fixed = TRUE)
  expect_match(workflow, "run_scale_benchmarks.R", fixed = TRUE)
  expect_match(benchmark, "maximum_allocated_bytes", fixed = TRUE)
  expect_match(benchmark, "bundle_bytes", fixed = TRUE)
})
