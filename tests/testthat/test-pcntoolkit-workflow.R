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
  expect_match(workflow, "generate_release_superiority.py", fixed = TRUE)
  expect_match(workflow, "--replicates 20", fixed = TRUE)
  expect_match(workflow, "run_release_superiority.R", fixed = TRUE)
  expect_match(workflow, "github.event_name == 'release'", fixed = TRUE)
  expect_match(
    workflow,
    "Release SHASH fitting and replicated superiority evidence",
    fixed = TRUE
  )
})
