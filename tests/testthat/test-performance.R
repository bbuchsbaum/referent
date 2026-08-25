test_that("frozen bundle is small and reproduces the fit after a round trip", {
  d <- perf_data()
  for (nm in names(perf_specs())) {
    fit <- ref_fit(perf_specs()[[nm]], d$ref, c("y", "marker_01"))
    fit <- ref_adapt(fit, d$new, by = site, parameters = c("location", "scale"))
    fit <- ref_calibrate(fit, d$new, by = site, uncertainty = "conditional")
    bundle <- ref_freeze(fit)
    path <- withr::local_tempfile(fileext = ".rds")
    # spec formulas made in the test helper env serialise with a package-env
    # warning under load_all; that is the harness, not the bundle.
    suppressWarnings(saveRDS(bundle, path))
    thawed <- readRDS(path)
    size <- function(x) length(suppressWarnings(serialize(x, NULL)))
    expect_lt(size(bundle), 0.8 * size(fit), label = paste(nm, "bundle size"))
    expect_equal(nrow(thawed$models$y$model$model), 0L)
    expect_null(thawed$models$y$model$family)
    expect_null(thawed$models$y$model$fitted.values)
    expect_equal(
      predict(thawed, d$new, uncertainty = "conditional"),
      predict(fit, d$new, uncertainty = "conditional"),
      tolerance = 1e-10, label = paste(nm, "conditional")
    )
    expect_equal(
      predict(thawed, d$new, uncertainty = "total", type = "distribution"),
      predict(fit, d$new, uncertainty = "total", type = "distribution"),
      tolerance = 1e-10, label = paste(nm, "total distribution")
    )
    expect_equal(ref_assess(thawed, d$new, uncertainty = "conditional",
                            allow_calibration_reuse = TRUE)$overall,
                 ref_assess(fit, d$new, uncertainty = "conditional",
                            allow_calibration_reuse = TRUE)$overall,
                 tolerance = 1e-10)
    expect_equal(ref_support(thawed, d$new), ref_support(fit, d$new))
    expect_equal(tidy(thawed)$n, tidy(fit)$n)
  }
})

test_that("one link prediction per outcome, one lp matrix only for draws", {
  d <- perf_data()
  fit <- ref_fit(perf_specs()$shash, d$ref, c("y", "marker_01"))
  calls <- character()
  local_mocked_bindings(
    predict_gam_quiet = function(model, newdata, ...) {
      calls <<- c(calls, list(...)$type %||% "response")
      stats::predict(model, newdata = newdata, ...)
    }
  )
  predict(fit, d$new, uncertainty = "conditional")
  expect_equal(calls, c("link", "link"))
  calls <- character()
  predict(fit, d$new, uncertainty = "total")
  expect_equal(calls, rep(c("link", "lpmatrix"), 2))
})
