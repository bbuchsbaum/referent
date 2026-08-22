perf_data <- function() {
  ref <- norm_simulate(300, sites = 4, scale = "age", seed = 11)
  new <- norm_simulate(80, sites = 4, scale = "age", seed = 12)
  list(ref = ref, new = new)
}

perf_specs <- function() {
  list(
    constant = norm_spec(norm_gaussian(), ~ s(age, k = 6) + sex + s(site, bs = "re")),
    gaulss = norm_spec(norm_gaussian(), ~ s(age, k = 6) + sex + s(site, bs = "re"),
                       scale = ~ s(age, k = 4)),
    shash = norm_spec(norm_shash(), ~ s(age, k = 6) + sex, scale = ~ s(age, k = 4))
  )
}

test_that("frozen bundle is small and reproduces the fit after a round trip", {
  d <- perf_data()
  for (nm in names(perf_specs())) {
    fit <- norm_fit(perf_specs()[[nm]], d$ref, c("y", "marker_01"))
    fit <- norm_adapt(fit, d$new, by = site, parameters = c("location", "scale"))
    fit <- norm_calibrate(fit, d$new, by = site)
    bundle <- norm_reference(fit)
    path <- withr::local_tempfile(fileext = ".rds")
    saveRDS(bundle, path)
    thawed <- readRDS(path)
    size <- function(x) length(serialize(x, NULL))
    expect_lt(size(bundle), 0.8 * size(fit), label = paste(nm, "bundle size"))
    expect_null(thawed$models$y$model$model)
    expect_null(thawed$models$y$model$family)
    expect_null(thawed$models$y$model$fitted.values)
    for (u in c("conditional", "total")) {
      expect_equal(
        predict(thawed, d$new, uncertainty = u),
        predict(fit, d$new, uncertainty = u),
        tolerance = 1e-10, label = paste(nm, u)
      )
    }
    expect_equal(norm_assess(thawed, d$new)$overall, norm_assess(fit, d$new)$overall,
                 tolerance = 1e-10)
    expect_equal(norm_support(thawed, d$new), norm_support(fit, d$new))
    expect_equal(tidy(thawed)$n, tidy(fit)$n)
  }
})

test_that("one link prediction per outcome, one lp matrix only for draws", {
  d <- perf_data()
  fit <- norm_fit(perf_specs()$shash, d$ref, c("y", "marker_01"))
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

test_that("augment matches the long score table", {
  d <- perf_data()
  fit <- norm_calibrate(norm_fit(perf_specs()$gaulss, d$ref, c("y", "marker_01")),
                        d$new, by = site)
  long <- predict(fit, d$new, uncertainty = "conditional")
  wide <- augment(fit, d$new, uncertainty = "conditional")
  for (nm in c("y", "marker_01")) {
    rows <- long[long$.outcome == nm, ]
    expect_equal(wide[[paste0(".z_", nm)]], rows$z)
    expect_equal(wide[[paste0(".centile_", nm)]], rows$centile)
  }
  expect_equal(wide$.support, long$support[long$.outcome == "y"])
})

test_that("draw parameters from the lp matrix match per-draw evaluation", {
  d <- perf_data()
  fit <- norm_fit(perf_specs()$shash, d$ref, "y")
  m <- fit$models$y
  lp <- predict(m$model, d$new[1:5, ], type = "lpmatrix")
  beta <- coef(m$model)
  draws <- cbind(beta, beta * 1.01, beta * 0.99)
  joint <- params_from_eta(eta_from_lp(lp, draws), m$model, "shash", m)
  for (j in 1:3) {
    single <- params_from_eta(eta_from_lp(lp, draws[, j]), m$model, "shash", m)
    for (p in c("location", "scale", "skew", "tail")) {
      expect_equal(joint[[p]][, j], single[[p]])
    }
  }
})

test_that("shifted log density agrees with the distribution methods", {
  d <- perf_data()
  fit <- norm_fit(perf_specs()$shash, d$ref, "y")
  dist <- predict(fit, d$new, type = "distribution", uncertainty = "conditional")$y
  y <- d$new$y
  u <- dist_unpack(dist)
  expect_equal(
    shifted_log_density(u, y, 0.3, 0.1),
    dist_log_density(shift_params(dist, 0.3, 0.1), y)
  )
  total <- predict(fit, d$new[1:10, ], type = "distribution", uncertainty = "total")$y
  expect_equal(
    shifted_log_density(dist_unpack(total), y[1:10], -0.2, 0.05),
    dist_log_density(shift_params(total, -0.2, 0.05), y[1:10])
  )
})
