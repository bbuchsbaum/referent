test_that("a frozen reference predicts identically after a saveRDS/readRDS round trip", {
  dat <- ref_simulate(200, seed = 100)
  fit <- ref_fit(simple_spec(), data = dat[1:140, ], outcomes = c("y", "marker_01"))
  fit <- ref_calibrate(fit, data = dat[141:170, ], uncertainty = "conditional")
  ref <- ref_freeze(fit, criteria = "healthy volunteers", units = "mm")
  expect_identical(ref$bundle_schema_version, "1.0.0")
  tmp <- withr::local_tempfile(fileext = ".rds")
  suppressWarnings(saveRDS(ref, tmp))
  ref2 <- readRDS(tmp)
  expect_s3_class(ref2, "ref_freeze")
  expect_s3_class(ref2, "ref_fit")
  new <- dat[171:200, ]
  a <- predict(fit, newdata = new, uncertainty = "conditional")
  b <- predict(ref2, newdata = new, uncertainty = "conditional")
  expect_equal(b$z, a$z, tolerance = 1e-12)
  expect_equal(b$centile, a$centile, tolerance = 1e-12)
  expect_true(all(b$calibrated))
  d <- predict(ref2, newdata = new, type = "distribution", uncertainty = "total")
  expect_named(d, c(".id", "y", "marker_01"))
  expect_equal(length(d$y), nrow(new))
  expect_equal(ref_support(ref2, new)$support, a$support[a$.outcome == "y"])
  expect_equal(ref2$criteria, "healthy volunteers")
  expect_equal(ref2$units, "mm")
  expect_equal(fit_statuses(ref2), c(y = "ok", marker_01 = "ok"))
  expect_equal(tidy(ref2)$outcome, c("y", "marker_01"))
})

test_that("versioned bundle IO migrates legacy bundles and rejects future schemas", {
  dat <- ref_simulate(120, seed = 104)
  fit <- ref_fit(simple_spec(), dat, "y")
  bundle <- ref_freeze(fit)
  path <- withr::local_tempfile(fileext = ".rds")
  expect_identical(ref_write(bundle, path), path)
  expect_error(ref_read(path), "Refusing to deserialize")
  loaded <- ref_read(path, trusted = TRUE)
  expect_identical(loaded$bundle_schema_version, "1.0.0")
  expect_equal(
    predict(loaded, dat[1:4, ], uncertainty = "conditional")$z,
    predict(bundle, dat[1:4, ], uncertainty = "conditional")$z
  )

  legacy <- bundle
  legacy$bundle_type <- NULL
  legacy$bundle_schema_version <- NULL
  legacy_path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(legacy, legacy_path)
  expect_error(predict(legacy, dat[1:2, ]), "unversioned")
  expect_warning(migrated <- ref_read(legacy_path, trusted = TRUE), "Migrating")
  expect_identical(migrated$migration$from, "unversioned")

  future <- bundle
  future$bundle_schema_version <- "2.0.0"
  expect_error(predict(future, dat[1:2, ]), "Unsupported.*2.0.0")
})

test_that("ref_read rejects an untrusted RDS before deserialization", {
  sentinel <- withr::local_tempfile()
  payload <- new.env(parent = emptyenv())
  makeActiveBinding(
    "bundle_schema_version",
    local({
      marker <- sentinel
      function(value) {
        file.create(marker)
        "1.0.0"
      }
    }),
    payload
  )
  class(payload) <- c("ref_freeze", "ref_fit")
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(payload, path)

  expect_error(ref_read(path), "Refusing to deserialize")
  expect_false(file.exists(sentinel))
})

test_that("a frozen reference does not serialise the frame its spec was built in", {
  build <- function() {
    ballast <- stats::rnorm(2e6)
    dat <- ref_simulate(120, seed = 101)
    spec <- ref_spec(ref_gaussian(), location = ~ s(age, k = 5) + sex)
    fit <- ref_fit(spec, data = dat, outcomes = "y")
    list(fit = fit, bundle = ref_freeze(fit), dat = dat, ballast = ballast)
  }
  built <- build()
  tmp <- withr::local_tempfile(fileext = ".rds")
  saveRDS(built$bundle, tmp)
  expect_lt(file.size(tmp), 2e6)
  thawed <- readRDS(tmp)
  expect_equal(
    predict(thawed, built$dat[1:10, ], uncertainty = "conditional")$z,
    predict(built$fit, built$dat[1:10, ], uncertainty = "conditional")$z
  )
  expect_s3_class(ref_fit(thawed$spec, data = built$dat, outcomes = "y"), "ref_fit")
})

test_that("the bam path is used for constant-scale Gaussian fits above bam_min_n", {
  dat <- ref_simulate(300, seed = 101)
  spec_bam <- ref_spec(ref_gaussian(), location = ~ s(age, k = 6) + sex, scale = ~1,
                        bam_min_n = 200)
  spec_gam <- ref_spec(ref_gaussian(), location = ~ s(age, k = 6) + sex, scale = ~1,
                        use_bam = FALSE)
  fit_bam <- ref_fit(spec_bam, data = dat, outcomes = "y")
  fit_gam <- ref_fit(spec_gam, data = dat, outcomes = "y")
  expect_s3_class(fit_bam$models$y$model, "bam")
  expect_false(inherits(fit_gam$models$y$model, "bam"))
  expect_equal(fit_bam$models$y$status, "ok")
  new <- dat[1:20, ]
  a <- predict(fit_bam, newdata = new, uncertainty = "conditional")
  b <- predict(fit_gam, newdata = new, uncertainty = "conditional")
  expect_equal(a$z, b$z, tolerance = 0.02)
  # total uncertainty is analytic for both
  da <- predict(fit_bam, newdata = new, type = "distribution", uncertainty = "total")$y
  expect_s3_class(vctrs::vec_data(da)[[1]], "dist_normal")
  expect_true(all(variance(da) > variance(
    predict(fit_bam, newdata = new, type = "distribution", uncertainty = "conditional")$y
  )))
  # bam is not used for a gaulss spec
  spec_ls <- ref_spec(ref_gaussian(), location = ~ s(age, k = 6) + sex,
                       scale = ~ s(age, k = 4), bam_min_n = 200)
  fit_ls <- ref_fit(spec_ls, data = dat, outcomes = "y")
  expect_false(inherits(fit_ls$models$y$model, "bam"))
  # the bam fit supports chart derivatives with standard errors
  dv <- ref_derivative(fit_bam, data.frame(age = c(40, 60), sex = factor("F", levels = c("F", "M"))),
                        with_respect_to = age, centiles = 0.5)
  expect_true(all(is.finite(dv$chart_velocity_se) & dv$chart_velocity_se > 0))
})

test_that("ref_derivative accepts a bare name or a string and rejects unknown columns", {
  dat <- ref_simulate(150, seed = 102)
  fit <- ref_fit(ref_spec(ref_gaussian(), location = ~ age + sex, scale = ~1),
                  data = dat, outcomes = "y")
  grid <- data.frame(age = c(30, 60), sex = factor("M", levels = c("F", "M")))
  a <- ref_derivative(fit, grid, with_respect_to = age, centiles = 0.5)
  b <- ref_derivative(fit, grid, with_respect_to = "age", centiles = 0.5)
  expect_equal(a, b)
  expect_error(ref_derivative(fit, grid, with_respect_to = years), "must contain")
})

test_that("a frozen reference keeps offset terms and factor levels", {
  dat <- ref_simulate(200, seed = 5)
  dat$log_icv <- stats::rnorm(200, 0, 0.3)
  dat$y <- dat$y + dat$log_icv
  fit <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 5) + offset(log_icv)), dat, "y")
  ref <- ref_freeze(fit)
  a <- predict(fit, dat[1:8, ], uncertainty = "conditional")
  b <- predict(ref, dat[1:8, ], uncertainty = "conditional")
  expect_equal(a$z, b$z, tolerance = 1e-10)
  expect_gt(stats::sd(a$z - predict(fit, transform(dat[1:8, ], log_icv = 0),
                                    uncertainty = "conditional")$z), 0)

  fit2 <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 5) + s(site, bs = "re")), dat, "y")
  ref2 <- ref_freeze(fit2)
  new <- dat[1:8, ]
  new$site <- as.character(new$site)
  new$site[1] <- "ZZ"
  a <- predict(fit2, new, uncertainty = "conditional", allow_extrapolation = TRUE)
  b <- predict(ref2, new, uncertainty = "conditional", allow_extrapolation = TRUE)
  expect_equal(sum(is.na(a$z)), 0L)
  expect_equal(a$z, b$z, tolerance = 1e-10)
  for (unc in c("conditional", "total")) {
    a <- predict(fit2, new, uncertainty = unc, type = "distribution")$y
    b <- predict(ref2, new, uncertainty = unc, type = "distribution")$y
    expect_equal(mean(a), mean(b), tolerance = 1e-10)
  }
})
