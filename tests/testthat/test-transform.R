# A response transform is a monotone reparameterisation of the outcome: the
# fitted model is the same object, the PIT is unchanged, and only the
# density (through the Jacobian) and the quantile scale move. These tests
# pin all four halves of that statement.

warp_data <- function(n = 300, seed = 5, sd = 0.35) {
  set.seed(seed)
  age <- stats::runif(n, 20, 80)
  sex <- factor(sample(c("F", "M"), n, replace = TRUE))
  y <- exp(stats::rnorm(n, 7 + 0.01 * age + 0.2 * (sex == "M"), sd))
  data.frame(age, sex, y)
}

warp_spec <- function(transform) {
  ref_spec(ref_shash(), location = ~ s(age, k = 5) + sex,
           scale = ~ s(age, k = 4), transform = transform)
}

test_that("a transform leaves the PIT alone and moves the density by the Jacobian", {
  dat <- warp_data()
  logged <- dat
  logged$y <- log(dat$y)
  warped <- ref_fit(warp_spec("log"), data = dat, outcomes = "y")
  plain <- ref_fit(warp_spec("identity"), data = logged, outcomes = "y")

  new <- dat[1:40, ]
  new_log <- logged[1:40, ]
  a <- predict(warped, newdata = new, uncertainty = "conditional",
               allow_extrapolation = TRUE)
  b <- predict(plain, newdata = new_log, uncertainty = "conditional",
               allow_extrapolation = TRUE)
  # same model, so the probability integral transform is identical
  expect_equal(a$centile, b$centile)
  expect_equal(a$z, b$z)
  # the density picks up log |dh/dy| = -log y
  expect_equal(a$log_density, b$log_density - log(new$y))
  # the median is back-transformed, not the exponential of nothing
  expect_equal(a$median, exp(b$median))
  expect_equal(a$residual, new$y - exp(b$median))
})

test_that("the warped predictive is a proper density with a consistent CDF", {
  dat <- warp_data()
  for (tf in list("log", "sqrt", 0.25)) {
    fit <- ref_fit(warp_spec(tf), data = dat, outcomes = "y")
    row <- dat[3, ]
    dens <- function(yy) {
      r <- row
      r$y <- yy
      exp(predict(fit, newdata = r, uncertainty = "conditional",
                  allow_extrapolation = TRUE)$log_density)
    }
    # (0, Inf) is the whole support: a power transform loses no mass to the
    # part of the fitted distribution that lies below h(0)
    total <- stats::integrate(Vectorize(dens), 0, Inf, subdivisions = 2000L)$value
    expect_equal(total, 1, tolerance = 1e-5)

    d <- predict(fit, newdata = row, type = "distribution",
                 uncertainty = "conditional")$y
    p <- c(0.001, 0.05, 0.5, 0.95, 0.999)
    q <- as.numeric(quantile(d, p)[[1L]])
    expect_true(all(diff(q) > 0))
    expect_true(all(q > 0))
    back <- vapply(q, function(qi) as.numeric(unlist(cdf(d, qi))), numeric(1))
    expect_equal(back, p, tolerance = 1e-8)
    # the density is the derivative of the CDF on the response scale
    h <- 1e-4 * q[[3L]]
    numeric_dens <- (unlist(cdf(d, q[[3L]] + h)) - unlist(cdf(d, q[[3L]] - h))) / (2 * h)
    expect_equal(unname(numeric_dens), unname(unlist(density(d, q[[3L]]))),
                 tolerance = 1e-5)
    # y outside the domain is outside the support
    off <- row
    off$y <- -1
    expect_equal(exp(predict(fit, newdata = off, uncertainty = "conditional",
                             allow_extrapolation = TRUE)$log_density), 0)
  }
})

test_that("a transform makes log densities and CRPS comparable across specs", {
  dat <- warp_data()
  ladder <- list(identity = warp_spec("identity"), log = warp_spec("log"))
  fits <- lapply(ladder, ref_fit, data = dat, outcomes = "y")
  new <- warp_data(120, seed = 6)
  sc <- lapply(fits, function(f) {
    predict(f, newdata = new, uncertainty = "conditional", allow_extrapolation = TRUE)
  })
  # both are densities of the same y in the same units, so the log-normal
  # generator must favour the log transform
  expect_gt(mean(sc$log$log_density), mean(sc$identity$log_density))
  # CRPS is on the response scale, so it is not the CRPS of log y
  crps <- vapply(fits, function(f) {
    d <- predict_dists(f, new, uncertainty = "conditional")$y
    mean(crps_from_dist(d, new$y))
  }, numeric(1))
  expect_true(all(crps > 100))
  expect_lt(abs(crps[["log"]] - crps[["identity"]]) / crps[["identity"]], 0.3)
})

test_that("ref_select can choose the transform, and out-of-domain rows drop out", {
  dat <- warp_data(500, seed = 7, sd = 0.6)
  set.seed(31)
  sel <- ref_select(
    list(identity = warp_spec("identity"), log = warp_spec("log")),
    data = dat, outcomes = "y", folds = 4
  )
  expect_equal(sel$selected_name, "log")
  expect_equal(spec_transform(sel$selected)$name, "log")

  bad <- dat
  bad$y[1:5] <- -bad$y[1:5]
  expect_message(
    fit <- ref_fit(warp_spec("log"), data = bad, outcomes = "y"),
    "outside the domain"
  )
  expect_equal(fit$models$y$n_obs, nrow(dat) - 5L)
  # the reference baseline stays on the response scale
  expect_equal(fit$reference_baseline$y$mean, mean(bad$y))
})

test_that("transform arguments are validated and the identity is a no-op", {
  expect_error(ref_spec(transform = "cuberoot"), "Box-Cox")
  expect_error(ref_spec(transform = -1), "Box-Cox")
  expect_error(ref_spec(transform = c(1, 2)), "Box-Cox")
  expect_equal(spec_transform(ref_spec())$name, "identity")
  expect_equal(spec_transform(ref_spec(transform = 0.5))$name, "sqrt")
  expect_equal(spec_transform(ref_spec(transform = 0.3))$name, "power 0.3")

  dat <- warp_data(150, seed = 8)
  a <- ref_fit(warp_spec("identity"), data = dat, outcomes = "y")
  b <- ref_fit(warp_spec(1), data = dat, outcomes = "y")
  expect_equal(
    predict(a, newdata = dat[1:5, ], uncertainty = "conditional"),
    predict(b, newdata = dat[1:5, ], uncertainty = "conditional")
  )
})
