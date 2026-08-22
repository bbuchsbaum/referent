test_that("empty newdata returns an empty distribution", {
  set.seed(40)
  dat <- norm_simulate(60, kind = "gaussian_location", seed = 40)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  d <- predict(fit, newdata = dat[0, ], type = "distribution")$y
  expect_equal(length(d), 0L)
  sc <- predict(fit, newdata = dat[0, ], uncertainty = "conditional")
  expect_equal(nrow(sc), 0L)
})

test_that("total uncertainty works for one row and gaulss", {
  set.seed(41)
  dat <- norm_simulate(80, kind = "gaussian_scale", seed = 41)
  fit <- norm_fit(
    norm_spec(
      family = norm_gaussian(),
      location = ~ s(age, k = 5) + sex,
      scale = ~ s(age, k = 4)
    ),
    data = dat,
    outcomes = "y"
  )
  d <- predict(fit, newdata = dat[1, ], type = "distribution", uncertainty = "total")$y
  expect_equal(length(d), 1L)
  expect_false(is.null(attr(d, "location_draws")))
  expect_equal(ncol(attr(d, "scale_draws")), ncol(attr(d, "location_draws")))
  expect_gt(field_or(d, "epistemic_sd"), 0)
})

test_that("unseen groups do not abort prediction", {
  set.seed(42)
  dat <- norm_simulate(50, kind = "gaussian_location", seed = 42)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  extra <- dat[1, ]
  extra$site <- factor("ZZ", levels = c(levels(dat$site), "ZZ"))
  sc <- predict(fit, newdata = extra, uncertainty = "conditional")
  expect_equal(nrow(sc), 1L)
  expect_equal(sc$support, "new_group")
})
