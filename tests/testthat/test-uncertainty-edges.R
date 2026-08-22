test_that("empty newdata returns an empty distribution", {
  set.seed(40)
  dat <- norm_simulate(60, seed = 40)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  d <- predict(fit, newdata = dat[0, ], type = "distribution")
  expect_equal(nrow(d), 0L)
  expect_equal(length(d$y), 0L)
  sc <- predict(fit, newdata = dat[0, ], uncertainty = "conditional")
  expect_equal(nrow(sc), 0L)
})

test_that("total uncertainty works for one row and gaulss", {
  set.seed(41)
  dat <- norm_simulate(80, kind = "gaussian", scale = "age", seed = 41)
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
  u <- dist_unpack(d)
  expect_s3_class(u, "dist_shash_mc")
  expect_equal(dim(u$sigma), dim(u$mu))
  expect_gt(stats::sd(u$mu), 0)
})

test_that("unseen groups do not abort prediction", {
  set.seed(42)
  dat <- norm_simulate(50, seed = 42)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  extra <- dat[1, ]
  extra$sex <- factor("X", levels = c(levels(dat$sex), "X"))
  sc <- predict(fit, newdata = extra, uncertainty = "conditional")
  expect_equal(nrow(sc), 1L)
  expect_equal(sc$support, "new_group")
})
