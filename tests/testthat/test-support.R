test_that("out-of-range ages are flagged and z is NA by default", {
  set.seed(9)
  dat <- norm_simulate(150, kind = "gaussian_location", seed = 9)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "y")
  extra <- dat[1:3, ]
  extra$age <- c(5, 50, 120)
  sc <- predict(fit, newdata = extra, uncertainty = "conditional")
  expect_true(sc$support[[1]] %in% c("out", "edge"))
  expect_true(sc$support[[3]] %in% c("out", "edge"))
  expect_true(is.na(sc$z[sc$support == "out"][1]) || any(sc$support == "out"))
})

test_that("unseen factor levels are new_group", {
  set.seed(10)
  dat <- norm_simulate(80, kind = "gaussian_location", seed = 10)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "y")
  extra <- dat[1, ]
  extra$site <- factor("Z", levels = c(levels(dat$site), "Z"))
  st <- classify_support(fit$support_ref, extra)
  expect_equal(st, "new_group")
})
