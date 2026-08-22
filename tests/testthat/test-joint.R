test_that("joint scores handle missing outcomes via a submatrix", {
  dat <- ref_simulate(100, seed = 26)
  fit <- ref_fit(simple_spec(), data = dat, outcomes = starts_with("marker_"))
  sc <- predict(fit, newdata = dat, uncertainty = "conditional", allow_extrapolation = TRUE)
  sc$z[sc$.outcome == "marker_03" & sc$.row <= 10] <- NA_real_
  jt <- ref_joint(sc, covariance = "shrinkage")
  expect_s3_class(jt, "ref_joint")
  ref <- jt$reference
  expect_true(all(ref$n_observed[1:10] == 2))
  expect_true(all(is.finite(ref$d2[1:10])))
  expect_true(all(ref$joint_centile > 0 & ref$joint_centile < 1))
  lam <- attr(jt$correlation, "lambda")
  expect_true(lam >= 0 && lam <= 1)
})

test_that("joint centiles are calibrated irrespective of how many outcomes are observed", {
  ref_dat <- ref_simulate(600, seed = 27)
  new_dat <- ref_simulate(3000, seed = 28)
  fit <- ref_fit(simple_spec(), data = ref_dat, outcomes = starts_with("marker_"))
  cf <- ref_crossfit(simple_spec(), data = ref_dat, outcomes = starts_with("marker_"), folds = 4)
  jt <- ref_joint(cf)
  sc <- predict(fit, newdata = new_dat, uncertainty = "conditional", allow_extrapolation = TRUE)
  # knock out one outcome in a third of the new subjects, two in another third
  drop1 <- sc$.outcome == "marker_03" & sc$.row %% 3 == 1
  drop2 <- sc$.outcome %in% c("marker_02", "marker_03") & sc$.row %% 3 == 2
  sc$z[drop1 | drop2] <- NA_real_
  pr <- predict(jt, sc)
  expect_equal(sort(unique(pr$n_observed)), 1:3)
  for (k in 1:3) {
    rate <- mean(pr$joint_centile[pr$n_observed == k] > 0.95)
    expect_lt(abs(rate - 0.05), 0.02)
  }
  expect_lt(abs(mean(pr$joint_centile > 0.95) - 0.05), 0.012)
  expect_lt(suppressWarnings(pit_ks(pr$joint_centile)), 0.05)
  # a grossly deviant subject lands in the far tail
  far <- sc[sc$.row == 1, ]
  far$z <- c(4, 4, 4)
  expect_gt(predict(jt, far)$joint_centile, 0.99)
})

test_that("ref_joint drops constant outcomes instead of producing a NaN correlation", {
  dat <- ref_simulate(100, seed = 3)
  dat$const <- 1
  fit <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 5)), dat,
                  c("y", "marker_01", "const"))
  aug <- augment(fit, dat, uncertainty = "conditional")
  expect_message(joint <- ref_joint(aug), "const")
  expect_equal(joint$outcomes, c("y", "marker_01"))
  expect_true(all(is.finite(joint$correlation)))
  expect_true(is.finite(attr(joint$correlation, "lambda")))
  pr <- predict(joint, aug)
  expect_true(all(is.finite(pr$joint_z)))
  expect_equal(pr$n_observed, rep(2L, nrow(aug)))
  expect_error(ref_joint(aug[, ".z_const", drop = FALSE]), "finite")
})
