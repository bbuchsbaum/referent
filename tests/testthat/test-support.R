test_that("out-of-range ages are out, near-boundary ages are edge, and z is masked", {
  dat <- norm_simulate(150, seed = 9)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "y")
  r <- fit$support_ref$numeric$age
  extra <- dat[1:4, ]
  extra$age <- c(5, 50, 120, r$edge_hi + 0.5 * (r$max - r$edge_hi))
  st <- norm_support(fit, extra)
  expect_equal(st$support, c("out", "in", "out", "edge"))
  sc <- predict(fit, newdata = extra, uncertainty = "conditional")
  expect_equal(sc$support, st$support)
  expect_true(all(is.na(sc$z[c(1, 3)])))
  expect_true(all(is.finite(sc$z[c(2, 4)])))
  # centiles are still computed; only z is masked for out-of-support rows
  expect_true(all(is.finite(sc$centile)))
  sc2 <- predict(fit, newdata = extra, uncertainty = "conditional", allow_extrapolation = TRUE)
  expect_true(all(is.finite(sc2$z)))
  expect_equal(sc2$z[c(2, 4)], sc$z[c(2, 4)])
})

test_that("unseen factor levels are new_group only for covariates in the spec", {
  dat <- norm_simulate(80, seed = 10)
  fit <- norm_fit(simple_spec(), data = dat, outcomes = "y")
  extra <- dat[1, ]
  extra$sex <- factor("X", levels = c(levels(dat$sex), "X"))
  st <- norm_support(fit, extra)
  expect_equal(st$support, "new_group")
  extra$site <- factor("Z", levels = c(levels(dat$site), "Z"))
  extra$sex <- dat$sex[[1]]
  expect_equal(norm_support(fit, extra)$support, "in")
})

test_that("joint support uses the Mahalanobis radius over several numeric covariates", {
  dat <- norm_simulate(300, seed = 11)
  dat$bmi <- 22 + 0.1 * (dat$age - 50) + stats::rnorm(nrow(dat), sd = 1)
  spec <- norm_spec(norm_gaussian(), location = ~ s(age, k = 5) + bmi + sex)
  fit <- norm_fit(spec, data = dat, outcomes = "y")
  expect_false(is.null(fit$support_ref$cov))
  new <- dat[1:3, ]
  # marginally in range on both covariates but jointly implausible
  new$age[2] <- 25
  new$bmi[2] <- 30
  st <- norm_support(fit, new)
  expect_true(all(is.finite(st$d2)))
  expect_gt(st$d2[[2]], stats::qchisq(0.99, df = 2))
  expect_true(st$support[[2]] %in% c("edge", "out"))
  expect_equal(st$support[[1]], "in")
})

test_that("new subject ids are not new_group support", {
  dat <- norm_simulate(50, seed = 75)
  dat$participant_id <- paste0("S", seq_len(nrow(dat)))
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y",
    id = participant_id
  )
  tgt <- dat[1:8, ]
  tgt$participant_id <- paste0("T", seq_len(nrow(tgt)))
  sc <- predict(fit, newdata = tgt, uncertainty = "conditional")
  expect_false(any(sc$support == "new_group"))
})
