test_that("autoplot returns a ggplot and norm_reference prints a card", {
  set.seed(30)
  dat <- norm_simulate(80, seed = 30)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  p <- ggplot2::autoplot(fit, type = "centiles")
  expect_s3_class(p, "ggplot")
  expect_s3_class(theme_referent(), "theme")
  card <- norm_reference(fit)
  expect_s3_class(card, "norm_reference")
  out <- paste(cli::cli_fmt(print(card)), collapse = "\n")
  expect_match(out, "model card")
})

test_that("centile fortify has paired ribbons and a median line", {
  dat <- norm_simulate(90, seed = 31)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  built <- fortify_centiles(fit, by = "sex")
  expect_true(all(c("x", "y", "centile", ".label", ".group") %in% names(built$lines)))
  expect_true(all(c("ymin", "ymax", "band") %in% names(built$ribbons)))
  expect_true("Median" %in% built$lines$.label)
  expect_true(all(c("F", "M") %in% built$lines$.group))
  p <- autoplot(fit, type = "centiles", by = sex, newdata = dat)
  expect_s3_class(p, "ggplot")
  expect_true(any(vapply(p$layers, function(ly) inherits(ly$geom, "GeomRibbon"), logical(1))))
})

test_that("visit and coverage plots consume tidy builders", {
  dat <- norm_simulate(60, kind = "longitudinal", seed = 32)
  tab <- fortify_visits(dat, participant_id)
  expect_true(all(tab$n_visits >= 1))
  expect_equal(sum(tab$n_subjects), length(unique(dat$participant_id)))
  cov <- fortify_coverage(dat, age, by = site)
  expect_true(all(c("time", ".group") %in% names(cov)))
  p1 <- autoplot(structure(list(outcomes = "y"), class = "norm_fit"),
                 type = "visits", data = dat, id = participant_id)
  expect_s3_class(p1, "ggplot")
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y",
    id = participant_id
  )
  p2 <- autoplot(fit, type = "coverage", data = dat, time = age, by = site)
  expect_s3_class(p2, "ggplot")
  p3 <- autoplot(fit, type = "trajectories", data = dat, id = participant_id, time = age)
  expect_s3_class(p3, "ggplot")
})

test_that("forecast fan uses history and named centiles", {
  skip_on_cran()
  dat <- norm_simulate(90, kind = "longitudinal", seed = 33)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age)
  hist <- dat[dat$participant_id == dat$participant_id[[1]], ]
  fc <- norm_forecast(dyn, history = hist, times = max(hist$age) + c(1, 2, 3))
  expect_true(!is.null(fc$history))
  p <- autoplot(fc, type = "fan")
  expect_s3_class(p, "ggplot")
  p2 <- autoplot(fc, type = "thrive")
  expect_s3_class(p2, "ggplot")
  pk <- autoplot(dyn, type = "kernel")
  expect_s3_class(pk, "ggplot")
  expect_true(all(fortify_kernel(dyn)$correlation <= 1 + 1e-8))
})

test_that("assessment and score plots stay ggplot objects", {
  dat <- norm_simulate(80, seed = 34)
  val <- norm_simulate(40, seed = 35)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  a <- norm_assess(fit, newdata = val)
  expect_s3_class(autoplot(a, type = "calibration"), "ggplot")
  expect_s3_class(autoplot(a, type = "worm"), "ggplot")
  sc <- predict(fit, newdata = val[1:8, ], uncertainty = "conditional")
  expect_s3_class(autoplot(sc, type = "profile"), "ggplot")
  expect_s3_class(autoplot(sc, type = "heatmap"), "ggplot")
})
