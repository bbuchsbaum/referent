test_that("autoplot returns a ggplot and norm_card prints", {
  set.seed(30)
  dat <- norm_simulate(80, kind = "gaussian_location", seed = 30)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat,
    outcomes = "y"
  )
  p <- ggplot2::autoplot(fit, type = "centiles")
  expect_s3_class(p, "ggplot")
  card <- norm_card(fit)
  expect_s3_class(card, "norm_card")
  out <- paste(utils::capture.output(print(card)), collapse = "\n")
  expect_match(out, "model card")
})
