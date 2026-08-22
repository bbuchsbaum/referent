# Chart-velocity plots (ref_derivative) and the predictive interval on the
# individual velocity plot (ref_transition). Ported from the gauntlet-2
# velocity lane.

vel_grid <- data.frame(age = seq(25, 75, by = 5), sex = factor("F", levels = c("F", "M")))

# A single smooth bump that a k = 4 basis can only just represent, so REML
# spends the whole basis on it. (A trajectory far beyond the basis is
# collapsed towards linear by the penalty instead; that case has low edf
# and is not what the saturation warning detects.)
wiggly_cross <- function(n, seed) {
  withr::with_seed(seed, {
    age <- stats::runif(n, 20, 80)
    sex <- factor(sample(c("F", "M"), n, TRUE), levels = c("F", "M"))
    mu <- 0.04 * age + 3 * sin((age - 20) / 12) + 0.3 * (sex == "M")
    data.frame(age = age, sex = sex, y = mu + stats::rnorm(n, sd = 0.4))
  })
}

test_that("autoplot(<ref_derivative>) draws a line and a delta-method ribbon", {
  d <- ref_simulate(600, seed = 99)
  fit <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 8) + sex, scale = ~ s(age, k = 5)), d, "y")
  dv <- ref_derivative(fit, vel_grid, with_respect_to = age, centiles = c(0.05, 0.5, 0.95))
  p <- autoplot(dv)
  expect_s3_class(p, "ggplot")
  expect_true(all(c("GeomRibbon", "GeomLine", "GeomHline") %in% layer_classes(p)))
  expect_builds(p)
  expect_match(p$labels$subtitle, "95% delta-method interval")
  expect_match(p$labels$subtitle, "does not cover smoothing bias")

  # the ribbon is exactly est +/- z * se at the requested level
  for (lev in c(0.5, 0.95, 0.99)) {
    q <- autoplot(dv, level = lev)
    rib <- ggplot2::ggplot_build(q)$data[[which(layer_classes(q) == "GeomRibbon")]]
    z <- stats::qnorm(1 - (1 - lev) / 2)
    ref <- dv[order(dv$.outcome, dv$centile, dv$time), ]
    got <- rib[order(rib$group, rib$x), ]
    expect_equal(got$ymin, ref$chart_velocity - z * ref$chart_velocity_se, tolerance = 1e-8)
    expect_equal(got$ymax, ref$chart_velocity + z * ref$chart_velocity_se, tolerance = 1e-8)
  }

  # bad level and an empty derivative are refused, not silently drawn
  expect_error(autoplot(dv, level = 0), "level")
  expect_error(autoplot(dv, level = 1.2), "level")
  empty <- dv
  empty$chart_velocity <- NA_real_
  expect_error(autoplot(empty), "No finite chart velocities")

  # a missing standard error is disclosed rather than dropped
  bad <- dv
  bad$chart_velocity_se[1:3] <- NA_real_
  q <- autoplot(bad)
  expect_match(q$labels$subtitle, "3 points without a standard error")
  expect_equal(nrow(ggplot2::ggplot_build(q)$data[[which(layer_classes(q) == "GeomRibbon")]]),
               nrow(dv) - 3L)

  # multiple outcomes are faceted
  dv2 <- ref_derivative(fit, vel_grid, with_respect_to = age, centiles = 0.5)
  dv2b <- dv2
  dv2b$.outcome <- "y2"
  both <- vctrs::vec_rbind(dv2, dv2b)
  class(both) <- class(dv2)
  expect_s3_class(autoplot(both)$facet, "FacetWrap")
})

test_that("the individual velocity plot carries its predictive interval", {
  ld <- ref_simulate(150, kind = "longitudinal", seed = 31)
  lf <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 6) + sex, scale = ~1), ld, "y")
  dyn <- ref_dynamics(lf, data = ld, id = participant_id, time = age, crossfit = 0)
  tr <- ref_transition(dyn, data = ld, id = participant_id, time = age)
  p <- autoplot(tr, type = "velocity")
  expect_true("GeomLinerange" %in% layer_classes(p))
  expect_match(p$labels$subtitle, "90% history-conditioned predictive interval")
  expect_equal(p$labels$x, "Elapsed time")
  b <- ggplot2::ggplot_build(p)$data[[which(layer_classes(p) == "GeomLinerange")]]
  ok <- is.finite(tr$velocity_lower) & is.finite(tr$velocity_upper) & is.finite(tr$observed_velocity)
  expect_equal(nrow(b), sum(ok))
  for (ty in c("innovation", "change")) {
    q <- autoplot(tr, type = ty)
    expect_false("GeomLinerange" %in% layer_classes(q))
    expect_equal(sum(layer_classes(q) == "GeomHline"), 2L)
    expect_null(q$labels$subtitle)
  }
})

test_that("ref_derivative warns when the time smooth has spent its basis", {
  d8 <- wiggly_cross(1500, seed = 77)
  f8 <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 4) + sex, scale = ~1), d8, "y")
  expect_warning(ref_derivative(f8, vel_grid, with_respect_to = age, centiles = 0.5),
                 "spent .* of its basis")
  f30 <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 30) + sex, scale = ~1), d8, "y")
  expect_silent(ref_derivative(f30, vel_grid, with_respect_to = age, centiles = 0.5))
  fp <- ref_fit(ref_spec(ref_gaussian(), ~ age + sex, scale = ~1), d8, "y")
  expect_silent(ref_derivative(fp, vel_grid, with_respect_to = age, centiles = 0.5))
})
