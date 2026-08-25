plot_fit <- function(seed = 30, n = 120) {
  dat <- ref_simulate(n, seed = seed)
  dat$participant_id <- paste0("P", seq_len(nrow(dat)))
  fit <- ref_fit(
    ref_spec(family = ref_gaussian(), location = ~ s(age, k = 5) + sex, scale = ~1),
    data = dat,
    outcomes = c("y", "marker_01"),
    id = participant_id
  )
  list(fit = fit, data = dat)
}


# The guide of the colour scale ("legend" or "none").
colour_guide <- function(p) {
  p$scales$get_scales("colour")$guide
}

test_that("centile chart: paired ribbons, median, facets, overlay, and stable appearance", {
  pf <- plot_fit(31)
  built <- fortify_centiles(pf$fit, by = "sex")
  expect_true(all(c("x", "y", "centile", ".label", ".group") %in% names(built$lines)))
  expect_equal(sort(unique(built$ribbons$band)), c("25th-75th", "5th-95th"))
  expect_true("Median" %in% built$lines$.label)
  expect_setequal(unique(built$lines$.group), c("F", "M"))
  # the 25-75 band is the darker (inner) band
  p <- autoplot(pf$fit, type = "centiles", by = sex, newdata = pf$data[1:30, ])
  b <- expect_builds(p)
  fills <- unique(b$data[[1]]$fill)
  alpha <- function(col) grDevices::col2rgb(col, alpha = TRUE)["alpha", ] / 255
  expect_equal(sort(alpha(fills)), c(0.12, 0.25), tolerance = 0.01)
  expect_true("GeomRibbon" %in% layer_classes(p))
  expect_true("GeomPoint" %in% layer_classes(p))
  expect_equal(length(unique(b$layout$layout$PANEL)), 2L)
  expect_error(fortify_centiles(pf$fit, centiles = c(0, 0.5)), "must lie in")
  expect_error(fortify_centiles(pf$fit, by = "site"), "factor covariate")
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("centiles-by-sex", p)
})

test_that("trajectories overlay subject paths on the centile chart", {
  dat <- ref_simulate(150, kind = "longitudinal", seed = 32)
  fit <- ref_fit(
    ref_spec(family = ref_gaussian(), location = ~ s(age, k = 5) + sex, scale = ~1),
    data = dat, outcomes = "y", id = participant_id
  )
  p <- autoplot(fit, type = "trajectories", data = dat[dat$participant_id <= 12, ],
                id = participant_id, time = age)
  b <- expect_builds(p)
  cls <- layer_classes(p)
  expect_true(all(c("GeomRibbon", "GeomLine", "GeomPath", "GeomPoint") %in% cls))
  path <- b$data[[which(cls == "GeomPath")]]
  expect_equal(length(unique(path$group)), 12L)
  expect_error(autoplot(fit, type = "trajectories"), "data")
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("trajectories", p)
})

test_that("support plot classifies newdata against the fit", {
  pf <- plot_fit(33)
  new <- pf$data[1:40, ]
  new$age[1:5] <- 120
  p <- autoplot(pf$fit, type = "support", newdata = new)
  b <- expect_builds(p)
  expect_true(all(c("in", "out") %in% b$plot$data$support))
  expect_false(is.null(plot_label(p, "fill")))
  # a single support level drops the legend
  r <- pf$fit$support_ref$numeric$age
  inside <- pf$data[pf$data$age > r$edge_lo & pf$data$age < r$edge_hi, ][1:20, ]
  p1 <- autoplot(pf$fit, type = "support", newdata = inside)
  expect_builds(p1)
  expect_identical(p1$theme$legend.position, "none")
  expect_error(autoplot(pf$fit, type = "support"), "newdata")
})

test_that("adaptation plot shows offsets with intervals per group", {
  pf <- plot_fit(35, n = 200)
  local <- ref_simulate(120, site_shift = c(0, 1, 0, 0), seed = 36)
  ad <- ref_adapt(pf$fit, data = local, by = site)
  p <- autoplot(ad, type = "adaptation")
  b <- expect_builds(p)
  cls <- layer_classes(p)
  expect_true(all(c("GeomLinerange", "GeomPoint") %in% cls))
  pts <- b$data[[which(cls == "GeomPoint")]]
  expect_equal(nrow(pts), 4L * 2L)
  expect_equal(plot_label(p, "x"), "site")
  expect_error(autoplot(pf$fit, type = "adaptation"), "ref_adapt")
  expect_identical(colour_guide(p), "legend")
  # pooled adaptation of a single outcome has no colour legend
  fit1 <- ref_fit(pf$fit$spec, data = pf$data, outcomes = "y")
  ad1 <- ref_adapt(fit1, data = local)
  p1 <- autoplot(ad1, type = "adaptation")
  expect_builds(p1)
  expect_identical(colour_guide(p1), "none")
})

test_that("assessment plots: calibration, qq, worm, conditional", {
  pf <- plot_fit(37, n = 150)
  val <- ref_simulate(80, seed = 38)
  a <- ref_assess(pf$fit, newdata = val)
  p <- autoplot(a, type = "calibration")
  b <- expect_builds(p)
  cls <- layer_classes(p)
  expect_equal(nrow(b$data[[which(cls == "GeomPoint")]]), 2L * 5L)
  expect_identical(colour_guide(p), "legend")
  # qq and worm share the envelope: the worm is the qq chart minus the expected line
  pq <- autoplot(a, type = "qq")
  pw <- autoplot(a, type = "worm")
  expect_true(grepl("\n", plot_label(pq, "subtitle"), fixed = TRUE))
  expect_true(grepl("\n", plot_label(pw, "subtitle"), fixed = TRUE))
  bq <- expect_builds(pq)
  bw <- expect_builds(pw)
  clq <- layer_classes(pq)
  expect_equal(sum(clq == "GeomRibbon"), 2L)
  expect_true(all(grepl("^simultaneous:", bq$data[[which(clq == "GeomText")]]$label)))
  expect_equal(nrow(bw$data[[which(layer_classes(pw) == "GeomPoint")]]), 2L * 80L)
  env_q <- pq$data
  env_w <- pw$data
  expect_equal(env_w$y, env_q$y - env_q$expected)
  expect_equal(env_w$sim_lo, env_q$sim_lo - env_q$expected)
  expect_equal(env_w$point_hi, env_q$point_hi - env_q$expected)
  expect_identical(env_w$outside, env_q$outside)
  # the pointwise band is the exact Beta order-statistic band
  n <- sum(env_q$.outcome == "y")
  one <- env_q[env_q$.outcome == "y", ]
  i <- seq_len(n)
  expect_equal(one$point_lo, stats::qnorm(stats::qbeta(0.025, i, n - i + 1)))
  expect_equal(one$point_hi, stats::qnorm(stats::qbeta(0.975, i, n - i + 1)))
  expect_true(all(one$sim_lo <= one$point_lo & one$sim_hi >= one$point_hi))
  expect_error(autoplot(a, type = "qq", level = 1.2), "probability")
  pc <- autoplot(a, type = "conditional")
  bc <- expect_builds(pc)
  expect_equal(nrow(bc$data[[which(layer_classes(pc) == "GeomCol")]]), 2L)
  # a single outcome leaks no legend
  a1 <- ref_assess(ref_fit(pf$fit$spec, pf$data, "y"), newdata = val)
  p1 <- autoplot(a1, type = "calibration")
  expect_builds(p1)
  expect_identical(colour_guide(p1), "none")
  expect_builds(autoplot(a1, type = "worm"))
})

test_that("score plots: profile as dot and segment with guides, heatmap with ordered ids", {
  pf <- plot_fit(39)
  new <- pf$data[1:8, ]
  new$participant_id <- paste0("S", 8:1)
  sc <- predict(pf$fit, newdata = new, uncertainty = "conditional")
  p <- autoplot(sc, type = "profile")
  b <- expect_builds(p)
  cls <- layer_classes(p)
  expect_true(all(c("GeomVline", "GeomSegment", "GeomPoint") %in% cls))
  guides <- b$data[[which(cls == "GeomVline")[2]]]$xintercept
  expect_equal(sort(guides), c(-2, 2))
  expect_match(plot_label(p, "subtitle"), "subject S8")
  p2 <- autoplot(sc, type = "profile", id = "S3")
  expect_match(plot_label(p2, "subtitle"), "subject S3")
  expect_error(autoplot(sc, type = "profile", id = "nope"), "No rows")
  ph <- autoplot(sc, type = "heatmap")
  bh <- expect_builds(ph)
  expect_equal(levels(ph$data$.id), paste0("S", 8:1))
  expect_equal(nrow(bh$data[[1]]), 16L)
  expect_equal(plot_label(ph, "fill"), "z")
})

test_that("dynamics plots: kernel, held-out calibration, transitions, and anchored fan", {
  skip_on_cran()
  dat <- ref_simulate(600, kind = "longitudinal", seed = 40)
  fit <- ref_fit(
    ref_spec(family = ref_gaussian(), location = ~ s(age, k = 5) + sex, scale = ~1),
    data = dat, outcomes = "y"
  )
  dyn <- ref_dynamics(fit, data = dat, id = participant_id, time = age, crossfit = 0)
  expect_true(dyn$processes$y$identified)
  pk <- autoplot(dyn, type = "kernel")
  bk <- expect_builds(pk)
  k <- fortify_kernel(dyn)
  expect_true(all(k$correlation <= 1 + 1e-8 & k$correlation >= 0))
  expect_equal(max(k$lag), 1.5 * dyn$lag_range[[2]])
  expect_true("GeomRect" %in% layer_classes(pk))
  held <- ref_simulate(300, kind = "longitudinal", seed = 41)
  pc <- autoplot(dyn, type = "calibration", data = held, id = participant_id, time = age)
  bc <- expect_builds(pc)
  expect_match(plot_label(pc, "subtitle"), "held-out innovation Z")
  expect_true(grepl("\n", plot_label(pc, "subtitle"), fixed = TRUE))
  expect_equal(plot_label(pc, "y"), "Innovation Z")
  tr <- ref_transition(dyn, data = held, id = participant_id, time = age)
  pts <- bc$data[[which(layer_classes(pc) == "GeomPoint")]]
  expect_equal(nrow(pts), sum(is.finite(tr$innovation_z)))
  expect_equal(sum(layer_classes(pc) == "GeomRibbon"), 2L)
  expect_error(autoplot(dyn, type = "calibration"), "data")
  expect_error(autoplot(dyn, type = "calibration", data = held), "id")
  for (ty in c("velocity", "innovation", "change")) {
    pt <- autoplot(tr, type = ty)
    expect_builds(pt)
    expect_equal(plot_label(pt, "x"), "Elapsed time")
  }
  hist <- dat[dat$participant_id == dat$participant_id[[1]], ]
  fc <- ref_forecast(dyn, history = hist, times = max(hist$age) + c(1, 2, 3))
  pf <- autoplot(fc, type = "fan")
  bf <- expect_builds(pf)
  cls <- layer_classes(pf)
  expect_true(all(c("GeomRibbon", "GeomLine", "GeomPath", "GeomPoint") %in% cls))
  # the fan is anchored at the last observed visit: every centile starts there
  ribbon <- bf$data[[which(cls == "GeomRibbon")[1]]]
  last <- hist[which.max(hist$age), ]
  first <- ribbon[ribbon$x == min(ribbon$x), ]
  expect_equal(unique(first$x), last$age)
  expect_equal(unique(first$ymin), last$y)
  expect_equal(unique(first$ymax), last$y)
  expect_error(autoplot(fc, type = "thrive"))
})

test_that("theme and model card", {
  expect_s3_class(theme_referent(), "theme")
  expect_s3_class(theme_referent(grid = "both"), "theme")
  expect_s3_class(theme_referent(grid = "none"), "theme")
  pf <- plot_fit(42, n = 80)
  card <- ref_freeze(pf$fit)
  expect_s3_class(card, "ref_freeze")
  out <- paste(cli::cli_fmt(print(card)), collapse = "\n")
  expect_match(out, "model card")
  expect_match(out, "age: \\[")
})

test_that("thrive lines are ref_forecast() on a grid of one-visit histories", {
  dat <- ref_simulate(300, kind = "longitudinal", seed = 71)
  fit <- ref_fit(
    ref_spec(family = ref_gaussian(), location = ~ age + sex, scale = ~1),
    data = dat, outcomes = "y"
  )
  dyn <- ref_dynamics(fit, data = dat, id = participant_id, time = age, crossfit = 0)
  th <- fortify_thrive(dyn, "y", anchors = c(0.05, 0.5, 0.95),
                       from = c(40, 55), horizon = 2, thrive = 0.025, level = 0.9)
  # the anchor sits exactly on the population centile curve and carries no
  # correlation uncertainty of its own
  anchor <- th[th$horizon == 0, ]
  expect_equal(anchor$centile, anchor$anchor_centile, tolerance = 1e-10)
  expect_equal(anchor$lower, anchor$value, tolerance = 1e-8)
  expect_equal(anchor$r, rep(1, nrow(anchor)), tolerance = 1e-12)
  # every forward point is ref_forecast()'s quantile at the thrive level
  tmpl <- centile_grid(fit, "age", n = 1)
  for (s in unique(th$.segment)) {
    seg <- th[th$.segment == s, ]
    hist <- tmpl
    hist$age <- seg$anchor_time[[1]]
    hist$y <- seg$value[seg$horizon == 0]
    fc <- ref_forecast(dyn, history = hist, times = seg$time[seg$horizon > 0],
                        outcome = "y")
    expect_equal(seg$value[seg$horizon > 0],
                 as.numeric(dist_quantile(fc$dist, 0.025)), tolerance = 1e-6)
  }
  # thrive = 0.5 is the conditional-median (regression-to-the-mean) path
  mid <- fortify_thrive(dyn, "y", anchors = 0.9, from = 50, horizon = 3, thrive = 0.5)
  fwd <- mid[mid$horizon == 3, ]
  expect_equal(fwd$z, fwd$r * stats::qnorm(0.9), tolerance = 1e-10)
  expect_true(fwd$z < stats::qnorm(0.9))
  # the correlation interval widens with the horizon
  wide <- fortify_thrive(dyn, "y", anchors = 0.5, from = 50,
                         horizon = c(1, 6), thrive = 0.025, level = 0.9)
  w <- wide$upper - wide$lower
  expect_true(w[wide$horizon == 6] > w[wide$horizon == 1])
  expect_true(all(w >= 0))
  # lags beyond the reference range are flagged, not silently drawn
  expect_true(any(wide$support == "extrapolated_lag"))
  expect_error(fortify_thrive(dyn, "y", anchors = 1.2), "must lie in")
  expect_error(fortify_thrive(dyn, "nope"), "not an outcome")
  p <- autoplot(dyn, type = "thrive", outcome = "y", from = c(40, 55))
  expect_builds(p)
  expect_true("GeomRibbon" %in% layer_classes(p))
})
