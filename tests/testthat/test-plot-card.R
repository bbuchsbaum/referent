plot_fit <- function(seed = 30, n = 120) {
  dat <- norm_simulate(n, seed = seed)
  dat$participant_id <- paste0("P", seq_len(nrow(dat)))
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ s(age, k = 5) + sex, scale = ~1),
    data = dat,
    outcomes = c("y", "marker_01"),
    id = participant_id
  )
  list(fit = fit, data = dat)
}

expect_builds <- function(p) {
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
  invisible(built)
}

layer_classes <- function(p) {
  vapply(p$layers, function(ly) class(ly$geom)[[1]], character(1))
}

plot_label <- function(p, aes) {
  ggplot2::get_labs(p)[[aes]]
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
  dat <- norm_simulate(150, kind = "longitudinal", seed = 32)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ s(age, k = 5) + sex, scale = ~1),
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

test_that("support, visits, coverage, and composition use the fit", {
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
  # visits come from the id declared at fit time
  pv <- autoplot(pf$fit, type = "visits")
  bv <- expect_builds(pv)
  expect_equal(sum(bv$data[[1]]$y), nrow(pf$data))
  long <- norm_simulate(120, kind = "longitudinal", seed = 34)
  pv2 <- autoplot(pf$fit, type = "visits", data = long, id = participant_id)
  bv2 <- expect_builds(pv2)
  expect_equal(sum(bv2$data[[1]]$y), length(unique(long$participant_id)))
  expect_equal(fortify_visits(long, participant_id)$n_subjects, bv2$data[[1]]$y)
  fit_noid <- norm_fit(pf$fit$spec, data = pf$data, outcomes = "y")
  expect_error(autoplot(fit_noid, type = "visits"), "identifier")
  # coverage draws the reference range from the fit
  pc <- autoplot(pf$fit, type = "coverage", data = long, by = site)
  bc <- expect_builds(pc)
  cls <- layer_classes(pc)
  expect_true("GeomVline" %in% cls)
  r <- pf$fit$support_ref$numeric$age
  expect_equal(sort(bc$data[[which(cls == "GeomVline")]]$xintercept), c(r$min, r$max))
  expect_equal(plot_label(pc, "x"), "age")
  expect_builds(autoplot(pf$fit, type = "coverage", data = long))
  expect_error(autoplot(pf$fit, type = "coverage"), "data")
  cov <- fortify_coverage(long, age, by = site)
  expect_true(all(c("time", ".group") %in% names(cov)))
  skip_if_not_installed("patchwork")
  pcomp <- autoplot(pf$fit, type = "composition", data = long, id = participant_id, by = site)
  expect_s3_class(pcomp, "patchwork")
})

test_that("adaptation plot shows offsets with intervals per group", {
  pf <- plot_fit(35, n = 200)
  local <- norm_simulate(120, site_shift = c(0, 1, 0, 0), seed = 36)
  ad <- norm_adapt(pf$fit, data = local, by = site)
  p <- autoplot(ad, type = "adaptation")
  b <- expect_builds(p)
  cls <- layer_classes(p)
  expect_true(all(c("GeomLinerange", "GeomPoint") %in% cls))
  pts <- b$data[[which(cls == "GeomPoint")]]
  expect_equal(nrow(pts), 4L * 2L)
  expect_equal(plot_label(p, "x"), "site")
  expect_error(autoplot(pf$fit, type = "adaptation"), "norm_adapt")
  # pooled adaptation of a single outcome has no colour legend
  fit1 <- norm_fit(pf$fit$spec, data = pf$data, outcomes = "y")
  ad1 <- norm_adapt(fit1, data = local)
  p1 <- autoplot(ad1, type = "adaptation")
  expect_builds(p1)
  expect_false("colour" %in% names(p1$mapping))
})

test_that("assessment plots: calibration, worm, conditional", {
  pf <- plot_fit(37, n = 150)
  val <- norm_simulate(80, seed = 38)
  a <- norm_assess(pf$fit, newdata = val)
  p <- autoplot(a, type = "calibration")
  b <- expect_builds(p)
  cls <- layer_classes(p)
  expect_equal(nrow(b$data[[which(cls == "GeomPoint")]]), 2L * 5L)
  expect_equal(plot_label(p, "colour"), "Outcome")
  pw <- autoplot(a, type = "worm")
  bw <- expect_builds(pw)
  expect_true("GeomSmooth" %in% layer_classes(pw))
  expect_equal(nrow(bw$data[[which(layer_classes(pw) == "GeomPoint")]]), 2L * 80L)
  pc <- autoplot(a, type = "conditional")
  bc <- expect_builds(pc)
  expect_equal(nrow(bc$data[[which(layer_classes(pc) == "GeomCol")]]), 2L)
  # a single outcome maps no colour and leaks no legend
  a1 <- norm_assess(norm_fit(pf$fit$spec, pf$data, "y"), newdata = val)
  p1 <- autoplot(a1, type = "calibration")
  expect_builds(p1)
  expect_null(plot_label(p1, "colour"))
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
  dat <- norm_simulate(600, kind = "longitudinal", seed = 40)
  fit <- norm_fit(
    norm_spec(family = norm_gaussian(), location = ~ s(age, k = 5) + sex, scale = ~1),
    data = dat, outcomes = "y"
  )
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age, crossfit = 0)
  expect_true(dyn$processes$y$identified)
  pk <- autoplot(dyn, type = "kernel")
  bk <- expect_builds(pk)
  k <- fortify_kernel(dyn)
  expect_true(all(k$correlation <= 1 + 1e-8 & k$correlation >= 0))
  expect_equal(max(k$lag), 1.5 * dyn$lag_range[[2]])
  expect_true("GeomRect" %in% layer_classes(pk))
  held <- norm_simulate(300, kind = "longitudinal", seed = 41)
  pc <- autoplot(dyn, type = "calibration", data = held, id = participant_id, time = age)
  bc <- expect_builds(pc)
  expect_match(plot_label(pc, "subtitle"), "held-out innovation Z")
  tr <- norm_transition(dyn, data = held, id = participant_id, time = age)
  expect_equal(nrow(bc$data[[2]]), sum(is.finite(tr$innovation_z)))
  expect_error(autoplot(dyn, type = "calibration"), "data")
  expect_error(autoplot(dyn, type = "calibration", data = held), "id")
  for (ty in c("velocity", "innovation", "change")) {
    pt <- autoplot(tr, type = ty)
    expect_builds(pt)
    expect_equal(plot_label(pt, "x"), "Elapsed time")
  }
  hist <- dat[dat$participant_id == dat$participant_id[[1]], ]
  fc <- norm_forecast(dyn, history = hist, times = max(hist$age) + c(1, 2, 3))
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
  card <- norm_reference(pf$fit)
  expect_s3_class(card, "norm_reference")
  out <- paste(cli::cli_fmt(print(card)), collapse = "\n")
  expect_match(out, "model card")
  expect_match(out, "age: \\[")
})
