# Snapshot tests of the print methods. Snapshots are skipped on CRAN; all
# inputs are seeded and small so the output is deterministic locally.

# Numbers -> "#" and whitespace runs collapsed, so pillar column widths (which
# depend on digit counts) cannot make snapshots brittle across platforms.
snap_transform <- function(x) {
  x <- gsub("-?[0-9]+\\.[0-9]+", "#", x)
  gsub("[ \t]+", " ", x)
}

snap_fit <- function() {
  dat <- norm_simulate(120, seed = 90)
  dat$marker_bad <- 1
  norm_fit(
    norm_spec(norm_gaussian(), location = ~ s(age, k = 5) + sex, scale = ~1),
    data = dat, outcomes = c("y", "marker_01", "marker_bad")
  )
}

test_that("print.norm_fit shows family, formulas, and status counts", {
  withr::local_options(cli.num_colors = 1L, width = 80)
  fit <- snap_fit()
  expect_snapshot(print(fit))
  expect_snapshot(print(fit$spec))
  expect_snapshot(print(norm_shash()))
})

test_that("print.norm_scores counts rows and warns about in-sample scores", {
  withr::local_options(cli.num_colors = 1L, width = 80)
  fit <- snap_fit()
  dat <- norm_simulate(120, seed = 90)
  sc <- predict(fit, newdata = dat[1:3, c("age", "sex", "y", "marker_01")],
                uncertainty = "conditional")
  expect_snapshot(print(sc[, c(".row", ".outcome", "z", "support", "status")]))
})

test_that("print.norm_reference prints the model card", {
  withr::local_options(cli.num_colors = 1L, width = 80)
  fit <- snap_fit()
  ref <- norm_reference(fit, units = "mm")
  ref$versions <- list(referent = "x.y.z", mgcv = "x.y", r = "x.y.z")
  expect_snapshot(print(ref))
})

test_that("print.norm_assessment shows the overall table", {
  withr::local_options(cli.num_colors = 1L, width = 80)
  fit <- snap_fit()
  val <- norm_simulate(60, seed = 91)
  a <- norm_assess(fit, newdata = val)
  expect_snapshot(print(a), transform = snap_transform)
})

test_that("print.norm_dynamics reports components, identifiability, and the kernel", {
  withr::local_options(cli.num_colors = 1L, width = 80)
  dat <- simulate_two_visit(200, r = 0.6, lag = 2, seed = 92)
  fit <- norm_fit(norm_spec(norm_gaussian(), location = ~ s(age, k = 5) + sex, scale = ~1),
                  data = dat, outcomes = "y")
  dyn <- norm_dynamics(fit, data = dat, id = participant_id, time = age, crossfit = 0)
  expect_snapshot(print(dyn), transform = snap_transform)
  single <- dat[dat$visit == 1, ]
  dyn0 <- norm_dynamics(fit, data = single, id = participant_id, time = age, crossfit = 0)
  expect_snapshot(print(dyn0), transform = snap_transform)
})

test_that("remaining print methods use cli and return invisibly", {
  withr::local_options(cli.num_colors = 1L, width = 80)
  dat <- norm_simulate(150, seed = 93)
  fit <- norm_fit(simple_spec(), data = dat[1:100, ], outcomes = "y")
  ad <- norm_adapt(fit, data = dat[101:130, ], by = site)
  expect_snapshot(print(ad$adaptation), transform = snap_transform)
  cal <- norm_calibrate(fit, data = dat[101:150, ], by = site)
  expect_snapshot(print(cal$calibration))
  sc <- predict(fit, newdata = dat[101:150, ], uncertainty = "conditional")
  jt <- norm_joint(sc)
  expect_snapshot(print(jt))
  msg <- cli::cli_fmt(res <- withVisible(print(jt)))
  expect_false(res$visible)
  dyn <- norm_dynamics(fit, data = norm_simulate(200, kind = "longitudinal", seed = 94),
                       id = participant_id, time = age, crossfit = 0)
  tr <- norm_transition(dyn, data = norm_simulate(30, kind = "longitudinal", seed = 95),
                        id = participant_id, time = age)
  body <- utils::capture.output(msg <- cli::cli_fmt(print(tr)))
  expect_match(paste(msg, collapse = "\n"), "norm_transition")
  expect_true(any(grepl("innovation_z", body)))
})
