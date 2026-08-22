ladder_specs <- function() {
  list(
    gaussian_const = ref_spec(ref_gaussian(), location = ~ s(age, k = 5) + sex),
    gaussian_scale = ref_spec(ref_gaussian(), location = ~ s(age, k = 5) + sex,
                               scale = ~ s(age, k = 4)),
    shash_const = ref_spec(ref_shash(), location = ~ s(age, k = 5) + sex,
                            scale = ~ s(age, k = 4)),
    shash_shape = ref_spec(ref_shash(), location = ~ s(age, k = 5) + sex,
                            scale = ~ s(age, k = 4), skew = ~ s(age, k = 4))
  )
}

test_that("the ladder picks a Gaussian level on Gaussian constant-scale data", {
  skip_on_cran()
  dat <- ref_simulate(400, seed = 61)
  withr::with_seed(61, {
    sel <- ref_select(ladder_specs(), data = dat, outcomes = "y", folds = 3)
  })
  expect_s3_class(sel, "ref_selection")
  expect_true(sel$selected_name %in% c("gaussian_const", "gaussian_scale"))
  tab <- sel$comparison
  expect_true(all(c("crps", "mean_z", "var_z", "cover_95", "tail_05", "shape") %in% names(tab)))
  expect_equal(tab$shape, c(FALSE, FALSE, FALSE, TRUE))
  expect_true(all(is.finite(tab$crps[tab$status == "ok"])))
  # gate values are out-of-fold: they differ from the in-sample assessment
  fit <- attr(sel$crossfit, "deployment")
  ins <- ref_assess(fit, newdata = dat)
  oof <- tab[tab$selected, ]
  expect_false(isTRUE(all.equal(oof$var_z, ins$marginal$var_z)))
  expect_equal(oof$var_z, stats::var(sel$crossfit$z, na.rm = TRUE))
  expect_true(all(!sel$crossfit$.in_sample))
  # the crossfit carries per-row CRPS and a global row index
  expect_equal(sort(sel$crossfit$.row), seq_len(nrow(dat)))
  expect_equal(mean(sel$crossfit$crps), oof$crps)
})

test_that("the ladder picks SHASH on strongly skewed data", {
  skip_on_cran()
  dat <- simulate_skewed(500, seed = 62, skew = 1.4, tail = 0.7)
  withr::with_seed(62, {
    sel <- ref_select(ladder_specs(), data = dat, outcomes = "y", folds = 3)
  })
  expect_true(grepl("^shash", sel$selected_name))
  tab <- sel$comparison
  expect_gt(
    max(tab$mean_log_score[grepl("shash", tab$model)]),
    max(tab$mean_log_score[grepl("gaussian", tab$model)])
  )
})

test_that("acceptable_calibration scales its tolerance with n", {
  make <- function(n, mean_z, var_z, cover_95, tail_05 = 0.05) {
    list(
      marginal = tibble::tibble(.outcome = "y", n = n, mean_z = mean_z, var_z = var_z,
                                cover_95 = cover_95),
      tail = tibble::tibble(.outcome = "y", tail_level = 0.05, observed = tail_05, n = n)
    )
  }
  expect_true(acceptable_calibration(make(1000, 0.05, 1.05, 0.95)))
  expect_false(acceptable_calibration(make(1000, 0.3, 1.05, 0.95)))
  expect_false(acceptable_calibration(make(1000, 0.0, 1.6, 0.95)))
  expect_false(acceptable_calibration(make(1000, 0.0, 1.0, 0.90)))
  expect_false(acceptable_calibration(make(1000, 0.0, 1.0, 0.95, tail_05 = 0.10)))
  # the same deviations are within noise at n = 30
  expect_true(acceptable_calibration(make(30, 0.3, 1.6, 0.90)))
})

test_that("the comparison table carries the paired standard errors the one-SE rule uses", {
  dat <- ref_simulate(200, seed = 15)
  specs <- list(
    linear = ref_spec(ref_gaussian(), ~ age + sex),
    smooth = ref_spec(ref_gaussian(), ~ s(age, k = 5) + sex)
  )
  set.seed(1)
  sel <- ref_select(specs, dat, "y", folds = 3)
  tab <- sel$comparison
  expect_false("se_log_score" %in% names(tab))
  expect_true(all(c("se_log_score_paired", "se_crps_paired") %in% names(tab)))
  best <- tab$model[which.max(tab$mean_log_score)]
  expect_equal(unname(tab$se_log_score_paired[tab$model == best]), 0)
  other <- tab[tab$model != best, ]
  cf_b <- sel$crossfit
  expect_true(all(is.finite(other$se_log_score_paired)))
  expect_true(all(other$se_log_score_paired > 0))
})
