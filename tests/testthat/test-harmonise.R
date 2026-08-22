# Harmonised outputs: conditional quantile mapping out of the group effect.

harmonise_check <- function(fit, dat, by = "site", tol = 1e-8) {
  h <- predict(fit, dat, type = "harmonised", uncertainty = "conditional")
  u0 <- predict(fit, dat, uncertainty = "conditional", allow_extrapolation = TRUE)$centile
  ok <- Filter(fit_ok, fit$models)
  excl <- unique(unlist(lapply(ok, function(m) smooths_using(m$model, by))))
  target <- predict_dists(fit, dat, uncertainty = "conditional", exclude = excl,
                          adapt = FALSE)$y
  u_h <- dist_cdf(target, h$y)
  expect_equal(u_h, u0, tolerance = tol)
  h
}

test_that("harmonised values keep their centile under the zero-site predictive", {
  dat <- ref_simulate(600, sites = 4, site_shift = c(0, 1.2, -0.8, 0.3), seed = 31)
  fit <- ref_fit(ref_spec(ref_gaussian(), ~ s(age, k = 6) + sex + s(site, bs = "re")),
                 dat, outcomes = "y")
  h <- harmonise_check(fit, dat)
  expect_equal(names(h), names(dat))
  expect_equal(h$site, dat$site)
  # site means of harmonised y agree within sampling error once the
  # age/sex structure is removed
  pop <- predict_dists(fit, dat, uncertainty = "conditional", exclude = "s(site)")$y
  r <- h$y - as.numeric(mean(pop))
  m <- tapply(r, dat$site, mean)
  se <- tapply(r, dat$site, function(x) stats::sd(x) / sqrt(length(x)))
  expect_true(all(abs(m - mean(r)) < 2.5 * se))
  raw <- tapply(dat$y - as.numeric(mean(pop)), dat$site, mean)
  expect_gt(diff(range(raw)), 1.5)
  expect_lt(diff(range(m)), 0.4)
  # mapping to a level is the same map plus that level's effect
  hB <- predict(fit, dat, type = "harmonised", to = "B", uncertainty = "conditional")
  expect_lt(stats::sd(hB$y - h$y), 1e-6)
  expect_gt(mean(hB$y - h$y), 0.5)
})

test_that("harmonisation works for SHASH and log-transformed fits and after adaptation", {
  dat <- ref_simulate(500, kind = "shash", sites = 3, site_shift = c(0, 1, -1), seed = 32)
  dat$y <- dat$y - min(dat$y) + 1
  spec <- ref_spec(ref_shash(), ~ s(age, k = 5) + s(site, bs = "re"), transform = "log")
  fit <- ref_fit(spec, dat, outcomes = "y")
  h <- harmonise_check(fit, dat, tol = 1e-6)
  expect_true(all(h$y > 0))
  # adaptation to a new site: harmonise the adapted predictive to the reference
  new <- ref_simulate(200, sites = 1, site_shift = 2, seed = 33)
  new$y <- new$y - min(new$y) + 3
  new$site <- factor("Z")
  ad <- ref_adapt(fit, new, by = site)
  expect_gt(abs(ad$adaptation$offsets$y$Z$location), 0.05)
  h_ad <- harmonise_check(ad, new, tol = 1e-6)
  h_raw <- predict(fit, new, type = "harmonised", uncertainty = "conditional")
  expect_false(isTRUE(all.equal(h_ad$y, h_raw$y)))
})

test_that("by is required when it cannot be inferred", {
  dat <- ref_simulate(100, seed = 34)
  fit <- ref_fit(ref_spec(ref_gaussian(), ~ age + sex), dat, outcomes = "y")
  expect_error(predict(fit, dat, type = "harmonised"), "by")
})

# Site-specific trajectories: factor smooths.

test_that("factor-smooth site trajectories fit, freeze, predict for unseen sites, and adapt", {
  dat <- ref_simulate(500, sites = 4, site_shift = c(0, 1.2, -0.8, 0.3), seed = 35)
  spec <- ref_spec(ref_gaussian(), ~ s(age, k = 6) + s(age, site, bs = "fs", k = 5))
  fit <- ref_fit(spec, dat, outcomes = "y")
  expect_equal(fit$models$y$status, "ok")
  new <- dat[1:6, ]
  new$site <- factor(c("A", "Z", "B", "Z", "C", "D"))
  sc <- predict(fit, new, uncertainty = "conditional", allow_extrapolation = TRUE)
  expect_true(all(is.finite(sc$z)))
  expect_equal(sc$support[2], "new_group")
  # an unseen site gets the population curve: the fs term excluded
  m <- fit$models$y$model
  pop <- suppressWarnings(predict(m, new[2, ], exclude = "s(age,site)"))
  expect_equal(sc$median[2], unname(as.numeric(pop)), tolerance = 1e-8)
  # seen sites keep their own curves
  expect_equal(sc$median[c(1, 3)],
               unname(as.numeric(predict(m, new[c(1, 3), ]))), tolerance = 1e-8)
  sct <- predict(fit, new, uncertainty = "total", allow_extrapolation = TRUE)
  expect_true(all(is.finite(sct$z)))
  frozen <- ref_freeze(fit)
  scf <- predict(frozen, new, uncertainty = "conditional", allow_extrapolation = TRUE)
  expect_equal(scf$z, sc$z, tolerance = 1e-10)
  local <- ref_simulate(150, sites = 1, site_shift = 1.5, seed = 36)
  local$site <- factor("Z")
  ad <- ref_adapt(fit, local, by = site)
  expect_gt(ad$adaptation$offsets$y$Z$location, 0.8)
  sca <- predict(ad, local, uncertainty = "conditional", allow_extrapolation = TRUE)
  expect_lt(abs(mean(sca$z)), 0.2)
  # harmonising out the factor smooth removes the site shifts
  h <- harmonise_check(fit, dat)
  shift <- tapply(h$y - dat$y, dat$site, mean)
  expect_lt(abs(shift[["B"]] + 1.2), 0.3)
  expect_lt(abs(shift[["C"]] - 0.8), 0.3)
})
