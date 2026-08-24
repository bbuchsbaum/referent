test_that("HBR draw averages and Referent mixtures are explicitly non-equivalent", {
  mu <- c(-1.2, 0.4, 1.7)
  sigma <- c(0.7, 1.1, 1.8)
  y <- 0.8
  p <- 0.9
  draw_z <- (y - mu) / sigma
  reported_z <- mean(draw_z)
  mixture_z <- stats::qnorm(mean(stats::pnorm(draw_z)))
  reported_quantile <- mean(stats::qnorm(p, mu, sigma))
  mixture_quantile <- stats::uniroot(
    function(q) mean(stats::pnorm(q, mu, sigma)) - p, c(-20, 20)
  )$root
  reported_logp <- mean(stats::dnorm(y, mu, sigma, log = TRUE))
  mixture_logp <- log(mean(stats::dnorm(y, mu, sigma)))

  expect_gt(abs(reported_z - mixture_z), 0.05)
  expect_gt(abs(reported_quantile - mixture_quantile), 0.05)
  expect_gt(mixture_logp, reported_logp) # Jensen: log(mean density) >= mean(log density)
})

test_that("site gates cannot hide a small-site failure in a pooled average", {
  set.seed(12)
  large_u <- (1:950 - 0.5) / 950
  small_u <- rep(0.995, 50)
  prediction <- data.frame(
    site = rep(c("large", "small"), c(950, 50)),
    observed = c(stats::qnorm(large_u), rep(3, 50)),
    centile = c(large_u, small_u),
    log_density = c(stats::dnorm(stats::qnorm(large_u), log = TRUE), rep(-6, 50)),
    q05 = stats::qnorm(0.05), q50 = 0, q95 = stats::qnorm(0.95)
  )
  summary <- pcn_site_summary(prediction, "negative_control", "transport")
  gated <- pcn_site_gate(summary)
  expect_true(gated$pass[gated$site == "large"])
  expect_false(gated$pass[gated$site == "small"])
  expect_false(all(gated$pass))
  expect_lt(abs(mean(prediction$centile) - 0.5), 0.03)
})

test_that("adaptation-prior selection uses every observed site and deterministic ties", {
  validation <- expand.grid(
    site = c("site-1", "site-2", "site-3"),
    location_prior_n = c(0, 5), scale_prior_n = c(10, 25),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  validation$mean_log_score <- with(
    validation,
    -1 - 0.02 * location_prior_n - 0.001 * abs(scale_prior_n - 25)
  )
  selected <- pcn_select_adaptation_candidate(validation)
  expect_equal(selected$selected$location_prior_n, 0)
  expect_equal(selected$selected$scale_prior_n, 25)
  expect_true(all(selected$summary$validated_sites == 3L))
  expect_identical(which(selected$summary$selected), 1L)
  expect_error(
    pcn_select_adaptation_candidate(validation[-1L, , drop = FALSE]),
    "incomplete"
  )
})

test_that("HBR release runner records mandatory convergence and aggregation fields", {
  script_path <- system.file(
    "pcntoolkit", "run_hbr_site_evidence.py", package = "referent"
  )
  script <- readLines(
    script_path, warn = FALSE
  )
  text <- paste(script, collapse = "\n")
  expect_match(text, '"nutpie==0.16.8"', fixed = TRUE)
  expect_match(text, '"max_rhat": 1.01', fixed = TRUE)
  expect_match(text, '"min_ess_bulk": 400', fixed = TRUE)
  expect_match(text, '"divergences": 0', fixed = TRUE)
  expect_match(text, '--target-accept", type=float, default=0.99', fixed = TRUE)
  expect_match(text, '"sampling_seeds"', fixed = TRUE)
  expect_match(text, "with_sampling_controls", fixed = TRUE)
  expect_match(text, '"reported_z": "mean of draw-specific z"', fixed = TRUE)
  expect_match(text, '"mixture_outputs"', fixed = TRUE)
  expect_match(text, '"hbr_public_random_seed_argument": False', fixed = TRUE)
})

test_that("retained HBR and site receipts satisfy every registered gate", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  evidence <- file.path(root, "docs", "evidence", "pcntoolkit", "v1.3.0")
  hbr_path <- file.path(evidence, "hbr_receipt.json")
  site_path <- file.path(evidence, "site_referent_receipt.json")
  skip_if_not(all(file.exists(hbr_path, site_path)), "release receipts are not bundled")
  hbr <- jsonlite::read_json(hbr_path, simplifyVector = TRUE)
  site <- jsonlite::read_json(site_path, simplifyVector = TRUE)
  expect_true(hbr$all_stages_pass)
  expect_true(all(hbr$stages$divergences == 0L))
  expect_equal(hbr$target_accept, 0.99)
  expect_true(site$all_sites_pass)
  expect_identical(site$adaptation_prior_selection$source_split, "reference_train")
  expect_equal(site$adaptation_prior_selection$selected$location_prior_n, 0)
})
