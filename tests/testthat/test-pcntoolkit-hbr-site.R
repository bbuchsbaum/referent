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
  expect_match(text, '"hbr_divergences.csv"', fixed = TRUE)
  expect_match(text, '"sha256"', fixed = TRUE)
  expect_match(text, '"reported_z": "mean of draw-specific z"', fixed = TRUE)
  expect_match(text, '"mixture_outputs"', fixed = TRUE)
  expect_match(text, '"hbr_public_random_seed_argument": False', fixed = TRUE)
})

test_that("retained HBR and site artifacts recompute every registered gate", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  evidence <- file.path(root, "docs", "evidence", "pcntoolkit", "v1.3.0")
  hbr_path <- file.path(evidence, "hbr_receipt.json")
  site_path <- file.path(evidence, "site_referent_receipt.json")
  skip_if_not(all(file.exists(hbr_path, site_path)), "release receipts are not bundled")
  validator_path <- system.file(
    "pcntoolkit", "validate_hbr_site_evidence.R", package = "referent"
  )
  validator <- new.env(parent = globalenv())
  sys.source(validator_path, envir = validator)
  validated <- validator$pcn_validate_site_evidence(
    evidence, evidence, "hbr_receipt.json", "site_referent_receipt.json"
  )
  expect_true(all(validated$hbr$stages$pass))
  expect_true(all(validated$hbr$stages$divergences == 0L))
  expect_true(all(validated$comparison$pass))
  expect_equal(
    validated$receipt$adaptation_prior_selection$selected$location_prior_n, 0
  )
})

test_that("HBR and site validators reject drift, false verdicts, and omissions", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  evidence <- file.path(root, "docs", "evidence", "pcntoolkit", "v1.3.0")
  validator_path <- system.file(
    "pcntoolkit", "validate_hbr_site_evidence.R", package = "referent"
  )
  validator <- new.env(parent = globalenv())
  sys.source(validator_path, envir = validator)
  files <- c(
    "site_data.csv", "pcntoolkit_predictions.csv", "pcntoolkit_site_summary.csv",
    "hbr_convergence.csv", "hbr_divergences.csv", "hbr_receipt.json",
    "referent_predictions.csv", "site_comparison.csv", "adaptation_selection.csv",
    "adaptation_selection_summary.csv", "site_referent_receipt.json"
  )
  skip_if_not(all(file.exists(file.path(evidence, files))), "raw receipts are not bundled")

  drift <- file.path(tempdir(), "referent-hbr-hash-drift")
  dir.create(drift, recursive = TRUE, showWarnings = FALSE)
  expect_true(all(file.copy(file.path(evidence, files), drift, overwrite = TRUE)))
  convergence <- utils::read.csv(file.path(drift, "hbr_convergence.csv"))
  convergence$r_hat[[1L]] <- convergence$r_hat[[1L]] + 0.001
  utils::write.csv(convergence, file.path(drift, "hbr_convergence.csv"), row.names = FALSE)
  expect_error(
    validator$pcn_validate_hbr_evidence(drift, "hbr_receipt.json"),
    "receipt mismatch"
  )

  hbr_verdict <- file.path(tempdir(), "referent-hbr-verdict-drift")
  dir.create(hbr_verdict, recursive = TRUE, showWarnings = FALSE)
  expect_true(all(file.copy(
    file.path(evidence, files), hbr_verdict, overwrite = TRUE
  )))
  hbr <- jsonlite::read_json(
    file.path(hbr_verdict, "hbr_receipt.json"), simplifyVector = FALSE
  )
  hbr$all_stages_pass <- FALSE
  jsonlite::write_json(
    hbr, file.path(hbr_verdict, "hbr_receipt.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null"
  )
  expect_error(
    validator$pcn_validate_hbr_evidence(hbr_verdict, "hbr_receipt.json"),
    "verdict does not match"
  )

  verdict <- file.path(tempdir(), "referent-site-verdict-drift")
  dir.create(verdict, recursive = TRUE, showWarnings = FALSE)
  expect_true(all(file.copy(file.path(evidence, files), verdict, overwrite = TRUE)))
  site <- jsonlite::read_json(
    file.path(verdict, "site_referent_receipt.json"), simplifyVector = FALSE
  )
  site$all_sites_pass <- FALSE
  jsonlite::write_json(
    site, file.path(verdict, "site_referent_receipt.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null"
  )
  expect_error(
    validator$pcn_validate_site_evidence(
      verdict, verdict, "hbr_receipt.json", "site_referent_receipt.json"
    ),
    "verdict does not match"
  )

  missing_draw <- file.path(tempdir(), "referent-hbr-resigned-missing-draw")
  dir.create(missing_draw, recursive = TRUE, showWarnings = FALSE)
  expect_true(all(file.copy(
    file.path(evidence, files), missing_draw, overwrite = TRUE
  )))
  divergence <- utils::read.csv(file.path(missing_draw, "hbr_divergences.csv"))
  divergence <- divergence[-nrow(divergence), , drop = FALSE]
  utils::write.csv(
    divergence, file.path(missing_draw, "hbr_divergences.csv"), row.names = FALSE
  )
  hbr <- jsonlite::read_json(
    file.path(missing_draw, "hbr_receipt.json"), simplifyVector = FALSE
  )
  hbr$files[["hbr_divergences.csv"]] <- validator$pcn_evidence_file_receipt(
    file.path(missing_draw, "hbr_divergences.csv")
  )
  jsonlite::write_json(
    hbr, file.path(missing_draw, "hbr_receipt.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null"
  )
  expect_error(
    validator$pcn_validate_hbr_evidence(missing_draw, "hbr_receipt.json"),
    "draw identities"
  )

  missing_prediction <- file.path(
    tempdir(), "referent-site-resigned-missing-prediction"
  )
  dir.create(missing_prediction, recursive = TRUE, showWarnings = FALSE)
  expect_true(all(file.copy(
    file.path(evidence, files), missing_prediction, overwrite = TRUE
  )))
  prediction <- utils::read.csv(
    file.path(missing_prediction, "referent_predictions.csv")
  )
  prediction <- prediction[-nrow(prediction), , drop = FALSE]
  utils::write.csv(
    prediction, file.path(missing_prediction, "referent_predictions.csv"),
    row.names = FALSE
  )
  site <- jsonlite::read_json(
    file.path(missing_prediction, "site_referent_receipt.json"),
    simplifyVector = FALSE
  )
  site$output_files[["referent_predictions.csv"]] <-
    validator$pcn_evidence_file_receipt(
      file.path(missing_prediction, "referent_predictions.csv")
    )
  jsonlite::write_json(
    site, file.path(missing_prediction, "site_referent_receipt.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null"
  )
  expect_error(
    validator$pcn_validate_site_evidence(
      missing_prediction, missing_prediction,
      "hbr_receipt.json", "site_referent_receipt.json"
    ),
    "cover every retained evaluation row"
  )
})
