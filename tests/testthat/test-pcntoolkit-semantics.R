pcn_semantic_root <- testthat::test_path("..", "fixtures", "pcntoolkit", "v1.3.0")

read_pcn_semantic <- function(filename) {
  utils::read.csv(file.path(pcn_semantic_root, filename), check.names = FALSE)
}

expect_scale_close <- function(actual, expected, atol, rtol = atol) {
  expect_true(
    all(abs(actual - expected) <= atol + rtol * abs(expected)),
    info = paste("maximum scaled error:", max(abs(actual - expected) / (atol + rtol * abs(expected))))
  )
}

test_that("fixed Gaussian semantics agree with PCNtoolkit and analytic R", {
  oracle <- read_pcn_semantic("semantic_gaussian.csv")
  dist <- distributional::dist_normal(oracle$mu, oracle$sigma)

  expect_scale_close(dist_log_density(dist, oracle$y), oracle$log_density, 1e-10)
  expect_scale_close(exp(dist_log_density(dist, oracle$y)), oracle$density, 1e-10)
  expect_scale_close(dist_cdf(dist, oracle$y), oracle$cdf, 1e-10)
  expect_scale_close(
    dist_eval(dist, log_tail, oracle$y, lower.tail = TRUE),
    oracle$log_cdf, 1e-10
  )
  expect_scale_close(
    exp(dist_eval(dist, log_tail, oracle$y, lower.tail = FALSE)),
    oracle$upper_tail, 1e-10
  )
  expect_scale_close(
    dist_eval(dist, log_tail, oracle$y, lower.tail = FALSE),
    oracle$log_upper_tail, 1e-10
  )
  expect_scale_close(dist_quantile(dist, oracle$p), oracle$quantile, 1e-10)
  expect_scale_close(dist_z(dist, oracle$y), oracle$z, 1e-10)

  # Independent analytic path, rather than fixture-to-implementation only.
  expect_equal(oracle$cdf, stats::pnorm(oracle$y, oracle$mu, oracle$sigma), tolerance = 1e-14)
  expect_equal(oracle$quantile, stats::qnorm(oracle$p, oracle$mu, oracle$sigma), tolerance = 1e-14)
  expect_true(all(is.finite(oracle$log_cdf)))
  expect_true(all(is.finite(oracle$log_upper_tail)))
})

test_that("PCNtoolkit SHASHb is the converted mgcv SHASH distribution", {
  oracle <- read_pcn_semantic("semantic_shashb.csv")
  variance_s <- oracle$moment_2 - oracle$moment_1^2
  expect_equal(
    oracle$sigma_r,
    oracle$sigma_b / (oracle$delta * sqrt(variance_s)),
    tolerance = 1e-14
  )
  expect_equal(
    oracle$mu_r,
    oracle$mu_b - oracle$sigma_b * oracle$moment_1 / sqrt(variance_s),
    tolerance = 1e-14
  )

  dist <- dist_shash(oracle$mu_r, oracle$sigma_r, oracle$epsilon, oracle$delta)
  expect_scale_close(dist_log_density(dist, oracle$y), oracle$log_density, 1e-7)
  expect_scale_close(exp(dist_log_density(dist, oracle$y)), oracle$density, 1e-7)
  expect_scale_close(dist_cdf(dist, oracle$y), oracle$cdf, 1e-7)
  expect_scale_close(
    dist_eval(dist, log_tail, oracle$y, lower.tail = TRUE),
    oracle$log_cdf, 1e-7
  )
  expect_scale_close(
    dist_eval(dist, log_tail, oracle$y, lower.tail = FALSE),
    oracle$log_upper_tail, 1e-7
  )
  expect_scale_close(dist_quantile(dist, oracle$p), oracle$quantile, 1e-7)
  expect_scale_close(dist_z(dist, oracle$y), oracle$z, 1e-7)
  expect_true(all(is.finite(oracle$log_density)))
  expect_true(all(is.finite(oracle$log_cdf)))
  expect_true(all(is.finite(oracle$log_upper_tail)))

  # The conversion preserves PCNtoolkit's declared mean and standard deviation.
  expect_scale_close(mean(dist), oracle$mu_b, 1e-7)
  expect_scale_close(sqrt(variance(dist)), oracle$sigma_b, 1e-7)
})

test_that("response transformation parity includes support and Jacobian", {
  mu <- c(-0.5, 0.2, 1.1)
  sigma <- c(0.4, 1.0, 1.8)
  y <- c(0.03, 1.5, 30)
  p <- c(1e-6, 0.5, 1 - 1e-6)
  warped <- dist_warped(distributional::dist_normal(mu, sigma), 0)

  expect_scale_close(dist_cdf(warped, y), stats::pnorm(log(y), mu, sigma), 1e-10)
  expect_scale_close(dist_quantile(warped, p), exp(stats::qnorm(p, mu, sigma)), 1e-10)
  expect_scale_close(
    dist_log_density(warped, y),
    stats::dnorm(log(y), mu, sigma, log = TRUE) - log(y),
    1e-10
  )
  expect_equal(dist_cdf(warped[1], -1), 0)
  expect_equal(dist_log_density(warped[1], -1), -Inf)
})

test_that("metric normalisations match or differ exactly as declared", {
  oracle <- read_pcn_semantic("semantic_metrics.csv")
  value <- setNames(oracle$value, oracle$metric)
  observed <- c(-2, -0.2, 0.1, 0.7, 1.4, 2.8, 5)
  predicted <- c(-1.5, -0.4, 0.3, 0.5, 1.8, 2.1, 4.2)
  z <- c(-1.8, -0.9, -0.2, 0.1, 0.5, 1.4, 2.7)
  residual <- observed - predicted
  pit <- stats::pnorm(z)

  expect_equal(sqrt(mean(residual^2)), value[["rmse"]], tolerance = 1e-12)
  expect_equal(
    mean(residual^2) / mean((observed - mean(observed))^2),
    value[["smse_population"]], tolerance = 1e-12
  )
  expect_equal(
    mean(residual^2) / stats::var(observed),
    value[["smse_population"]] * (length(observed) - 1) / length(observed),
    tolerance = 1e-12
  )
  expect_equal(1 - stats::var(residual) / stats::var(observed), value[["ev"]], tolerance = 1e-12)
  expect_equal(stats::cor(observed, predicted, method = "spearman"), value[["rho_spearman"]], tolerance = 1e-12)
  expect_equal(stats::cor(observed, predicted), value[["pearson"]], tolerance = 1e-12)
  expect_gt(abs(value[["rho_spearman"]] - value[["pearson"]]), 1e-3)
  expect_equal(mace(pit), value[["mace_no_batch"]], tolerance = 1e-12)
  expect_equal(std_moment(z, 3), value[["skew_referent_sample_sd"]], tolerance = 1e-12)
  expect_equal(std_moment(z, 4) - 3, value[["excess_kurtosis_referent_sample_sd"]], tolerance = 1e-12)
  expect_gt(abs(value[["skew_bias_corrected"]] - value[["skew_referent_sample_sd"]]), 0.05)
  expect_gt(abs(value[["excess_kurtosis_bias_corrected"]] - value[["excess_kurtosis_referent_sample_sd"]]), 0.1)
  expect_gt(abs(value[["mace_equal_batch_average"]] - value[["mace_pooled_rows"]]), 0.01)
})
