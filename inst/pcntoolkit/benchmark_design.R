# Preregistered synthetic truths for the PCNtoolkit utility benchmark.

pcn_benchmark_registry <- function(root = NULL) {
  installed <- system.file("pcntoolkit", "benchmark_scenarios.csv", package = "referent")
  path <- if (is.null(root) && nzchar(installed)) {
    installed
  } else {
    root <- if (is.null(root)) "." else root
    source_path <- file.path(root, "tools", "pcntoolkit", "benchmark_scenarios.csv")
    if (file.exists(source_path)) source_path else
      file.path(root, "inst", "pcntoolkit", "benchmark_scenarios.csv")
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

pcn_simulate_benchmark <- function(scenario, seed, root = NULL) {
  registry <- pcn_benchmark_registry(root)
  design <- registry[registry$scenario == scenario, , drop = FALSE]
  if (nrow(design) != 1L) stop("unknown benchmark scenario: ", scenario, call. = FALSE)
  n_part <- unlist(design[c("n_train", "n_validation", "n_test")], use.names = FALSE)
  n <- sum(n_part)
  set.seed(seed)
  split <- rep(c("train", "validation", "test"), n_part)
  x <- if (scenario == "covariate_shift") {
    c(stats::runif(n_part[[1L]], -2, 0.75),
      stats::runif(n_part[[2L]] + n_part[[3L]], 0.25, 2))
  } else {
    stats::runif(n, -2, 2)
  }
  site_levels <- paste0("site-", 1:5)
  site <- if (scenario == "unequal_sites") {
    sample(site_levels, n, replace = TRUE, prob = c(0.6, 0.2, 0.1, 0.07, 0.03))
  } else if (scenario == "unseen_site") {
    c(sample(site_levels[1:4], n_part[[1L]] + n_part[[2L]], replace = TRUE),
      sample(site_levels, n_part[[3L]], replace = TRUE, prob = c(rep(0.125, 4), 0.5)))
  } else {
    sample(site_levels, n, replace = TRUE)
  }
  site_index <- match(site, site_levels)
  mu <- 0.3 + 0.9 * x
  sigma <- rep(0.75, n)
  epsilon <- rep(0, n)
  delta <- rep(1, n)

  if (scenario %in% c("nonlinear_mean", "small_n", "large_n", "skew_heavy")) {
    mu <- 0.4 + sin(1.35 * x) + 0.15 * x^2
  }
  if (scenario %in% c("heteroskedastic", "small_n", "large_n", "skew_heavy")) {
    sigma <- exp(-0.35 + 0.22 * x)
  }
  if (scenario %in% c("skew_heavy", "small_n", "large_n")) {
    epsilon <- rep(0.9, n)
    delta <- rep(0.7, n)
  }
  if (scenario %in% c("unequal_sites", "unseen_site")) {
    mu <- mu + c(-0.8, -0.25, 0.2, 0.65, 1.1)[site_index]
    sigma <- sigma * c(0.7, 0.85, 1, 1.15, 1.35)[site_index]
  }
  if (scenario == "site_slopes") {
    mu <- mu + c(-0.4, -0.2, 0, 0.2, 0.4)[site_index] * x
  }
  latent <- stats::rnorm(n)
  y <- mu + sigma * delta * sinh((asinh(latent) + epsilon) / delta)
  data.frame(
    scenario = scenario,
    row_id = sprintf("%s-%s-%05d", scenario, split, ave(seq_len(n), split, FUN = seq_along)),
    split = split,
    cluster = if (design$site_structure == "none") seq_len(n) else site,
    x = x,
    site = factor(site, levels = site_levels),
    y = y,
    truth_mu = mu,
    truth_sigma = sigma,
    truth_epsilon = epsilon,
    truth_delta = delta,
    stringsAsFactors = FALSE
  )
}
