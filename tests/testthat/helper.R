simple_spec <- function(family = "gaussian", scale = FALSE) {
  loc <- ~ s(age, k = 5) + sex
  sc <- if (isTRUE(scale)) ~ s(age, k = 4) else ~ 1
  if (identical(family, "shash")) {
    ref_spec(
      family = ref_shash(),
      location = loc,
      scale = sc,
      skew = ~1,
      tail = ~1
    )
  } else {
    ref_spec(family = ref_gaussian(), location = loc, scale = sc)
  }
}

expect_near <- function(x, y, tol = 1e-6) {
  expect_true(all(abs(x - y) < tol, na.rm = TRUE))
}

simulate_skewed <- function(n, seed = NULL, skew = 1.4, tail = 0.7) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  age <- stats::runif(n, 20, 80)
  sex <- factor(sample(c("F", "M"), n, replace = TRUE))
  site <- factor(sample(LETTERS[1:4], n, replace = TRUE))
  mu <- 10 + 0.08 * (age - 50) - 0.001 * (age - 50)^2 + 0.4 * (sex == "M")
  sigma <- 1.2 + 0.02 * pmax(age - 40, 0)
  y <- shash_quantile(stats::runif(n), mu, sigma, skew, tail)
  data.frame(age, sex, site, y)
}

# Two-visit longitudinal data with a fixed lag and correlation `r` between
# the visits' normal scores (stable rank + nugget, no Matern component).
simulate_two_visit <- function(n_id, r, lag = 2, seed = NULL) {
  ref_simulate(
    2L * n_id, kind = "longitudinal", lag = lag,
    tau_b = sqrt(r), tau_g = 0, sigma_e = sqrt(1 - r), seed = seed
  )
}

# Vectorised CDF and log density of a distribution vector.
dist_cdf <- function(d, q) {
  exp(dist_eval(d, log_tail, q))
}

dist_log_density <- function(d, y) {
  dist_eval(d, log_dens, y)
}

worm_rmse <- function(z) {
  z <- sort(z[is.finite(z)])
  expected <- stats::qnorm(stats::ppoints(length(z)))
  sqrt(mean((z - expected)^2))
}

pit_ks <- function(u) {
  u <- u[is.finite(u) & u > 0 & u < 1]
  as.numeric(stats::ks.test(u, "punif")$statistic)
}

false_positive_rate <- function(z, threshold = stats::qnorm(0.975)) {
  z <- z[is.finite(z)]
  mean(abs(z) > threshold)
}

# Everything a print method emits: cli messages and standard output.
print_text <- function(x) {
  body <- utils::capture.output(msg <- cli::cli_fmt(print(x)))
  paste(c(msg, body), collapse = "\n")
}

# Multi-outcome fits used by the prediction and performance tests.
perf_data <- function() {
  ref <- ref_simulate(300, sites = 4, scale = "age", seed = 11)
  new <- ref_simulate(80, sites = 4, scale = "age", seed = 12)
  list(ref = ref, new = new)
}

perf_specs <- function() {
  list(
    constant = ref_spec(ref_gaussian(), ~ s(age, k = 6) + sex + s(site, bs = "re")),
    gaulss = ref_spec(ref_gaussian(), ~ s(age, k = 6) + sex + s(site, bs = "re"),
                       scale = ~ s(age, k = 4)),
    shash = ref_spec(ref_shash(), ~ s(age, k = 6) + sex, scale = ~ s(age, k = 4))
  )
}


# ggplot helpers shared by the plot tests.
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
