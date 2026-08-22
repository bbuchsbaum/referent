simple_spec <- function(family = "gaussian", scale = FALSE) {
  loc <- ~ s(age, k = 5) + sex
  sc <- if (isTRUE(scale)) ~ s(age, k = 4) else ~ 1
  if (identical(family, "shash")) {
    norm_spec(
      family = norm_shash(),
      location = loc,
      scale = sc,
      skew = ~1,
      tail = ~1
    )
  } else {
    norm_spec(family = norm_gaussian(), location = loc, scale = sc)
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
  d <- norm_dist("shash", location = mu, scale = sigma, skew = skew, tail = tail)
  y <- as.numeric(quantile(d, stats::runif(n)))
  data.frame(age, sex, site, y)
}

simulate_two_visit <- function(n_id, r, lag = 2, sigma_e = 0, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  r <- max(min(r, 0.98), -0.98)
  id <- rep(seq_len(n_id), each = 2L)
  age1 <- stats::runif(n_id, 30, 70)
  age <- as.numeric(rbind(age1, age1 + lag))
  sex <- factor(rep(sample(c("F", "M"), n_id, replace = TRUE), each = 2L))
  site <- factor(rep("A", 2L * n_id))
  z1 <- stats::rnorm(n_id)
  z2 <- r * z1 + sqrt(1 - r^2) * stats::rnorm(n_id)
  if (sigma_e > 0) {
    z1 <- z1 + stats::rnorm(n_id, sd = sigma_e)
    z2 <- z2 + stats::rnorm(n_id, sd = sigma_e)
  }
  z <- as.numeric(rbind(z1, z2))
  mu <- 10 + 0.08 * (age - 50)
  y <- mu + 1.3 * z
  data.frame(
    participant_id = id,
    age = age,
    sex = sex,
    site = site,
    y = y,
    visit = rep(1:2, n_id)
  )
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
