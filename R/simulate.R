#' Simulate reference-model data for law tests
#'
#' @param n Number of rows.
#' @param kind Generator name.
#' @param seed Optional seed.
#' @return A data frame with covariates and one or more outcomes.
#' @export
norm_simulate <- function(n = 400,
                          kind = c(
                            "gaussian_location",
                            "gaussian_scale",
                            "shash",
                            "panel",
                            "shift",
                            "longitudinal",
                            "ordinal"
                          ),
                          seed = NULL) {
  kind <- match.arg(kind)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  age <- stats::runif(n, 20, 80)
  sex <- factor(sample(c("F", "M"), n, replace = TRUE))
  site <- factor(sample(LETTERS[1:4], n, replace = TRUE))
  mu <- 10 + 0.08 * (age - 50) - 0.001 * (age - 50)^2 + 0.4 * (sex == "M")
  sigma <- 1.2 + 0.02 * pmax(age - 40, 0)
  switch(
    kind,
    gaussian_location = {
      y <- mu + stats::rnorm(n, sd = 1.3)
      data.frame(age, sex, site, y, marker_01 = y, marker_02 = y + stats::rnorm(n, sd = 0.4))
    },
    gaussian_scale = {
      y <- mu + stats::rnorm(n, sd = sigma)
      data.frame(age, sex, site, y, marker_01 = y)
    },
    shash = {
      eps <- 0.6
      delta <- 0.85
      u <- stats::runif(n)
      y <- shash_quantile(u, mu, sigma, eps, delta)
      data.frame(age, sex, site, y, marker_01 = y)
    },
    panel = {
      y1 <- mu + stats::rnorm(n, sd = 1.1)
      y2 <- 0.6 * scale(y1)[, 1] + stats::rnorm(n, sd = 0.9)
      y3 <- -0.2 * (age - 50) / 10 + stats::rnorm(n, sd = 1)
      data.frame(age, sex, site, marker_01 = y1, marker_02 = as.numeric(y2), marker_03 = y3)
    },
    shift = {
      delta_mu <- c(A = 0, B = 1.2, C = -0.8, D = 0.3)
      delta_log_s <- c(A = 0, B = 0.15, C = 0, D = -0.1)
      y <- mu + delta_mu[as.character(site)] +
        stats::rnorm(n, sd = sigma * exp(delta_log_s[as.character(site)]))
      data.frame(age, sex, site, y, marker_01 = y)
    },
    longitudinal = {
      id <- rep(seq_len(ceiling(n / 3)), each = 3, length.out = n)
      t <- ave(age, id, FUN = function(z) sort(z))
      if (length(t) != n) {
        t <- age
      }
      z0 <- stats::rnorm(length(unique(id)))
      names(z0) <- unique(id)
      rho <- exp(-abs(0.15))
      z <- z0[as.character(id)] + stats::rnorm(n, sd = sqrt(1 - 0.6))
      z <- as.numeric(scale(z)[, 1])
      y <- mu + 1.2 * z
      data.frame(
        participant_id = id,
        age = t,
        sex,
        site,
        y,
        marker_01 = y,
        visit = ave(seq_len(n), id, FUN = seq_along)
      )
    },
    ordinal = {
      z <- mu / 3 + stats::rnorm(n)
      y <- cut(z, breaks = c(-Inf, -0.5, 0.5, Inf), labels = c("low", "mid", "high"))
      data.frame(age, sex, site, y, marker_01 = as.numeric(y))
    }
  )
}
