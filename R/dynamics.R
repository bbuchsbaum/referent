#' Longitudinal dependence process
#'
#' Attaches a Gaussian-copula process on marginal normal scores
#' \eqn{Z=\Phi^{-1}(F(y\mid x))}. The default kernel is stable rank plus
#' a Matern-3/2 process plus a measurement nugget, which is positive
#' semidefinite for every irregular visit schedule.
#'
#' Because the marginal scores have unit variance, the three variance
#' components are estimated on the simplex
#' \eqn{\tau_b^2 + \tau_g^2 + \sigma_e^2 = 1} and the Matern length-scale
#' is bounded by the observed lags. When all within-subject lags are
#' (nearly) identical the length-scale is not identified and the model is
#' automatically reduced to the stable-rank-plus-nugget kernel
#' (`process = "stable"`), which has a single parameter: the correlation
#' at the common lag. Optimisation is multi-start L-BFGS-B; when no start
#' converges the process is reported as `identified = FALSE` and every
#' history-conditioned quantity downstream is `NA`.
#'
#' @param reference A [ref_fit].
#' @param data Longitudinal reference data.
#' @param id Subject identifier.
#' @param time Time variable (typically age).
#' @param process A kernel from [ref_process()].
#' @param crossfit Number of subject-level folds used to obtain
#'   out-of-fold Z scores via [ref_crossfit()] (whole subjects stay in one
#'   fold). `0` or `NULL` uses in-sample scores from `reference`.
#' @param uncertainty Marginal predictive estimand used to create every Z
#'   score and every downstream transition or forecast. Defaults to
#'   `"total"`, matching [predict.ref_fit()].
#' @param n_draw Coefficient draws for total marginal uncertainty. The value
#'   is stored and reused downstream.
#' @return An object of class `ref_dynamics` with, per outcome, a fitted
#'   process (`$processes`), and a `components` table of normalised
#'   variance fractions with approximate standard errors.
#' @export
ref_dynamics <- function(reference,
                          data,
                          id,
                          time,
                          process = ref_process(),
                          crossfit = 5,
                          uncertainty = c("total", "conditional"),
                          n_draw = NULL) {
  uncertainty <- match.arg(uncertainty)
  data <- tibble::as_tibble(data)
  id_quo <- rlang::enquo(id)
  id_vec <- pull_column(data, id_quo)
  id_name <- as_col_name(id_quo, "id")
  time_quo <- rlang::enquo(time)
  time_vec <- pull_column(data, time_quo)
  time_name <- as_col_name(time_quo, "time")
  if (!inherits(process, "ref_process")) {
    cli::cli_abort("{.arg process} must be created by {.fn ref_process}.")
  }
  if (!is.null(reference$calibration)) {
    cli::cli_abort(
      "Longitudinal conditioning does not yet support a PIT-calibrated reference; use the uncalibrated fit and validate its marginal scores independently."
    )
  }
  crossfit <- as.integer(crossfit %||% 0L)
  if (crossfit >= 2L && !is.null(reference$adaptation)) {
    cli::cli_abort(
      "Site adaptation cannot be reproduced inside subject-level cross-fitting; use an unadapted reference or set {.arg crossfit = 0} for a labelled in-sample process fit."
    )
  }
  scores <- if (crossfit >= 2L) {
    # Subject-level out-of-fold Z; `.row` indexes rows of `data`.
    data$.dyn_id <- as.character(id_vec)
    ref_crossfit(reference$spec, data = data, outcomes = reference$outcomes,
                  folds = crossfit, cluster = !!rlang::sym(".dyn_id"),
                  uncertainty = uncertainty, n_draw = n_draw)
  } else {
    predict(reference, newdata = data, type = "scores",
            uncertainty = uncertainty, n_draw = n_draw,
            allow_extrapolation = TRUE)
  }
  ident <- process_identifiability(id_vec, time_vec)
  processes <- lapply(reference$outcomes, function(nm) {
    sc <- scores[scores$.outcome == nm, , drop = FALSE]
    fit_process(process, z = sc$z, id = id_vec[sc$.row], time = time_vec[sc$.row],
                ident = ident)
  })
  names(processes) <- reference$outcomes
  structure(
    list(
      reference = reference,
      process = process,
      processes = processes,
      outcomes = reference$outcomes,
      time_name = time_name,
      time_range = range(time_vec, na.rm = TRUE),
      lag_range = ident$lag_range,
      n_subject = length(unique(id_vec)),
      identifiability = ident,
      z_source = if (crossfit >= 2L) "out_of_fold" else "in_sample",
      uncertainty = uncertainty,
      n_draw = if (identical(uncertainty, "total")) {
        as.integer(n_draw %||% reference$spec$control$n_draw %||% 200L)
      } else {
        NA_integer_
      },
      kernel_uncertainty = "plug_in",
      data_hash = digest_data(data[, unique(c(id_name, time_name, reference$covariates,
                                               reference$outcomes)), drop = FALSE]),
      data_hash_version = 2L,
      data_columns = unique(c(id_name, time_name, reference$covariates,
                              reference$outcomes)),
      data_n = nrow(data),
      components = process_components(processes)
    ),
    class = "ref_dynamics"
  )
}

#' Longitudinal process kernel
#'
#' `ref_process()` builds the dependence kernel for [ref_dynamics()].
#' `"matern32"` is stable rank + Matern-3/2 + nugget; `"stable"` is stable
#' rank + nugget only (a single correlation shared by every positive lag),
#' the right model when visits share one fixed lag. `ell = Inf` also
#' selects the stable kernel.
#'
#' @param kernel `"matern32"` or `"stable"`.
#' @param stable_rank Include a subject-level intercept.
#' @param ell Initial length-scale (`Inf` selects the stable kernel).
#' @return A `ref_process` object.
#' @export
ref_process <- function(kernel = c("matern32", "stable"),
                         stable_rank = TRUE,
                         ell = 5) {
  kernel <- match.arg(kernel)
  if (is.infinite(ell)) {
    kernel <- "stable"
  }
  structure(
    list(
      name = kernel,
      stable_rank = isTRUE(stable_rank),
      ell = if (identical(kernel, "stable")) Inf else ell
    ),
    class = "ref_process"
  )
}

matern32_kernel <- function(d, ell) {
  if (is.infinite(ell)) {
    return(array(1, dim = dim(d) %||% length(d)))
  }
  a <- sqrt(3) * abs(d) / ell
  (1 + a) * exp(-a)
}

# Process correlation at a lag (any array shape): the shared variance
# fraction at that lag, plus the nugget at lag zero.
process_kernel <- function(lag, psi, process) {
  num <- 0
  if (isTRUE(process$stable_rank)) {
    num <- num + psi$tau_b^2
  }
  if (identical(process$name, "matern32")) {
    num <- num + psi$tau_g^2 * matern32_kernel(lag, psi$ell)
  }
  tot <- (if (isTRUE(process$stable_rank)) psi$tau_b^2 else 0) + psi$tau_g^2 + psi$sigma_e^2
  num / tot + (lag == 0) * psi$sigma_e^2 / tot
}

process_correlation <- function(times, psi, process) {
  process_kernel(abs(outer(times, times, `-`)), psi, process)
}

# ---- estimator ------------------------------------------------------------

# theta -> psi. Matern: theta = (logit_b, logit_g, log_ell) with
# (tau_b^2, tau_g^2, sigma_e^2) = softmax(logit_b, logit_g, 0).
# Stable: theta = (logit_b) with tau_b^2 = plogis(logit_b).
process_par <- function(theta, process) {
  theta <- as.numeric(theta)
  if (identical(process$name, "stable")) {
    v_b <- stats::plogis(theta[[1]])
    return(list(tau_b = sqrt(v_b), tau_g = 0, sigma_e = sqrt(1 - v_b), ell = Inf))
  }
  w <- exp(c(theta[1:2], 0) - max(theta[1:2], 0))
  w <- w / sum(w)
  list(tau_b = sqrt(w[[1]]), tau_g = sqrt(w[[2]]), sigma_e = sqrt(w[[3]]),
       ell = exp(theta[[3]]))
}

# Group subjects by visit count once: each group holds an m x n matrix of
# scores and an m x n x n array of lags, so the likelihood is evaluated
# with a Cholesky factorisation vectorised across the m subjects.
process_data <- function(z, id, time) {
  ok <- is.finite(z) & is.finite(time)
  z <- z[ok]
  id <- as.character(id[ok])
  time <- as.numeric(time[ok])
  groups <- split(seq_along(z), id)
  # duplicate visit times would make the correlation matrix singular
  groups <- lapply(groups, function(i) i[!duplicated(time[i])])
  groups <- groups[lengths(groups) >= 2L]
  by_n <- split(groups, lengths(groups))
  batches <- lapply(by_n, function(gs) {
    n <- length(gs[[1]])
    zm <- t(vapply(gs, function(i) z[i], numeric(n)))
    tm <- t(vapply(gs, function(i) time[i], numeric(n)))
    lag <- array(0, c(length(gs), n, n))
    for (a in seq_len(n)) {
      for (b in seq_len(n)) {
        lag[, a, b] <- abs(tm[, a] - tm[, b])
      }
    }
    list(z = zm, lag = lag, n = n, m = length(gs))
  })
  list(batches = batches, n_subject = length(groups), n_obs = sum(lengths(groups)))
}

process_nll <- function(theta, process, pd) {
  psi <- process_par(theta, process)
  ll <- 0
  for (b in pd$batches) {
    ll <- ll + batched_loglik(b$z, process_kernel(b$lag, psi, process))
  }
  if (!is.finite(ll)) {
    return(1e10)
  }
  -ll
}

# Gaussian log-likelihood of m subjects with n observations each, with
# correlation array r (m x n x n), by a Cholesky factorisation vectorised
# across subjects. Returns -Inf when any subject's matrix is not PD.
batched_loglik <- function(z, r) {
  m <- nrow(z)
  n <- ncol(z)
  l <- array(0, c(m, n, n))
  w <- matrix(0, m, n)
  ldet <- numeric(m)
  for (j in seq_len(n)) {
    djj <- r[, j, j]
    if (j > 1L) {
      djj <- djj - rowSums(matrix(l[, j, seq_len(j - 1L)], m)^2)
    }
    if (any(!is.finite(djj)) || any(djj <= 1e-12)) {
      return(-Inf)
    }
    ljj <- sqrt(djj)
    l[, j, j] <- ljj
    ldet <- ldet + log(ljj)
    if (j < n) {
      for (i in (j + 1L):n) {
        s <- r[, i, j]
        if (j > 1L) {
          s <- s - rowSums(matrix(l[, i, seq_len(j - 1L)], m) * matrix(l[, j, seq_len(j - 1L)], m))
        }
        l[, i, j] <- s / ljj
      }
    }
    wj <- z[, j]
    if (j > 1L) {
      wj <- wj - rowSums(matrix(l[, j, seq_len(j - 1L)], m) * w[, seq_len(j - 1L), drop = FALSE])
    }
    w[, j] <- wj / ljj
  }
  sum(-0.5 * (n * log(2 * pi) + 2 * ldet + rowSums(w^2)))
}

# Lags are "fixed" when their spread is small relative to their size; the
# Matern length-scale is then unidentified and the stable kernel is used.
lags_are_fixed <- function(lag_range, tol = 0.05) {
  lag_range[[2]] <= 0 || diff(lag_range) <= tol * stats::median(lag_range)
}

fit_process <- function(process, z, id, time, ident = NULL) {
  ident <- ident %||% process_identifiability(id, time)
  unfit <- function(reason) {
    list(
      process = process, psi = NULL, theta = NULL, se = NULL, vcov = NULL,
      nll = NA_real_, identified = FALSE, ell_identified = FALSE,
      r_median = NA_real_, r_median_se = NA_real_, median_lag = NA_real_,
      reason = reason, components = c(stable = NA_real_, dynamic = NA_real_,
                                      measurement = NA_real_)
    )
  }
  if (!isTRUE(ident$change)) {
    return(unfit(ident$reason))
  }
  pd <- process_data(z, id, time)
  if (pd$n_subject < 5L) {
    return(unfit("too few subjects with finite repeated scores"))
  }
  lag_range <- ident$lag_range
  reduced <- identical(process$name, "matern32") && lags_are_fixed(lag_range)
  if (reduced) {
    process <- ref_process("stable", stable_rank = process$stable_rank)
  }
  stable <- identical(process$name, "stable")
  med_lag <- ident$median_lag
  if (stable) {
    starts <- list(stats::qlogis(0.3), stats::qlogis(0.6), stats::qlogis(0.85))
    lower <- -12
    upper <- 12
  } else {
    min_pos <- max(ident$min_positive_lag, 1e-6)
    ell_lo <- log(min_pos / 2)
    ell_hi <- log(3 * lag_range[[2]])
    ell0 <- min(max(log(med_lag), ell_lo), ell_hi)
    ell1 <- min(max(log(2 * med_lag), ell_lo), ell_hi)
    starts <- list(
      c(0, 0, ell0),
      c(log(0.7 / 0.15), 0, ell0),
      c(log(0.2 / 0.1), log(0.7 / 0.1), ell1)
    )
    lower <- c(-12, -12, ell_lo)
    upper <- c(12, 12, ell_hi)
  }
  fits <- lapply(starts, function(s) {
    tryCatch(
      stats::optim(s, process_nll, process = process, pd = pd,
                   method = "L-BFGS-B", lower = lower, upper = upper,
                   control = list(maxit = 500)),
      error = function(e) NULL
    )
  })
  ok <- vapply(fits, function(f) {
    !is.null(f) && isTRUE(f$convergence == 0) && is.finite(f$value) && f$value < 1e9
  }, logical(1))
  if (!any(ok)) {
    msgs <- unique(unlist(lapply(fits, function(f) f$message)))
    return(unfit(paste0("optimiser did not converge",
                        if (length(msgs)) paste0(" (", paste(msgs, collapse = "; "), ")"))))
  }
  best <- fits[ok][[which.min(vapply(fits[ok], `[[`, numeric(1), "value"))]]
  theta <- best$par
  # An optimum on the box is not an interior maximum: the curvature there
  # says nothing about the sampling variability, so no standard errors
  # are reported and the fit is flagged.
  at_boundary <- any(abs(theta - lower) < 1e-6 | abs(theta - upper) < 1e-6)
  hess <- if (at_boundary) {
    NULL
  } else {
    tryCatch(
      stats::optimHess(theta, process_nll, process = process, pd = pd),
      error = function(e) NULL
    )
  }
  vc <- NULL
  if (!is.null(hess) && all(is.finite(hess))) {
    vc <- tryCatch(solve(hess), error = function(e) NULL)
    if (!is.null(vc) && any(diag(vc) < 0)) {
      vc <- NULL
    }
  }
  se <- if (is.null(vc)) rep(NA_real_, length(theta)) else sqrt(diag(vc))
  psi <- process_par(theta, process)
  r_fun <- function(th) process_kernel(med_lag, process_par(th, process), process)
  r_med <- r_fun(theta)
  r_se <- NA_real_
  if (!is.null(vc)) {
    grad <- numeric_gradient(r_fun, theta)
    r_se <- sqrt(max(drop(t(grad) %*% vc %*% grad), 0))
  }
  nm <- if (stable) "logit_stable" else c("logit_stable", "logit_dynamic", "log_ell")
  names(theta) <- names(se) <- nm
  list(
    process = process,
    psi = psi,
    theta = theta,
    se = se,
    vcov = vc,
    nll = best$value,
    identified = TRUE,
    ell_identified = !stable && !at_boundary,
    at_boundary = at_boundary,
    reduced = reduced,
    r_median = r_med,
    r_median_se = r_se,
    median_lag = med_lag,
    n_subject = pd$n_subject,
    reason = if (at_boundary) {
      "optimum at the parameter boundary; estimates are limits, no standard errors"
    } else if (reduced) {
      "lags do not vary; Matern length-scale not identified, stable kernel used"
    } else {
      NULL
    },
    components = c(stable = psi$tau_b^2, dynamic = psi$tau_g^2, measurement = psi$sigma_e^2)
  )
}

numeric_gradient <- function(f, x, h = 1e-4) {
  vapply(seq_along(x), function(j) {
    e <- replace(numeric(length(x)), j, h)
    (f(x + e) - f(x - e)) / (2 * h)
  }, numeric(1))
}

process_components <- function(processes) {
  dplyr_bind(lapply(names(processes), function(nm) {
    pr <- processes[[nm]]
    tibble::tibble(
      .outcome = nm,
      process = pr$process$name,
      identified = isTRUE(pr$identified),
      at_boundary = isTRUE(pr$at_boundary),
      stable = unname(pr$components[["stable"]]),
      dynamic = unname(pr$components[["dynamic"]]),
      measurement = unname(pr$components[["measurement"]]),
      ell = if (isTRUE(pr$ell_identified)) pr$psi$ell else NA_real_,
      ell_identified = isTRUE(pr$ell_identified),
      median_lag = pr$median_lag,
      r_median_lag = pr$r_median,
      r_median_lag_se = pr$r_median_se
    )
  }))
}

process_identifiability <- function(id, time) {
  id <- as.character(id)
  time <- as.numeric(time)
  n_per <- table(id[is.finite(time)])
  n_repeat <- sum(n_per >= 2L)
  lags <- unlist(lapply(split(time, id), function(t) {
    t <- t[is.finite(t)]
    if (length(t) < 2L) {
      return(numeric())
    }
    as.numeric(stats::dist(t))
  }), use.names = FALSE)
  pos <- lags[lags > 0]
  lag_rng <- if (length(pos)) range(pos) else c(0, 0)
  change_ok <- n_repeat >= 5L && length(lags) >= 8L
  meas_ok <- any(n_per >= 3L) ||
    (length(pos) >= 8L && stats::median(pos) > 0 &&
       min(pos) <= 0.25 * stats::median(pos))
  reason <- NULL
  if (!change_ok) {
    reason <- "too few repeated observations to identify within-person dependence"
  } else if (!meas_ok) {
    reason <- "no short-interval repeats; measurement noise is not separated from rank dynamics"
  }
  list(
    change = change_ok,
    measurement = meas_ok,
    n_repeat_subject = as.integer(n_repeat),
    n_lag = length(lags),
    lag_range = lag_rng,
    median_lag = if (length(pos)) stats::median(pos) else NA_real_,
    min_positive_lag = if (length(pos)) min(pos) else NA_real_,
    fixed_lag = length(pos) > 0 && lags_are_fixed(lag_rng),
    reason = reason
  )
}

# ---- conditioning and support --------------------------------------------

# Conditional mean and sd of the normal score at each `t_new` given the
# history scores `z_hist` at `t_hist`, under a fitted process (a list with
# `psi` and `process`). Without history, or without a fitted process, the
# marginal N(0, 1) is returned.
condition_history <- function(z_hist, t_hist, t_new, fitted) {
  n_h <- length(t_hist)
  n_new <- length(t_new)
  if (!n_h || is.null(fitted$psi)) {
    return(list(m = rep(0, n_new), s = rep(1, n_new)))
  }
  psi <- fitted$psi
  process <- fitted$process
  r_hh <- process_correlation(t_hist, psi, process) + diag(1e-8, n_h)
  inv <- tryCatch(solve(r_hh), error = function(e) NULL)
  if (is.null(inv)) {
    return(list(m = rep(NA_real_, n_new), s = rep(NA_real_, n_new)))
  }
  r_hs <- process_kernel(abs(outer(t_hist, t_new, `-`)), psi, process)
  w <- inv %*% r_hs
  list(
    m = drop(crossprod(w, z_hist)),
    s = sqrt(pmax(1 - colSums(r_hs * w), 1e-8))
  )
}

# Lag support is always checked against the reference lag range: a lag-10
# forecast under a fixed-lag-2 reference is extrapolated, not "in".
classify_temporal_support <- function(dyn, t_from, t_to, history_n) {
  status <- rep("in", length(t_from))
  lag <- abs(t_to - t_from)
  rng <- dyn$lag_range
  tol <- 0.05 * max(rng[[2]], .Machine$double.eps)
  status[lag > rng[[2]] + tol | lag < rng[[1]] - tol] <- "extrapolated_lag"
  status[t_from < dyn$time_range[1] | t_to > dyn$time_range[2]] <- "unsupported_age"
  status[history_n < 1] <- "insufficient_history"
  status
}

#' @export
print.ref_dynamics <- function(x, ...) {
  cli::cli_text(
    "{.cls ref_dynamics} {length(x$outcomes)} outcome{?s}; requested kernel: {x$process$name}; Z: {x$z_source}, {x$uncertainty} uncertainty"
  )
  cli::cli_text("kernel uncertainty: plug-in (not propagated)")
  cli::cli_text(
    "subjects: {x$n_subject}; time range [{signif(x$time_range[1], 4)}, {signif(x$time_range[2], 4)}]; lag range [{signif(x$lag_range[1], 3)}, {signif(x$lag_range[2], 3)}]"
  )
  ident <- x$identifiability
  if (!isTRUE(ident$change)) {
    cli::cli_alert_warning("change not identified: {ident$reason}")
  } else if (!isTRUE(ident$measurement)) {
    cli::cli_alert_warning("measurement not separated: {ident$reason}")
  }
  print(x$components)
  invisible(x)
}
