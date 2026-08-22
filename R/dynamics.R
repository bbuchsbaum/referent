#' Longitudinal dependence process
#'
#' Attaches a Gaussian-copula process on calibrated normal scores
#' \eqn{Z=\Phi^{-1}(F(y\mid x))}. The default kernel is stable rank plus
#' a Matern-3/2 process plus a measurement nugget, which is positive
#' semidefinite for every irregular visit schedule.
#'
#' @param reference A [norm_fit].
#' @param data Longitudinal reference data.
#' @param id Subject identifier.
#' @param time Time variable (typically age).
#' @param process A process constructor such as [norm_matern32()].
#' @param crossfit Number of subject-level folds used to obtain Z scores.
#' @return An object of class `norm_dynamics`.
#' @export
norm_dynamics <- function(reference,
                          data,
                          id,
                          time,
                          process = norm_matern32(),
                          crossfit = 5) {
  data <- tibble::as_tibble(data)
  id_vec <- pull_column(data, rlang::enquo(id))
  time_vec <- pull_column(data, rlang::enquo(time))
  time_name <- tryCatch(rlang::as_name(rlang::enquo(time)), error = function(e) "time")
  scores <- predict(reference, newdata = data, type = "scores",
                    uncertainty = "conditional", allow_extrapolation = TRUE)
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
      components = c(
        stable_subject_variance = process$stable_rank,
        dynamic_process_variance = ident$change,
        measurement_variance = ident$measurement
      )
    ),
    class = "norm_dynamics"
  )
}

#' Matern-3/2 longitudinal process
#'
#' @param stable_rank Include a subject-level intercept.
#' @param measurement Optional formula for measurement-noise covariates
#'   (currently recorded, not fully expanded).
#' @param ell Initial length-scale.
#' @export
norm_matern32 <- function(stable_rank = TRUE,
                          measurement = ~1,
                          ell = 5) {
  structure(
    list(
      name = "matern32",
      stable_rank = isTRUE(stable_rank),
      measurement = measurement,
      ell = ell
    ),
    class = c("norm_process_matern32", "norm_process")
  )
}

matern32_kernel <- function(d, ell) {
  a <- sqrt(3) * abs(d) / ell
  (1 + a) * exp(-a)
}

process_covariance <- function(times, psi, process) {
  n <- length(times)
  d <- abs(outer(times, times, `-`))
  k <- matrix(0, n, n)
  if (isTRUE(process$stable_rank)) {
    k <- k + psi$tau_b^2
  }
  if (identical(process$name, "matern32")) {
    k <- k + psi$tau_g^2 * matern32_kernel(d, psi$ell)
  }
  diag(k) <- diag(k) + psi$sigma_e^2
  k
}

process_correlation <- function(times, psi, process) {
  k <- process_covariance(times, psi, process)
  s <- sqrt(pmax(diag(k), .Machine$double.eps))
  k / outer(s, s)
}

fit_process <- function(process, z, id, time, ident = NULL) {
  ok <- is.finite(z) & is.finite(time)
  z <- z[ok]
  id <- as.character(id[ok])
  time <- as.numeric(time[ok])
  ident <- ident %||% process_identifiability(id, time)
  start <- c(tau_b = 0.6, tau_g = 0.4, sigma_e = 0.3, ell = process$ell %||% 5)
  if (!isTRUE(ident$change)) {
    return(list(
      process = process,
      psi = process_par(log(start), process),
      nll = NA_real_,
      identified = FALSE,
      reason = ident$reason
    ))
  }
  nll <- function(par) {
    psi <- process_par(par, process)
    ll <- 0
    for (s in unique(id)) {
      idx <- which(id == s)
      if (length(idx) < 2L) {
        next
      }
      r <- process_correlation(time[idx], psi, process)
      r <- r + diag(1e-6, length(idx))
      ll <- ll + gaussian_loglik(z[idx], r)
    }
    -ll
  }
  opt <- tryCatch(
    stats::optim(log(pmax(start, 1e-3)), nll, method = "BFGS"),
    error = function(e) NULL
  )
  psi <- if (is.null(opt)) {
    process_par(log(start), process)
  } else {
    process_par(opt$par, process)
  }
  list(
    process = process,
    psi = psi,
    nll = if (is.null(opt)) NA_real_ else opt$value,
    identified = TRUE,
    reason = NULL
  )
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
  lag_rng <- if (length(lags)) range(lags) else c(0, 0)
  change_ok <- n_repeat >= 5L && length(lags) >= 8L
  meas_ok <- any(n_per >= 3L) ||
    (length(lags) >= 8L && stats::median(lags) > 0 &&
       min(lags) <= 0.25 * stats::median(lags))
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
    reason = reason
  )
}

process_par <- function(par, process) {
  par <- exp(as.numeric(par))
  list(tau_b = par[[1]], tau_g = par[[2]], sigma_e = par[[3]], ell = par[[4]])
}

gaussian_loglik <- function(z, r) {
  ev <- eigen(r, symmetric = TRUE)
  ev$values <- pmax(ev$values, 1e-8)
  inv <- ev$vectors %*% diag(1 / ev$values) %*% t(ev$vectors)
  ldet <- sum(log(ev$values))
  -0.5 * (length(z) * log(2 * pi) + ldet + drop(t(z) %*% inv %*% z))
}


condition_z <- function(z_hist, t_hist, t_new, psi, process) {
  t_all <- c(t_hist, t_new)
  r <- process_correlation(t_all, psi, process)
  n_h <- length(t_hist)
  if (!n_h) {
    return(list(m = 0, s = 1, r = 0))
  }
  r_hh <- r[seq_len(n_h), seq_len(n_h), drop = FALSE]
  r_hs <- r[seq_len(n_h), n_h + 1L, drop = FALSE]
  r_hh <- r_hh + diag(1e-8, n_h)
  inv <- tryCatch(solve(r_hh), error = function(e) NULL)
  if (is.null(inv)) {
    return(list(m = 0, s = 1, r = r[1, 2]))
  }
  m <- drop(t(r_hs) %*% inv %*% z_hist)
  s2 <- pmax(1 - drop(t(r_hs) %*% inv %*% r_hs), 1e-8)
  list(m = m, s = sqrt(s2), r = if (n_h == 1L) r[1, 2] else NA_real_)
}

classify_temporal_support <- function(dyn, t_from, t_to, history_n) {
  status <- rep("in", length(t_from))
  lag <- abs(t_to - t_from)
  status[t_from < dyn$time_range[1] | t_to > dyn$time_range[2]] <- "unsupported_age"
  if (diff(dyn$lag_range) > 0) {
    status[lag > dyn$lag_range[2] * 1.05] <- "unsupported_lag"
  }
  status[history_n < 1] <- "insufficient_history"
  status
}

#' @export
print.norm_dynamics <- function(x, ...) {
  cli::cli_text("{.cls norm_dynamics} {x$process$name} for {length(x$outcomes)} outcome{?s}")
  cli::cli_text(
    "subjects: {x$n_subject}; time range [{signif(x$time_range[1], 4)}, {signif(x$time_range[2], 4)}]"
  )
  ident <- x$identifiability
  if (!is.null(ident) && !isTRUE(ident$change)) {
    cli::cli_alert_warning("change not identified: {ident$reason %||% 'insufficient repeats'}")
  } else if (!is.null(ident) && !isTRUE(ident$measurement)) {
    cli::cli_alert_warning("measurement not separated: {ident$reason %||% ''}")
  }
  invisible(x)
}
