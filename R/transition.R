#' Transitions, change, and velocity
#'
#' Distinguishes the longitudinal quantities of the velocity design memo:
#' observed velocity, expected velocity, velocity centile, innovation Z,
#' and change Z.
#'
#' * `innovation_z` is the history-conditioned score
#'   \eqn{(Z_2 - m)/s}; with one previous visit it equals the Z-gain
#'   \eqn{(Z_2 - rZ_1)/\sqrt{1-r^2}}.
#' * `velocity_centile` is the conditional centile of the end value given
#'   the history, i.e. `pnorm(innovation_z)` after the copula map.
#' * `change_z` is the unconditional standardised difference
#'   \eqn{(Z_2 - Z_1)/\sqrt{2(1-r)}} and `change_centile = pnorm(change_z)`.
#'
#' When the process for an outcome is not identified
#' (`dynamic$processes[[outcome]]$identified` is `FALSE`), every
#' history-conditioned quantity is `NA` and `support` is `"unidentified"`.
#' Duplicate visit times within a subject give a warning and `NA`
#' velocities for the zero-length transition.
#'
#' @param dynamic A [norm_dynamics] object.
#' @param data Subject visits.
#' @param id Subject identifier.
#' @param time Time variable. Every transition is conditioned on the
#'   subject's full prior history.
#' @return A `norm_transition` tibble, one row per consecutive visit pair
#'   and outcome (zero rows, full column set, when no subject has two
#'   visits).
#' @export
norm_transition <- function(dynamic, data, id, time) {
  data <- tibble::as_tibble(data)
  id_vec <- pull_column(data, rlang::enquo(id))
  time_vec <- pull_column(data, rlang::enquo(time))
  dup <- stats::ave(time_vec, id_vec, FUN = function(t) duplicated(t))
  if (any(dup > 0, na.rm = TRUE)) {
    cli::cli_warn(
      "{sum(dup > 0, na.rm = TRUE)} duplicate visit time{?s} within subject; velocities for zero-length transitions are NA."
    )
  }
  rows <- lapply(dynamic$outcomes, function(nm) {
    transition_outcome(dynamic, data, nm, id_vec, time_vec)
  })
  rows <- Filter(Negate(is.null), rows)
  out <- if (length(rows)) dplyr_bind(rows) else transition_template()
  structure(out, class = c("norm_transition", class(out)))
}

# Typed zero-row table for the case without any transition.
transition_template <- function() {
  tibble::tibble(
    .id = character(), .outcome = character(), .from_time = numeric(),
    .to_time = numeric(), .dt = numeric(), start_value = numeric(),
    end_value = numeric(), start_centile = numeric(), end_centile = numeric(),
    observed_change = numeric(), expected_change = numeric(),
    change_centile = numeric(), change_z = numeric(),
    observed_velocity = numeric(), expected_velocity = numeric(),
    velocity_lower = numeric(), velocity_upper = numeric(),
    velocity_centile = numeric(), innovation_z = numeric(),
    history_n = integer(), predictive_sd = numeric(), measurement_sd = numeric(),
    support = character(), calibrated = logical()
  )
}

transition_outcome <- function(dynamic, data, outcome, id_vec, time_vec) {
  if (!outcome %in% names(data)) {
    return(NULL)
  }
  pr <- dynamic$processes[[outcome]]
  identified <- isTRUE(pr$identified)
  marg <- predict_dists(dynamic$reference, data, uncertainty = "conditional")[[outcome]]
  z <- dist_z(marg, data[[outcome]])
  sd_marg <- sqrt(variance(marg))
  meas_sd <- if (identified) pr$psi$sigma_e else NA_real_
  pieces <- lapply(unique(id_vec), function(s) {
    idx <- which(id_vec == s)
    idx <- idx[order(time_vec[idx])]
    if (length(idx) < 2L) {
      return(NULL)
    }
    lapply(seq_along(idx)[-1], function(j) {
      i1 <- idx[[j - 1L]]
      i2 <- idx[[j]]
      hist <- idx[seq_len(j - 1L)]
      y1 <- data[[outcome]][i1]
      y2 <- data[[outcome]][i2]
      dt <- time_vec[i2] - time_vec[i1]
      d_end <- marg[i2]
      exp_end <- NA_real_
      v_lo <- v_hi <- NA_real_
      change_z <- innov <- vel_u <- NA_real_
      dt_ok <- is.finite(dt) && dt > 0
      div <- if (dt_ok) dt else NA_real_
      if (identified && dt_ok) {
        cond <- condition_history(z[hist], time_vec[hist], time_vec[i2], pr)
        r12 <- process_kernel(abs(dt), pr$psi, pr$process)
        if (is.finite(cond$s)) {
          d_cond <- dist_conditioned(d_end, cond$m, cond$s)
          innov <- (z[i2] - cond$m) / cond$s
          vel_u <- stats::pnorm(innov)
          exp_end <- dist_quantile(d_cond, 0.5)
          v_lo <- dist_quantile(d_cond, 0.05)
          v_hi <- dist_quantile(d_cond, 0.95)
        }
        change_z <- (z[i2] - z[i1]) / sqrt(pmax(2 * (1 - r12), 1e-8))
      }
      support <- if (!identified) {
        "unidentified"
      } else {
        classify_temporal_support(dynamic, time_vec[i1], time_vec[i2], length(hist))
      }
      tibble::tibble(
        .id = s,
        .outcome = outcome,
        .from_time = time_vec[i1],
        .to_time = time_vec[i2],
        .dt = dt,
        start_value = y1,
        end_value = y2,
        start_centile = exp(dist_eval(marg[i1], log_tail, y1)),
        end_centile = exp(dist_eval(d_end, log_tail, y2)),
        observed_change = y2 - y1,
        expected_change = exp_end - y1,
        change_centile = stats::pnorm(change_z),
        change_z = change_z,
        observed_velocity = (y2 - y1) / div,
        expected_velocity = (exp_end - y1) / div,
        velocity_lower = (v_lo - y1) / div,
        velocity_upper = (v_hi - y1) / div,
        velocity_centile = vel_u,
        innovation_z = innov,
        history_n = length(hist),
        predictive_sd = sd_marg[[i2]],
        measurement_sd = meas_sd,
        support = support,
        calibrated = identified && isTRUE(dynamic$identifiability$measurement)
      )
    })
  })
  pieces <- Filter(Negate(is.null), unlist(pieces, recursive = FALSE))
  if (!length(pieces)) {
    return(NULL)
  }
  dplyr_bind(pieces)
}

#' @export
print.norm_transition <- function(x, ...) {
  cli::cli_text("{.cls norm_transition} {nrow(x)} transition{?s}")
  NextMethod()
}
