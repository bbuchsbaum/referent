#' Transitions, change, and velocity
#'
#' Distinguishes six longitudinal quantities: chart velocity, observed
#' velocity, expected velocity, velocity centile, innovation Z, and
#' change Z. `z_gain` is an interoperability alias for `innovation_z`.
#'
#' @param dynamic A [norm_dynamics] object.
#' @param data Subject visits.
#' @param id Subject identifier.
#' @param time Time variable.
#' @param conditioning `"all"` uses the full prior history.
#' @return A `norm_transition` table.
#' @export
norm_transition <- function(dynamic,
                            data,
                            id,
                            time,
                            conditioning = c("all", "last")) {
  conditioning <- match.arg(conditioning)
  data <- tibble::as_tibble(data)
  id_vec <- pull_column(data, rlang::enquo(id))
  time_vec <- pull_column(data, rlang::enquo(time))
  rows <- lapply(dynamic$outcomes, function(nm) {
    transition_outcome(dynamic, data, nm, id_vec, time_vec, conditioning)
  })
  out <- dplyr_bind(Filter(Negate(is.null), rows))
  structure(out, class = c("norm_transition", class(out)))
}

transition_outcome <- function(dynamic, data, outcome, id_vec, time_vec, conditioning) {
  if (!outcome %in% names(data)) {
    return(NULL)
  }
  pr <- dynamic$processes[[outcome]]
  marg <- predict(dynamic$reference, newdata = data, type = "distribution",
                  uncertainty = "conditional")[[outcome]]
  z <- stats::qnorm(clamp_prob(cdf(marg, data[[outcome]])))
  pieces <- lapply(unique(id_vec), function(s) {
    idx <- which(id_vec == s)
    idx <- idx[order(time_vec[idx])]
    if (length(idx) < 2L) {
      return(NULL)
    }
    lapply(seq_along(idx)[-1], function(j) {
      i1 <- idx[[j - 1L]]
      i2 <- idx[[j]]
      hist <- if (identical(conditioning, "all")) idx[seq_len(j - 1L)] else i1
      cond <- condition_z(z[hist], time_vec[hist], time_vec[i2], pr$psi, pr$process)
      r12 <- process_correlation(time_vec[c(i1, i2)], pr$psi, pr$process)[1, 2]
      y1 <- data[[outcome]][i1]
      y2 <- data[[outcome]][i2]
      dt <- time_vec[i2] - time_vec[i1]
      d_end <- vctrs::vec_slice(marg, i2)
      d_cond <- condition_norm_dist(d_end, cond$m, cond$s)
      change_u <- cdf(d_cond, y2)
      # unconditional change: (Z2 - Z1) / sqrt(2(1-r))
      change_z <- (z[i2] - z[i1]) / sqrt(pmax(2 * (1 - r12), 1e-8))
      innov <- (z[i2] - cond$m) / cond$s
      obs_v <- (y2 - y1) / dt
      exp_end <- center(d_cond)
      exp_v <- (exp_end - y1) / dt
      v_lo <- (quantile(d_cond, 0.05) - y1) / dt
      v_hi <- (quantile(d_cond, 0.95) - y1) / dt
      vel_u <- cdf(d_cond, y2)
      identified <- if (is.null(dynamic$identifiability)) {
        TRUE
      } else {
        isTRUE(dynamic$identifiability$change)
      }
      if (!identified) {
        change_z <- NA_real_
        innov <- NA_real_
      }
      tibble::tibble(
        .id = s,
        .outcome = outcome,
        .from_time = time_vec[i1],
        .to_time = time_vec[i2],
        .dt = dt,
        start_value = y1,
        end_value = y2,
        start_centile = cdf(vctrs::vec_slice(marg, i1), y1),
        end_centile = cdf(d_end, y2),
        observed_change = y2 - y1,
        expected_change = exp_end - y1,
        change_centile = change_u,
        change_z = change_z,
        observed_velocity = obs_v,
        expected_velocity = exp_v,
        velocity_lower = v_lo,
        velocity_upper = v_hi,
        velocity_centile = vel_u,
        innovation_z = innov,
        z_gain = innov,
        history_n = length(hist),
        aleatoric_sd = field_or(d_end, "aleatoric_sd"),
        measurement_sd = pr$psi$sigma_e,
        epistemic_sd = field_or(d_end, "epistemic_sd"),
        support = if (!identified) {
          "insufficient_history"
        } else {
          classify_temporal_support(dynamic, time_vec[i1], time_vec[i2], length(hist))
        },
        calibrated = identified && isTRUE(dynamic$identifiability$measurement)
      )
    })
  })
  dplyr_bind(unlist(pieces, recursive = FALSE))
}

#' Velocity summaries from a transition table
#'
#' @param transition A [norm_transition] object.
#' @param time_unit Label only; values are already in the time variable's unit.
#' @param scale `"response"` or `"centile"`.
#' @export
norm_velocity <- function(transition,
                          time_unit = "year",
                          scale = c("response", "centile")) {
  scale <- match.arg(scale)
  out <- transition
  attr(out, "time_unit") <- time_unit
  attr(out, "scale") <- scale
  if (identical(scale, "centile")) {
    out$observed_velocity <- (out$end_centile - out$start_centile) / out$.dt
  }
  class(out) <- unique(c("norm_velocity", class(out)))
  out
}

#' @export
print.norm_transition <- function(x, ...) {
  cli::cli_text("{.cls norm_transition} {nrow(x)} transition{?s}")
  NextMethod()
}
