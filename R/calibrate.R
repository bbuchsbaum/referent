#' Held-out probability recalibration
#'
#' Maps PIT values through an empirical CDF estimated on a calibration
#' sample: \eqn{F^\star(y\mid x)=G(F_0(y\mid x))}. This targets global
#' (or group-specific) calibration. It does not create full conditional
#' calibration.
#'
#' @param fit A [norm_fit].
#' @param data Calibration reference data, not reused for evaluation.
#' @param method `"rank"` empirical CDF.
#' @param by Optional grouping column (e.g. site).
#' @return A `norm_fit` with a calibration map attached.
#' @export
norm_calibrate <- function(fit, data, method = c("rank"), by = NULL) {
  method <- match.arg(method)
  data <- tibble::as_tibble(data)
  by_vec <- pull_column(data, rlang::enquo(by), default = NULL)
  scores <- predict(fit, newdata = data, type = "scores",
                    uncertainty = "conditional", allow_extrapolation = TRUE)
  pre <- norm_assess(fit, newdata = data)
  maps <- lapply(split(scores, scores$.outcome), function(sc) {
    if (is.null(by_vec)) {
      list(.global = pit_map(sc$centile))
    } else {
      groups <- split(sc$centile, by_vec[sc$.row])
      lapply(groups, pit_map)
    }
  })
  fit$calibration <- list(
    method = method,
    by = if (is.null(by_vec)) NULL else as.character(rlang::as_name(rlang::enquo(by))),
    maps = maps,
    pre = pre
  )
  class(fit) <- unique(c("norm_calibrated", class(fit)))
  fit
}

pit_map <- function(u) {
  u <- clamp_prob(u[is.finite(u)])
  stats::ecdf(u)
}

apply_calibration <- function(fit, scores, newdata = NULL) {
  if (is.null(fit$calibration)) {
    return(scores)
  }
  maps <- fit$calibration$maps
  by_nm <- fit$calibration$by
  by_vec <- if (!is.null(by_nm) && !is.null(newdata) && by_nm %in% names(newdata)) {
    as.character(newdata[[by_nm]])
  } else {
    NULL
  }
  for (i in seq_len(nrow(scores))) {
    nm <- scores$.outcome[[i]]
    mp <- maps[[nm]]
    if (is.null(mp)) {
      next
    }
    g <- if (is.null(by_vec)) {
      mp$.global
    } else {
      mp[[by_vec[scores$.row[[i]]]]] %||% mp$.global
    }
    if (is.null(g) || !is.finite(scores$centile[[i]])) {
      next
    }
    scores$centile[[i]] <- as.numeric(g(scores$centile[[i]]))
  }
  scores$z <- stats::qnorm(clamp_prob(scores$centile))
  scores$tail_prob <- 2 * pmin(scores$centile, 1 - scores$centile)
  scores$tail_surprisal <- -safe_log(scores$tail_prob)
  scores$calibrated <- TRUE
  scores
}

#' @export
predict.norm_calibrated <- function(object, newdata, ...) {
  sc <- NextMethod()
  if (inherits(sc, "norm_scores")) {
    sc <- apply_calibration(object, sc, newdata)
  }
  sc
}
