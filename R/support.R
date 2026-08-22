support_reference <- function(data, covariate_names) {
  covs <- data[, covariate_names, drop = FALSE]
  numeric_names <- names(covs)[vapply(covs, is.numeric, logical(1))]
  factor_names <- names(covs)[vapply(covs, function(z) is.factor(z) || is.character(z), logical(1))]
  ranges <- lapply(numeric_names, function(nm) {
    z <- covs[[nm]]
    z <- z[is.finite(z)]
    qs <- stats::quantile(z, probs = c(0, 0.02, 0.98, 1), na.rm = TRUE)
    list(min = qs[[1]], edge_lo = qs[[2]], edge_hi = qs[[3]], max = qs[[4]],
         mean = mean(z), sd = stats::sd(z) %||% 1)
  })
  names(ranges) <- numeric_names
  levels <- lapply(factor_names, function(nm) unique(as.character(covs[[nm]])))
  names(levels) <- factor_names
  X <- NULL
  if (length(numeric_names)) {
    X <- scale(as.matrix(covs[, numeric_names, drop = FALSE]))
    X[!is.finite(X)] <- 0
  }
  list(
    numeric = ranges,
    factor_levels = levels,
    numeric_names = numeric_names,
    factor_names = factor_names,
    center = if (!is.null(X)) colMeans(X, na.rm = TRUE) else NULL,
    cov = if (!is.null(X) && nrow(X) > ncol(X) + 2) {
      stats::cov(X, use = "pairwise.complete.obs")
    } else {
      NULL
    }
  )
}

#' Classify whether covariates lie in the reference support
#'
#' Statuses: `in`, `edge`, `out`, `new_group`.
#'
#' @param ref A support reference from a [norm_fit].
#' @param newdata Target data.
#' @export
classify_support <- function(ref, newdata) {
  n <- nrow(newdata)
  status <- rep("in", n)
  for (nm in ref$numeric_names) {
    if (!nm %in% names(newdata)) {
      next
    }
    z <- newdata[[nm]]
    r <- ref$numeric[[nm]]
    status[is.finite(z) & (z < r$edge_lo | z > r$edge_hi) & status == "in"] <- "edge"
    status[is.finite(z) & (z < r$min | z > r$max)] <- "out"
  }
  for (nm in ref$factor_names) {
    if (!nm %in% names(newdata)) {
      next
    }
    unseen <- !as.character(newdata[[nm]]) %in% ref$factor_levels[[nm]]
    status[unseen] <- "new_group"
  }
  if (!is.null(ref$cov) && length(ref$numeric_names)) {
    X <- scale_with_ref(newdata[, ref$numeric_names, drop = FALSE], ref)
    d2 <- mahalanobis_safe(X, ref$cov)
    q <- stats::qchisq(0.99, df = max(ncol(X), 1))
    status[is.finite(d2) & d2 > q & status == "in"] <- "edge"
    status[is.finite(d2) & d2 > q * 2] <- "out"
  }
  status
}

scale_with_ref <- function(data, ref) {
  X <- as.matrix(data)
  for (j in seq_along(ref$numeric_names)) {
    nm <- ref$numeric_names[[j]]
    s <- ref$numeric[[nm]]$sd
    if (!is.finite(s) || s == 0) {
      s <- 1
    }
    X[, j] <- (X[, j] - ref$numeric[[nm]]$mean) / s
  }
  X[!is.finite(X)] <- 0
  X
}

mahalanobis_safe <- function(X, cov) {
  cov <- as.matrix(cov)
  cov <- cov + diag(1e-6, nrow(cov))
  inv <- tryCatch(solve(cov), error = function(e) NULL)
  if (is.null(inv)) {
    return(rep(0, nrow(X)))
  }
  rowSums((X %*% inv) * X)
}
