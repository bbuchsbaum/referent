fit_engine <- function(spec, data, outcome, ...) {
  switch(
    spec$engine,
    mgcv = fit_engine_mgcv(spec, data, outcome, ...),
    cli::cli_abort("Unknown engine {.val {spec$engine}}.")
  )
}

predict_engine_dist <- function(fit_one, newdata, uncertainty = "conditional",
                                n_draw = 200L, seed = 1L) {
  switch(
    fit_one$engine,
    mgcv = predict_mgcv_dist(fit_one, newdata, uncertainty, n_draw, seed),
    cli::cli_abort("Unknown engine {.val {fit_one$engine}}.")
  )
}

fit_engine_mgcv <- function(spec, data, outcome, ...) {
  fam_name <- spec$family$name
  n <- nrow(data)
  use_bam <- isTRUE(spec$use_bam) && n >= spec$bam_min_n &&
    fam_name == "gaussian" && formula_is_intercept_only(spec$scale)
  fitter0 <- if (use_bam) mgcv::bam else mgcv::gam
  # mgcv step-failure warnings are recorded on the fit rather than raised;
  # `status` reflects convergence.
  warnings <- character()
  fitter <- function(...) {
    withCallingHandlers(
      fitter0(...),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  }
  # The outcome is copied to a syntactic temporary column so that names
  # such as "brain volume" fit; `y_name` keeps the real name.
  data$.referent_y <- data[[outcome]]
  rhs_loc <- spec$location[[length(spec$location)]]
  loc_f <- stats::as.formula(eval(bquote(.referent_y ~ .(rhs_loc))))
  model <- tryCatch(
    {
      if (fam_name == "gaussian" && formula_is_intercept_only(spec$scale)) {
        fitter(loc_f, data = data, method = spec$method, ...)
      } else if (fam_name == "gaussian") {
        fitter(
          list(loc_f, spec$scale),
          data = data,
          family = mgcv::gaulss(b = spec$family$min_scale %||% 0.01),
          method = spec$method,
          ...
        )
      } else if (fam_name == "shash") {
        fitter(
          list(loc_f, spec$scale, spec$skew, spec$tail),
          data = data,
          family = mgcv::shash(b = spec$family$min_scale %||% 1e-2),
          method = spec$method,
          ...
        )
      } else {
        stop("unsupported_family", call. = FALSE)
      }
    },
    error = function(e) e
  )
  if (inherits(model, "error")) {
    msg <- conditionMessage(model)
    status <- if (grepl("unsupported_family", msg)) {
      "unsupported_family"
    } else {
      "nonconverged"
    }
    return(list(
      engine = "mgcv",
      family = spec$family,
      outcome = outcome,
      status = status,
      message = msg,
      model = NULL,
      spec = spec
    ))
  }
  conv <- mgcv_converged(model)
  status <- if (conv) "ok" else "nonconverged"
  list(
    engine = "mgcv",
    family = spec$family,
    outcome = outcome,
    status = status,
    message = if (length(warnings)) paste(unique(warnings), collapse = "; ") else NULL,
    model = model,
    spec = spec,
    y_name = outcome
  )
}

mgcv_converged <- function(model) {
  if (is.null(model)) {
    return(FALSE)
  }
  conv <- model$converged
  if (is.null(conv)) {
    return(TRUE)
  }
  isTRUE(conv)
}

predict_mgcv_dist <- function(fit_one, newdata, uncertainty = c("conditional", "total"),
                              n_draw = 200L, seed = 1L) {
  uncertainty <- match.arg(uncertainty)
  model <- fit_one$model
  fam <- fit_one$family$name
  if (is.null(model)) {
    cli::cli_abort("Cannot predict from a failed fit.")
  }
  if (is.null(newdata) || !nrow(newdata)) {
    return(distributional::dist_normal(numeric(), numeric()))
  }
  pars <- mgcv_parameters(model, newdata, fam, fit_one)
  if (identical(uncertainty, "total") && !is.null(model$Vp)) {
    if (is_plain_gaussian(model, fam)) {
      # identity location, constant scale: the total predictive is exactly
      # N(mu, sigma^2 + se^2), so no draws are needed.
      return(distributional::dist_normal(
        pars$location, sqrt(pars$scale^2 + pars$epistemic_sd^2)
      ))
    }
    return(mgcv_dist_total(model, newdata, fam, fit_one, n_draw, seed))
  }
  make_dist(fam, pars)
}

# Distribution vector from a list of parameter vectors.
make_dist <- function(fam, p) {
  if (identical(fam, "gaussian")) {
    return(distributional::dist_normal(p$location, p$scale))
  }
  dist_shash(p$location, p$scale, p$skew, p$tail)
}

# predict.gam warns about unseen factor levels; `support == "new_group"`
# already carries that information, so the warning is muffled here.
predict_gam_quiet <- function(model, newdata, ...) {
  withCallingHandlers(
    stats::predict(model, newdata = newdata, ...),
    warning = function(w) {
      if (grepl("not in original fit", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

is_plain_gaussian <- function(model, fam) {
  fam == "gaussian" && is.null(model$family$n.theta) &&
    !inherits(model$family, "general.family")
}

mgcv_parameters <- function(model, newdata, fam, fit_one) {
  n <- nrow(newdata)
  pr <- tryCatch(
    predict_gam_quiet(model, newdata, type = "response"),
    error = function(e) NULL
  )
  if (is.null(pr)) {
    return(list(
      location = rep(NA_real_, n),
      scale = rep(NA_real_, n),
      skew = 0,
      tail = 1,
      epistemic_sd = 0
    ))
  }
  se <- tryCatch(
    predict_gam_quiet(model, newdata, type = "link", se.fit = TRUE),
    error = function(e) NULL
  )
  epistemic <- 0
  if (fam == "gaussian" && is.null(model$family$n.theta) &&
      !inherits(model$family, "general.family")) {
    mu <- as.numeric(pr)
    sigma <- residual_scale(model)
    if (!is.null(se) && !is.null(se$se.fit)) {
      epistemic <- as.numeric(se$se.fit)
    }
    return(list(location = mu, scale = rep(sigma, length(mu)), skew = 0, tail = 1,
                epistemic_sd = epistemic))
  }
  pr <- as.matrix(pr)
  if (fam == "gaussian") {
    mu <- pr[, 1]
    prec <- pr[, 2]
    sigma <- 1 / pmax(prec, .Machine$double.eps)
    return(list(location = mu, scale = sigma, skew = 0, tail = 1,
                epistemic_sd = epistemic_from_se(se)))
  }
  if (fam == "shash") {
    eta <- tryCatch(
      as.matrix(predict_gam_quiet(model, newdata, type = "link")),
      error = function(e) as.matrix(pr)
    )
    mu <- eta[, 1]
    linfo <- model$family$linfo
    if (is.list(linfo) && length(linfo) >= 2 && !is.null(linfo[[2]]$linkinv)) {
      tau <- linfo[[2]]$linkinv(eta[, 2])
      sigma <- exp(tau)
    } else {
      sigma <- exp(pr[, 2])
    }
    eps <- if (ncol(eta) >= 3) eta[, 3] else 0
    phi <- if (ncol(eta) >= 4) eta[, 4] else 0
    delta <- exp(phi)
    return(list(
      location = as.numeric(mu),
      scale = pmax(as.numeric(sigma), 1e-6),
      skew = as.numeric(eps),
      tail = pmax(as.numeric(delta), 1e-3),
      epistemic_sd = epistemic_from_se(se)
    ))
  }
  mu <- as.numeric(pr)
  list(
    location = mu,
    scale = rep(residual_scale(model), length(mu)),
    skew = 0,
    tail = 1,
    epistemic_sd = epistemic
  )
}

epistemic_from_se <- function(se) {
  if (is.null(se) || is.null(se$se.fit)) {
    return(0)
  }
  sf <- se$se.fit
  if (is.matrix(sf)) {
    return(as.numeric(sf[, 1]))
  }
  as.numeric(sf)
}

residual_scale <- function(model) {
  s <- model$sig2
  if (is.null(s) || !is.finite(s) || s <= 0) {
    s <- mean(model$residuals^2, na.rm = TRUE)
  }
  sqrt(max(s, .Machine$double.eps))
}

# Equal-weight mixture over coefficient draws from N(beta_hat, Vp).
mgcv_dist_total <- function(model, newdata, fam, fit_one, n_draw, seed = 1L) {
  cond <- mgcv_parameters(model, newdata, fam, fit_one)
  n <- length(cond$location)
  lp <- tryCatch(
    predict_gam_quiet(model, newdata, type = "lpmatrix"),
    error = function(e) NULL
  )
  if (is.null(lp) || is.null(model$Vp)) {
    return(make_dist(fam, cond))
  }
  beta_hat <- stats::coef(model)
  vp <- as.matrix(model$Vp)
  keep <- is.finite(beta_hat)
  if (ncol(vp) == length(beta_hat)) {
    keep <- keep & colSums(is.finite(vp)) > 0
  }
  if (!all(keep)) {
    beta_hat <- beta_hat[keep]
    vp <- vp[keep, keep, drop = FALSE]
    lp <- lp[, keep, drop = FALSE]
    lpi <- attr(lp, "lpi")
    if (!is.null(lpi)) {
      idx <- which(keep)
      attr(lp, "lpi") <- lapply(lpi, function(cols) {
        match(cols[cols %in% idx], idx)
      })
    }
  }
  draws <- tryCatch(
    withr::with_seed(seed, mvtnorm_draw(n_draw, beta_hat, vp)),
    error = function(e) {
      matrix(beta_hat, nrow = n_draw, ncol = length(beta_hat), byrow = TRUE)
    }
  )
  loc_mat <- matrix(cond$location, n, n_draw)
  scale_mat <- matrix(cond$scale, n, n_draw)
  skew_mat <- matrix(cond$skew, n, n_draw)
  tail_mat <- matrix(cond$tail, n, n_draw)
  for (j in seq_len(n_draw)) {
    etas <- eta_from_lp(lp, draws[j, ])
    par_j <- params_from_eta(etas, model, fam, fit_one, cond)
    loc_mat[, j] <- par_j$location
    scale_mat[, j] <- par_j$scale
    skew_mat[, j] <- par_j$skew
    tail_mat[, j] <- par_j$tail
  }
  dist_shash_mc(loc_mat, scale_mat, skew_mat, tail_mat)
}

eta_from_lp <- function(lp, beta) {
  lpi <- attr(lp, "lpi")
  if (is.null(lpi)) {
    return(list(as.numeric(as.matrix(lp) %*% beta)))
  }
  lapply(lpi, function(cols) {
    cols <- cols[is.finite(cols) & cols >= 1L & cols <= length(beta)]
    if (!length(cols)) {
      return(rep(NA_real_, nrow(lp)))
    }
    as.numeric(lp[, cols, drop = FALSE] %*% beta[cols])
  })
}

params_from_eta <- function(etas, model, fam, fit_one, fallback) {
  n <- length(fallback$location)
  loc <- rep_len(etas[[1]] %||% fallback$location, n)
  if (fam == "gaussian" && length(etas) == 1L) {
    return(list(
      location = loc,
      scale = rep_len(fallback$scale, n),
      skew = rep(0, n),
      tail = rep(1, n)
    ))
  }
  if (fam == "gaussian" && length(etas) >= 2L) {
    min_s <- fit_one$family$min_scale %||% 0.01
    sigma <- min_s + exp(etas[[2]])
    return(list(
      location = loc,
      scale = pmax(sigma, 1e-6),
      skew = rep(0, n),
      tail = rep(1, n)
    ))
  }
  if (fam == "shash") {
    linfo <- model$family$linfo
    if (is.list(linfo) && length(linfo) >= 2 && !is.null(linfo[[2]]$linkinv)) {
      tau <- linfo[[2]]$linkinv(etas[[2]])
      sigma <- exp(tau)
    } else {
      sigma <- exp(etas[[2]])
    }
    eps <- if (length(etas) >= 3L) etas[[3]] else fallback$skew
    phi <- if (length(etas) >= 4L) etas[[4]] else 0
    return(list(
      location = loc,
      scale = pmax(as.numeric(sigma), 1e-6),
      skew = rep_len(eps, n),
      tail = pmax(rep_len(exp(phi), n), 1e-3)
    ))
  }
  fallback
}

mvtnorm_draw <- function(n, mean, sigma) {
  ev <- eigen(sigma, symmetric = TRUE)
  ev$values <- pmax(ev$values, 0)
  a <- ev$vectors %*% diag(sqrt(ev$values), nrow = length(ev$values))
  z <- matrix(stats::rnorm(n * length(mean)), n, length(mean))
  sweep(z %*% t(a), 2, mean, "+")
}
