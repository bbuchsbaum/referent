# The mgcv family is chosen from the referent family and the scale
# formula: a Gaussian with constant scale is a plain `gam()` (or `bam()`
# for large n, see below); a Gaussian with a modelled scale is `gaulss()`;
# SHASH is `shash()`. `mgcv_family()` rebuilds the family object from this
# name and the minimum scale, so a frozen bundle need not carry mgcv's
# family closures.
mgcv_family_name <- function(spec) {
  if (spec$family$name == "shash") {
    "shash"
  } else if (formula_is_intercept_only(spec$scale)) {
    "gaussian"
  } else {
    "gaulss"
  }
}

mgcv_family <- function(fit_one) {
  b <- fit_one$family$min_scale %||% 0.01
  switch(
    fit_one$mgcv_family,
    gaussian = stats::gaussian(),
    gaulss = mgcv::gaulss(b = b),
    shash = mgcv::shash(b = b)
  )
}

fit_engine_mgcv <- function(spec, data, outcome, ...) {
  mgcv_fam <- mgcv_family_name(spec)
  n <- nrow(data)
  # `bam()` is used only for the constant-scale Gaussian (a single linear
  # predictor, no general family) when n >= spec$bam_min_n; every other
  # model goes through `gam()`.
  use_bam <- isTRUE(spec$use_bam) && n >= spec$bam_min_n && mgcv_fam == "gaussian"
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
  # such as "brain volume" fit; `y_name` keeps the real name. The formula
  # keeps the spec's environment: a formula created here would capture
  # this frame (and with it `data`) and serialise it with every gam.
  data$.referent_y <- data[[outcome]]
  rhs_loc <- spec$location[[length(spec$location)]]
  loc_f <- eval(bquote(.referent_y ~ .(rhs_loc)))
  environment(loc_f) <- environment(spec$location) %||% globalenv()
  stub <- list(
    engine = "mgcv", family = spec$family, outcome = outcome, spec = spec,
    mgcv_family = mgcv_fam, y_name = outcome, n_obs = n
  )
  family <- mgcv_family(stub)
  model <- tryCatch(
    switch(
      mgcv_fam,
      gaussian = fitter(loc_f, data = data, method = spec$method, ...),
      gaulss = fitter(list(loc_f, spec$scale), data = data, family = family,
                      method = spec$method, ...),
      shash = fitter(list(loc_f, spec$scale, spec$skew, spec$tail), data = data,
                     family = family, method = spec$method, ...)
    ),
    error = function(e) e
  )
  if (inherits(model, "error")) {
    return(c(stub, list(status = "nonconverged", message = conditionMessage(model),
                        model = NULL)))
  }
  # Random-effect smooths keep a model formula whose environment is
  # mgcv's construction frame; give it the spec's environment as well.
  model$smooth <- lapply(model$smooth, function(sm) {
    if (!is.null(sm$form)) {
      environment(sm$form) <- environment(loc_f)
    }
    sm
  })
  c(stub, list(
    status = if (isTRUE(model$converged %||% TRUE)) "ok" else "nonconverged",
    message = if (length(warnings)) paste(unique(warnings), collapse = "; ") else NULL,
    model = model
  ))
}

# Frozen bundles drop the family from each gam (see `strip_gam()`); the
# model must carry one for `predict.gam()` to take the right branch.
thaw_model <- function(fit_one) {
  model <- fit_one$model
  if (is.null(model$family)) {
    model$family <- mgcv_family(fit_one)
  }
  model
}

predict_mgcv_dist <- function(fit_one, newdata, uncertainty = "conditional",
                              n_draw = 200L, seed = 1L, exclude = NULL) {
  if (is.null(fit_one$model)) {
    cli::cli_abort("Cannot predict from a failed fit.")
  }
  if (is.null(newdata) || !nrow(newdata)) {
    return(distributional::dist_normal(numeric(), numeric()))
  }
  model <- thaw_model(fit_one)
  fam <- fit_one$family$name
  pars <- mgcv_parameters(model, newdata, fam, fit_one, exclude = exclude)
  if (identical(uncertainty, "total") && !is.null(model$Vp)) {
    if (identical(fit_one$mgcv_family, "gaussian")) {
      # identity location, constant scale: the total predictive is exactly
      # N(mu, sigma^2 + se^2), so no draws are needed.
      return(distributional::dist_normal(
        pars$location, sqrt(pars$scale^2 + pars$epistemic_sd^2)
      ))
    }
    return(mgcv_dist_total(model, newdata, fam, fit_one, pars, n_draw, seed, exclude))
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

# Rows of `newdata` whose parametric factor level was not seen in the fit.
# predict.gam cannot build the parametric model matrix for such rows and
# fails for the whole call, so they are held out of prediction and given
# NA parameters.
unseen_level_rows <- function(model, newdata) {
  bad <- rep(FALSE, nrow(newdata))
  for (nm in intersect(names(model$xlevels), names(newdata))) {
    v <- as.character(newdata[[nm]])
    bad <- bad | (!is.na(v) & !v %in% model$xlevels[[nm]])
  }
  bad
}

# Labels of the smooths whose terms include any of `vars` (random-effect
# and factor-smooth terms of a grouping covariate, for instance).
smooths_using <- function(model, vars) {
  sm <- model$smooth
  uses <- vapply(sm, function(s) any(c(s$term, s$fterm) %in% vars), logical(1))
  vapply(sm[uses], `[[`, "", "label")
}

# Rows whose level of a factor-smooth (`bs = "fs"`) factor was not seen
# in the fit. mgcv predicts NA for them, whereas an unseen level of a
# `bs = "re"` smooth is given a zero effect; the two are reconciled by
# excluding the factor smooth for such rows (the population curve).
unseen_fs_rows <- function(model, newdata) {
  bad <- rep(FALSE, nrow(newdata))
  for (s in model$smooth) {
    if (!inherits(s, "fs.interaction") || !s$fterm %in% names(newdata)) {
      next
    }
    v <- as.character(newdata[[s$fterm]])
    bad <- bad | (!is.na(v) & !v %in% s$flev)
  }
  bad
}

# predict.gam on the rows without unseen parametric levels; NULL when
# prediction fails. `pad` rebuilds a full-length result from the kept
# rows' values. Rows with an unseen factor-smooth level are predicted
# with that smooth excluded.
predict_gam_rows <- function(model, newdata, pad, ..., exclude = NULL) {
  bad <- unseen_level_rows(model, newdata)
  if (all(bad)) {
    return(pad(NULL, bad))
  }
  fs_new <- unseen_fs_rows(model, newdata) & !bad
  fs_labels <- vapply(Filter(function(s) inherits(s, "fs.interaction"), model$smooth),
                      `[[`, "", "label")
  one <- function(rows, excl) {
    tryCatch(
      predict_gam_quiet(model, newdata[rows, , drop = FALSE], ..., exclude = excl),
      error = function(e) NULL
    )
  }
  pr <- if (!any(fs_new)) {
    one(!bad, exclude)
  } else {
    pieces <- list(one(!bad & !fs_new, exclude), one(fs_new, union(exclude, fs_labels)))
    if (any(vapply(pieces, is.null, logical(1)))) NULL else bind_pred(pieces, !bad, fs_new)
  }
  if (is.null(pr)) {
    return(NULL)
  }
  pad(pr, bad)
}

# Reassemble two predict.gam results (plain or `se.fit` lists, vectors or
# matrices) into one in the order of the kept rows.
bind_pred <- function(pieces, kept, second) {
  ord <- order(c(which(kept & !second), which(second)))
  join <- function(a, b) {
    if (is.list(a)) {
      return(Map(join, a, b))
    }
    out <- if (is.matrix(a)) rbind(a, b)[ord, , drop = FALSE] else c(a, b)[ord]
    attr(out, "lpi") <- attr(a, "lpi")
    out
  }
  join(pieces[[1L]], pieces[[2L]])
}

# One `predict(type = "link", se.fit = TRUE)` call: a list of linear
# predictors (one vector per lp) and the standard error of the first.
mgcv_link <- function(model, newdata, exclude = NULL) {
  predict_gam_rows(model, newdata, type = "link", se.fit = TRUE, exclude = exclude,
                   pad = function(pr, bad) {
    n <- length(bad)
    if (is.null(pr)) {
      return(list(eta = list(rep(NA_real_, n)), se = rep(NA_real_, n)))
    }
    eta <- as.matrix(pr$fit)
    se <- as.matrix(pr$se.fit)
    full <- function(v) replace(rep(NA_real_, n), !bad, v)
    list(eta = lapply(seq_len(ncol(eta)), function(k) full(eta[, k])),
         se = full(as.numeric(se[, 1L])))
  })
}

# The lp matrix with NA rows for unseen parametric levels; NULL when
# prediction fails or no row can be predicted.
mgcv_lpmatrix <- function(model, newdata, exclude = NULL) {
  predict_gam_rows(model, newdata, type = "lpmatrix", exclude = exclude, pad = function(lp, bad) {
    if (is.null(lp)) {
      return(NULL)
    }
    out <- matrix(NA_real_, length(bad), ncol(lp))
    out[!bad, ] <- lp
    attr(out, "lpi") <- attr(lp, "lpi")
    out
  })
}

mgcv_parameters <- function(model, newdata, fam, fit_one, exclude = NULL) {
  n <- nrow(newdata)
  link <- mgcv_link(model, newdata, exclude)
  if (is.null(link)) {
    return(list(location = rep(NA_real_, n), scale = rep(NA_real_, n),
                skew = rep(0, n), tail = rep(1, n), epistemic_sd = rep(0, n)))
  }
  c(params_from_eta(link$eta, model, fam, fit_one), list(epistemic_sd = link$se))
}

# Distribution parameters from the linear predictors, through the mgcv
# link functions. Every element of `etas` may be a vector or an n x K
# matrix of draws; the result has the same shape. `fallback` supplies
# the scale when a Gaussian has a single linear predictor.
params_from_eta <- function(etas, model, fam, fit_one, fallback = NULL) {
  loc <- etas[[1L]]
  full <- function(v) if (is.matrix(loc)) array(v, dim(loc)) else rep_len(v, length(loc))
  lp <- function(k) if (length(etas) >= k) etas[[k]] else NULL
  if (fam == "gaussian" && length(etas) == 1L) {
    scale <- if (is.null(fallback)) residual_scale(model) else fallback$scale
    return(list(location = loc, scale = full(scale), skew = full(0), tail = full(1)))
  }
  linfo <- (model$family %||% mgcv_family(fit_one))$linfo
  if (fam == "gaussian") {
    # gaulss: the second link inverse returns the precision 1 / sigma.
    sigma <- 1 / pmax(linfo[[2L]]$linkinv(lp(2L)), .Machine$double.eps)
    return(list(location = loc, scale = sigma, skew = full(0), tail = full(1)))
  }
  # shash: the second link inverse returns log sigma; the fourth lp is
  # log delta.
  sigma <- exp(linfo[[2L]]$linkinv(lp(2L)))
  list(
    location = loc,
    scale = pmax(sigma, 1e-6),
    skew = lp(3L) %||% full(0),
    tail = pmax(exp(lp(4L) %||% full(0)), 1e-3)
  )
}

residual_scale <- function(model) {
  s <- model$sig2
  if (is.null(s) || !is.finite(s) || s <= 0) {
    s <- mean(model$residuals^2, na.rm = TRUE)
  }
  sqrt(max(s, .Machine$double.eps))
}

# Equal-weight mixture over coefficient draws from N(beta_hat, Vp). Each
# linear predictor is one product of the lp matrix with the draw matrix.
mgcv_dist_total <- function(model, newdata, fam, fit_one, cond, n_draw, seed = 1L,
                            exclude = NULL) {
  lp <- mgcv_lpmatrix(model, newdata, exclude)
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
  par <- params_from_eta(eta_from_lp(lp, t(draws)), model, fam, fit_one, cond)
  dist_shash_draws(par$location, par$scale, par$skew, par$tail)
}

# Linear predictors for a coefficient vector, or an n x K matrix of them
# for a p x K matrix of coefficient draws.
eta_from_lp <- function(lp, beta) {
  beta <- as.matrix(beta)
  shape <- function(m) if (ncol(beta) == 1L) as.numeric(m) else unname(m)
  lpi <- attr(lp, "lpi")
  if (is.null(lpi)) {
    return(list(shape(as.matrix(lp) %*% beta)))
  }
  lapply(lpi, function(cols) {
    cols <- cols[is.finite(cols) & cols >= 1L & cols <= nrow(beta)]
    if (!length(cols)) {
      return(shape(array(NA_real_, c(nrow(lp), ncol(beta)))))
    }
    shape(lp[, cols, drop = FALSE] %*% beta[cols, , drop = FALSE])
  })
}

mvtnorm_draw <- function(n, mean, sigma) {
  ev <- eigen(sigma, symmetric = TRUE)
  ev$values <- pmax(ev$values, 0)
  a <- ev$vectors %*% diag(sqrt(ev$values), nrow = length(ev$values))
  z <- matrix(stats::rnorm(n * length(mean)), n, length(mean))
  sweep(z %*% t(a), 2, mean, "+")
}

# Reduce a fitted gam to what prediction needs: coefficients and their
# covariance, the smooth constructions, the formula and terms objects,
# factor levels, and a few scalars. The model frame is kept with zero
# rows: predict.gam reads its column classes and factor levels (to
# coerce newdata) and its terms attribute (to find offset terms). Fitted
# values, residuals, weights, and the family (rebuilt by `thaw_model()`)
# are dropped, and the formula environments are replaced by an empty
# child of the global environment so no enclosing frame is serialised.
strip_gam <- function(model) {
  keep <- c(
    "coefficients", "Vp", "nsdf", "smooth", "formula", "pred.formula", "terms",
    "pterms", "xlevels", "contrasts", "cmX", "assign", "paraPen", "rank", "sp",
    "sig2", "edf", "converged", "df.residual", "var.summary", "lpi"
  )
  out <- model[intersect(names(model), keep)]
  if (!is.null(model$model)) {
    mf <- model$model[0, , drop = FALSE]
    attr(mf, "terms") <- attr(model$model, "terms")
    out$model <- mf
  }
  env <- new.env(parent = globalenv())
  reset_env <- function(f) {
    if (is.list(f)) {
      f <- lapply(f, reset_env)
    } else if (!is.null(f) && (inherits(f, "formula") || inherits(f, "terms"))) {
      environment(f) <- env
    }
    f
  }
  if (!is.null(attr(out$model, "terms"))) {
    attr(out$model, "terms") <- reset_env(attr(out$model, "terms"))
  }
  for (nm in c("formula", "pred.formula", "terms", "pterms")) {
    if (!is.null(out[[nm]])) {
      attrs <- attributes(out[[nm]])
      out[[nm]] <- reset_env(out[[nm]])
      if (is.list(out[[nm]])) {
        attributes(out[[nm]]) <- attrs
      }
    }
  }
  out$smooth <- lapply(out$smooth, function(s) {
    if (!is.null(s$form)) {
      environment(s$form) <- env
    }
    s
  })
  class(out) <- class(model)
  out
}
