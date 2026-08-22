fit_engine_gamlss <- function(spec, data, outcome, ...) {
  if (!has_pkg("gamlss")) {
    return(list(
      engine = "gamlss",
      family = spec$family,
      outcome = outcome,
      status = "unsupported_family",
      message = "Package 'gamlss' is not installed.",
      model = NULL,
      spec = spec
    ))
  }
  rhs <- spec$location[[length(spec$location)]]
  form <- stats::as.formula(eval(bquote(.(as.name(outcome)) ~ .(rhs))))
  fam <- if (spec$family$name == "shash") {
    gamlss::SHASH()
  } else {
    gamlss::NO()
  }
  model <- tryCatch(
    gamlss::gamlss(
      form,
      sigma.formula = spec$scale,
      nu.formula = spec$skew,
      tau.formula = spec$tail,
      data = data,
      family = fam,
      ...
    ),
    error = function(e) e
  )
  if (inherits(model, "error")) {
    return(list(
      engine = "gamlss",
      family = spec$family,
      outcome = outcome,
      status = "nonconverged",
      message = conditionMessage(model),
      model = NULL,
      spec = spec
    ))
  }
  list(
    engine = "gamlss",
    family = spec$family,
    outcome = outcome,
    status = "ok",
    message = NULL,
    model = model,
    spec = spec
  )
}

fit_engine_gamlss2 <- function(spec, data, outcome, ...) {
  if (!has_pkg("gamlss2")) {
    return(list(
      engine = "gamlss2",
      family = spec$family,
      outcome = outcome,
      status = "unsupported_family",
      message = "Package 'gamlss2' is not installed.",
      model = NULL,
      spec = spec
    ))
  }
  fit_engine_gamlss(spec, data, outcome, ...)
}

predict_gamlss_dist <- function(fit_one, newdata, uncertainty = "conditional",
                                n_draw = 64L) {
  if (is.null(fit_one$model)) {
    cli::cli_abort("Cannot predict from a failed gamlss fit.")
  }
  mu <- as.numeric(stats::predict(fit_one$model, newdata = newdata, what = "mu"))
  sigma <- tryCatch(
    as.numeric(stats::predict(fit_one$model, newdata = newdata, what = "sigma")),
    error = function(e) rep(residual_scale(fit_one$model), length(mu))
  )
  eps <- tryCatch(
    as.numeric(stats::predict(fit_one$model, newdata = newdata, what = "nu")),
    error = function(e) 0
  )
  delta <- tryCatch(
    as.numeric(stats::predict(fit_one$model, newdata = newdata, what = "tau")),
    error = function(e) 1
  )
  fam <- if (fit_one$family$name == "shash") "shash" else "gaussian"
  norm_dist(
    family = fam,
    location = mu,
    scale = pmax(sigma, 1e-6),
    skew = eps,
    tail = pmax(delta, 1e-3),
    aleatoric_sd = pmax(sigma, 1e-6),
    epistemic_sd = 0
  )
}
