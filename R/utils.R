`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

clamp_prob <- function(p, eps = 1e-12) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

recycle_to <- function(x, n) {
  if (is.null(x)) {
    return(NULL)
  }
  x <- as.numeric(x)
  if (length(x) == 1L) {
    return(rep(x, n))
  }
  if (length(x) != n) {
    cli::cli_abort("Parameter length {length(x)} is incompatible with n = {n}.")
  }
  x
}

# Recycle a distribution and a vector so they share a common length.
# A scalar distribution may be evaluated at many y/p values.
recycle_pair <- function(distribution, y) {
  y <- as.numeric(y)
  n_d <- length(distribution)
  n_y <- length(y)
  if (n_d == n_y) {
    return(list(distribution = distribution, y = y))
  }
  if (n_d == 1L) {
    extra <- attributes(distribution)
    extra <- extra[setdiff(names(extra), c("names", "class", "row.names"))]
    d2 <- vctrs::vec_rep(distribution, n_y)
    if (inherits(distribution, "norm_dist_conditional")) {
      class(d2) <- unique(c("norm_dist_conditional", class(d2)))
    }
    for (nm in setdiff(names(extra), names(attributes(d2)))) {
      attr(d2, nm) <- extra[[nm]]
    }
    return(list(distribution = d2, y = y))
  }
  if (n_y == 1L) {
    return(list(distribution = distribution, y = rep(y, n_d)))
  }
  cli::cli_abort(
    "Cannot recycle distribution of length {n_d} with a vector of length {n_y}."
  )
}




select_outcomes <- function(quo, data) {
  expr <- rlang::quo_get_expr(quo)
  nms <- if (is.character(expr)) {
    expr
  } else {
    val <- tryCatch(rlang::eval_tidy(quo), error = function(e) NULL)
    if (is.character(val)) {
      val
    } else {
      names(tidyselect::eval_select(quo, data = data))
    }
  }
  missing <- setdiff(nms, names(data))
  if (length(missing)) {
    cli::cli_abort("Unknown outcome{?s}: {.field {missing}}.")
  }
  if (!length(nms)) {
    cli::cli_abort("No outcomes selected.")
  }
  nms
}


safe_log <- function(x) {
  log(pmax(as.numeric(x), .Machine$double.xmin))
}


se_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return(NA_real_)
  }
  stats::sd(x) / sqrt(length(x))
}


warn_in_sample <- function(scores, used_for = "comparison") {
  if (isTRUE(attr(scores, "in_sample")) ||
      (is.data.frame(scores) && isTRUE(scores$.in_sample[[1L]]) &&
       all(scores$.in_sample))) {
    cli::cli_warn(
      "These scores are in-sample ({.field .in_sample} = TRUE) and should not be used for {used_for}."
    )
  }
  invisible(scores)
}

has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}
