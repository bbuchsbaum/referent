`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
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
