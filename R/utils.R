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
  if (is.data.frame(scores) && nrow(scores) && all(scores$.in_sample)) {
    cli::cli_warn(
      "These scores are in-sample ({.field .in_sample} = TRUE) and should not be used for {used_for}."
    )
  }
  invisible(scores)
}

has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

# Column name named by a quosure, or `default` when it is NULL, missing, or
# not a name.
as_col_name <- function(quo, default = NULL) {
  if (rlang::quo_is_null(quo) || rlang::quo_is_missing(quo)) {
    return(default)
  }
  tryCatch(rlang::as_name(quo), error = function(e) default)
}

# Default x-axis covariate: `age` when the fit uses it, else the first
# numeric covariate, else the first covariate.
default_x <- function(fit) {
  if ("age" %in% fit$covariates) {
    return("age")
  }
  names(fit$support_ref$numeric)[1] %||% fit$covariates[[1]]
}

# Pivot a long score table to a rows x outcomes matrix of `value`, or read
# the `.<value>_<outcome>` columns of an augment() table.
scores_matrix <- function(scores, value = "z") {
  wide <- grep(paste0("^\\.", value, "_"), names(scores), value = TRUE)
  if (!".outcome" %in% names(scores) && length(wide)) {
    mat <- as.matrix(scores[, wide, drop = FALSE])
    colnames(mat) <- sub(paste0("^\\.", value, "_"), "", wide)
    rows <- seq_len(nrow(mat))
    return(list(matrix = mat, .row = rows, .id = scores[[".id"]] %||% rows))
  }
  if (!value %in% names(scores)) {
    cli::cli_abort("Unknown score column {.field {value}}.")
  }
  row_col <- if (".row" %in% names(scores)) ".row" else ".id"
  rows <- unique(scores[[row_col]])
  outs <- unique(scores$.outcome)
  mat <- matrix(NA_real_, length(rows), length(outs), dimnames = list(NULL, outs))
  mat[cbind(match(scores[[row_col]], rows), match(scores$.outcome, outs))] <-
    as.numeric(scores[[value]])
  ids <- if (".id" %in% names(scores)) {
    scores$.id[match(rows, scores[[row_col]])]
  } else {
    rows
  }
  list(matrix = mat, .row = rows, .id = ids)
}
