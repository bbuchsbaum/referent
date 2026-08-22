#' referent: Distributional reference models
#'
#' A small distributional reference-modeling layer: given a reference
#' population, estimate the conditional distribution of one or more
#' outcomes; validate that distribution; place new observations within it;
#' and distinguish trustworthy deviation from extrapolation, site shift,
#' measurement noise, and model uncertainty.
#'
#' The package is CDF-first. A Z-score is only one representation of a
#' centile. The term "abnormality score" is not part of the core API.
#'
#' @keywords internal
#' @importFrom rlang .data :=
#' @importFrom ggplot2 autoplot
#' @importFrom stats density quantile predict
"_PACKAGE"
