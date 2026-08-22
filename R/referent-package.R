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
#' @importFrom stats pnorm qnorm dnorm rnorm ppoints median sd var density
#'   quantile predict coef vcov model.matrix complete.cases setNames
#'   as.formula terms delete.response na.omit pchisq p.adjust
#'   ave dist
#' @importFrom utils packageVersion head
"_PACKAGE"
