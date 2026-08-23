# Classification and uncertainty helpers for paired held-out comparisons.

pcn_benchmark_mace <- function(u, grid = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  mean(vapply(grid, function(q) abs(q - mean(u <= q)), numeric(1)))
}

pcn_cluster_bootstrap <- function(delta, cluster, B = 1999L, seed = 20260823L) {
  stopifnot(length(delta) == length(cluster), length(delta) > 1L, B >= 99L)
  keep <- is.finite(delta) & !is.na(cluster)
  parts <- split(delta[keep], as.character(cluster[keep]))
  if (!length(parts)) stop("no finite paired score differences", call. = FALSE)
  draws <- withr::with_seed(seed, replicate(B, {
    selected <- sample.int(length(parts), length(parts), replace = TRUE)
    mean(unlist(parts[selected], use.names = FALSE))
  }))
  c(
    estimate = mean(unlist(parts, use.names = FALSE)),
    lower = unname(stats::quantile(draws, 0.025, names = FALSE)),
    upper = unname(stats::quantile(draws, 0.975, names = FALSE)),
    clusters = length(parts)
  )
}

pcn_calibration_summary <- function(prediction) {
  required <- c("observed", "q05", "q95", "centile")
  if (!all(required %in% names(prediction))) {
    stop("prediction is missing calibration columns", call. = FALSE)
  }
  u <- prediction$centile
  tail_prob <- 2 * pmin(u, 1 - u)
  c(
    coverage90 = mean(prediction$observed >= prediction$q05 & prediction$observed <= prediction$q95),
    mace = pcn_benchmark_mace(u),
    tail05 = mean(tail_prob <= 0.05)
  )
}

pcn_calibration_gate <- function(referent, comparator,
                                 coverage_regret = 0.01,
                                 mace_regret = 0.01,
                                 tail_regret = 0.01) {
  ref <- pcn_calibration_summary(referent)
  pcn <- pcn_calibration_summary(comparator)
  regret <- c(
    coverage90 = abs(ref[["coverage90"]] - 0.90) - abs(pcn[["coverage90"]] - 0.90),
    mace = ref[["mace"]] - pcn[["mace"]],
    tail05 = abs(ref[["tail05"]] - 0.05) - abs(pcn[["tail05"]] - 0.05)
  )
  limits <- c(coverage90 = coverage_regret, mace = mace_regret, tail05 = tail_regret)
  list(pass = all(regret <= limits), referent = ref, comparator = pcn,
       regret = regret, limits = limits)
}

pcn_classify_benchmark <- function(estimate, lower, upper, calibration_pass,
                                   conformance_pass = FALSE,
                                   equivalence_margin = 0.02,
                                   noninferiority_margin = 0.02) {
  if (isTRUE(conformance_pass)) return("conformant")
  if (!isTRUE(calibration_pass)) {
    return(if (estimate > 0) "tradeoff" else "inferior")
  }
  if (lower > 0) return("superior")
  if (upper < -noninferiority_margin) return("inferior")
  if (lower >= -equivalence_margin && upper <= equivalence_margin) return("equivalent")
  if (lower >= -noninferiority_margin) return("non_inferior")
  "inconclusive"
}

pcn_compare_predictions <- function(referent, comparator, cluster,
                                    scenario = NA_character_, B = 1999L,
                                    seed = 20260823L) {
  if (!identical(referent$row_id, comparator$row_id)) {
    stop("paired prediction row identities differ", call. = FALSE)
  }
  interval <- pcn_cluster_bootstrap(
    referent$log_density - comparator$log_density, cluster, B = B, seed = seed
  )
  calibration <- pcn_calibration_gate(referent, comparator)
  data.frame(
    scenario = scenario,
    n_test = nrow(referent),
    clusters = unname(interval[["clusters"]]),
    log_score_difference = unname(interval[["estimate"]]),
    ci_lower = unname(interval[["lower"]]),
    ci_upper = unname(interval[["upper"]]),
    referent_coverage90 = calibration$referent[["coverage90"]],
    pcntoolkit_coverage90 = calibration$comparator[["coverage90"]],
    referent_mace = calibration$referent[["mace"]],
    pcntoolkit_mace = calibration$comparator[["mace"]],
    referent_tail05 = calibration$referent[["tail05"]],
    pcntoolkit_tail05 = calibration$comparator[["tail05"]],
    calibration_pass = calibration$pass,
    classification = pcn_classify_benchmark(
      interval[["estimate"]], interval[["lower"]], interval[["upper"]],
      calibration$pass
    ),
    stringsAsFactors = FALSE
  )
}

pcn_release_superiority_summary <- function(results, B = 4999L,
                                            seed = 20260823L) {
  required <- c(
    "replicate", "log_score_difference", "referent_coverage90",
    "pcntoolkit_coverage90", "referent_mace", "pcntoolkit_mace",
    "referent_tail05", "pcntoolkit_tail05", "fit_status", "fit_converged"
  )
  if (!all(required %in% names(results))) {
    stop("release results are missing required fields", call. = FALSE)
  }
  if (nrow(results) < 5L || anyDuplicated(results$replicate)) {
    stop("release superiority requires at least five unique replicates", call. = FALSE)
  }
  numeric_evidence <- results[c(
    "log_score_difference", "referent_coverage90", "pcntoolkit_coverage90",
    "referent_mace", "pcntoolkit_mace", "referent_tail05", "pcntoolkit_tail05"
  )]
  if (any(!is.finite(as.matrix(numeric_evidence))) || anyNA(results$fit_converged)) {
    stop("release results contain incomplete evidence", call. = FALSE)
  }
  interval <- pcn_cluster_bootstrap(
    results$log_score_difference, results$replicate, B = B, seed = seed
  )
  fit_pass <- all(results$fit_status == "ok" & results$fit_converged)
  calibration_regret <- cbind(
    coverage90 = abs(results$referent_coverage90 - 0.90) -
      abs(results$pcntoolkit_coverage90 - 0.90),
    mace = results$referent_mace - results$pcntoolkit_mace,
    tail05 = abs(results$referent_tail05 - 0.05) -
      abs(results$pcntoolkit_tail05 - 0.05)
  )
  calibration_limits <- c(coverage90 = 0.01, mace = 0.01, tail05 = 0.01)
  critical_limits <- c(coverage90 = 0.05, mace = 0.05, tail05 = 0.05)
  calibration_upper <- vapply(seq_len(ncol(calibration_regret)), function(j) {
    pcn_cluster_bootstrap(
      calibration_regret[, j], results$replicate, B = B,
      seed = seed + j
    )[["upper"]]
  }, numeric(1))
  names(calibration_upper) <- colnames(calibration_regret)
  mean_regret <- colMeans(calibration_regret)
  max_regret <- apply(calibration_regret, 2L, max)
  calibration_pass <- all(calibration_upper <= calibration_limits)
  critical_pass <- all(max_regret <= critical_limits)
  classification <- if (fit_pass) {
    pcn_classify_benchmark(
      interval[["estimate"]], interval[["lower"]], interval[["upper"]],
      calibration_pass = calibration_pass && critical_pass
    )
  } else {
    "fit_failure"
  }
  data.frame(
    replicates = nrow(results),
    log_score_difference = unname(interval[["estimate"]]),
    ci_lower = unname(interval[["lower"]]),
    ci_upper = unname(interval[["upper"]]),
    all_shash_fits_converged = fit_pass,
    coverage_regret = mean_regret[["coverage90"]],
    coverage_regret_upper = calibration_upper[["coverage90"]],
    mace_regret = mean_regret[["mace"]],
    mace_regret_upper = calibration_upper[["mace"]],
    tail_regret = mean_regret[["tail05"]],
    tail_regret_upper = calibration_upper[["tail05"]],
    maximum_replicate_regret = max(max_regret),
    calibration_noninferior = calibration_pass,
    no_critical_replicate_regression = critical_pass,
    classification = classification,
    pass = identical(classification, "superior"),
    stringsAsFactors = FALSE
  )
}
