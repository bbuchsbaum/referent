# Render checked-in, machine-readable evidence from pinned fixtures.

source("tools/pcntoolkit/compare_results.R", local = TRUE)
source("tools/pcntoolkit/run_concordance.R", local = TRUE)
source("tools/pcntoolkit/run_superiority.R", local = TRUE)

render_pcntoolkit_evidence <- function(output_dir, fixture_root = pcn_fixture_root()) {
  root <- fixture_root
  manifest <- pcn_validate_fixture(root)
  concordance <- pcn_run_concordance(root)
  superiority <- pcn_run_superiority(root)
  gaussian <- utils::read.csv(file.path(root, "semantic_gaussian.csv"))
  gaussian_error <- max(abs(
    stats::pnorm(gaussian$y, gaussian$mu, gaussian$sigma) - gaussian$cdf
  ))
  shash <- utils::read.csv(file.path(root, "semantic_shashb.csv"))
  shash_dist <- referent::dist_shash(shash$mu_r, shash$sigma_r, shash$epsilon, shash$delta)
  shash_error <- max(abs(
    referent:::dist_eval(shash_dist, referent:::log_dens, shash$y) - shash$log_density
  ))
  matched <- concordance$evidence_class == "matched_estimator"
  evidence <- data.frame(
    lane = c(
      "Gaussian distribution semantics", "SHASHb parameter conversion",
      "matched fitted estimators", "skew-heavy held-out utility",
      "HBR posterior aggregation", "new-site transport"
    ),
    evidence_class = c(
      "exact", "parameter_converted", "matched_estimator", "same_estimand",
      "non_equivalent", "non_equivalent"
    ),
    status = c(
      if (gaussian_error <= 1e-10) "pass" else "fail",
      if (shash_error <= 1e-7) "pass" else "fail",
      if (all(pcn_concordance_pass(concordance))) "pass" else "fail",
      superiority$classification,
      "release_evidence_required", "release_evidence_required"
    ),
    primary_quantity = c(
      "maximum absolute CDF error", "maximum absolute log-density error",
      "matched scenarios passing all margins", "mean paired log-score difference",
      "draw aggregation semantics", "site-stratified calibration"
    ),
    value = c(
      gaussian_error, shash_error, sum(pcn_concordance_pass(concordance)),
      superiority$log_score_difference, NA_real_, NA_real_
    ),
    threshold_or_interval = c(
      "<= 1e-10", "<= 1e-7",
      paste0(sum(matched), " of ", sum(matched)),
      sprintf("95%% CI %.6f to %.6f", superiority$ci_lower, superiority$ci_upper),
      "R-hat <= 1.01; ESS bulk/tail >= 400; zero divergences",
      "every site must pass; pooled metrics are insufficient"
    ),
    comparator_version = manifest$comparator$version,
    stringsAsFactors = FALSE
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(concordance, file.path(output_dir, "concordance.csv"), row.names = FALSE)
  utils::write.csv(superiority, file.path(output_dir, "superiority.csv"), row.names = FALSE)
  utils::write.csv(evidence, file.path(output_dir, "evidence_summary.csv"), row.names = FALSE)
  invisible(evidence)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 3L ||
      (length(args) >= 2L && !identical(args[[2L]], "--fixture-root")) ||
      length(args) == 2L) {
    stop("usage: render_evidence.R [OUTPUT_DIR [--fixture-root FIXTURE_DIR]]")
  }
  output <- if (length(args)) args[[1L]] else {
    file.path("docs", "evidence", "pcntoolkit", "v1.3.0")
  }
  fixture_root <- if (length(args) == 3L) args[[3L]] else pcn_fixture_root()
  print(render_pcntoolkit_evidence(output, fixture_root), row.names = FALSE)
}
