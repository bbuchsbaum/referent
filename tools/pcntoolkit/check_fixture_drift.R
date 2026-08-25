# Compare a freshly generated fixture with the reviewed baseline. Exact bytes
# are not required across BLAS/libm platforms; schemas, environments, row
# identities, and numerical values are checked at much tighter tolerances than
# any Referent-to-PCNtoolkit fitted concordance claim.

source("tools/pcntoolkit/compare_results.R")

pcn_compare_fixture_file <- function(baseline, candidate, filename,
                                     atol, rtol = atol) {
  a <- utils::read.csv(file.path(baseline, filename), check.names = FALSE,
                       stringsAsFactors = FALSE)
  b <- utils::read.csv(file.path(candidate, filename), check.names = FALSE,
                       stringsAsFactors = FALSE)
  if (!identical(names(a), names(b)) || nrow(a) != nrow(b)) {
    stop("fixture schema/row drift: ", filename, call. = FALSE)
  }
  numeric <- vapply(a, is.numeric, logical(1))
  if (!identical(numeric, vapply(b, is.numeric, logical(1)))) {
    stop("fixture column-type drift: ", filename, call. = FALSE)
  }
  if (any(!vapply(names(a)[!numeric], function(nm) identical(a[[nm]], b[[nm]]), logical(1)))) {
    stop("fixture identity drift: ", filename, call. = FALSE)
  }
  aa <- as.matrix(a[numeric])
  bb <- as.matrix(b[numeric])
  scale <- atol + rtol * abs(aa)
  scaled_error <- abs(aa - bb) / scale
  if (any(!is.finite(scaled_error)) || any(scaled_error > 1)) {
    stop(
      "fixture numerical drift exceeds tolerance in ", filename,
      "; maximum scaled error = ", max(scaled_error, na.rm = TRUE),
      call. = FALSE
    )
  }
  data.frame(file = filename, max_scaled_error = max(scaled_error),
             atol = atol, rtol = rtol)
}

pcn_check_fixture_drift <- function(baseline, candidate) {
  a <- pcn_validate_fixture(baseline)
  b <- pcn_validate_fixture(candidate)
  stable_manifest_fields <- c(
    "schema_version", "comparator", "python", "dependencies", "seeds",
    "centiles", "scenario_classes", "normalisations"
  )
  for (field in stable_manifest_fields) {
    if (!identical(a[[field]], b[[field]])) {
      stop("fixture manifest drift: ", field, call. = FALSE)
    }
  }
  do.call(rbind, list(
    pcn_compare_fixture_file(baseline, candidate, "semantic_gaussian.csv", 1e-12),
    pcn_compare_fixture_file(baseline, candidate, "semantic_shashb.csv", 1e-10),
    pcn_compare_fixture_file(baseline, candidate, "semantic_metrics.csv", 1e-10),
    pcn_compare_fixture_file(baseline, candidate, "fitted_inputs.csv", 1e-14),
    pcn_compare_fixture_file(baseline, candidate, "fitted_predictions.csv", 5e-6)
  ))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop("usage: check_fixture_drift.R BASELINE_DIR CANDIDATE_DIR")
  }
  print(pcn_check_fixture_drift(args[[1L]], args[[2L]]), row.names = FALSE)
}
