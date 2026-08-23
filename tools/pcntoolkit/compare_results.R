# Utilities for consuming the checked-in PCNtoolkit comparison fixture.
# This file is deliberately outside R/ because it is validation infrastructure,
# not part of referent's user-facing API.

pcn_fixture_root <- function(version = "1.3.0", package_root = ".") {
  file.path(package_root, "tests", "fixtures", "pcntoolkit", paste0("v", version))
}

pcn_read_manifest <- function(root) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to validate PCNtoolkit fixtures", call. = FALSE)
  }
  path <- file.path(root, "manifest.json")
  if (!file.exists(path)) stop("missing PCNtoolkit manifest: ", path, call. = FALSE)
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

pcn_validate_fixture <- function(root, expected_version = "1.3.0") {
  manifest <- pcn_read_manifest(root)
  if (!identical(manifest$schema_version, "1.0.0")) {
    stop("unsupported PCNtoolkit fixture schema: ", manifest$schema_version, call. = FALSE)
  }
  if (!identical(manifest$comparator$package, "pcntoolkit") ||
      !identical(manifest$comparator$version, expected_version)) {
    stop("fixture comparator/version does not match the declared target", call. = FALSE)
  }
  required <- c(
    "semantic_gaussian.csv", "semantic_shashb.csv", "semantic_metrics.csv",
    "fitted_inputs.csv", "fitted_predictions.csv"
  )
  if (!setequal(names(manifest$files), required)) {
    stop("fixture manifest has an unexpected file set", call. = FALSE)
  }
  for (filename in required) {
    receipt <- manifest$files[[filename]]
    path <- file.path(root, filename)
    if (!file.exists(path)) stop("missing fixture file: ", filename, call. = FALSE)
    actual_md5 <- unname(tools::md5sum(path))
    if (!identical(actual_md5, receipt$md5)) {
      stop("fixture hash mismatch: ", filename, call. = FALSE)
    }
    dat <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!identical(nrow(dat), as.integer(receipt$rows))) {
      stop("fixture row-count mismatch: ", filename, call. = FALSE)
    }
    if (!identical(names(dat), unlist(receipt$columns, use.names = FALSE))) {
      stop("fixture column mismatch: ", filename, call. = FALSE)
    }
  }
  invisible(manifest)
}

pcn_read_fixture <- function(root, filename) {
  pcn_validate_fixture(root)
  utils::read.csv(file.path(root, filename), check.names = FALSE,
                  stringsAsFactors = FALSE)
}

pcn_concordance <- function(referent, comparator, outcome_scale) {
  stopifnot(nrow(referent) == nrow(comparator), outcome_scale > 0)
  delta_median <- referent$median - comparator$median
  delta_sd <- referent$predictive_sd - comparator$predictive_sd
  delta_z <- referent$z - comparator$z
  c(
    median_nrmse = sqrt(mean(delta_median^2)) / outcome_scale,
    predictive_sd_nrmse = sqrt(mean(delta_sd^2)) / outcome_scale,
    z_rmse = sqrt(mean(delta_z^2)),
    z_abs_q95 = unname(stats::quantile(abs(delta_z), 0.95, names = FALSE)),
    coverage90_difference = abs(
      mean(referent$observed >= referent$q05 & referent$observed <= referent$q95) -
        mean(comparator$observed >= comparator$q05 & comparator$observed <= comparator$q95)
    )
  )
}
