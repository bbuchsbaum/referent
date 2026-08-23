# Run independent SHASH fits against independently generated PCNtoolkit
# predictions. This is a release/manual gate, not an ordinary package test.

source("tools/pcntoolkit/compare_results.R", local = TRUE)
source("tools/pcntoolkit/run_concordance.R", local = TRUE)
source("tools/pcntoolkit/benchmark_helpers.R", local = TRUE)

pcn_validate_release_generation <- function(input_dir) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required for release evidence", call. = FALSE)
  }
  receipt_path <- file.path(input_dir, "release_generation_receipt.json")
  if (!file.exists(receipt_path)) {
    stop("missing release generation receipt", call. = FALSE)
  }
  receipt <- jsonlite::read_json(receipt_path, simplifyVector = FALSE)
  if (!identical(receipt$schema_version, "1.0.0") ||
      !identical(receipt$comparator$package, "pcntoolkit") ||
      !identical(receipt$comparator$version, "1.3.0") ||
      !identical(receipt$scenario, "skew_heavy") ||
      !identical(receipt$selection_rule,
                 "prespecified SHASH after family validation")) {
    stop("release generation receipt does not match the declared contract", call. = FALSE)
  }
  required <- c("release_inputs.csv", "release_pcntoolkit_predictions.csv")
  if (!setequal(names(receipt$files), required)) {
    stop("release generation receipt has an unexpected file set", call. = FALSE)
  }
  for (filename in required) {
    path <- file.path(input_dir, filename)
    if (!file.exists(path)) stop("missing release evidence file: ", filename, call. = FALSE)
    recorded <- receipt$files[[filename]]
    if (!identical(unname(tools::md5sum(path)), recorded$md5)) {
      stop("release evidence hash mismatch: ", filename, call. = FALSE)
    }
    dat <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!identical(nrow(dat), as.integer(recorded$rows)) ||
        !identical(names(dat), unlist(recorded$columns, use.names = FALSE))) {
      stop("release evidence schema mismatch: ", filename, call. = FALSE)
    }
  }
  invisible(receipt)
}

pcn_run_release_superiority <- function(input_dir, output_dir,
                                        within_B = 499L,
                                        aggregate_B = 4999L,
                                        n_draw = 2000L,
                                        seed = 20260823L,
                                        expected_replicates = 20L) {
  receipt <- pcn_validate_release_generation(input_dir)
  inputs <- utils::read.csv(
    file.path(input_dir, "release_inputs.csv"), stringsAsFactors = FALSE
  )
  comparator <- utils::read.csv(
    file.path(input_dir, "release_pcntoolkit_predictions.csv"),
    stringsAsFactors = FALSE
  )
  required_inputs <- c("replicate", "seed", "scenario", "row_id", "split", "site", "y")
  required_predictions <- c(
    "replicate", "seed", "scenario", "row_id", "observed", "log_density",
    "centile", "q05", "q95"
  )
  if (!all(required_inputs %in% names(inputs)) ||
      !all(required_predictions %in% names(comparator))) {
    stop("release benchmark tables are missing required columns", call. = FALSE)
  }
  replicates <- sort(unique(inputs$replicate))
  receipt_seeds <- sort(as.integer(unlist(receipt$seeds, use.names = FALSE)))
  observed_seeds <- sort(unique(inputs$seed))
  declared_splits <- c(train = 500L, validation = 100L, test = 500L)
  receipt_splits <- vapply(names(declared_splits), function(split) {
    as.integer(receipt$split_sizes[[split]])
  }, integer(1))
  split_counts <- table(inputs$replicate, inputs$split)
  if (length(replicates) != expected_replicates ||
      as.integer(receipt$replicates) != expected_replicates ||
      !identical(receipt_seeds, observed_seeds) ||
      !identical(receipt_splits, declared_splits) ||
      !all(names(declared_splits) %in% colnames(split_counts)) ||
      any(vapply(names(declared_splits), function(split) {
        any(split_counts[, split] != declared_splits[[split]])
      }, logical(1))) ||
      !identical(replicates, sort(unique(comparator$replicate))) ||
      anyDuplicated(inputs$row_id) || anyDuplicated(comparator$row_id)) {
    stop("release benchmark replicate or row identities are invalid", call. = FALSE)
  }

  rows <- lapply(replicates, function(replicate_id) {
    replicate_inputs <- inputs[inputs$replicate == replicate_id, , drop = FALSE]
    replicate_seed <- unique(replicate_inputs$seed)
    if (length(replicate_seed) != 1L) {
      stop("release replicate has multiple generation seeds: ", replicate_id,
           call. = FALSE)
    }
    fitted <- pcn_referent_predictions(
      replicate_inputs, "skew_heavy", n_draw = n_draw, return_fit = TRUE
    )
    referent <- fitted$predictions
    pcn <- comparator[comparator$replicate == replicate_id, , drop = FALSE]
    pcn <- pcn[match(referent$row_id, pcn$row_id), , drop = FALSE]
    if (anyNA(pcn$row_id) || !identical(referent$row_id, pcn$row_id)) {
      stop("release prediction row identity mismatch in replicate ", replicate_id,
           call. = FALSE)
    }
    if (any(!is.finite(pcn$observed)) ||
        max(abs(referent$observed - pcn$observed)) > 1e-12) {
      stop("release observed outcomes differ in replicate ", replicate_id,
           call. = FALSE)
    }
    comparison <- pcn_compare_predictions(
      referent, pcn, cluster = referent$row_id, scenario = "skew_heavy",
      B = within_B, seed = seed + as.integer(replicate_id)
    )
    fit_one <- fitted$fit$models$y
    comparison$replicate <- replicate_id
    comparison$seed <- replicate_seed
    comparison$fit_status <- fit_one$status
    comparison$fit_converged <- identical(fit_one$status, "ok") &&
      !is.null(fit_one$model) && !identical(fit_one$model$converged, FALSE)
    comparison
  })
  replicate_results <- do.call(rbind, rows)
  summary <- pcn_release_superiority_summary(
    replicate_results, B = aggregate_B, seed = seed
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  replicate_path <- file.path(output_dir, "release_superiority_replicates.csv")
  summary_path <- file.path(output_dir, "release_superiority_summary.csv")
  utils::write.csv(replicate_results, replicate_path, row.names = FALSE)
  utils::write.csv(summary, summary_path, row.names = FALSE)
  release_receipt <- list(
    schema_version = "1.0.0",
    comparator = receipt$comparator,
    scenario = "skew_heavy",
    generated_replicates = length(replicates),
    expected_replicates = expected_replicates,
    within_replicate_bootstrap = within_B,
    aggregate_replicate_bootstrap = aggregate_B,
    referent_prediction_draws = n_draw,
    requested_seed = seed,
    primary_contrast = "mean paired response-scale log score by simulation replicate",
    calibration_rule = paste(
      "replicate-bootstrap upper regret <= 0.01 for coverage, MACE, and tail;",
      "no single-replicate regret > 0.05"
    ),
    all_shash_fits_converged = summary$all_shash_fits_converged,
    calibration_noninferior = summary$calibration_noninferior,
    no_critical_replicate_regression = summary$no_critical_replicate_regression,
    classification = summary$classification,
    pass = summary$pass,
    files = list(
      release_superiority_replicates.csv = unname(tools::md5sum(replicate_path)),
      release_superiority_summary.csv = unname(tools::md5sum(summary_path))
    )
  )
  jsonlite::write_json(
    release_receipt,
    file.path(output_dir, "release_superiority_receipt.json"),
    auto_unbox = TRUE, pretty = TRUE
  )
  invisible(list(replicates = replicate_results, summary = summary))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L || length(args) > 6L) {
    stop(paste(
      "usage: run_release_superiority.R INPUT_DIR OUTPUT_DIR",
      "[WITHIN_BOOTSTRAP] [AGGREGATE_BOOTSTRAP] [N_DRAW] [EXPECTED_REPLICATES]"
    ))
  }
  within_B <- if (length(args) >= 3L) as.integer(args[[3L]]) else 499L
  aggregate_B <- if (length(args) >= 4L) as.integer(args[[4L]]) else 4999L
  n_draw <- if (length(args) >= 5L) as.integer(args[[5L]]) else 2000L
  expected_replicates <- if (length(args) >= 6L) as.integer(args[[6L]]) else 20L
  if (anyNA(c(within_B, aggregate_B, n_draw, expected_replicates)) ||
      within_B < 99L || aggregate_B < 99L || n_draw < 100L ||
      expected_replicates < 5L) {
    stop("bootstrap counts and prediction draws are below the release minimum")
  }
  result <- pcn_run_release_superiority(
    args[[1L]], args[[2L]], within_B = within_B,
    aggregate_B = aggregate_B, n_draw = n_draw,
    expected_replicates = expected_replicates
  )
  print(result$summary, row.names = FALSE)
  if (!isTRUE(result$summary$pass)) quit(status = 1L)
}
