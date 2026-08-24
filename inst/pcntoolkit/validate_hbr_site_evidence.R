# Independently validate retained HBR and site-transfer evidence. Receipt
# booleans are never accepted as verdicts: hashes, schemas, row identities,
# convergence summaries, site metrics, prior selection, and gates are all
# recomputed from the retained tables.

pcn_evidence_sha256 <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest is required to validate evidence", call. = FALSE)
  }
  digest::digest(file = path, algo = "sha256")
}

pcn_evidence_read_csv <- function(path) {
  if (!file.exists(path)) {
    stop("missing evidence file: ", path, call. = FALSE)
  }
  utils::read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("NA", "")
  )
}

pcn_evidence_file_receipt <- function(path) {
  dat <- pcn_evidence_read_csv(path)
  list(
    sha256 = pcn_evidence_sha256(path),
    rows = nrow(dat),
    columns = names(dat)
  )
}

pcn_evidence_blob_receipt <- function(path) {
  if (!file.exists(path)) {
    stop("missing evidence file: ", path, call. = FALSE)
  }
  list(
    sha256 = pcn_evidence_sha256(path),
    bytes = unname(file.info(path)$size)
  )
}

pcn_evidence_validate_files <- function(root, records, expected,
                                        aliases = character()) {
  if (!setequal(names(records), expected)) {
    stop("evidence receipt has an unexpected file set", call. = FALSE)
  }
  for (name in expected) {
    filename <- if (name %in% names(aliases)) aliases[[name]] else name
    path <- file.path(root, filename)
    recorded <- records[[name]]
    if (!is.null(recorded$bytes)) {
      actual <- pcn_evidence_blob_receipt(path)
      valid <- identical(actual$sha256, recorded$sha256) &&
        identical(actual$bytes, as.numeric(recorded$bytes))
    } else {
      actual <- pcn_evidence_file_receipt(path)
      valid <- identical(actual$sha256, recorded$sha256) &&
        identical(actual$rows, as.integer(recorded$rows)) &&
        identical(actual$columns, unlist(recorded$columns, use.names = FALSE))
    }
    if (!valid) {
      stop("evidence file receipt mismatch: ", name, call. = FALSE)
    }
  }
  invisible(TRUE)
}

pcn_evidence_assert_columns <- function(dat, required, label) {
  if (!all(required %in% names(dat))) {
    stop(label, " is missing required columns", call. = FALSE)
  }
}

pcn_evidence_assert_exact_columns <- function(dat, expected, label) {
  if (!identical(names(dat), expected)) {
    stop(label, " does not match the registered schema", call. = FALSE)
  }
}

pcn_evidence_assert_finite <- function(dat, columns, label) {
  if (any(!vapply(dat[columns], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(dat[columns])))) {
    stop(label, " contains non-finite numeric evidence", call. = FALSE)
  }
}

pcn_evidence_assert_frame <- function(actual, expected, keys, label,
                                      tolerance = 1e-11) {
  if (!setequal(names(actual), names(expected))) {
    stop(label, " schema does not match recomputed evidence", call. = FALSE)
  }
  actual <- actual[, names(expected), drop = FALSE]
  order_rows <- function(dat) {
    values <- unname(lapply(dat[keys], as.character))
    do.call(order, c(values, list(na.last = TRUE)))
  }
  actual <- actual[order_rows(actual), , drop = FALSE]
  expected <- expected[order_rows(expected), , drop = FALSE]
  rownames(actual) <- rownames(expected) <- NULL
  if (nrow(actual) != nrow(expected)) {
    stop(label, " row count does not match recomputed evidence", call. = FALSE)
  }
  for (name in names(expected)) {
    a <- actual[[name]]
    b <- expected[[name]]
    if (!identical(is.na(a), is.na(b))) {
      stop(label, " NA pattern differs in ", name, call. = FALSE)
    }
    keep <- !is.na(a)
    if (is.numeric(a) && is.numeric(b)) {
      scale <- pmax(1, abs(a[keep]), abs(b[keep]))
      if (any(abs(a[keep] - b[keep]) > tolerance * scale)) {
        stop(label, " differs from recomputed values in ", name, call. = FALSE)
      }
    } else if (!identical(as.character(a[keep]), as.character(b[keep]))) {
      stop(label, " differs from recomputed values in ", name, call. = FALSE)
    }
  }
  invisible(TRUE)
}

pcn_evidence_mace <- function(centile) {
  grid <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  mean(vapply(grid, function(p) abs(p - mean(centile <= p)), numeric(1)))
}

pcn_evidence_site_summary <- function(prediction, method) {
  required <- c(
    "lane", "site", "observed", "centile", "log_density", "q05", "q50", "q95"
  )
  pcn_evidence_assert_columns(prediction, required, "prediction table")
  numeric <- c("observed", "centile", "log_density", "q05", "q50", "q95")
  pcn_evidence_assert_finite(prediction, numeric, "prediction table")
  if (any(prediction$centile <= 0 | prediction$centile >= 1)) {
    stop("prediction centiles must lie strictly inside (0, 1)", call. = FALSE)
  }
  groups <- split(
    seq_len(nrow(prediction)),
    interaction(prediction$lane, prediction$site, drop = TRUE, lex.order = TRUE)
  )
  rows <- lapply(groups, function(index) {
    d <- prediction[index, , drop = FALSE]
    z <- stats::qnorm(d$centile)
    data.frame(
      method = method,
      lane = unique(d$lane),
      site = unique(d$site),
      n = nrow(d),
      mean_log_score = mean(d$log_density),
      median_rmse = sqrt(mean((d$observed - d$q50)^2)),
      coverage90 = mean(d$observed >= d$q05 & d$observed <= d$q95),
      mace = pcn_evidence_mace(d$centile),
      mean_z = mean(z),
      var_z = stats::var(z),
      coverage90_se = sqrt(0.9 * 0.1 / nrow(d)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

pcn_evidence_site_gate <- function(summary) {
  summary$pass <- summary$coverage90 >= 0.80 & summary$mace <= 0.08 &
    abs(summary$mean_z) <= 0.35
  summary
}

pcn_evidence_as_logical <- function(x, label) {
  if (is.logical(x)) return(x)
  if (is.numeric(x) && all(x %in% c(0, 1))) return(as.logical(x))
  value <- tolower(as.character(x))
  if (!all(value %in% c("true", "false"))) {
    stop(label, " is not a valid logical column", call. = FALSE)
  }
  value == "true"
}

pcn_validate_hbr_evidence <- function(root, receipt_file = "receipt.json") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to validate evidence", call. = FALSE)
  }
  receipt_path <- file.path(root, receipt_file)
  if (!file.exists(receipt_path)) {
    stop("missing HBR receipt: ", receipt_path, call. = FALSE)
  }
  receipt <- jsonlite::read_json(receipt_path, simplifyVector = FALSE)
  thresholds <- receipt$convergence_thresholds
  undefined_rhat_policy <- paste(
    "allowed only for retained posterior-summary rows with sd exactly zero"
  )
  aggregation <- receipt$aggregation
  if (!identical(receipt$schema_version, "1.1.0") ||
      !identical(receipt$comparator, "pcntoolkit") ||
      !identical(receipt$version, "1.3.0") ||
      !identical(as.integer(receipt$draws), 1500L) ||
      !identical(as.integer(receipt$tune), 1000L) ||
      !identical(as.integer(receipt$chains), 4L) ||
      !identical(as.integer(receipt$cores), 4L) ||
      !identical(receipt$nuts_sampler, "nutpie") ||
      !isTRUE(all.equal(receipt$target_accept, 0.99)) ||
      !identical(as.integer(receipt$requested_seed), 20260823L) ||
      !identical(as.integer(receipt$sampling_seeds$base), 20260823L) ||
      !identical(as.integer(receipt$sampling_seeds$transfer), 20260824L) ||
      !identical(receipt$hbr_public_random_seed_argument, FALSE) ||
      !isTRUE(all.equal(thresholds$max_rhat, 1.01)) ||
      !identical(as.integer(thresholds$min_ess_bulk), 400L) ||
      !identical(as.integer(thresholds$min_ess_tail), 400L) ||
      !identical(as.integer(thresholds$divergences), 0L) ||
      !identical(
        receipt$convergence_policy$undefined_rhat, undefined_rhat_policy
      ) ||
      !setequal(
        names(aggregation),
        c("reported_z", "reported_quantile", "reported_logp", "mixture_outputs")
      ) ||
      !identical(aggregation$reported_z, "mean of draw-specific z") ||
      !identical(
        aggregation$reported_quantile, "mean of draw-specific quantiles"
      ) ||
      !identical(
        aggregation$reported_logp, "mean of draw-specific log densities"
      ) ||
      !identical(
        aggregation$mixture_outputs,
        "CDF and density averaged over posterior draws, then transformed"
      )) {
    stop("HBR receipt does not match the registered contract", call. = FALSE)
  }
  expected_files <- c(
    "site_data.csv", "pcntoolkit_predictions.csv",
    "pcntoolkit_site_summary.csv", "hbr_convergence.csv",
    "hbr_divergences.csv"
  )
  pcn_evidence_validate_files(root, receipt$files, expected_files)

  convergence <- pcn_evidence_read_csv(file.path(root, "hbr_convergence.csv"))
  divergences <- pcn_evidence_read_csv(file.path(root, "hbr_divergences.csv"))
  pcn_evidence_assert_exact_columns(
    convergence,
    c(
      "stage", "parameter", "mean", "sd", "hdi_3%", "hdi_97%",
      "mcse_mean", "mcse_sd", "ess_bulk", "ess_tail", "r_hat"
    ),
    "HBR convergence table"
  )
  pcn_evidence_assert_exact_columns(
    divergences, c("stage", "chain", "draw", "diverging"),
    "HBR divergence table"
  )
  pcn_evidence_assert_finite(
    convergence, c("sd", "ess_bulk", "ess_tail"), "HBR convergence table"
  )
  pcn_evidence_assert_finite(
    divergences, c("chain", "draw"), "HBR divergence table"
  )
  if (anyDuplicated(convergence[c("stage", "parameter")]) ||
      anyDuplicated(divergences[c("stage", "chain", "draw")])) {
    stop("HBR diagnostic row identities are duplicated", call. = FALSE)
  }
  undefined_rhat <- is.na(convergence$r_hat)
  if (any(is.infinite(convergence$r_hat)) ||
      any(undefined_rhat & convergence$sd != 0)) {
    stop("HBR R-hat is undefined outside an exactly constant row",
         call. = FALSE)
  }
  divergences$diverging <- pcn_evidence_as_logical(
    divergences$diverging, "HBR diverging"
  )
  stages <- c("base", "transfer")
  if (!setequal(unique(convergence$stage), stages) ||
      !setequal(unique(divergences$stage), stages)) {
    stop("HBR diagnostic stages are incomplete", call. = FALSE)
  }
  expected_draws <- expand.grid(
    stage = stages,
    chain = 0:3,
    draw = 0:1499,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pcn_evidence_assert_frame(
    divergences[c("stage", "chain", "draw")], expected_draws,
    c("stage", "chain", "draw"), "HBR divergence draw identities"
  )
  recomputed_stages <- do.call(rbind, lapply(stages, function(stage) {
    conv <- convergence[convergence$stage == stage, , drop = FALSE]
    div <- divergences[divergences$stage == stage, , drop = FALSE]
    if (!any(is.finite(conv$r_hat))) {
      stop("HBR stage has no finite R-hat values: ", stage, call. = FALSE)
    }
    values <- data.frame(
      stage = stage,
      divergences = sum(div$diverging),
      max_rhat = max(conv$r_hat, na.rm = TRUE),
      min_ess_bulk = min(conv$ess_bulk),
      min_ess_tail = min(conv$ess_tail),
      undefined_rhat_constant_rows = sum(is.na(conv$r_hat) & conv$sd == 0),
      stringsAsFactors = FALSE
    )
    values$pass <- values$divergences == 0L && values$max_rhat <= 1.01 &&
      values$min_ess_bulk >= 400 && values$min_ess_tail >= 400
    values
  }))
  receipt_stages <- do.call(rbind, lapply(receipt$stages, function(stage) {
    data.frame(
      stage = stage$stage,
      divergences = as.integer(stage$divergences),
      max_rhat = stage$max_rhat,
      min_ess_bulk = stage$min_ess_bulk,
      min_ess_tail = stage$min_ess_tail,
      undefined_rhat_constant_rows = as.integer(
        stage$undefined_rhat_constant_rows
      ),
      pass = isTRUE(stage$pass),
      stringsAsFactors = FALSE
    )
  }))
  pcn_evidence_assert_frame(
    receipt_stages, recomputed_stages, "stage", "HBR stage receipt"
  )
  if (!identical(isTRUE(receipt$all_stages_pass), all(recomputed_stages$pass))) {
    stop("HBR receipt verdict does not match recomputed diagnostics", call. = FALSE)
  }

  data <- pcn_evidence_read_csv(file.path(root, "site_data.csv"))
  prediction <- pcn_evidence_read_csv(file.path(root, "pcntoolkit_predictions.csv"))
  pcn_evidence_assert_exact_columns(
    data, c("row_id", "site", "split", "x", "y"), "site data"
  )
  pcn_evidence_assert_exact_columns(
    prediction,
    c(
      "lane", "row_id", "site", "observed",
      "reported_mean_draw_z", "reported_mean_draw_log_density",
      "reported_mean_draw_q05", "reported_mean_draw_q50",
      "reported_mean_draw_q95", "mixture_centile", "mixture_z",
      "mixture_log_density", "mixture_q05", "mixture_q50", "mixture_q95"
    ),
    "PCNtoolkit prediction table"
  )
  pcn_evidence_assert_finite(data, c("x", "y"), "site data")
  pcn_evidence_assert_finite(
    prediction,
    c(
      "observed", "reported_mean_draw_z", "reported_mean_draw_log_density",
      "reported_mean_draw_q05", "reported_mean_draw_q50",
      "reported_mean_draw_q95", "mixture_centile", "mixture_z",
      "mixture_log_density", "mixture_q05", "mixture_q50", "mixture_q95"
    ),
    "PCNtoolkit prediction table"
  )
  if (anyDuplicated(data$row_id) || anyDuplicated(prediction$row_id)) {
    stop("HBR evidence has duplicate row identities", call. = FALSE)
  }
  expected_structure <- data.frame(
    site = c("site-1", "site-1", "site-2", "site-2", "site-3", "site-3",
             "site-4", "site-4", "site-5", "site-5"),
    split = rep(
      c("reference_train", "observed_site_test"), 4L
    ) |> c("adapt", "transport_test"),
    n = c(rep(c(50L, 30L), 4L), 25L, 75L),
    stringsAsFactors = FALSE
  )
  actual_structure <- stats::aggregate(
    list(n = data$row_id), list(site = data$site, split = data$split), length
  )
  pcn_evidence_assert_frame(
    actual_structure, expected_structure, c("site", "split"),
    "site data split structure"
  )
  if (any(!grepl(
    "^site-[1-5]-(reference_train|observed_site_test|adapt|transport_test)-[0-9]{4}$",
    data$row_id
  ))) {
    stop("site data row identities do not match the registered design",
         call. = FALSE)
  }
  evaluation <- data[data$split %in% c("observed_site_test", "transport_test"), ]
  evaluation$lane <- ifelse(
    evaluation$split == "observed_site_test", "observed_site", "transferred_site"
  )
  if (nrow(prediction) != nrow(evaluation) ||
      !setequal(prediction$row_id, evaluation$row_id)) {
    stop("HBR predictions do not cover every retained evaluation row",
         call. = FALSE)
  }
  evaluation <- evaluation[match(prediction$row_id, evaluation$row_id), ]
  if (anyNA(evaluation$row_id) ||
      !identical(prediction$row_id, evaluation$row_id) ||
      !identical(prediction$site, evaluation$site) ||
      !identical(prediction$lane, evaluation$lane) ||
      any(abs(prediction$observed - evaluation$y) > 1e-12)) {
    stop("HBR predictions are not bound to the retained site rows", call. = FALSE)
  }
  if (any(prediction$mixture_centile <= 0 | prediction$mixture_centile >= 1) ||
      any(prediction$mixture_q05 > prediction$mixture_q50) ||
      any(prediction$mixture_q50 > prediction$mixture_q95) ||
      any(abs(
        prediction$mixture_z - stats::qnorm(prediction$mixture_centile)
      ) > 1e-11 * pmax(1, abs(prediction$mixture_z)))) {
    stop("PCNtoolkit mixture predictions violate distribution invariants",
         call. = FALSE)
  }
  standardized <- data.frame(
    lane = prediction$lane,
    site = prediction$site,
    observed = prediction$observed,
    centile = prediction$mixture_centile,
    log_density = prediction$mixture_log_density,
    q05 = prediction$mixture_q05,
    q50 = prediction$mixture_q50,
    q95 = prediction$mixture_q95,
    stringsAsFactors = FALSE
  )
  recomputed_summary <- pcn_evidence_site_summary(
    standardized, "pcntoolkit_hbr_mixture"
  )
  recorded_summary <- pcn_evidence_read_csv(
    file.path(root, "pcntoolkit_site_summary.csv")
  )
  pcn_evidence_assert_frame(
    recorded_summary, recomputed_summary, c("method", "lane", "site"),
    "PCNtoolkit site summary"
  )
  invisible(list(
    receipt = receipt,
    stages = recomputed_stages,
    site_summary = recomputed_summary,
    data = data
  ))
}

pcn_evidence_selection_summary <- function(validation) {
  required <- c("site", "location_prior_n", "scale_prior_n", "mean_log_score")
  pcn_evidence_assert_columns(validation, required, "adaptation selection table")
  pcn_evidence_assert_finite(
    validation, c("location_prior_n", "scale_prior_n", "mean_log_score"),
    "adaptation selection table"
  )
  if (anyDuplicated(validation[c("site", "location_prior_n", "scale_prior_n")])) {
    stop("adaptation selection has duplicate site-candidate rows", call. = FALSE)
  }
  counts <- stats::aggregate(
    site ~ location_prior_n + scale_prior_n, validation,
    function(x) length(unique(x))
  )
  names(counts)[names(counts) == "site"] <- "validated_sites"
  expected_sites <- length(unique(validation$site))
  if (any(counts$validated_sites != expected_sites)) {
    stop("adaptation selection omits an observed site", call. = FALSE)
  }
  scores <- stats::aggregate(
    mean_log_score ~ location_prior_n + scale_prior_n, validation, mean
  )
  out <- merge(
    scores, counts, by = c("location_prior_n", "scale_prior_n"), sort = FALSE
  )
  out <- out[order(-out$mean_log_score, out$location_prior_n, out$scale_prior_n), ]
  rownames(out) <- NULL
  out$selected <- seq_len(nrow(out)) == 1L
  out
}

pcn_validate_site_evidence <- function(input_dir, output_dir,
                                       hbr_receipt_file = "receipt.json",
                                       site_receipt_file = "referent_receipt.json") {
  hbr <- pcn_validate_hbr_evidence(input_dir, hbr_receipt_file)
  receipt_path <- file.path(output_dir, site_receipt_file)
  if (!file.exists(receipt_path)) {
    stop("missing site receipt: ", receipt_path, call. = FALSE)
  }
  receipt <- jsonlite::read_json(receipt_path, simplifyVector = FALSE)
  thresholds <- receipt$gate_thresholds
  if (!identical(receipt$schema_version, "1.1.0") ||
      !identical(
        receipt$estimand,
        "separate adaptation and transfer policies; no rowwise parity claim"
      ) ||
      !identical(receipt$referent_version, "0.1.0") ||
      !isTRUE(all.equal(thresholds$coverage90_floor, 0.80)) ||
      !isTRUE(all.equal(thresholds$mace_ceiling, 0.08)) ||
      !isTRUE(all.equal(thresholds$mean_z_absolute_ceiling, 0.35))) {
    stop("site receipt does not match the registered contract", call. = FALSE)
  }
  hbr_files <- c(
    "site_data.csv", "pcntoolkit_predictions.csv",
    "pcntoolkit_site_summary.csv", "hbr_convergence.csv",
    "hbr_divergences.csv", "receipt.json"
  )
  pcn_evidence_validate_files(
    input_dir, receipt$input_files, hbr_files,
    aliases = stats::setNames(hbr_receipt_file, "receipt.json")
  )
  output_files <- c(
    "referent_predictions.csv", "site_comparison.csv",
    "adaptation_selection.csv", "adaptation_selection_summary.csv"
  )
  pcn_evidence_validate_files(output_dir, receipt$output_files, output_files)

  prediction <- pcn_evidence_read_csv(
    file.path(output_dir, "referent_predictions.csv")
  )
  pcn_evidence_assert_exact_columns(
    prediction,
    c("method", "lane", "row_id", "site", "observed", "centile",
      "log_density", "q05", "q50", "q95"),
    "Referent prediction table"
  )
  pcn_evidence_assert_finite(
    prediction,
    c("observed", "centile", "log_density", "q05", "q50", "q95"),
    "Referent prediction table"
  )
  if (anyDuplicated(prediction$row_id) ||
      any(prediction$method != "referent_adaptation")) {
    stop("Referent prediction identities are invalid", call. = FALSE)
  }
  evaluation <- hbr$data[
    hbr$data$split %in% c("observed_site_test", "transport_test"),
  ]
  evaluation$lane <- ifelse(
    evaluation$split == "observed_site_test", "observed_site", "transferred_site"
  )
  if (nrow(prediction) != nrow(evaluation) ||
      !setequal(prediction$row_id, evaluation$row_id)) {
    stop("Referent predictions do not cover every retained evaluation row",
         call. = FALSE)
  }
  evaluation <- evaluation[match(prediction$row_id, evaluation$row_id), ]
  if (anyNA(evaluation$row_id) ||
      !identical(prediction$row_id, evaluation$row_id) ||
      !identical(prediction$site, evaluation$site) ||
      !identical(prediction$lane, evaluation$lane) ||
      any(abs(prediction$observed - evaluation$y) > 1e-12)) {
    stop("Referent predictions are not bound to the retained site rows", call. = FALSE)
  }
  if (any(prediction$centile <= 0 | prediction$centile >= 1) ||
      any(prediction$q05 > prediction$q50) ||
      any(prediction$q50 > prediction$q95)) {
    stop("Referent predictions violate distribution invariants", call. = FALSE)
  }
  ref_summary <- pcn_evidence_site_summary(
    prediction, "referent_adaptation"
  )
  recomputed_comparison <- pcn_evidence_site_gate(rbind(
    ref_summary[, names(hbr$site_summary), drop = FALSE],
    hbr$site_summary
  ))
  recorded_comparison <- pcn_evidence_read_csv(
    file.path(output_dir, "site_comparison.csv")
  )
  pcn_evidence_assert_frame(
    recorded_comparison, recomputed_comparison,
    c("method", "lane", "site"), "site comparison"
  )

  validation <- pcn_evidence_read_csv(
    file.path(output_dir, "adaptation_selection.csv")
  )
  pcn_evidence_assert_exact_columns(
    validation,
    c(
      "site", "adaptation_n", "validation_n", "location_prior_n",
      "scale_prior_n", "mean_log_score"
    ),
    "adaptation selection table"
  )
  selection <- pcn_evidence_selection_summary(validation)
  recorded_selection <- pcn_evidence_read_csv(
    file.path(output_dir, "adaptation_selection_summary.csv")
  )
  pcn_evidence_assert_frame(
    recorded_selection, selection, c("location_prior_n", "scale_prior_n"),
    "adaptation selection summary"
  )
  reference_sites <- sort(unique(hbr$data$site[hbr$data$split == "reference_train"]))
  expected_candidates <- expand.grid(
    site = reference_sites,
    location_prior_n = c(0, 2, 5, 10),
    scale_prior_n = c(0, 5, 10, 25),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pcn_evidence_assert_frame(
    validation[c("site", "location_prior_n", "scale_prior_n")],
    expected_candidates, c("site", "location_prior_n", "scale_prior_n"),
    "adaptation candidate grid"
  )
  if (any(validation$adaptation_n != 25L) ||
      any(validation$validation_n != 25L) ||
      !identical(
        receipt$adaptation_prior_selection$rule,
        "leave-one-observed-site-out mean conditional log score"
      )) {
    stop("adaptation selection does not match the registered folds",
         call. = FALSE)
  }
  if (!identical(sort(unique(validation$site)), reference_sites) ||
      !identical(
        sort(unlist(receipt$adaptation_prior_selection$source_sites, use.names = FALSE)),
        reference_sites
      ) ||
      !identical(receipt$adaptation_prior_selection$source_split, "reference_train")) {
    stop("adaptation selection is not restricted to observed reference sites",
         call. = FALSE)
  }
  for (site in reference_sites) {
    fold <- receipt$adaptation_prior_selection$fold_sizes[[site]]
    if (!identical(as.integer(fold$adaptation), 25L) ||
        !identical(as.integer(fold$validation), 25L)) {
      stop("adaptation fold-size receipt does not match retained evidence",
           call. = FALSE)
    }
  }
  selected <- selection[selection$selected, c("location_prior_n", "scale_prior_n")]
  recorded_selected <- receipt$adaptation_prior_selection$selected
  if (nrow(selected) != 1L ||
      !isTRUE(all.equal(
        selected$location_prior_n[[1L]], recorded_selected$location_prior_n
      )) ||
      !isTRUE(all.equal(
        selected$scale_prior_n[[1L]], recorded_selected$scale_prior_n
      ))) {
    stop("adaptation-prior receipt does not match recomputed selection",
         call. = FALSE)
  }
  if (!identical(
    isTRUE(receipt$all_sites_pass), all(recomputed_comparison$pass)
  )) {
    stop("site receipt verdict does not match recomputed gates", call. = FALSE)
  }
  expected_failed <- split(
    as.character(recomputed_comparison$site[!recomputed_comparison$pass]),
    as.character(recomputed_comparison$method[!recomputed_comparison$pass])
  )
  recorded_failed <- lapply(
    receipt$failed_sites,
    function(x) sort(unlist(x, use.names = FALSE))
  )
  expected_failed <- lapply(expected_failed, sort)
  if (!setequal(names(recorded_failed), names(expected_failed)) ||
      any(vapply(names(expected_failed), function(name) {
        !identical(recorded_failed[[name]], expected_failed[[name]])
      }, logical(1)))) {
    stop("site receipt failed-site list does not match recomputed gates",
         call. = FALSE)
  }
  invisible(list(
    hbr = hbr,
    receipt = receipt,
    comparison = recomputed_comparison,
    selection = selection
  ))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L || length(args) > 4L) {
    stop(paste(
      "usage: validate_hbr_site_evidence.R INPUT_DIR OUTPUT_DIR",
      "[HBR_RECEIPT_FILE] [SITE_RECEIPT_FILE]"
    ), call. = FALSE)
  }
  hbr_receipt <- if (length(args) >= 3L) args[[3L]] else "receipt.json"
  site_receipt <- if (length(args) >= 4L) args[[4L]] else "referent_receipt.json"
  validated <- pcn_validate_site_evidence(
    args[[1L]], args[[2L]], hbr_receipt, site_receipt
  )
  cat(
    "validated HBR stages:", nrow(validated$hbr$stages),
    "; validated site gates:", nrow(validated$comparison), "\n"
  )
}
