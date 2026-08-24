# Compare Referent adaptation with the separately defined PCNtoolkit HBR
# transfer evidence. Usage:
# Rscript tools/pcntoolkit/run_site_transport_evidence.R INPUT_DIR OUTPUT_DIR

source("tools/pcntoolkit/benchmark_helpers.R")
source("tools/pcntoolkit/site_evidence.R")

pcn_fit_site_model <- function(data) {
  data$site <- factor(data$site)
  referent::ref_fit(
    referent::ref_spec(
      referent::ref_gaussian(),
      location = ~ s(x, k = 6) + s(site, bs = "re"),
      scale = ~ s(site, bs = "re")
    ),
    data, outcomes = "y"
  )
}

pcn_validate_adaptation_priors <- function(reference, grid = pcn_adaptation_grid()) {
  sites <- sort(unique(as.character(reference$site)))
  if (length(sites) < 3L) {
    stop("adaptation-prior validation requires at least three observed sites",
         call. = FALSE)
  }
  rows <- list()
  fold_sizes <- list()
  for (held_site in sites) {
    pool <- reference[as.character(reference$site) != held_site, , drop = FALSE]
    held <- reference[as.character(reference$site) == held_site, , drop = FALSE]
    held <- held[order(held$row_id), , drop = FALSE]
    n_adapt <- floor(nrow(held) / 2L)
    if (n_adapt < 5L || nrow(held) - n_adapt < 5L) {
      stop("each adaptation-prior validation site needs at least ten rows",
           call. = FALSE)
    }
    local <- held[seq_len(n_adapt), , drop = FALSE]
    validation <- held[seq.int(n_adapt + 1L, nrow(held)), , drop = FALSE]
    local$site <- factor(local$site)
    validation$site <- factor(validation$site)
    fit <- pcn_fit_site_model(pool)
    fold_sizes[[held_site]] <- list(
      adaptation = nrow(local), validation = nrow(validation)
    )
    for (candidate in seq_len(nrow(grid))) {
      adapted <- referent::ref_adapt(
        fit, local, by = site, parameters = c("location", "scale"),
        location_prior_n = grid$location_prior_n[[candidate]],
        scale_prior_n = grid$scale_prior_n[[candidate]]
      )
      score <- stats::predict(
        adapted, validation, uncertainty = "conditional",
        allow_extrapolation = TRUE
      )
      rows[[length(rows) + 1L]] <- data.frame(
        site = held_site,
        adaptation_n = nrow(local),
        validation_n = nrow(validation),
        location_prior_n = grid$location_prior_n[[candidate]],
        scale_prior_n = grid$scale_prior_n[[candidate]],
        mean_log_score = mean(score$log_density),
        stringsAsFactors = FALSE
      )
    }
  }
  validation <- do.call(rbind, rows)
  selection <- pcn_select_adaptation_candidate(validation)
  list(
    selected = selection$selected,
    summary = selection$summary,
    validation = validation,
    source_sites = sites,
    fold_sizes = fold_sizes
  )
}

run_referent_site_evidence <- function(input_dir, output_dir) {
  data <- utils::read.csv(file.path(input_dir, "site_data.csv"), stringsAsFactors = FALSE)
  train <- data[data$split == "reference_train", , drop = FALSE]
  observed <- data[data$split == "observed_site_test", , drop = FALSE]
  adaptation <- data[data$split == "adapt", , drop = FALSE]
  transported <- data[data$split == "transport_test", , drop = FALSE]
  train$site <- factor(train$site)
  observed$site <- factor(observed$site, levels = levels(train$site))
  adaptation$site <- factor(adaptation$site)
  transported$site <- factor(transported$site)

  selection <- pcn_validate_adaptation_priors(train)
  fit <- pcn_fit_site_model(train)
  adapted <- referent::ref_adapt(
    fit, adaptation, by = site, parameters = c("location", "scale"),
    location_prior_n = selection$selected$location_prior_n,
    scale_prior_n = selection$selected$scale_prior_n
  )

  score_rows <- function(model, frame, lane) {
    score <- stats::predict(
      model, frame, uncertainty = "total", allow_extrapolation = TRUE
    )
    dist <- stats::predict(
      model, frame, type = "distribution", uncertainty = "total",
      allow_extrapolation = TRUE
    )$y
    data.frame(
      method = "referent_adaptation", lane = lane, row_id = frame$row_id,
      site = as.character(frame$site), observed = frame$y,
      centile = score$centile, log_density = score$log_density,
      q05 = referent:::dist_quantile(dist, 0.05),
      q50 = referent:::dist_quantile(dist, 0.5),
      q95 = referent:::dist_quantile(dist, 0.95),
      stringsAsFactors = FALSE
    )
  }
  predictions <- rbind(
    score_rows(fit, observed, "observed_site"),
    score_rows(adapted, transported, "transferred_site")
  )
  ref_summary <- do.call(rbind, lapply(split(predictions, predictions$lane), function(d) {
    pcn_site_summary(d, "referent_adaptation", unique(d$lane))
  }))
  pcn_summary <- utils::read.csv(
    file.path(input_dir, "pcntoolkit_site_summary.csv"), stringsAsFactors = FALSE
  )
  if (!"coverage90_se" %in% names(pcn_summary)) {
    pcn_summary$coverage90_se <- sqrt(0.9 * 0.1 / pcn_summary$n)
  }
  summary <- rbind(
    ref_summary[, names(pcn_summary), drop = FALSE],
    pcn_summary
  )
  gated <- pcn_site_gate(summary)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(predictions, file.path(output_dir, "referent_predictions.csv"), row.names = FALSE)
  utils::write.csv(gated, file.path(output_dir, "site_comparison.csv"), row.names = FALSE)
  selection_validation_path <- file.path(output_dir, "adaptation_selection.csv")
  selection_summary_path <- file.path(output_dir, "adaptation_selection_summary.csv")
  utils::write.csv(selection$validation, selection_validation_path, row.names = FALSE)
  utils::write.csv(selection$summary, selection_summary_path, row.names = FALSE)
  receipt <- list(
    estimand = "separate adaptation and transfer policies; no rowwise parity claim",
    referent_version = as.character(utils::packageVersion("referent")),
    r_version = R.version.string,
    all_sites_pass = all(gated$pass),
    failed_sites = split(gated$site[!gated$pass], gated$method[!gated$pass]),
    adaptation_prior_selection = list(
      rule = "leave-one-observed-site-out mean conditional log score",
      source_split = "reference_train",
      source_sites = selection$source_sites,
      fold_sizes = selection$fold_sizes,
      selected = unclass(selection$selected),
      files = list(
        adaptation_selection.csv = unname(tools::md5sum(selection_validation_path)),
        adaptation_selection_summary.csv = unname(tools::md5sum(selection_summary_path))
      )
    ),
    referent_adaptation = unclass(adapted$adaptation)
  )
  jsonlite::write_json(
    receipt, file.path(output_dir, "referent_receipt.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null"
  )
  invisible(gated)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("usage: run_site_transport_evidence.R INPUT_DIR OUTPUT_DIR")
  run_referent_site_evidence(args[[1L]], args[[2L]])
}
