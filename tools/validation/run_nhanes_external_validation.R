# Independent real-cohort validation using public CDC NHANES files.
#
# Train: 2015-2016 adults. External evaluation: 2017-2018 adults.
# This lane is deliberately unweighted: it validates transport between two
# independently sampled cohorts, not survey-weighted U.S. prevalence.

nhanes_files <- function() {
  c(
    demo_train = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/DEMO_I.XPT",
    body_train = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/BMX_I.XPT",
    demo_test = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.XPT",
    body_test = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BMX_J.XPT"
  )
}

nhanes_download <- function(cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  urls <- nhanes_files()
  paths <- file.path(cache_dir, basename(urls))
  for (i in seq_along(urls)) {
    if (!file.exists(paths[[i]]) || file.info(paths[[i]])$size == 0) {
      utils::download.file(urls[[i]], paths[[i]], mode = "wb", quiet = TRUE)
    }
  }
  stats::setNames(paths, names(urls))
}

nhanes_cohort <- function(demo_path, body_path, cohort) {
  demo <- foreign::read.xport(demo_path)
  body <- foreign::read.xport(body_path)
  keep_demo <- c("SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3")
  keep_body <- c("SEQN", "BMXHT", "BMXWT", "BMXBMI", "BMXWAIST")
  missing <- c(setdiff(keep_demo, names(demo)), setdiff(keep_body, names(body)))
  if (length(missing)) {
    stop("NHANES schema is missing: ", paste(unique(missing), collapse = ", "),
         call. = FALSE)
  }
  out <- merge(demo[keep_demo], body[keep_body], by = "SEQN", all = FALSE)
  out <- out[is.finite(out$RIDAGEYR) & out$RIDAGEYR >= 18 & out$RIDAGEYR <= 79, ]
  data.frame(
    participant_id = sprintf("%s-%d", cohort, as.integer(out$SEQN)),
    cohort = cohort,
    age = as.numeric(out$RIDAGEYR),
    sex = factor(out$RIAGENDR, levels = c(1, 2), labels = c("male", "female")),
    race_ethnicity = factor(out$RIDRETH3),
    height_cm = as.numeric(out$BMXHT),
    weight_kg = as.numeric(out$BMXWT),
    bmi = as.numeric(out$BMXBMI),
    waist_cm = as.numeric(out$BMXWAIST),
    stringsAsFactors = FALSE
  )
}

nhanes_external_gate <- function(fit, assessment) {
  statuses <- vapply(fit$models, function(x) x$status, character(1))
  marginal <- assessment$marginal
  calibration_pass <- referent:::acceptable_calibration(assessment)
  data.frame(
    all_fits_valid = all(statuses == "ok"),
    held_out = !isTRUE(assessment$in_sample),
    minimum_outcome_n = min(marginal$n),
    complete_evidence = all(is.finite(unlist(marginal[c(
      "mean_z", "var_z", "mace", "cover_95"
    )]))),
    calibration_pass = calibration_pass,
    pass = all(statuses == "ok") && !isTRUE(assessment$in_sample) &&
      min(marginal$n) >= 500L && calibration_pass,
    stringsAsFactors = FALSE
  )
}

nhanes_external_assess <- function(fit, data, n_draw) {
  scores <- stats::predict(
    fit, data, type = "scores", uncertainty = "total", n_draw = n_draw,
    allow_extrapolation = FALSE
  )
  by_outcome <- split(scores, scores$.outcome)
  overall <- do.call(rbind, lapply(names(by_outcome), function(outcome) {
    sc <- by_outcome[[outcome]]
    ok <- is.finite(sc$observed) & is.finite(sc$log_density)
    data.frame(
      .outcome = outcome,
      mean_log_score = mean(sc$log_density[ok]),
      mae = mean(abs(sc$residual[ok])),
      rmse = sqrt(mean(sc$residual[ok]^2)),
      stringsAsFactors = FALSE
    )
  }))
  structure(
    list(
      overall = overall,
      marginal = referent:::assess_marginal(scores),
      conditional = referent:::assess_conditional(scores, data, fit),
      tail = referent:::assess_tail(scores),
      scores = scores,
      n = nrow(data),
      in_sample = nrow(scores) > 0L && all(scores$.in_sample)
    ),
    class = "ref_assessment"
  )
}

run_nhanes_external_validation <- function(output_dir, cache_dir) {
  if (!requireNamespace("referent", quietly = TRUE) ||
      !requireNamespace("foreign", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE) ||
      !requireNamespace("digest", quietly = TRUE)) {
    stop("referent, foreign, jsonlite, and digest are required", call. = FALSE)
  }
  paths <- nhanes_download(cache_dir)
  train <- nhanes_cohort(paths[["demo_train"]], paths[["body_train"]], "2015-2016")
  test_population <- nhanes_cohort(
    paths[["demo_test"]], paths[["body_test"]], "2017-2018"
  )
  set.seed(20260824L)
  evaluation_n <- min(2000L, nrow(test_population))
  test <- test_population[sample.int(nrow(test_population), evaluation_n), , drop = FALSE]
  outcomes <- c("height_cm", "weight_kg", "bmi", "waist_cm")
  spec <- referent::ref_spec(
    referent::ref_shash(),
    location = ~ s(age, k = 8) + sex + race_ethnicity,
    scale = ~ s(age, k = 5) + sex + race_ethnicity,
    skew = ~ 1,
    tail = ~ 1,
    control = list(seed = 20260824L, n_draw = 1000L)
  )
  fit <- referent::ref_fit(spec, train, outcomes = outcomes, id = participant_id)
  assessment <- nhanes_external_assess(fit, test, n_draw = 1000L)
  gate <- nhanes_external_gate(fit, assessment)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  products <- list(
    overall = assessment$overall,
    marginal = assessment$marginal,
    conditional = assessment$conditional,
    tail = assessment$tail,
    gate = gate
  )
  output_paths <- vapply(names(products), function(name) {
    path <- file.path(output_dir, paste0(name, ".csv"))
    utils::write.csv(products[[name]], path, row.names = FALSE)
    path
  }, character(1))
  input_receipts <- lapply(seq_along(paths), function(i) {
    list(
      url = unname(nhanes_files()[[i]]),
      sha256 = digest::digest(file = paths[[i]], algo = "sha256"),
      bytes = unname(file.info(paths[[i]])$size)
    )
  })
  names(input_receipts) <- names(paths)
  output_receipts <- lapply(output_paths, function(path) {
    list(
      sha256 = digest::digest(file = path, algo = "sha256"),
      rows = nrow(utils::read.csv(path, check.names = FALSE))
    )
  })
  receipt <- list(
    schema_version = "1.0.0",
    estimand = "unweighted adult cohort transport, conditional predictive distribution",
    train_cohort = "NHANES 2015-2016",
    evaluation_cohort = "NHANES 2017-2018",
    age_range = c(18L, 79L),
    outcomes = outcomes,
    train_rows = nrow(train),
    evaluation_population_rows = nrow(test_population),
    evaluation_rows = nrow(test),
    evaluation_sampling = "fixed-seed simple random sample without replacement (seed 20260824)",
    uncertainty = "total",
    n_draw = 1000L,
    package_version = as.character(utils::packageVersion("referent")),
    r_version = R.version.string,
    input_files = input_receipts,
    output_files = output_receipts,
    all_fits_valid = isTRUE(gate$all_fits_valid),
    calibration_pass = isTRUE(gate$calibration_pass),
    pass = isTRUE(gate$pass)
  )
  jsonlite::write_json(
    receipt, file.path(output_dir, "receipt.json"), auto_unbox = TRUE,
    pretty = TRUE, null = "null"
  )
  invisible(list(fit = fit, assessment = assessment, gate = gate))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1L || length(args) > 2L) {
    stop("usage: run_nhanes_external_validation.R OUTPUT_DIR [CACHE_DIR]", call. = FALSE)
  }
  cache <- if (length(args) == 2L) args[[2L]] else file.path(tempdir(), "nhanes-cache")
  result <- run_nhanes_external_validation(args[[1L]], cache)
  print(result$gate, row.names = FALSE)
  if (!isTRUE(result$gate$pass)) quit(status = 1L)
}
