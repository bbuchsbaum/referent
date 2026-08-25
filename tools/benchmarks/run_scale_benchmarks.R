# Enforce broad release budgets for fit, prediction, allocation, and bundle size.

profiled <- function(expr) {
  path <- tempfile("referent-rprofmem-")
  gc()
  utils::Rprofmem(path)
  on.exit(utils::Rprofmem(NULL), add = TRUE)
  timing <- system.time(value <- force(expr))
  utils::Rprofmem(NULL)
  lines <- readLines(path, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  list(value = value, seconds = unname(timing[["elapsed"]]),
       allocated_bytes = sum(bytes, na.rm = TRUE))
}

run_scale_benchmarks <- function(output_dir, mode = c("smoke", "release")) {
  mode <- match.arg(mode)
  if (!requireNamespace("referent", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("referent and jsonlite are required", call. = FALSE)
  }
  size <- if (mode == "release") {
    list(gaussian_train = 10000L, gaussian_test = 3000L,
         shash_train = 1500L, shash_test = 750L)
  } else {
    list(gaussian_train = 2000L, gaussian_test = 750L,
         shash_train = 500L, shash_test = 250L)
  }
  limits <- if (mode == "release") {
    c(gaussian_fit_seconds = 60, gaussian_predict_seconds = 30,
      shash_fit_seconds = 60, shash_predict_seconds = 120,
      maximum_allocated_bytes = 8e9, bundle_bytes = 1e7)
  } else {
    c(gaussian_fit_seconds = 20, gaussian_predict_seconds = 10,
      shash_fit_seconds = 30, shash_predict_seconds = 60,
      maximum_allocated_bytes = 3e9, bundle_bytes = 1e7)
  }

  gaussian_train <- referent::ref_simulate(size$gaussian_train, seed = 20260824)
  gaussian_test <- referent::ref_simulate(size$gaussian_test, seed = 20260825)
  outcomes <- c("y", "marker_01", "marker_02", "marker_03")
  gaussian_spec <- referent::ref_spec(
    referent::ref_gaussian(), location = ~ s(age, k = 8) + sex,
    scale = ~ 1, bam_min_n = 1000L
  )
  gaussian_fit <- profiled(
    referent::ref_fit(gaussian_spec, gaussian_train, outcomes = outcomes)
  )
  gaussian_predict <- profiled(stats::predict(
    gaussian_fit$value, gaussian_test, uncertainty = "total"
  ))
  bundle <- referent::ref_freeze(gaussian_fit$value)
  bundle_bytes <- length(serialize(bundle, NULL, version = 3))

  shash_train <- referent::ref_simulate(size$shash_train, kind = "shash",
                                        seed = 20260826)
  shash_test <- referent::ref_simulate(size$shash_test, kind = "shash",
                                       seed = 20260827)
  shash_spec <- referent::ref_spec(
    referent::ref_shash(), location = ~ s(age, k = 6) + sex,
    scale = ~ s(age, k = 5), skew = ~ 1, tail = ~ 1,
    control = list(seed = 20260824L, n_draw = 200L)
  )
  shash_fit <- profiled(referent::ref_fit(shash_spec, shash_train, outcomes = "y"))
  shash_predict <- profiled(stats::predict(
    shash_fit$value, shash_test, uncertainty = "total", n_draw = 200L
  ))

  observed <- c(
    gaussian_fit_seconds = gaussian_fit$seconds,
    gaussian_predict_seconds = gaussian_predict$seconds,
    shash_fit_seconds = shash_fit$seconds,
    shash_predict_seconds = shash_predict$seconds,
    maximum_allocated_bytes = max(
      gaussian_fit$allocated_bytes, gaussian_predict$allocated_bytes,
      shash_fit$allocated_bytes, shash_predict$allocated_bytes
    ),
    bundle_bytes = bundle_bytes
  )
  results <- data.frame(
    metric = names(observed), observed = unname(observed),
    limit = unname(limits[names(observed)]),
    pass = unname(observed <= limits[names(observed)]),
    stringsAsFactors = FALSE
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(output_dir, "scale_benchmarks.csv")
  utils::write.csv(results, result_path, row.names = FALSE)
  receipt <- list(
    schema_version = "1.0.0",
    mode = mode,
    sizes = size,
    package_version = as.character(utils::packageVersion("referent")),
    r_version = R.version.string,
    platform = R.version$platform,
    all_pass = all(results$pass),
    result_md5 = unname(tools::md5sum(result_path))
  )
  jsonlite::write_json(receipt, file.path(output_dir, "receipt.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  invisible(results)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L || !args[[2L]] %in% c("smoke", "release")) {
    stop("usage: run_scale_benchmarks.R OUTPUT_DIR smoke|release", call. = FALSE)
  }
  results <- run_scale_benchmarks(args[[1L]], args[[2L]])
  print(results, row.names = FALSE)
  if (!all(results$pass)) quit(status = 1L)
}
