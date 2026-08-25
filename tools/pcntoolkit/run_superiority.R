# Run the locked skew/heavy-tail superiority benchmark.

source("tools/pcntoolkit/run_concordance.R")
source("tools/pcntoolkit/benchmark_helpers.R")

pcn_run_superiority <- function(root, B = 4999L, n_draw = 2000L,
                               seed = 91L) {
  inputs <- utils::read.csv(file.path(root, "fitted_inputs.csv"), stringsAsFactors = FALSE)
  comparator <- utils::read.csv(
    file.path(root, "fitted_predictions.csv"), stringsAsFactors = FALSE
  )
  referent <- pcn_referent_predictions(inputs, "skew_heavy", n_draw = n_draw)
  comparator <- comparator[comparator$scenario == "skew_heavy", , drop = FALSE]
  comparator <- comparator[match(referent$row_id, comparator$row_id), , drop = FALSE]
  pcn_compare_predictions(
    referent, comparator, cluster = referent$row_id,
    scenario = "skew_heavy", B = B, seed = seed
  )
}

if (sys.nframe() == 0L) {
  root <- file.path("tests", "fixtures", "pcntoolkit", "v1.3.0")
  result <- pcn_run_superiority(root)
  print(result, row.names = FALSE)
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 2L && identical(args[[1L]], "--output")) {
    dir.create(dirname(args[[2L]]), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(result, args[[2L]], row.names = FALSE)
  } else if (length(args)) {
    stop("usage: run_superiority.R [--output FILE]")
  }
  if (!identical(result$classification, "superior")) quit(status = 1L)
}
