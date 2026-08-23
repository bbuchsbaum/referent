pcn_site_summary <- function(prediction, method, lane) {
  required <- c("site", "observed", "centile", "log_density", "q05", "q50", "q95")
  if (!all(required %in% names(prediction))) {
    stop("site prediction table is missing required columns", call. = FALSE)
  }
  parts <- split(seq_len(nrow(prediction)), prediction$site)
  rows <- lapply(names(parts), function(site) {
    d <- prediction[parts[[site]], , drop = FALSE]
    z <- stats::qnorm(d$centile)
    data.frame(
      method = method, lane = lane, site = site, n = nrow(d),
      mean_log_score = mean(d$log_density),
      median_rmse = sqrt(mean((d$observed - d$q50)^2)),
      coverage90 = mean(d$observed >= d$q05 & d$observed <= d$q95),
      coverage90_se = sqrt(0.9 * 0.1 / nrow(d)),
      mace = pcn_benchmark_mace(d$centile),
      mean_z = mean(z), var_z = stats::var(z),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

pcn_site_gate <- function(summary, coverage_floor = 0.80, mace_ceiling = 0.08,
                          mean_z_ceiling = 0.35) {
  per_site <- summary$coverage90 >= coverage_floor &
    summary$mace <= mace_ceiling & abs(summary$mean_z) <= mean_z_ceiling
  data.frame(summary, pass = per_site, stringsAsFactors = FALSE)
}
