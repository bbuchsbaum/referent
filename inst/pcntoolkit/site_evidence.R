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

pcn_adaptation_grid <- function() {
  expand.grid(
    location_prior_n = c(0, 2, 5, 10),
    scale_prior_n = c(0, 5, 10, 25),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

pcn_select_adaptation_candidate <- function(validation) {
  required <- c("site", "location_prior_n", "scale_prior_n", "mean_log_score")
  if (!all(required %in% names(validation)) || !nrow(validation) ||
      any(!is.finite(validation$mean_log_score))) {
    stop("adaptation validation evidence is incomplete", call. = FALSE)
  }
  candidate_site <- validation[
    c("site", "location_prior_n", "scale_prior_n")
  ]
  if (anyDuplicated(candidate_site)) {
    stop("adaptation validation evidence contains duplicate site results",
         call. = FALSE)
  }
  expected_sites <- length(unique(validation$site))
  site_counts <- stats::aggregate(
    site ~ location_prior_n + scale_prior_n, validation,
    function(x) length(unique(x))
  )
  names(site_counts)[names(site_counts) == "site"] <- "validated_sites"
  if (any(site_counts$validated_sites != expected_sites)) {
    stop("adaptation validation evidence is incomplete", call. = FALSE)
  }
  summary <- stats::aggregate(
    mean_log_score ~ location_prior_n + scale_prior_n,
    validation, mean
  )
  summary <- merge(
    summary, site_counts,
    by = c("location_prior_n", "scale_prior_n"), sort = FALSE
  )
  summary <- summary[order(
    -summary$mean_log_score, summary$location_prior_n, summary$scale_prior_n
  ), , drop = FALSE]
  rownames(summary) <- NULL
  summary$selected <- seq_len(nrow(summary)) == 1L
  list(
    selected = summary[1L, c("location_prior_n", "scale_prior_n"), drop = FALSE],
    summary = summary
  )
}
