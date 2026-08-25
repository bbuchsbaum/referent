# Validation helpers must also be available when tests run against an installed
# package. Top-level tools/ files are in the source tarball but are not installed.
pcn_helper_files <- c(
  "compare_results.R", "run_concordance.R", "benchmark_design.R",
  "benchmark_helpers.R", "site_evidence.R"
)
for (pcn_helper_file in pcn_helper_files) {
  pcn_helper_path <- system.file("pcntoolkit", pcn_helper_file, package = "referent")
  if (!nzchar(pcn_helper_path)) {
    pcn_helper_path <- testthat::test_path(
      "..", "..", "inst", "pcntoolkit", pcn_helper_file
    )
  }
  source(pcn_helper_path, local = environment())
}
rm(pcn_helper_file, pcn_helper_files, pcn_helper_path)
