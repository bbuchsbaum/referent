# Scale and bundle evidence

`scale_benchmarks.csv` is the retained macOS arm64 release-budget run for
Referent 0.1.0. It covers a 10,000-row, four-outcome Gaussian/BAM fit and
prediction; a 1,500-row SHASH fit and 750-row, 200-draw prediction; allocated
bytes; and the serialized four-outcome frozen-bundle size.

Every declared budget passed. These timings are guardrails against order-of-
magnitude regressions, not portable performance promises: hosted CI reruns the
smoke budget on Linux, and release/manual CI reruns the larger budget with its
own platform receipt.
