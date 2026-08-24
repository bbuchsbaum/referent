# referent 0.1.0

Initial scoped release of the CDF-first normative-modelling API.

- Gaussian and SHASH distributional reference fits, total coefficient
  uncertainty, response transforms, support checks, site adaptation, PIT
  recalibration, harmonisation, joint scores, and FDR-aware flags.
- Fail-closed model selection now requires valid fits and complete marginal,
  tail, shape, and conditional calibration; uncalibrated exploratory selection
  requires an explicit opt-in.
- Training, calibration, and evaluation roles carry full SHA-256 provenance.
  Calibration is bound to its uncertainty estimand and coefficient-draw count.
- Longitudinal dynamics, transitions, and forecasts share one recorded
  marginal-uncertainty estimand and disclose plug-in kernel uncertainty.
- Frozen references use bundle schema 1.0.0 with validated `ref_write()` and
  `ref_read()` migration for unversioned development bundles. Because RDS is
  executable serialization, `ref_read()` requires an explicit `trusted = TRUE`
  acknowledgement before deserializing a bundle from an authenticated source.
- Numerical validation includes offline PCNtoolkit 1.3.0 oracles, a replicated
  five-scenario release matrix, HBR/site receipts, independent NHANES cohort
  validation, cross-platform R CMD check, and scale/allocation budgets.

The 0.1.0 release-candidate evidence now passes its registered gates without
relaxing them. Four replicated PCNtoolkit scenarios are equivalent and the
independently seeded skew-heavy confirmation scenario is superior with a valid
comparator; deterministic high-acceptance HBR sampling passes convergence,
leave-one-site-out-selected adaptation passes every site, and the independent
NHANES lane passes marginal, tail, shape, and conditional calibration. Claims
remain limited to the named evidence classes, scenarios, and cohorts.
