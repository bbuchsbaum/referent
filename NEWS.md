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
  `ref_read()` migration for unversioned development bundles.
- Numerical validation includes offline PCNtoolkit 1.3.0 oracles, a replicated
  five-scenario release matrix, HBR/site receipts, independent NHANES cohort
  validation, cross-platform R CMD check, and scale/allocation budgets.

The 0.1.0 release-candidate evidence is intentionally mixed: four replicated
PCNtoolkit scenarios are equivalent, while the skew-heavy comparator emitted
critical optimizer warnings; HBR/site gates and conditional NHANES transport
also identify failures. These are retained evidence, not passing claims.
