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
  five-scenario release matrix, HBR/site receipts, NHANES cohort diagnostics,
  cross-platform R CMD check, and scale/allocation budgets.

The 0.1.0 release-candidate evidence passes its registered computational gates
without relaxing them. Four replicated PCNtoolkit scenarios are equivalent and
the fresh-seed skew-heavy matrix is superior with a valid comparator;
deterministic high-acceptance HBR sampling passes convergence and
leave-one-site-out-selected adaptation passes every site. The NHANES 2017-2018
rerun passes computationally but is explicitly post-hoc model-development
evidence because its earlier failure informed the scale model. A separately
preregistered, untouched NHANES 2013-2014 confirmation passes the unchanged
full calibration contract and is the evidence used for the scoped
conditional-transport claim.
