# PCNtoolkit 1.3.0 evidence

This directory is generated from the pinned fixture and contains the current
machine-readable comparison result. It supports scoped claims, not blanket
equivalence.

- `evidence_summary.csv` is the claim-level ledger.
- `concordance.csv` contains every fitted comparison scenario and discrepancy.
- `superiority.csv` contains the paired held-out proper-score interval and
  calibration guardrails for the locked skew/heavy-tail case.
- `release_superiority_summary.csv` and `release_superiority_replicates.csv` contain the
  20-replicate, five-scenario 0.1.0 release-candidate result.
- `release_generation_receipt.json` and `release_matrix_receipt.json` bind the
  comparator warnings, dependencies, seeds, row/file identities, package
  version, and release verdict.
- `hbr_receipt.json`, `hbr_convergence.csv`, `site_comparison.csv`, and
  `site_referent_receipt.json` retain the HBR/transport result.

The exact and converted distribution lanes pass their strict semantic gates.
Five deliberately matched fitted scenarios pass all prediction, scale, Z, and
coverage margins. Estimator-divergent scenarios remain diagnostics.

The single checked-in skew/heavy-tail fixture remains a diagnostic result. In
the stronger replicated release matrix, linear Gaussian, nonlinear
heteroskedastic, unequal-site, and covariate-shift scenarios are all classified
`equivalent`. The skew-heavy log-score contrast favours Referent by 0.397
(95% replicate-bootstrap interval 0.377 to 0.416), with passing calibration
guardrails, but the claim is classified `comparator_failure`: PCNtoolkit
emitted critical ill-conditioning, non-positive-definite, and optimizer-retry
warnings in one or more replicates. A numerically favourable contrast against
an invalid comparator fit is not called superiority.

HBR and new-site transfer stay marked `non_equivalent`. The 0.1.0
release-candidate HBR run failed with 10 divergences in the base stage and 5 in
transfer, despite acceptable R-hat and effective sample sizes. Referent passed
the four observed-site gates but failed the independently evaluated site-5
adaptation gate. A pooled average cannot override either failure.

Regenerate these tables with:

```sh
Rscript tools/pcntoolkit/render_evidence.R
```
