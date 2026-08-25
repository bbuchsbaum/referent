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
- `site_data.csv`, `pcntoolkit_predictions.csv`,
  `pcntoolkit_site_summary.csv`, `hbr_convergence.csv`, and
  `hbr_divergences.csv` retain the raw HBR evidence bound by
  `hbr_receipt.json`.
- `referent_predictions.csv`, `site_comparison.csv`,
  `adaptation_selection.csv`, and `adaptation_selection_summary.csv` retain the
  site evidence bound by `site_referent_receipt.json`.

The exact and converted distribution lanes pass their strict semantic gates.
Five deliberately matched fitted scenarios pass all prediction, scale, Z, and
coverage margins. Estimator-divergent scenarios remain diagnostics.

The single checked-in skew/heavy-tail fixture remains a diagnostic result. In
the stronger fresh-seed matrix, linear Gaussian,
nonlinear heteroskedastic, unequal-site, and covariate-shift scenarios are all
classified `equivalent`. The skew-heavy log-score contrast favours Referent by
0.3961 (95% replicate-bootstrap interval 0.3785 to 0.4120), with passing
calibration guardrails and valid fits in every replicate, so this named
scenario is classified `superior`. The PCNtoolkit receipt records the 0.01
L-BFGS-B finite-difference step that removed the earlier invalid objective
probes. The seed block and result first entered version control together, so
this is not described as independently preregistered.

HBR and new-site transfer stay marked `non_equivalent`: PCNtoolkit draw
averages and Referent predictive mixtures are different estimands. The 0.1.0
release run nevertheless passes its operational gates with zero divergences in
both base and transfer stages, maximum R-hat 1.00, and minimum bulk/tail ESS
above 1,700. Every site passes separately. Referent's location/scale priors
were chosen by leave-one-observed-site-out conditional log score using only
`reference_train` rows; the independent site-5 result has 90% coverage 0.907,
MACE 0.0267, mean Z 0.0563, and variance Z 1.051. A pooled average is never used
to override a site result.

The receipts do not establish their own verdicts. Every retained CSV is bound
by SHA-256, row count, and exact column schema; the validator independently
recomputes row identities, per-draw divergence counts, R-hat/ESS summaries,
prediction summaries, adaptation-prior selection, and all ten per-site gates.
The two transfer-stage parameters with undefined R-hat are exactly constant
posterior rows (`sd = 0`), are retained and counted explicitly, and do not enter
the maximum finite R-hat calculation.

Regenerate the fixture-derived tables with:

```sh
Rscript tools/pcntoolkit/render_evidence.R
```

The release matrix, HBR convergence, and site-transfer receipts are generated
by their dedicated pinned workflows in `tools/pcntoolkit`; their JSON receipts
record the exact commands, environments, seeds, and file hashes.

Validate the retained HBR/site bundle with:

```sh
Rscript tools/pcntoolkit/validate_hbr_site_evidence.R \
  docs/evidence/pcntoolkit/v1.3.0 \
  docs/evidence/pcntoolkit/v1.3.0 \
  hbr_receipt.json site_referent_receipt.json
```
