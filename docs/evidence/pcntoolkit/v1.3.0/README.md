# PCNtoolkit 1.3.0 evidence

This directory is generated from the pinned fixture and contains the current
machine-readable comparison result. It supports scoped claims, not blanket
equivalence.

- `evidence_summary.csv` is the claim-level ledger.
- `concordance.csv` contains every fitted comparison scenario and discrepancy.
- `superiority.csv` contains the paired held-out proper-score interval and
  calibration guardrails for the locked skew/heavy-tail case.

The exact and converted distribution lanes pass their strict semantic gates.
Five deliberately matched fitted scenarios pass all prediction, scale, Z, and
coverage margins. Estimator-divergent scenarios remain diagnostics.

In the locked skew/heavy-tail case, Referent's prespecified SHASH model is
called superior only because the paired response-scale log-score interval is
above zero and the coverage, MACE, and tail-calibration non-inferiority gates
all pass. This result is scenario-specific; it is not a general claim that
Referent dominates PCNtoolkit.

HBR and new-site transfer stay marked `non_equivalent`. Release evidence must
include convergence receipts and every site's calibration result. The release
job fails if any HBR stage misses its R-hat, effective-sample-size, or
divergence gate; a pooled site average cannot override a failing small site.
The separate release superiority job also generates 20 independent comparison
replicates, requires every SHASH fit to converge, gates calibration on
replicate-bootstrap non-inferiority with a per-replicate fail-closed ceiling,
and gates the claim on the replicate-level paired log-score interval rather
than the single checked-in fixture.

Regenerate these tables with:

```sh
Rscript tools/pcntoolkit/render_evidence.R
```
