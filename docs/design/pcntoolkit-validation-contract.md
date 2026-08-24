# PCNtoolkit validation contract

Status: normative for the `referent` comparison suite

Comparator: PCNtoolkit 1.3.0

Legacy comparator: PCNtoolkit 0.35 is out of scope and must be named explicitly

## Purpose

`referent` and PCNtoolkit both estimate a conditional predictive distribution,
but their estimators, parameterisations, uncertainty calculations, and site
models are not generally identical. The validation target is therefore not
universal bit-for-bit parity. It has two independently reported parts:

1. **Conformance:** shared mathematical quantities and deliberately matched
   fitted models agree within a declared tolerance.
2. **Held-out utility:** when estimators differ, `referent` is non-inferior or
   superior on locked test data without sacrificing calibration.

Passing the second part never excuses failure of the first. Conversely, a
conformant implementation is not called superior without held-out evidence.

## Evidence classes

Every comparison must carry exactly one of these labels.

| Class | Claim permitted | Required evidence |
|---|---|---|
| `exact` | The same formula and parameterisation were evaluated. | Analytic R oracle and pinned PCNtoolkit fixture; absolute plus relative tolerance. |
| `parameter_converted` | The same distribution is represented after a documented, invertible conversion. | Forward and inverse conversion tests plus density, CDF, and quantile agreement. |
| `matched_estimator` | Fitted models target the same likelihood and design closely enough for numerical concordance. | Locked data, declared settings, parameter/prediction discrepancies, and predictive checks. |
| `same_estimand` | Methods estimate the same broad conditional predictive distribution but use different estimators. | Held-out proper scores and calibration, not parameter equality. |
| `non_equivalent` | The quantities or inferential procedures differ materially. | An explicit counterexample or semantic explanation; no parity threshold. |

A failure is recorded, not resolved by relabelling a fixture or loosening a
tolerance. A classification may be changed only by changing the contract and
adding an explanatory regression test.

## Comparison unit

Each result records the package versions, fixture schema, scenario and seed,
training/validation/test row identifiers, outcome transformation, design
matrix, family, uncertainty mode, site handling, and missing-row policy.
Comparisons are made on the same finite test rows and outcome scale.

The prediction record contains, where available:

- predictive median and standard deviation;
- CDF/centile, Gaussianised deviation score, and selected quantiles;
- pointwise log density and aggregate log score;
- interval coverage and calibration error.

Conditional predictions use plug-in distributional parameters. Total
predictions also include fitted-curve uncertainty. These are distinct
estimands and are never pooled in one parity statistic.

## Shared semantics

### Gaussian

For fixed `mu` and positive `sigma`, density, log density, CDF, quantiles,
median, variance, and `z = qnorm(CDF(y))` are `exact`. The gate is
`abs(error) <= 1e-10 + 1e-10 * abs(oracle)`, except log-tail cases where the
comparison is performed directly in log space.

### SHASH

PCNtoolkit's `SHASHb` standardises a sinh-arcsinh random variable to mean zero
and variance one. `mgcv::shash`, used by `referent`, exposes a location-scale
form. With PCNtoolkit parameters `(mu_b, sigma_b, epsilon, delta)`, and `m1`
and `v` the mean and variance of the unstandardised sinh-arcsinh variable,
the equivalent `referent` parameters are

```
sigma_r = sigma_b / (delta * sqrt(v))
mu_r    = mu_b - sigma_b * m1 / sqrt(v)
epsilon_r = epsilon
delta_r   = delta
```

This lane is `parameter_converted`. Density, CDF, and quantiles must satisfy
`abs(error) <= 1e-7 + 1e-7 * abs(oracle)`, including asymmetric and tail
cases. This proves distributional identity only; it does not claim identical
SHASH fitting procedures.

### Transformations

For monotone response transforms, forward/inverse round trips, transformed
CDFs and quantiles, and log-density Jacobians are `exact`. Support violations
must fail closed. The gate is `1e-10` absolute plus relative tolerance unless
a fixture documents a tighter numerical limit.

### Metrics

Metrics are compared only after normalising sign, variance convention,
baseline, centile grid, batch aggregation, and finite-row filtering. RMSE and
EV can then be exact. SMSE requires conversion between population and sample
variance. PCNtoolkit's Spearman `Rho` is not `referent`'s Pearson `cor`.
With standardised outputs, PCNtoolkit 1.3.0 returns `logp` on the standardised
outcome scale even though predictions and centiles are mapped back; the
comparison fixture subtracts the corresponding log-scale Jacobian before any
proper-score contrast.
PCNtoolkit 1.3.0 uses bias-corrected skewness and excess kurtosis, whereas
`referent` reports unadjusted moment estimates. MACE is exact without batch
effects; with batch combinations, PCNtoolkit averages batch-level values
equally while `referent` currently pools rows. The latter cases are
`non_equivalent` unless explicitly recomputed to one definition.

## Fitted concordance lane

Pinned BLR fixtures cover intercept-only, linear, multiple/categorical,
heteroskedastic/nonlinear, transformed-outcome, and balanced-site designs.
Only scenarios with a genuinely matched likelihood/design are labelled
`matched_estimator`; the rest are `same_estimand` diagnostics.

The provisional matched-estimator gates, frozen before routine use, are:

| Quantity | Gate |
|---|---:|
| predictive median normalised RMSE | at most 0.01 |
| predictive standard-deviation normalised RMSE | at most 0.02 |
| deviation-score RMSE | at most 0.03 |
| 95th percentile of absolute deviation-score difference | at most 0.08 |
| absolute 90% interval-coverage difference | at most 0.01 |

Exact semantic gates are not relaxed to these fitted-model gates. Any
platform-specific exception requires a stored result, a mechanistic reason,
and a narrower invariant that still detects scientific drift.

## Superiority and non-inferiority lane

Estimator-divergent scenarios use immutable train/validation/test splits.
They include nonlinear means, heteroskedasticity, skew/heavy tails, unequal
site sizes, site shifts/slopes, covariate shift, unseen sites, and sample-size
variation. Hyperparameters are selected without test outcomes.

The primary contrast is paired test-row log score (`referent` minus
PCNtoolkit), aggregated within simulation replicate. A cluster-aware bootstrap
is used when observations share a site or subject. CRPS, median RMSE, interval
coverage, PIT calibration, and tail calibration are confirmatory.

Results are classified as follows:

- `superior`: the 95% interval for the primary contrast is above zero and all
  prespecified calibration gates pass;
- `equivalent`: the interval lies within a prespecified practical-equivalence
  margin and calibration gates pass;
- `tradeoff`: point accuracy or score improves but a calibration gate fails;
- `inferior`: the interval is below the non-inferiority margin or a critical
  calibration gate fails;
- `inconclusive`: none of the above.

The suite includes negative controls: identical predictions cannot generate
superiority; deliberately overconfident predictions may improve RMSE but must
fail calibration; and controlled perturbations must worsen the relevant
metric monotonically. We report scenarios individually as well as in aggregate
so a favourable average cannot hide a failure regime.

## HBR, SHASH fitting, and site transport

These are separate evidence lanes. PCNtoolkit HBR averages draw-specific
statistics, while `referent` total uncertainty constructs a predictive mixture
and derives its statistics. They are `non_equivalent`; convergence diagnostics,
Monte Carlo error, chain/draw counts, and draw seeds accompany HBR results.

Site evidence distinguishes observed-site fit, small-sample adaptation,
leave-one-site-out transport, and genuinely unseen-site prediction. Balanced
and unequal site sizes are both required. Site-pooled performance never
substitutes for site-stratified calibration.

## Reproducibility and CI

The R test suite consumes checked-in fixtures and never installs Python or
uses the network. Fixtures live under
`tests/fixtures/pcntoolkit/v1.3.0/`, contain hashes and complete provenance,
and are regenerated only by the pinned script under `tools/pcntoolkit/`.

- Pull requests: validate hashes/schema and run semantic plus small fitted
  oracle tests offline.
- Weekly/manual: create a clean Python environment, regenerate PCNtoolkit
  fixtures, fail on schema, environment, identity, or tightly bounded
  numerical drift, and run the BLR scenario matrix. Byte differences alone
  are diagnostic because BLAS and elementary-math libraries vary by platform.
- Release/manual: run stochastic HBR, SHASH-fitting, multi-site transport, and
  the replicated superiority benchmark with convergence receipts.

The release matrix uses 20 independently seeded train/validation/test
simulations for linear Gaussian, nonlinear heteroskedastic, skew-heavy,
unequal-site, and covariate-shift cases. PCNtoolkit predictions are generated
before Referent is fit. Every Referent fit must be valid, and any critical
optimizer, conditioning, covariance, overflow, or convergence warning
invalidates the corresponding comparator replicate. For coverage, MACE, and
tail calibration, the 95% replicate-bootstrap upper bound on mean regret must
be at most 0.01, and no single-replicate regret may exceed 0.05. Matched and
same-estimand scenarios must be equivalent, non-inferior, or superior;
skew-heavy must meet the stronger superiority rule. Scenarios are reported
separately, so one favourable result cannot hide another scenario's failure.
The generation receipt records seeds, dependency versions, row identities,
warning verdicts, and file hashes.

Machine-readable results are retained alongside a human-readable evidence
table. Documentation may say “matches PCNtoolkit” only for a named evidence
class and scenario, and may say “outperforms” only for a locked held-out
contrast that meets the superiority definition above.
