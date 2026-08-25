# referent

Distributional reference models and individual deviation scores.

`referent` estimates the conditional predictive distribution of one or more
outcomes given a declared reference population, then places new observations
within that distribution. It reports centiles, normal-score representations,
tail probabilities, calibration, site adaptation, joint deviation across
outcomes, and history-conditioned longitudinal forecasts. A large deviation is
a threshold exceedance with a known chance rate, not an abnormality.

The package is CDF-first: a Z-score is one representation of a centile, not a
raw standardized residual.

> **Development status:** `referent` 0.1.0 is an initial, deliberately scoped,
> source-only release candidate in [PR #4](https://github.com/bbuchsbaum/referent/pull/4).
> It is not yet on the default branch or available as a tagged release. The
> current fitting surface supports numeric outcomes with
> Gaussian or sinh-arcsinh (SHASH) predictive families; categorical outcomes
> are reported as unsupported instead of being silently coerced. Exact and
> matched-estimator PCNtoolkit lanes pass. The replicated release matrix finds
> four named scenarios equivalent and the locked skew-heavy scenario superior;
> HBR convergence and every-site transport also pass under their recorded
> controls. The NHANES 2017-2018 rerun is post-hoc model-development evidence,
> while a separately preregistered, untouched 2013-2014 confirmation passes the
> same frozen conditional-calibration contract. These are scoped results, not a
> blanket equivalence or superiority claim.

For the candidate documentation, use `browseVignettes("referent")` after
installing it. The [published documentation site](https://bbuchsbaum.github.io/referent/)
tracks the default branch and may lag an open candidate. Start with [Getting started](https://bbuchsbaum.github.io/referent/articles/getting-started/), and keep
[Troubleshooting reference-model workflows](https://bbuchsbaum.github.io/referent/articles/troubleshooting/)
nearby for support, transport, calibration, panel-fit, and longitudinal
failure modes. The complete article map is below.

## Installation

```r
# install.packages("pak")
# Install the current source-only 0.1.0 candidate from PR #4:
pak::pak("bbuchsbaum/referent#4")
```

The PR-qualified reference is intentional while 0.1.0 remains source-only;
`bbuchsbaum/referent` without `#4` installs the repository's default branch.

## A first reference model

```r
library(referent)

ref <- ref_simulate(300, kind = "gaussian", scale = "age", seed = 1)
spec <- ref_spec(
  family = ref_gaussian(),
  location = ~ s(age, k = 8) + sex,
  scale = ~ s(age, k = 5)
)
fit <- ref_fit(spec, data = ref, outcomes = "y")
fit

target <- ref_simulate(100, kind = "gaussian", scale = "age", seed = 2)
scores <- predict(fit, newdata = target, type = "scores")
scores[1:3, c(".id", "observed", "median", "centile", "z", "tail_prob", "support")]
table(scores$support)

assessment <- ref_assess(fit, newdata = target)
assessment$marginal[, c("n", "mean_z", "var_z", "cover_95")]
autoplot(fit, type = "centiles", by = sex, newdata = target)
```

Support is part of the scoring result: the table makes boundary masking
visible, while `n` states the number of rows that actually entered the held-out
assessment.

Predictive distributions are
[distributional](https://pkg.mitchelloharawild.com/distributional/)
vectors, so the usual generics apply (`cdf()`, `quantile()`, `density()`,
`generate()`, `mean()`, `variance()`, `hilo()`), they sit in tibble columns,
and `ggdist` can draw them. `predict(type = "distribution")` returns one
`<distribution>` column per outcome; `tidy()`, `glance()`, and `augment()`
give per-outcome, per-fit, and wide per-observation summaries.

```r
dists <- predict(fit, newdata = target[1:3, ], type = "distribution")
dists$y
hilo(dists$y, 90)
as_scores(dists$y, target$y[1:3])

tidy(fit)
augment(fit, target[1:3, ])   # target + .z_y, .centile_y, .support
```

Gaussian and SHASH fits, coefficient uncertainty, and history-conditioned
forecasts all share this public `<distribution>` contract. Code should use the
generics above rather than depend on an internal representation.

## Workflow

| Step | Function |
|---|---|
| Specify a family, formulas, and an optional response transform | `ref_spec()`, `ref_gaussian()`, `ref_shash()` |
| Fit one model per outcome (in parallel under a `future` plan) | `ref_fit()` |
| Score observations or get distributions | `predict()`, `as_scores()`, `augment()` |
| Choose among candidate distributions out of sample | `ref_select()`, `ref_crossfit()` |
| Check calibration on held-out data | `ref_assess()` |
| Transport to a new site | `ref_adapt()`, `ref_calibrate()`, `ref_support()` |
| Remove site effects from the data themselves | `predict(type = "harmonised")` |
| Joint deviation across outcomes | `ref_joint()` |
| Threshold exceedances with FDR | `ref_flag()` |
| Freeze, persist, and validate a trusted versioned reference | `ref_freeze()`, `ref_write()`, `ref_read()` |
| Longitudinal change | `ref_dynamics()`, `ref_transition()`, `ref_forecast()`, `ref_derivative()` |
| Graphics | `autoplot()` methods, `fortify_centiles()`, `theme_referent()` |

Longitudinal change is a dynamic extension of the reference model, not a
standalone calculation: `ref_fit()` gives the marginal distributions,
`ref_dynamics()` estimates within-person dependence on normal scores, and
`ref_transition()` / `ref_forecast()` read innovation Z, change Z,
velocity centiles, and conditional forecasts off the same process.

## Articles

- [Getting started](https://bbuchsbaum.github.io/referent/articles/getting-started/): fit, model ladder,
  scoring, support, predictive distributions, and interpretation.
- [Validate and choose a reference model](https://bbuchsbaum.github.io/referent/articles/validate-reference/):
  out-of-fold model selection, held-out calibration, and conditional drift.
- [Brain charts across sites](https://bbuchsbaum.github.io/referent/articles/brain-charts-across-sites/):
  random site effects, out-of-fold reference scores, scaling to many
  outcomes, scoring a new site (population curve vs adapt vs calibrate vs
  refit), joint deviation, FDR flags, shipping a frozen reference.
- [Longitudinal change](https://bbuchsbaum.github.io/referent/articles/longitudinal-change/): dynamics,
  identifiability and the stable model, six longitudinal quantities,
  forecasts.
- [Coming from PCNtoolkit](https://bbuchsbaum.github.io/referent/articles/pcntoolkit/): concept and metric
  mapping (warped BLR, HBR batch effects and random slopes, harmonised
  outputs, MSLL/SMSE/EV/MACE), numerical-validation policy, and what is
  deliberately missing. The scoped, machine-readable evidence is in
  [docs/evidence/pcntoolkit/v1.3.0](docs/evidence/pcntoolkit/v1.3.0/README.md).
- [NHANES cohort evidence](docs/evidence/nhanes/README.md): retained post-hoc
  2015-2016 to 2017-2018 development diagnostics plus a preregistered,
  genuinely untouched 2013-2014 cycle confirmation.
- [Troubleshooting reference-model workflows](https://bbuchsbaum.github.io/referent/articles/troubleshooting/):
  missing and unsupported covariates, transported unseen sites, small local
  samples, partial panel fits, honest assessment, and unidentified dynamics.

Contributor design notes:

- [A design for an R normative-modeling library](docs/design/distributional-reference-models.md)
- [Velocity should be a core consequence of the model, not a bolt-on](docs/design/longitudinal-velocity.md)
- [PCNtoolkit validation contract](docs/design/pcntoolkit-validation-contract.md)
