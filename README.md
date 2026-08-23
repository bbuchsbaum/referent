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

> **Development status:** `referent` is an early development package
> (`0.0.0.9000`). The current fitting surface supports numeric outcomes with
> Gaussian or sinh-arcsinh (SHASH) predictive families; categorical outcomes
> are reported as unsupported instead of being silently coerced.

Browse the [documentation site](https://bbuchsbaum.github.io/referent/)
(or `browseVignettes("referent")` after installing), and start with [Getting started](https://bbuchsbaum.github.io/referent/articles/getting-started/). Keep
[Troubleshooting reference-model workflows](https://bbuchsbaum.github.io/referent/articles/troubleshooting/)
nearby for support, transport, calibration, panel-fit, and longitudinal
failure modes. The complete article map is below.

## Installation

```r
# install.packages("pak")
pak::pak("bbuchsbaum/referent")
```

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

ref_assess(fit, newdata = target)$marginal[, c("mean_z", "var_z", "cover_95")]
autoplot(fit, type = "centiles", by = sex, newdata = target)
```

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
| Freeze and document a reference | `ref_freeze()` |
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
- [Troubleshooting reference-model workflows](https://bbuchsbaum.github.io/referent/articles/troubleshooting/):
  missing and unsupported covariates, transported unseen sites, small local
  samples, partial panel fits, honest assessment, and unidentified dynamics.

Contributor design notes:

- [A design for an R normative-modeling library](docs/design/distributional-reference-models.md)
- [Velocity should be a core consequence of the model, not a bolt-on](docs/design/longitudinal-velocity.md)
- [PCNtoolkit validation contract](docs/design/pcntoolkit-validation-contract.md)
