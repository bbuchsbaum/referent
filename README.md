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

Gaussian fits return `dist_normal()`; SHASH fits return `dist_shash()`;
`uncertainty = "total"` returns an equal-weight mixture over coefficient
draws (an internal `dist_shash_draws` class); forecasts return
`dist_conditioned()`.

## Workflow

| Step | Function |
|---|---|
| Specify a family and formulas | `ref_spec()`, `ref_gaussian()`, `ref_shash()` |
| Fit one model per outcome | `ref_fit()` |
| Score observations or get distributions | `predict()`, `as_scores()`, `augment()` |
| Choose among candidate distributions out of sample | `ref_select()`, `ref_crossfit()` |
| Check calibration on held-out data | `ref_assess()` |
| Transport to a new site | `ref_adapt()`, `ref_calibrate()`, `ref_support()` |
| Joint deviation across outcomes | `ref_joint()` |
| Threshold exceedances with FDR | `ref_flag()` |
| Freeze and document a reference | `ref_reference()` |
| Longitudinal change | `ref_dynamics()`, `ref_transition()`, `ref_forecast()`, `ref_derivative()` |
| Graphics | `autoplot()` methods, `fortify_centiles()`, `theme_referent()` |

Longitudinal change is a dynamic extension of the reference model, not a
standalone calculation: `ref_fit()` gives the marginal distributions,
`ref_dynamics()` estimates within-person dependence on normal scores, and
`ref_transition()` / `ref_forecast()` read innovation Z, change Z,
velocity centiles, and conditional forecasts off the same process.

## Articles

- [Getting started](vignettes/getting-started.Rmd): fit, model ladder,
  calibration, scoring, tidiers.
- [Brain charts across sites](vignettes/brain-charts-across-sites.Rmd):
  random site effects, out-of-fold reference scores, scoring a new site
  (population curve vs adapt vs calibrate vs refit), joint deviation, FDR
  flags, what to report.
- [Longitudinal change](vignettes/longitudinal-change.Rmd): dynamics,
  identifiability and the stable model, six longitudinal quantities,
  forecasts.
- [Coming from PCNtoolkit](vignettes/pcntoolkit.Rmd): concept and metric
  mapping.

Design source of truth:

- [A design for an R normative-modeling library](docs/design/distributional-reference-models.md)
- [Velocity should be a core consequence of the model, not a bolt-on](docs/design/longitudinal-velocity.md)
