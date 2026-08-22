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

ref <- norm_simulate(300, kind = "gaussian", scale = "age", seed = 1)
spec <- norm_spec(
  family = norm_gaussian(),
  location = ~ s(age, k = 8) + sex,
  scale = ~ s(age, k = 5)
)
fit <- norm_fit(spec, data = ref, outcomes = "y")
fit

target <- norm_simulate(100, kind = "gaussian", scale = "age", seed = 2)
scores <- predict(fit, newdata = target, type = "scores")
scores[1:3, c(".id", "observed", "median", "centile", "z", "tail_prob", "support")]

norm_assess(fit, newdata = target)$marginal[, c("mean_z", "var_z", "cover_95")]
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
draws (`dist_shash_mc()`); forecasts return `dist_conditioned()`.

## Workflow

| Step | Function |
|---|---|
| Specify a family and formulas | `norm_spec()`, `norm_gaussian()`, `norm_shash()` |
| Fit one model per outcome | `norm_fit()` |
| Score observations or get distributions | `predict()`, `as_scores()`, `augment()` |
| Choose among candidate distributions out of sample | `norm_select()`, `norm_crossfit()` |
| Check calibration on held-out data | `norm_assess()`, `acceptable_calibration()` |
| Transport to a new site | `norm_adapt()`, `norm_calibrate()`, `norm_support()` |
| Joint deviation across outcomes | `norm_joint()` |
| Threshold exceedances with FDR | `norm_flag()` |
| Freeze and document a reference | `norm_reference()` |
| Longitudinal change | `norm_dynamics()`, `norm_transition()`, `norm_forecast()`, `norm_derivative()` |
| Graphics | `autoplot()` methods, `fortify_centiles()`, `theme_referent()` |

Longitudinal change is a dynamic extension of the reference model, not a
standalone calculation: `norm_fit()` gives the marginal distributions,
`norm_dynamics()` estimates within-person dependence on normal scores, and
`norm_transition()` / `norm_forecast()` read innovation Z, change Z,
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
