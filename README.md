# referent

Distributional reference models and individual deviation scores.

`referent` estimates the conditional predictive distribution of one or more
outcomes given a declared reference population, then places new observations
within that distribution. It reports centiles, normal-score representations,
tail probabilities, calibration, site adaptation, and history-conditioned
forecasts. A large deviation is not treated as clinically or morally abnormal.

The package is CDF-first: a Z-score is one representation of a centile, not a
raw standardized residual.

```r
library(referent)

ref <- norm_simulate(300, kind = "gaussian", scale = "age", seed = 1)
spec <- norm_spec(
  family = norm_gaussian(),
  location = ~ s(age, k = 8) + sex,
  scale = ~ s(age, k = 5)
)
fit <- norm_fit(spec, data = ref, outcomes = "y")

target <- norm_simulate(20, kind = "gaussian", scale = "age", seed = 2)
scores <- predict(fit, newdata = target, type = "scores")
norm_assess(fit, newdata = target)$overall
```

Predictive distributions are
[distributional](https://pkg.mitchelloharawild.com/distributional/) vectors,
so the usual generics apply (`cdf()`, `quantile()`, `density()`,
`generate()`, `mean()`, `variance()`, `hilo()`), they sit in tibble columns,
and `ggdist` can draw them. `predict(type = "distribution")` returns one
`<distribution>` column per outcome; `tidy()`, `glance()`, and `augment()`
give per-outcome, per-fit, and wide per-observation summaries.

```r
dists <- predict(fit, newdata = target, type = "distribution")
dists$y
hilo(dists$y, 90)
as_scores(dists$y, target$y)

tidy(fit)
augment(fit, target)   # target + .z_y, .centile_y, .support
```

Gaussian fits return `dist_normal()`; SHASH fits return `dist_shash()`;
`uncertainty = "total"` returns an equal-weight mixture over coefficient
draws (`dist_shash_mc()`); forecasts return `dist_conditioned()`.

Longitudinal change is a first-class dynamic extension, not a standalone
`norm_change()` calculation:

`norm_fit` → `norm_dynamics` → `norm_forecast` / `norm_transition` / `norm_derivative`

Design source of truth:

- [A design for an R normative-modeling library](docs/design/distributional-reference-models.md)
- [Velocity should be a core consequence of the model, not a bolt-on](docs/design/longitudinal-velocity.md)
