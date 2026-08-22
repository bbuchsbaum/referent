# A design for an R normative-modeling library

Yes—there is a genuinely valuable package here. The right product is **not** a port of PCNtoolkit, a thin wrapper around `mgcv::gam()`, or a collection of every method ever described as normative modeling. It should be a small **distributional reference-modeling layer**:

> Given a reference population, estimate the conditional distribution of one or more outcomes; validate that distribution; place new observations within it; and distinguish trustworthy deviation from extrapolation, site shift, measurement noise, and model uncertainty.

The attached work already points in this direction. It treats normative modeling as a general regression framework rather than a particular algorithm, and emphasizes individual centiles, predictive uncertainty, and deviation scores. The GAMLSS paper adds the essential ability to model nonlinear location, heteroskedastic scale, non-Gaussian shape, hierarchical effects, and transfer to new sites.

There are already R packages for continuous psychometric norms—`NormData`, `cNORM`, and the new GAMLSS-based `normref`—but their center of gravity is psychological or educational test norming rather than arbitrary feature panels, out-of-sample deviation estimation, transport, longitudinal change, and joint deviation profiles. The current PCNtoolkit abstraction is broader, but still handles multiple responses by copying a separate regression model for each response. That leaves a distinct and worthwhile gap.

A good working package name is **`referent`**, with the subtitle *Distributional reference models and individual deviation scores*. It says what the model actually is without implying that the reference population is morally or clinically “normal.” The attached work itself stresses that “normative” is a statistical, reference-dependent term, not a moral judgment or an assertion that a large deviation is unhealthy.

---

## 1. The central abstraction

For outcome or feature \(j\), the model is

\[
\mathcal M_j:\;x \longmapsto F_j(\,\cdot\mid x,\mathcal R),
\]

where:

- \(x\) is the subject’s covariate vector;
- \(\mathcal R\) is the declared reference population;
- \(F_j\) is the complete conditional predictive distribution of outcome \(Y_j\).

For an observed value \(y_{ij}\), the primary derived quantities are

\[
u_{ij}=F_j(y_{ij}\mid x_i),
\]

\[
z_{ij}=\Phi^{-1}(u_{ij}),
\]

\[
p^{\text{tail}}_{ij}=2\min(u_{ij},1-u_{ij}),
\]

along with predictive quantiles, median, mean where defined, and uncertainty summaries.

The package should therefore be **CDF-first, not Z-score-first**. A Z-score is only one representation of a centile. This distinction matters because:

1. It works identically for Gaussian, skewed, heavy-tailed, bounded, count, and ordinal outcomes.
2. It separates fitting a distribution from choosing how to summarize extremity.
3. It permits proper probability calibration after fitting.
4. It makes the meaning of a deviation independent of the response’s units.
5. It prevents the common mistake of interpreting a standardized raw residual as if it were necessarily a calibrated normal score.

For discrete distributions, the package would use randomized or interval probability-integral-transform residuals rather than pretending that a continuous Z-score exists uniquely.

The user-facing score object should distinguish:

| Quantity | Interpretation |
|---|---|
| `centile` | Conditional percentile \(u\) |
| `z` | Normal-score representation \(\Phi^{-1}(u)\) |
| `tail_prob` | Two-sided reference-tail probability |
| `tail_surprisal` | \(-\log p^{\text{tail}}\), scale-free extremity |
| `residual` | Observed minus predicted center, on the original scale |
| `log_density` | Predictive log density; useful for comparing models within an outcome |
| `aleatoric_sd` | Irreducible conditional variation |
| `epistemic_sd` | Uncertainty due to estimating the model |
| `support` | Whether the covariates lie inside the reference domain |

The term **“abnormality score” should not be part of the core API**. The package reports deviation and extremity; the scientific or clinical meaning belongs to the application.

---

## 2. Why `mgcv` should be the first engine

R is unusually well suited to this package because `mgcv` already supplies most of the difficult numerical machinery:

- penalized nonlinear smooths with automatic smoothness selection;
- random-effect and structured smooth terms;
- location–scale–shape families;
- approximate Bayesian coefficient covariance for uncertainty propagation;
- `bam()` for reduced-memory large-\(n\) fitting where the selected family is supported.

In particular, `mgcv::shash()` models location, scale, skewness, and tail weight through separate additive predictors. Its own documentation sensibly warns that such flexibility must be justified by the available data. That warning should be embodied in our package’s model-selection policy rather than left to users to discover after an unstable fit.

The first release should therefore use:

- **`mgcv` as the required engine**;
- classic `gamlss` as a later optional adapter;
- `gamlss2` as a promising extension target, not a required foundation. It is currently a young, modular reimplementation at version 0.1-0.
- `brms` or another full-Bayesian engine only through an extension interface.

The library’s value is not reimplementing splines or likelihood optimization. It is supplying the missing **normative-model semantics, validation, transport, scoring, and graphics** around those engines.

---

## 3. The proposed interface

The common case should be concise, while the underlying specification remains explicit.

```r
library(referent)

spec <- norm_spec(
  family = norm_shash(),

  location = ~ s(age, k = 10) +
               sex +
               s(site, bs = "re"),

  scale = ~ s(age, k = 6) +
            s(site, bs = "re"),

  skew = ~ 1,
  tail = ~ 1
)

fit <- norm_fit(
  spec,
  data = reference,
  outcomes = starts_with("marker_"),
  id = participant_id
)

assessment <- norm_assess(
  fit,
  newdata = validation
)

autoplot(assessment)
```

Prediction should use the ordinary R generic:

```r
scores <- predict(
  fit,
  newdata = target,
  type = "scores",
  uncertainty = "total"
)
```

The result would be a long-form `norm_scores` object:

```text
.row  .id   .outcome   observed   median   centile      z
1     P101  marker_01     12.8      10.2     0.964    1.80
1     P101  marker_02      4.1       6.7     0.071   -1.47
...

tail_prob  aleatoric_sd  epistemic_sd  support  calibrated
0.072          1.31          0.18       in        TRUE
0.142          1.04          0.11       edge      TRUE
```

Useful transformations would be explicit:

```r
as_wide(scores, value = "z")
as_wide(scores, value = "centile")
filter_scores(scores, abs(z) > 2)
```

Model adaptation, cross-fitting, and longitudinal scoring should feel like natural extensions:

```r
local_fit <- norm_adapt(
  fit,
  data = local_reference,
  by = site,
  parameters = c("location", "scale")
)

oof_scores <- norm_crossfit(
  spec,
  data = reference,
  outcomes = starts_with("marker_"),
  folds = 5,
  strata = site,
  cluster = participant_id
)

change <- norm_change(
  local_fit,
  data = repeated_observations,
  id = participant_id,
  time = visit_date,
  residual_correlation = "estimate"
)
```

---

## 4. A small object model

S3 is preferable to R6 here. The objects should be inspectable, serializable, and compatible with familiar R generics.

| Class | Responsibility |
|---|---|
| `norm_spec` | Family, parameter formulas, engine, fitting controls |
| `norm_fit` | A common schema plus named per-outcome fitted models |
| `norm_dist` | Vectorized predictive distribution returned by an engine |
| `norm_scores` | Subject-by-outcome centiles, deviations, tails, uncertainty |
| `norm_assessment` | Proper scores, calibration, support and convergence results |
| `norm_adaptation` | Documented target-domain adjustment |
| `norm_reference` | Frozen, shareable model bundle with provenance |

The crucial internal contract is not “return coefficients.” Every engine must be able to return a distribution object supporting:

```text
cdf(distribution, y)
quantile(distribution, p)
log_density(distribution, y)
draw(distribution, n)
center(distribution)
variance(distribution)
```

That contract is what allows the rest of the package—deviation scores, calibration, plots, longitudinal simulation, and joint summaries—to be completely engine-agnostic.

A failed fit for one outcome should not corrupt the panel. Each outcome should carry a status such as:

```text
ok
nonconverged
insufficient_variation
miscalibrated
unsupported_family
```

`print(fit)` should summarize failures and warnings rather than dumping hundreds of `gam` objects.

---

## 5. The default modeling strategy should be a ladder, not a zoo

The package should not default to the most flexible available model. It should consider a short, ordered ladder:

| Level | Conditional model |
|---|---|
| 1 | Gaussian, nonlinear location, constant scale |
| 2 | Gaussian, nonlinear location and scale |
| 3 | SHASH, nonlinear location and scale, constant skew/tail |
| 4 | SHASH with selected covariate-dependent shape parameters |

The Dinga work found substantial gains from non-Gaussian models, but it also recommends adding complexity in a coherent order and preferring the simpler model when fitted predictions are practically indistinguishable.

Model selection should be **out-of-sample and lexicographic**:

1. Exclude failed or numerically unstable fits.
2. Exclude models with unacceptable calibration.
3. Rank the survivors by out-of-sample logarithmic score.
4. Use CRPS as a secondary, less tail-dominated score.
5. Select the simplest model within one standard error of the best candidate.
6. Require a stronger evidential threshold before allowing covariate-dependent skewness or tail weight.

This is better than either “always use SHASH” or unrestricted automated family search.

Mean-oriented metrics such as \(R^2\), correlation, MAE, and RMSE should remain available, but they are secondary. The attached GAMLSS paper correctly argues that a normative model must be assessed with proper scoring rules because mean prediction alone ignores scale and shape.

---

## 6. The methodological advances worth building

### 6.1 Out-of-sample deviation scores by construction

A surprisingly important issue is whether a person’s deviation was computed using a model that included that same person.

In-sample reference residuals are systematically too well behaved. This becomes especially problematic when:

- reference subjects and target subjects are subsequently compared;
- deviation scores are used as predictors;
- model complexity was selected on the same observations;
- controls contribute to the normative fit but cases do not.

The package should therefore have two clearly separated workflows:

#### Deployment

Fit one final model to all suitable reference observations, then score genuinely new observations.

#### Scientific evaluation

Use `norm_crossfit()` so that every reference observation receives an out-of-fold score. The package then refits the final deployment model separately.

Training-data prediction should be visibly marked:

```text
.in_sample = TRUE
```

and should trigger a warning if used for a group comparison or downstream predictive model.

This is a meaningful improvement over many existing workflows and costs very little conceptually.

---

### 6.2 Calibration should be a first-class operation

Calibration and predictive accuracy are not interchangeable. The 2026 paper gives concrete examples in which outcomes with mediocre SMSE, MSLL, or correlation nevertheless had reasonable centile calibration, while other outcomes had acceptable mean-fit metrics but poor Z-score calibration.

`norm_assess()` should report at least four distinct dimensions:

#### Overall probabilistic fit

- mean log score;
- standardized log score;
- CRPS;
- tail-weighted score when tail behavior is scientifically important.

#### Marginal calibration

- PIT histogram;
- Z-score Q–Q and worm plots;
- mean and variance of \(z\);
- skewness and excess kurtosis;
- observed coverage of selected centiles.

#### Conditional calibration

Calibration should be examined against every important covariate and grouping variable:

\[
E[z\mid x]\approx 0,
\qquad
E[z^2\mid x]\approx 1.
\]

The package can fit light diagnostic GAMs to \(z\) and \(z^2\), showing where location or scale remains misspecified. This is more informative than one global Shapiro–Wilk statistic.

#### Tail calibration

For nominal tail levels such as 0.5%, 1%, 2.5%, and 5%, report observed versus expected exceedance rates with uncertainty. Normative modeling is often used precisely in these sparse regions, where a globally adequate model may still be poor.

#### A backend-independent recalibration layer

A particularly coherent innovation is to allow held-out probability calibration. Suppose

\[
u_i=F_0(y_i\mid x_i)
\]

are PIT values on a held-out reference sample, with empirical CDF \(G\). Define

\[
F^\star(y\mid x)=G\!\left(F_0(y\mid x)\right).
\]

Because \(G\) is monotone, \(F^\star\) remains a valid CDF. This can correct moderate global or group-specific probability miscalibration without refitting the underlying trajectory.

```r
calibrated_fit <- norm_calibrate(
  fit,
  data = calibration_reference,
  method = "rank",
  by = site
)
```

The package must state the limitation clearly: global recalibration targets global calibration; it does not magically create full conditional calibration. Pre- and post-calibration diagnostics should always be retained.

---

### 6.3 Transport and adaptation should be explicit

The attached work correctly treats transfer to an unseen site or acquisition context as a distinct problem. Modeling both location and scale site effects can materially improve predictions and reduce residual site structure.

For a new group \(g\), the default adaptation should freeze the shared biological or behavioral trajectory and estimate only penalized offsets:

\[
\eta_{\mu,g}^{\star}(x)
  =\eta_\mu(x)+\delta_{\mu,g},
\]

\[
\eta_{\sigma,g}^{\star}(x)
  =\eta_\sigma(x)+\delta_{\sigma,g}.
\]

Important design rules:

- location-only adaptation is allowed when the local reference sample is small;
- scale adaptation requires stronger evidence;
- offsets use shrinkage toward zero;
- calibration observations are never reused for evaluation;
- adaptation uncertainty is propagated into predictions;
- the resulting object records exactly which parameters changed.

A recent application used more than 25 local controls per new site and reported good performance, but that should not become a magical hard-coded threshold. The package should instead report effective sample size, parameter uncertainty, and expected calibration precision.

It should also distinguish:

- **adaptation**: changing location or scale parameters for a target domain;
- **recalibration**: correcting probabilities through a monotone PIT mapping;
- **refitting**: changing the shared trajectory itself.

These are scientifically different operations.

---

### 6.4 Longitudinal deviation must model change, not subtract Z-scores

The longitudinal paper makes the decisive conceptual point: a cross-sectional chart is a **trajectory of distributions**, not a **distribution over individual trajectories**. People may cross centiles without pathology.

Consequently,

\[
z(t_2)-z(t_1)
\]

is generally not a calibrated score of unusual change. It double-counts some uncertainty and ignores the covariance between repeated observations. The proposed z-diff score instead scales the observed change using the predictive uncertainty of the difference and the within-person measurement-noise component.

Our package should generalize that idea.

For Gaussian location–scale models:

\[
z_{\Delta}
=
\frac{
(y_2-y_1)-[\mu(x_2)-\mu(x_1)]
}{
\sqrt{
v_{\text{model},\Delta}
+\sigma_1^2+\sigma_2^2
-2\rho_{12}\sigma_1\sigma_2
}
}.
\]

For non-Gaussian models, the package should not force an approximate Gaussian formula. It should:

1. draw parameter values from the fitted uncertainty distribution;
2. draw paired residuals using an estimated or supplied residual copula;
3. construct the predictive distribution of \(Y_2-Y_1\);
4. locate the observed change within that distribution;
5. return `change_centile` and `z_change`.

The within-person correlation can be:

- estimated from repeated reference observations;
- modeled as a function of time separation;
- supplied from test–retest reliability;
- estimated from a held-out local longitudinal control set.

If none of those are available, the function should decline to label change as normatively unusual rather than silently assume independence.

The output should preserve both:

```text
baseline_z
followup_z
change_z
```

because position and change answer different questions.

---

### 6.5 Marginal models first, with one restrained joint extension

Fitting one model per feature is computationally attractive and easy to interpret, but it ignores cross-feature dependence. The attached protocol explicitly acknowledges that separate univariate models do not model correlations among outcomes.

The package should not respond by adding VAEs, Gaussian processes, factor models, neural processes, and every possible multivariate method. Instead, it should provide one transparent second-stage option:

#### Gaussian-copula deviation model

1. Obtain out-of-fold marginal normal scores \(z_{ij}\).
2. Estimate a shrinkage correlation matrix \(R\).
3. For a new subject, calculate a global deviation measure such as

\[
D_i^2=z_i^\top R^{-1}z_i.
\]

4. Calibrate \(D_i\) empirically on held-out reference subjects rather than relying blindly on a \(\chi^2\) approximation.
5. Handle missing outcomes using the corresponding covariance submatrix.
6. Return feature-level contribution summaries.

```r
joint <- norm_joint(
  oof_scores,
  method = "gaussian_copula",
  covariance = "shrinkage"
)
```

This captures the most important dependence structure while keeping the core package intelligible. More elaborate joint engines can live elsewhere.

---

### 6.6 Reference support and extrapolation must be visible

A precise centile is meaningless if a 90-year-old target is scored using a reference sample spanning ages 20–45.

Every prediction should carry a support classification based on:

- numeric range;
- unseen factor levels;
- leverage in the model matrix;
- robust distance from the training covariate distribution;
- local reference density.

Suggested statuses:

```text
in
edge
out
new_group
```

By default, severe extrapolation should produce `z = NA` unless the user explicitly opts in.

A shareable reference bundle should include:

- reference inclusion and exclusion criteria;
- outcome definitions and units;
- covariate names, coding and ranges;
- factor levels;
- transformations;
- family and formulas;
- calibration data provenance;
- missing-data policy;
- package and engine versions;
- feature-level fit and calibration status.

This should be generated automatically as a model card:

```r
norm_card(fit)
```

---

### 6.7 Do not turn \(|z|>2\) into a universal decision rule

With \(p\) well-calibrated independent outcomes, approximately \(0.0455p\) will exceed \(|z|>2\) by chance. With 500 outcomes, seeing several such values is expected.

The package should therefore expose:

- raw two-sided tail probabilities;
- FDR-adjusted feature flags;
- observed versus expected extreme counts;
- empirical subject-level deviation burden;
- joint/global scores when appropriate.

It can still offer:

```r
flag(scores, threshold = 2)
```

but the output should call these **threshold exceedances**, not abnormalities.

---

### 6.8 Aleatoric and epistemic uncertainty should remain separate

The existing framework emphasizes the benefit of separating irreducible outcome variation from uncertainty in model estimation.

For `mgcv`, the package can propagate coefficient uncertainty using the fitted covariance matrix and the linear-predictor matrix. The predictive CDF can then be obtained by Monte Carlo integration:

\[
F(y\mid x,\mathcal D)
=
\int F(y\mid x,\theta)\,
p(\theta\mid\mathcal D)\,d\theta.
\]

The user should be able to choose:

```r
uncertainty = "conditional"  # fitted distribution, parameters fixed
uncertainty = "total"        # integrate model uncertainty
```

`"total"` should be the default for final individual scoring. The faster conditional calculation can be used during exploratory diagnostics and large candidate searches.

---

## 7. Graphics are part of the method, not decoration

Every plotting function should return an ordinary `ggplot` object.

| Plot | Scientific question |
|---|---|
| `type = "centiles"` | What is the fitted conditional reference distribution? |
| `type = "calibration"` | Are nominal centiles attained? |
| `type = "worm"` | Where is distributional shape misspecified? |
| `type = "conditional"` | Does calibration drift with age, score, or site? |
| `type = "support"` | Where is the model interpolating versus extrapolating? |
| `type = "profile"` | How does one individual deviate across outcomes? |
| `type = "heatmap"` | What are the subject-by-feature deviation patterns? |
| `type = "change"` | Is longitudinal change unusual given within-person covariance? |
| `type = "adaptation"` | What changed when the model was transported? |

The default centile plot should show the median and a modest set of reference bands, not twenty equally prominent curves. Individual observations should be overlays, not part of the fitted geometry.

The calibration dashboard should place predictive fit and calibration side by side, because neither substitutes for the other.

---

## 8. Package boundaries

A compact first-class dependency set would be:

```text
Imports:
  mgcv
  ggplot2
  rlang
  vctrs
  tibble
  cli
  Matrix
  stats

Suggests:
  scoringRules
  future
  future.apply
  gamlss
  gamlss2
  testthat
  vdiffr
```

`scoringRules` already provides CRPS, logarithmic scores, and weighted scoring rules for many distributions, so it is a sensible optional dependency rather than something to duplicate.

The package should avoid:

- a full tidyverse dependency;
- Rcpp unless profiling shows a genuine bottleneck;
- neuroimaging formats or visualization;
- downstream classifiers and clustering algorithms;
- imputation frameworks;
- giant family registries;
- hidden preprocessing.

Per-feature fitting is naturally parallel. The core should iterate without first pivoting a large outcome matrix into an enormous long table, and optional `future` support can parallelize outcomes with deterministic seeds.

A frozen model bundle should omit raw training observations by default. It should preserve only what is required for prediction, support checks, calibration, and provenance.

---

## 9. Validation plan and release gates

The tests should establish statistical validity, not merely code coverage.

| Gate | Required demonstration |
|---|---|
| Distribution laws | \(F(Q(p))\approx p\); simulated draws reproduce declared moments and quantiles |
| Gaussian special case | Package Z-scores agree with \((y-\mu)/\sigma\) to numerical tolerance |
| Marginal calibration | Simulated 50%, 80%, 90%, 95%, and 99% coverage lies within Monte Carlo uncertainty |
| Conditional calibration | No systematic residual location or scale drift over modeled covariates |
| Model ladder | Simple generators select simple models; complex generators gain predictive score without inflated tail errors |
| Cross-fitting | Out-of-fold reference scores have the nominal variance and tail rate; in-sample shrinkage is detected |
| Adaptation | Known simulated location and scale shifts are recovered, with calibrated held-out target scores |
| Longitudinal inference | Nominal false-positive rate is maintained across residual correlations, time gaps, and measurement-noise levels |
| Missingness | Feature-specific missing outcomes do not change unaffected fits; predictor missingness is never silently imputed |
| Extrapolation | Out-of-support observations are reliably flagged |
| Serialization | Save/read round trips reproduce predictions and scores within numerical tolerance |
| Parallelism | Sequential and parallel fits are statistically and numerically reproducible |

The benchmark suite should contain:

### Distributional recovery

Simulated nonlinear location, heteroskedastic scale, skewness, light and heavy tails, bounded outcomes, and contamination.

### Reproduction of the GAMLSS results

Recreate the central comparison between Gaussian and non-Gaussian normative models: out-of-sample log score, PIT calibration, skewness, kurtosis, and worm plots.

### Reproduction of the z-diff simulation

Vary true disruption and residual autocorrelation, and verify that `norm_change()` maintains the nominal false-positive rate while naïve Z-score subtraction does not.

### Multi-group transfer

Train on several groups, leave one group out, adapt with local reference samples of varying size, and assess held-out calibration.

### Cross-language check

On shared simulated data, compare Gaussian and SHASH predictions against direct `mgcv` calls and, where model definitions coincide, PCNtoolkit.

### Scale tests

Benchmark at least:

- \(n=1{,}000,\ p=10\);
- \(n=10{,}000,\ p=100\);
- \(n=100{,}000,\ p=100\) for supported `bam` models;
- \(p=1{,}000\) with streaming or directory-backed storage.

---

## 10. Development sequence

### Phase 1: the scientific kernel

Implement and freeze:

1. the predictive-distribution contract;
2. Gaussian location–scale and SHASH `mgcv` engines;
3. centile, Z, tail and uncertainty calculations;
4. out-of-sample assessment;
5. cross-fitting;
6. support detection;
7. simulation laws.

There should be almost no polished visualization until these pass.

### Phase 2: the usable reference-model package

Add:

- centile and calibration graphics;
- multi-outcome ergonomics;
- location/scale adaptation;
- PIT recalibration;
- model bundles and model cards;
- parallel execution;
- three excellent vignettes.

The vignettes should use a synthetic growth measure, a cognitive or psychometric score, and a multi-site laboratory or sensor panel—not brain data.

### Phase 3: restrained extensions

Add:

- longitudinal `norm_change()`;
- discrete and ordinal outcomes with randomized PIT;
- the Gaussian-copula joint layer;
- optional `gamlss`/`gamlss2` engines;
- domain-specific extensions in separate packages.

---

## Bottom line

The package is worth building. Its paper-worthy contribution would not be “GAMLSS in a convenient wrapper.” It would be the coherent combination of:

1. **a distribution-first contract** in which centiles and Z-scores are derived from a full conditional predictive distribution;
2. **out-of-sample and calibration-first inference**, including cross-fitted reference scores and conditional calibration diagnostics;
3. **explicit transport and temporal semantics**, distinguishing adaptation, probability recalibration, baseline deviation, and unusual change;
4. **a restrained multivariate layer** that captures dependence without turning the package into a method zoo.

That would make `referent` both simpler and methodologically sharper than the attached workflows. I would begin by formalizing the distribution contract and its laws, then build the simulation harness, and only after that settle the final public API.
