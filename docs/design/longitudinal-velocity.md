# Velocity should be a core consequence of the model, not a bolt-on

Agreed. I would promote longitudinal change from the earlier `norm_change()` idea into a **first-class dynamic extension of the reference model**.

The central abstraction becomes:

\[
\boxed{
\text{marginal reference distribution}
+
\text{within-person dependence process}
=
\text{longitudinal reference distribution}
}
\]

The static model tells us what values are expected at a given time and covariate configuration. The dynamic component tells us how observations from the **same person** are expected to move through those distributions. Forecasts, change scores, velocity centiles, “thrive lines,” regression-to-the-mean adjustment, and history-conditioned anomaly scores all then emerge from the same conditional predictive distribution.

The attached paper gets the conceptual diagnosis right: ordinary population centiles describe the distribution at successive ages, but not how particular individuals move between those distributions. It correctly emphasizes centile crossing, developmental asynchrony, regression to the mean, and conditioning forecasts on a person’s history. The diagram on page 5 is especially useful because it separates these failure modes visually and shows why the same raw slope can have different significance depending on age, starting centile, and prior history.

But its method should be treated as a good proof of concept, not as the final architecture for our package.

---

## 1. What the paper’s “Z-gain” actually is

On pages 12–14, the paper transforms each observation through a cross-sectional normative model into a marginal standard-normal score \(Z_t\). For two observations with correlation \(r\), it defines

\[
Z_{\mathrm{gain}}
=
\frac{Z_2-rZ_1}{\sqrt{1-r^2}}.
\]

This is the standardized residual from

\[
Z_2\mid Z_1
\sim
N(rZ_1,\,1-r^2).
\]

It is therefore not literally a velocity in units per year. It is a **conditional predictive innovation**: how surprising the follow-up position is, given the baseline position and the expected longitudinal correlation. Elapsed time enters indirectly through \(r(t_1,t_2)\), rather than through division by \(t_2-t_1\).

The multiple-time-point extension is the corresponding multivariate Gaussian calculation:

\[
Z_*
\mid
Z_{1:m}
\sim
N\!\left(
\boldsymbol\beta^\top Z_{1:m},
\operatorname{SD}_*^2
\right),
\]

with regression weights obtained from the correlation matrix using a sweep operation.

That is a sensible foundation, but the terminology hides several distinct quantities that our package should keep separate.

### The key distinction

For two standard-normal observations with correlation \(r\),

\[
I
=
\frac{Z_2-rZ_1}{\sqrt{1-r^2}}
\]

asks:

> Is the follow-up position surprising after conditioning on baseline?

By contrast,

\[
D
=
\frac{Z_2-Z_1}{\sqrt{2(1-r)}}
\]

asks:

> Is the amount of change unusual in the population distribution of changes?

The first is a **conditional innovation**. The second is a **standardized difference**. The earlier z-diff work focuses on the latter and correctly shows that naïvely subtracting two cross-sectional Z-scores mis-scales the common person-level component and modeling uncertainty. Its generalized denominator explicitly contains the within-person autocorrelation through \(2\sigma^2(1-\rho)\).

Neither is universally “the” velocity statistic. They answer different questions, and our model can produce both.

---

## 2. The natural dynamic extension of our distribution-first model

For outcome \(j\), the static model already gives

\[
Y_{ij}(t)
\sim
F_j\!\left(
\,\cdot\mid x_{ij}(t)
\right),
\]

where \(F_j\) may be Gaussian location–scale, SHASH, or another supported distribution. GAMLSS-style modeling lets location, scale, and shape vary with time and other covariates, which is precisely what we need for a credible marginal reference chart.

Define the probability integral transform

\[
U_{ij}(t)
=
F_j\!\left(
Y_{ij}(t)\mid x_{ij}(t)
\right),
\]

and the latent normal score

\[
Z_{ij}(t)
=
\Phi^{-1}\!\left(U_{ij}(t)\right).
\]

If the marginal model is calibrated, each \(Z_{ij}(t)\) is marginally standard normal. The longitudinal extension models their **joint dependence**:

\[
\mathbf Z_{ij}
=
\left[
Z_{ij}(t_1),\ldots,Z_{ij}(t_m)
\right]^\top
\sim
N\!\left(
\mathbf 0,
R_j(\mathbf t_i;\psi_j)
\right).
\]

This is a continuous-time Gaussian-copula process attached to the existing marginal distribution model. It is not neuroimaging-specific and does not require a separate normative-modeling algorithm.

### Conditioning on a person’s history

Suppose \(Z_H\) contains the person’s previous latent scores, and \(Z_*\) is the score at a future time. Partition the correlation matrix as

\[
R =
\begin{bmatrix}
R_{HH} & r_{H*}\\
r_{*H} & 1
\end{bmatrix}.
\]

Then

\[
Z_*\mid Z_H
\sim
N(m_*,s_*^2),
\]

where

\[
m_*
=
r_{*H}R_{HH}^{-1}Z_H,
\]

and

\[
s_*^2
=
1-r_{*H}R_{HH}^{-1}r_{H*}.
\]

The standardized innovation is

\[
z_{\mathrm{innovation}}
=
\frac{Z_*-m_*}{s_*}.
\]

With one previous observation, this reduces exactly to the paper’s Z-gain. With several previous observations, it gives the paper’s history-conditioned forecast, but directly from a coherent stochastic process rather than from a separately completed correlation table.

Adding valid history cannot increase \(s_*^2\) at fixed model parameters. The paper’s example in which the score changes from 1.78 using only the latest observation to 4.44 using the full history is precisely this conditioning effect.

---

## 3. The dynamic forecast is still a `norm_dist`

The especially attractive feature is that the dynamic result preserves our original predictive-distribution contract.

Let

\[
z(y)
=
\Phi^{-1}\!\left(F_*(y)\right)
\]

for the static marginal distribution at the forecast time. Since

\[
Z_*\mid H
\sim
N(m_*,s_*^2),
\]

the history-conditioned response CDF is

\[
\boxed{
F_*(y\mid H)
=
\Phi\!\left(
\frac{z(y)-m_*}{s_*}
\right).
}
\]

Its quantile function is

\[
\boxed{
Q_*(p\mid H)
=
F_*^{-1}
\left[
\Phi\left(
m_*+s_*\Phi^{-1}(p)
\right)
\right].
}
\]

And, when the marginal density \(f_*(y)\) exists, the conditional density is

\[
f_*(y\mid H)
=
f_*(y)
\frac{
\phi\!\left((z(y)-m_*)/s_*\right)
}{
s_*\phi(z(y))
}.
\]

Consequently, `norm_forecast()` can return exactly the same type of distribution object as `norm_fit()`:

```text
cdf()
quantile()
log_density()
draw()
center()
variance()
```

No special “velocity score engine” is needed. The dynamic process modifies the predictive distribution, after which all scores and graphics follow from the existing package laws.

---

## 4. Actual velocity is an affine transform of the forecast

For an observed baseline \(y_0\) at \(t_0\) and a future time \(t_1\), define the genuine finite-interval velocity

\[
V
=
\frac{Y(t_1)-y_0}{\Delta t},
\qquad
\Delta t=t_1-t_0.
\]

Conditional on the person’s history,

\[
G_V(v\mid H)
=
P(V\le v\mid H)
=
F_{t_1}\!\left(
y_0+\Delta t\,v
\mid H
\right).
\]

Its quantile function is simply

\[
Q_V(p\mid H)
=
\frac{
Q_{t_1}(p\mid H)-y_0
}{\Delta t}.
\]

Therefore the package can return:

- expected or median change in the original units;
- expected or median velocity in units per day, month, or year;
- an asymmetric predictive interval for velocity;
- the velocity centile;
- the normal-score representation of that centile.

For a fixed observed baseline, the normal-score velocity centile is the same conditional innovation as Z-gain. But the user also receives the actual rate and its predictive distribution, which the term “velocity” ordinarily implies.

---

## 5. We should distinguish six longitudinal quantities

| Quantity | Scientific question | Definition |
|---|---|---|
| **Chart velocity** | How does a reference centile itself move with time? | \(\partial Q_t(p)/\partial t\) |
| **Observed velocity** | At what raw rate did this person change? | \((y_2-y_1)/\Delta t\) |
| **Expected velocity** | What rate was expected given this person’s history? | Center of \(V\mid H\) |
| **Velocity centile** | Where does the observed rate fall in the conditional reference distribution? | \(G_V(v_{\rm obs}\mid H)\) |
| **Innovation Z** | How surprising is the follow-up after conditioning on history? | \(\Phi^{-1}\{G_V(v_{\rm obs}\mid H)\}\) |
| **Change Z** | How unusual is the change in the unconditional distribution of changes? | Standardized \(\Delta Y\), z-diff-like |

This distinction should be visible in the API and documentation.

The attached paper’s `Z-gain` maps principally onto **innovation Z**. Its “thrive lines” are conditional forecast quantiles. The earlier z-diff maps onto **change Z**. Our package should not force users to choose one method globally when the estimands differ.

---

## 6. Population velocity and individual rank velocity

Our distributional formulation also permits a useful decomposition that the paper only gestures toward.

Write the response as

\[
Y_i(t)
=
Q_t\!\left(
\Phi(Z_i(t))
\right),
\]

where \(Q_t\) is the marginal quantile function and \(Z_i(t)\) describes the person’s position in the reference distribution.

When derivatives are supported,

\[
\frac{dY_i(t)}{dt}
=
\underbrace{
\frac{\partial Q_t(u)}{\partial t}
}_{\text{reference-chart drift}}
+
\underbrace{
\frac{\partial Q_t(u)}{\partial u}
\phi(Z_i(t))
\frac{dZ_i(t)}{dt}
}_{\text{movement across centiles}}.
\]

This cleanly separates:

1. **Reference drift:** the outcome changes because the population chart itself changes while the person tracks a fixed centile.
2. **Rank drift:** the outcome changes because the person moves upward or downward relative to the reference population.

That is a potentially important methodological contribution. Two people can have the same raw velocity while differing sharply in these components.

For example, a decline of one unit per year may be:

- expected because the median curve is declining rapidly;
- unexpectedly rapid because the chart is nearly flat;
- a fixed-centile decline for one person;
- a major downward centile crossing for another.

The package should calculate finite-interval versions by default. Instantaneous derivatives should be available only when the process is differentiable and the longitudinal reference data are sufficiently informative.

---

## 7. A better longitudinal dependence model

The paper estimates Pearson correlations in one-year age bins, accepts an age-pair estimate when at least five individuals contribute both observations, Fisher-transforms those correlations, fits a six-parameter surface, and completes a correlation matrix within a ten-year band. That is pragmatic, but it should not be our computational foundation.

A fitted surface of pairwise correlations does not, by itself, guarantee a valid positive-semidefinite covariance matrix for every arbitrary set of visit times. One-year binning also discards exact timing, and correlations based on as few as five pairs will be extremely uncertain.

Instead, `referent` should make covariance validity a construction invariant.

### A useful default process

For outcome \(j\), define an unnormalized covariance

\[
K_j(t,s)
=
\tau_{b,j}^2
+
\tau_{g,j}^2
k_{3/2}\!\left(|t-s|;\ell_j\right)
+
\mathbf 1(t=s)\,
\sigma_{e,j}^2(w_t),
\]

where:

- \(\tau_b^2\) is stable person-specific rank;
- \(k_{3/2}\) is a Matérn-\(3/2\) temporal kernel;
- \(\tau_g^2\) is smooth within-person rank movement;
- \(\sigma_e^2(w_t)\) is measurement noise, potentially depending on site, device, quality, or protocol.

Normalize it to a correlation:

\[
R_j(t,s)
=
\frac{K_j(t,s)}
{\sqrt{K_j(t,t)K_j(s,s)}}.
\]

This construction is automatically positive semidefinite, supports irregular observation times, separates stable rank from gradual evolution and measurement noise, and permits derivative inference because the Matérn-\(3/2\) process has a mean-square derivative.

For each subject’s small vector of repeat observations, parameters can be estimated through the full Gaussian likelihood

\[
\ell_j(\psi_j)
=
\sum_i
\log
N\!\left(
Z_{ij};
0,
R_{ij}(\psi_j)
\right).
\]

Most longitudinal normative datasets have few observations per subject, so these subject-level Cholesky factorizations are inexpensive even with thousands of people.

### Nonstationary dynamics

When stability changes materially across the lifespan, offer a low-rank spline process:

\[
K_j(t,s)
=
B(t)^\top
L_jL_j^\top
B(s)
+
\mathbf 1(t=s)\sigma_{e,j}^2(t).
\]

Here \(B(t)\) is a small smooth basis and \(L_jL_j^\top\) guarantees positive semidefiniteness. This gives age-dependent correlation without binning and without constructing pairwise correlations separately.

It should be the more flexible rung in a short model ladder, not the universal default.

---

## 8. Measurement noise must be distinct from biological dynamics

A strong velocity system must distinguish:

\[
\text{observed change}
=
\text{latent within-person change}
+
\text{measurement error}.
\]

The earlier z-diff work makes this point especially clearly: common person-level variance cancels differently from visit-specific noise, which is why simple Z-score subtraction is not properly standardized.

Our dynamic kernel should therefore expose at least:

```text
stable_subject_variance
dynamic_process_variance
measurement_variance
```

The measurement component may depend on:

```r
measurement = ~ site + device + quality
```

At a new site, adaptation should be separable:

```r
local <- norm_adapt(
  dyn,
  data = local_reference,
  components = c("location", "scale", "measurement")
)
```

The shared temporal biology remains fixed unless the local longitudinal reference sample is sufficiently large. A small set of local repeat controls can recalibrate the measurement-noise component without pretending to re-estimate the entire lifespan process.

There is an identifiability limit: without short-interval repeats, replicate measurements, or external reliability information, rapid biological variability and measurement noise may not be distinguishable. The package should report that limitation rather than silently assign all short-range variance to one component.

---

## 9. Marginal normality is not joint Gaussianity

Mapping calibrated observations into Z-space guarantees marginal standard normality. It does **not** guarantee that the repeated observations have a Gaussian copula.

The source paper itself notes that reliable inference requires the Z-gain distribution to be Gaussian and free from residual site effects. Our package should treat this as an empirical requirement.

Dynamic assessment should test one-step-ahead innovations for:

- mean zero and variance one;
- calibrated tail rates;
- normal Q–Q and worm plots;
- residual dependence after conditioning;
- calibration across age, interval length, baseline centile, site, and history length;
- predictive interval coverage;
- conditional log score and CRPS.

Calibration of static pseudo-Z scores is already central to the GAMLSS framework. The same principle should be extended to **conditional innovations**.

A useful optional safeguard is subject-level split-conformal recalibration of forecast intervals. It would provide empirical marginal coverage under subject-level exchangeability even when the Gaussian copula is somewhat misspecified. The package should label model-based and conformally calibrated intervals separately.

A more flexible non-Gaussian transition engine can be considered later, but it should not complicate the initial core. The Gaussian process plus held-out conditional recalibration is a strong, transparent starting point.

---

## 10. Support should be temporal as well as cross-sectional

The dynamic object needs a support assessment distinct from ordinary covariate support:

```text
in
edge
unsupported_age
unsupported_lag
insufficient_history
new_measurement_context
```

The paper limits forecasts to ten years because most repeat observations fell within that range, and it cautions that thrive lines should match the actual follow-up period and that sparsely sampled ages require caution.

Our package should estimate and display:

- effective numbers of reference transitions near the requested age and lag;
- observed lag range;
- history-length support;
- uncertainty in the fitted dependence;
- whether instantaneous derivatives are identifiable;
- whether the requested visit schedule is interpolation or extrapolation.

A user should never receive a polished five-year velocity centile from a reference dataset containing almost exclusively one-year intervals without an explicit warning.

---

## 11. Proposed public API

The dynamic layer should introduce only a few public verbs.

```r
reference <- norm_fit(
  spec,
  data = reference_data,
  outcomes = starts_with("marker_")
)

dynamic <- norm_dynamics(
  reference,
  data = longitudinal_reference,
  id = participant_id,
  time = age,
  process = norm_matern32(
    stable_rank = TRUE,
    measurement = ~ site + quality
  ),
  crossfit = 5
)
```

### Reference-chart derivatives

```r
chart_velocity <- norm_derivative(
  reference,
  newdata = age_grid,
  with_respect_to = age,
  centiles = c(.05, .25, .50, .75, .95),
  type = "chart"
)
```

This estimates derivatives of the modeled quantile curves using the `mgcv` prediction matrix and coefficient uncertainty. It is explicitly a population-chart derivative, not an individual longitudinal estimate.

### History-conditioned forecast

```r
forecast <- norm_forecast(
  dynamic,
  history = subject_history,
  times = c(70.5, 71, 72)
)
```

This returns a `norm_forecast` containing predictive `norm_dist` objects at each requested time.

### Transition and velocity analysis

```r
transition <- norm_transition(
  dynamic,
  data = subject_visits,
  id = participant_id,
  time = age,
  conditioning = "all"
)

velocity <- norm_velocity(
  transition,
  time_unit = "year",
  scale = "response"
)
```

A row of the resulting `norm_transition` object would contain:

```text
.id
.outcome
.from_time
.to_time
.dt

.start_value
.end_value
.start_centile
.end_centile

.observed_change
.expected_change
.change_centile
.change_z

.observed_velocity
.expected_velocity
.velocity_lower
.velocity_upper
.velocity_centile

.innovation_z
.history_n

.aleatoric_sd
.measurement_sd
.epistemic_sd

.support
.calibrated
```

`z_gain` may be offered as an interoperability alias for `.innovation_z`, but it should not be the package’s primary internal name.

### Plotting

```r
autoplot(forecast, type = "fan")
autoplot(transition, type = "velocity")
autoplot(transition, type = "innovation")
autoplot(dynamic, type = "calibration")
```

A “thrive line” should simply be a plotting representation of a conditional quantile path:

\[
h\mapsto Q_{Y(t+h)\mid H}(p).
\]

The underlying object remains a forecast distribution rather than a collection of line segments tied to one arbitrary interval.

---

## 12. A restrained dynamic model ladder

The package should fit only as much temporal structure as the data support.

| Level | Process | Appropriate data | Permitted inference |
|---|---|---|---|
| 1 | Stable rank + nugget | Mostly two visits, narrow lag range | Two-point change and conditional innovation |
| 2 | Stable rank + continuous-time Matérn + nugget | Varied intervals, substantial repeat sample | Arbitrary-interval velocity and full-history forecasts |
| 3 | Low-rank nonstationary spline process | Broad lifespan coverage with repeated observations | Age-varying stability and richer trajectory forecasts |
| 4 | Differentiable individual process | Enough subjects with three or more informative visits | Instantaneous latent velocity and possibly acceleration |

The package should decline to estimate an instantaneous individual derivative when the data only identify a distribution of two-point changes.

For multi-outcome panels, the temporal parameters can optionally be stabilized through empirical-Bayes pooling across outcomes. This is especially useful when hundreds of features share broadly similar reliability and temporal scales, but each feature has sparse repeats. The outcome-specific marginal distributions remain separate.

The existing `norm_joint()` concept can then be applied to transition innovations:

```r
joint_change <- norm_joint(
  transition,
  value = "innovation_z",
  method = "gaussian_copula"
)
```

No separate multivariate velocity subsystem is required.

---

## 13. Essential mathematical and computational laws

The dynamic implementation should pass several exact or simulation-based laws.

### Paper-reduction law

With Gaussian marginals, one previous observation, and a fixed correlation \(r\),

\[
\texttt{innovation\_z}
=
\frac{Z_2-rZ_1}{\sqrt{1-r^2}}.
\]

It must agree numerically with the paper’s Z-gain.

### Difference-reduction law

Without conditioning on baseline,

\[
\texttt{change\_z}
=
\frac{Z_2-Z_1}{\sqrt{2(1-r)}}
\]

when model uncertainty is negligible. With parameter and measurement uncertainty, it must reduce to the appropriate generalized z-diff expression.

### History law

For a valid covariance model, adding an informative history point cannot increase the conditional variance at fixed parameters.

### Positive-definiteness law

Every covariance matrix generated for every valid irregular visit schedule must be positive semidefinite by construction.

### Time-unit law

Changing years to months rescales raw velocity by \(1/12\), while change centiles, velocity centiles, and innovation Z remain unchanged.

### Transformation law

Under a monotone transformation of the response, marginal centiles and innovation Z remain invariant, while raw change and raw velocity transform appropriately.

### Calibration law

Held-out reference trajectories must produce:

\[
z_{\mathrm{innovation}}
\sim N(0,1)
\]

both marginally and, within tolerances, conditional on age, lag, starting centile, site, and history length.

### Forecast-density law

For every dynamic forecast distribution,

\[
F\{Q(p)\}=p
\]

and numerical integration of the density must equal one.

---

## 14. The package-level design decision

I would make the following architectural change to the earlier proposal:

\[
\boxed{
\texttt{norm\_fit}
\longrightarrow
\texttt{norm\_dynamics}
\longrightarrow
\begin{cases}
\texttt{norm\_forecast}\\
\texttt{norm\_transition}\\
\texttt{norm\_velocity}\\
\texttt{norm\_derivative}
\end{cases}
}
\]

`norm_change()` would no longer be an isolated calculation. It would either become a convenience wrapper around `norm_transition()` or disappear from the primary API.

The package’s substantive methodological contribution would then be broader than “R GAMLSS for normative modeling”:

> **A calibrated distributional reference model whose static centiles, temporal dependence, actual rates of change, conditional innovations, and history-conditioned forecasts all obey one predictive-distribution interface.**

That is cleaner than the paper’s combination of marginal SHASH models, age-binned correlation estimates, matrix completion, and separate Z-gain formulas. It preserves the paper’s strongest insight—model the movement of an individual through a reference distribution—but gives it a valid continuous-time process, original-scale velocity distributions, explicit measurement error, full uncertainty propagation, honest support diagnostics, and a precise separation between chart drift, centile crossing, change, and predictive surprise.
