# NHANES external-cohort evidence

This release lane fits adult body-measure reference distributions on the
public CDC NHANES 2015-2016 cohort and evaluates them once on the independent
2017-2018 cohort. The outcomes are height, weight, body-mass index, and waist
circumference; age, sex, and race/ethnicity are declared covariates.

The lane is a cohort-transport check, not a U.S. prevalence estimate. It does
not apply NHANES survey weights, so results must not be described as nationally
representative. Raw XPT files are downloaded from CDC at runtime and are not
committed. The receipt records their URLs, SHA-256 hashes, cohort sizes,
package/R versions, and hashes of every result table.

The gate requires valid fits, at least 500 evaluated observations per outcome,
held-out provenance, and the package's full marginal, tail, shape, and
conditional calibration contract. Artifacts are uploaded even when the gate
fails, so cohort drift remains visible rather than disappearing from CI.

The retained 0.1.0 run used all 5,392 eligible training rows and a fixed-seed
sample of 2,000 of 5,151 eligible evaluation rows. All fits were valid.
Marginal performance was strong across the four outcomes: MACE ranged from
0.005 to 0.015 and 95% coverage from 0.942 to 0.951. The full gate nevertheless
failed because conditional scale drift by race/ethnicity exceeded the declared
practical thresholds. The result supports marginal cohort transport for this
sample, but not a blanket conditional-calibration claim.

Run with:

```sh
Rscript tools/validation/run_nhanes_external_validation.R \
  evidence/nhanes-external /tmp/nhanes-cache
```
