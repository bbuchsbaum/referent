# NHANES cohort evidence

The retained root-level artifacts fit adult body-measure reference
distributions on public CDC NHANES 2015-2016 data and evaluate them on
2017-2018 data. They are **post-hoc model-development evidence**, not untouched
external confirmation: the earlier 2017-2018 failure informed the addition of
sex and race/ethnicity to the scale model, after which the same fixed sample
was re-evaluated.

The lane is a cohort-transport check, not a U.S. prevalence estimate. It does
not apply NHANES survey weights, so results must not be described as nationally
representative. Raw XPT files are downloaded from CDC at runtime and are not
committed. The receipt records their URLs, SHA-256 hashes, cohort sizes,
package/R versions, and hashes of every result table.

The gate requires valid fits, at least 500 evaluated observations per outcome,
held-out provenance, and the package's full marginal, tail, shape, and
conditional calibration contract. Artifacts are uploaded even when the gate
fails, so cohort drift remains visible rather than disappearing from CI.

The retained development run used all 5,392 eligible training rows and a fixed-seed
sample of 2,000 of 5,151 eligible evaluation rows. All fits were valid. The
scale model conditions on age, sex, and race/ethnicity, matching the declared
conditional diagnostic contract. MACE ranged from 0.0048 to 0.0176 and 95%
coverage from 0.947 to 0.950 across the four outcomes. Marginal, tail, shape,
and conditional gates all passed. That result demonstrates that the revised
model addresses the observed diagnostic failure on the development cohort; it
does not support an independent conditional-transport claim.

Before any NHANES 2013-2014 outcomes were inspected, Mote issue
`referent-32n` and published commit `b1301a4` froze a genuinely untouched
confirmation: the already-fixed 2015-2016 model, 2013-2014 evaluation cycle,
eligibility rules, outcomes, sample size and seed, uncertainty mode, and
unchanged calibration gates.

That confirmation was then executed once. It used all 5,392 eligible training
rows and a fixed-seed sample of 2,000 of 5,592 eligible 2013-2014 participants;
the minimum outcome count was 1,892. All fits were valid and the marginal,
tail, shape, and conditional gates passed. Across outcomes, MACE ranged from
0.0049 to 0.0142 and 95% coverage from 0.944 to 0.955. The retained artifacts
are in `confirmation-2013-2014/`, and their receipt points back to the prior
Mote/Git registration. This supports the named, unweighted cross-cycle
conditional-transport claim; it is not a national prevalence claim or evidence
for arbitrary cohorts.

Run with:

```sh
Rscript tools/validation/run_nhanes_external_validation.R \
  evidence/nhanes-development /tmp/nhanes-cache development

Rscript tools/validation/run_nhanes_external_validation.R \
  evidence/nhanes-confirmation /tmp/nhanes-cache confirmation
```
