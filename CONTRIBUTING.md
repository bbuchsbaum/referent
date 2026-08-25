# Contributing to referent

Bug reports should include a minimal data-generating example, the intended
predictive estimand (`conditional` or `total`), package versions, and whether
the data are training, calibration, or held-out evaluation rows. Never include
identifiable participant data.

For code changes:

1. Add a regression test that states the statistical contract or invariant.
2. Run `devtools::test()` and `R CMD check --as-cran` under a UTF-8 locale.
3. For numerical changes, include an independent oracle, metamorphic law, or
   adversarial counterexample; do not update a tolerance only to make a result
   pass.
4. Keep PCNtoolkit evidence classes and uncertainty estimands explicit. A
   better held-out score never excuses a conformance or calibration failure.
5. Do not commit raw external cohort data. Commit only small result tables and
   receipts with source URLs and hashes.

Pull requests that alter public functions, bundle schemas, comparison gates,
or evidence classifications should update the relevant vignette and NEWS.
