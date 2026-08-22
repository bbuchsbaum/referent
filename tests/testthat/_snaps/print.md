# print.ref_fit shows family, formulas, and status counts

    Code
      print(fit)
    Message
      <ref_fit> gaussian via mgcv
      3 outcomes, n = 120
      covariates: age and sex
      status: ok=2, insufficient_variation=1
      ! Failed or flagged outcome: marker_bad

---

    Code
      print(fit$spec)
    Message
      <ref_spec> gaussian via mgcv
      location: ~s(age, k = 5) + sex
      scale: ~1

---

    Code
      print(ref_shash())
    Message
      <ref_family> shash (4 parameters)

# print.ref_scores counts rows and warns about in-sample scores

    Code
      print(sc[, c(".row", ".outcome", "z", "support", "status")])
    Message
      <ref_scores> 9 rows
    Output
      # A tibble: 9 x 5
         .row .outcome         z support status                
        <int> <chr>        <dbl> <chr>   <chr>                 
      1     1 y           1.31   in      ok                    
      2     2 y          -1.30   in      ok                    
      3     3 y          -0.0248 in      ok                    
      4     1 marker_01   0.807  in      ok                    
      5     2 marker_01  -1.62   in      ok                    
      6     3 marker_01   0.0934 in      ok                    
      7     1 marker_bad NA      in      insufficient_variation
      8     2 marker_bad NA      in      insufficient_variation
      9     3 marker_bad NA      in      insufficient_variation

# print.ref_freeze prints the model card

    Code
      print(ref)
    Message
      
      -- referent model card ---------------------------------------------------------
      family: gaussian
      engine: mgcv
      outcomes: y, marker_01, and marker_bad
      covariates: age and sex
      n: 120
      missing-data policy: complete-case per outcome; predictors are never imputed
      package x.y.z, mgcv x.y, R x.y.z
      statuses: y=ok, marker_01=ok, marker_bad=insufficient_variation
      Covariate ranges
      age: [20.21, 79.85]

# print.ref_assessment shows the overall table

    Code
      print(a)
    Message
      <ref_assessment> n = 60
    Output
      # A tibble: 3 x 9
       .outcome mean_log_score standardized_log_score crps mae rmse smse
       <chr> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
      1 marker_01 # # # # # #
      2 marker_bad NaN NA NA NaN NaN NA 
      3 y # # # # # #
      # i 2 more variables: ev <dbl>, cor <dbl>

# print.ref_dynamics reports components, identifiability, and the kernel

    Code
      print(dyn)
    Message
      <ref_dynamics> 1 outcome; requested kernel: matern32; Z: in_sample
      subjects: 200; time range [#, #]; lag range [2, 2]
      ! measurement not separated: no short-interval repeats; measurement noise is not separated from rank dynamics
    Output
      # A tibble: 1 x 12
       .outcome process identified at_boundary stable dynamic measurement ell
       <chr> <chr> <lgl> <lgl> <dbl> <dbl> <dbl> <dbl>
      1 y stable TRUE FALSE # 0 # NA
      # i 4 more variables: ell_identified <lgl>, median_lag <dbl>,
      # r_median_lag <dbl>, r_median_lag_se <dbl>

---

    Code
      print(dyn0)
    Message
      <ref_dynamics> 1 outcome; requested kernel: matern32; Z: in_sample
      subjects: 200; time range [#, #]; lag range [0, 0]
      ! change not identified: too few repeated observations to identify within-person dependence
    Output
      # A tibble: 1 x 12
       .outcome process identified at_boundary stable dynamic measurement ell
       <chr> <chr> <lgl> <lgl> <dbl> <dbl> <dbl> <dbl>
      1 y matern32 FALSE FALSE NA NA NA NA
      # i 4 more variables: ell_identified <lgl>, median_lag <dbl>,
      # r_median_lag <dbl>, r_median_lag_se <dbl>

# remaining print methods use cli and return invisibly

    Code
      print(ad$adaptation)
    Message
      <ref_adaptation> parameters: location
      local n = 30
      y / A: location #, scale x1 (n = 10)
      y / B: location #, scale x1 (n = 8)
      y / C: location #, scale x1 (n = 4)
      y / D: location #, scale x1 (n = 8)

---

    Code
      print(cal$calibration)
    Message
      <ref_calibration> n = 50, by site

---

    Code
      print(jt)
    Message
      <ref_joint> gaussian copula, 1 outcome, n = 49 reference subjects

