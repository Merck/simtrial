# Simulate group sequential designs with fixed sample size

This function uses the option "stop" for the error-handling behavior of
the foreach loop. This will cause the entire function to stop when
errors are encountered and return the first error encountered instead of
returning errors for each individual simulation.

## Usage

``` r
sim_gs_n(
  n_sim = 1000,
  sample_size = 500,
  stratum = data.frame(stratum = "All", p = 1),
  enroll_rate = data.frame(duration = c(2, 2, 10), rate = c(3, 6, 9)),
  fail_rate = data.frame(stratum = "All", duration = c(3, 100), fail_rate = log(2)/c(9,
    18), hr = c(0.9, 0.6), dropout_rate = rep(0.001, 2)),
  block = rep(c("experimental", "control"), 2),
  test = wlr,
  cut = NULL,
  original_design = NULL,
  ia_alpha_spending = c("min_planned_actual", "actual"),
  fa_alpha_spending = c("full_alpha", "info_frac"),
  ...
)
```

## Arguments

- n_sim:

  Number of simulations to perform.

- sample_size:

  Total sample size per simulation.

- stratum:

  A data frame with stratum specified in `stratum`, probability
  (incidence) of each stratum in `p`.

- enroll_rate:

  Piecewise constant enrollment rates by time period. Note that these
  are overall population enrollment rates and the `stratum` argument
  controls the random distribution between stratum.

- fail_rate:

  Piecewise constant control group failure rates, hazard ratio for
  experimental vs. control, and dropout rates by stratum and time
  period.

- block:

  As in
  [`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md).
  Vector of treatments to be included in each block.

- test:

  One or more test functions such as
  [`wlr()`](https://merck.github.io/simtrial/reference/wlr.md),
  [`rmst()`](https://merck.github.io/simtrial/reference/rmst.md), or
  [`milestone()`](https://merck.github.io/simtrial/reference/milestone.md)
  ([`maxcombo()`](https://merck.github.io/simtrial/reference/maxcombo.md)
  can only be applied by itself). If a single test function is provided,
  it will be applied at each cut. Alternatively a list of functions
  created by
  [`create_test()`](https://merck.github.io/simtrial/reference/create_test.md).
  The list form is experimental and currently limited. It only accepts
  one test per cutting (in the future multiple tests may be accepted),
  and all the tests must consistently return the same exact results
  (again this may be more flexible in the future). Importantly, note
  that the simulated data set is always passed as the first positional
  argument to each test function provided.

- cut:

  A list of cutting functions created by
  [`create_cut()`](https://merck.github.io/simtrial/reference/create_cut.md),
  see examples. Alternatively, if cut is `NULL` (the default) and a
  design object is provided via the argument `original_design`, the cut
  functions are automatically created from this object.

- original_design:

  A design object from the gsDesign2 package, which is required when
  users want to calculate updated bounds. The default is NULL leaving
  the updated bounds uncalculated.

- ia_alpha_spending:

  Spend alpha at interim analysis based on

  - `"min_planned_actual"`: the minimal of planned and actual alpha
    spending.

  - `"actual"`: the actual alpha spending.

- fa_alpha_spending:

  If targeted final event count is not achieved (under-running at final
  analysis), specify how to do final spending. Generally, this should be
  specified in analysis plan.

  - `"info_frac"` = spend final alpha according to final information
    fraction

  - `"full_alpha"` = spend full alpha at final analysis.

- ...:

  Arguments passed to the test function(s) provided by the argument
  `test`.

## Value

A data frame summarizing the simulation ID, analysis date, z statistics
or p-values.

## Details

WARNING: This experimental function is a work-in-progress. The function
arguments will change as we add additional features.

## Examples

``` r
library(gsDesign2)

# Parameters for enrollment
enroll_rampup_duration <- 4 # Duration for enrollment ramp up
enroll_duration <- 16 # Total enrollment duration
enroll_rate <- define_enroll_rate(
  duration = c(
    enroll_rampup_duration,
    enroll_duration - enroll_rampup_duration
  ),
  rate = c(10, 30)
)

# Parameters for treatment effect
delay_effect_duration <- 3 # Delay treatment effect in months
median_ctrl <- 9 # Survival median of the control arm
median_exp <- c(9, 14) # Survival median of the experimental arm
dropout_rate <- 0.001
fail_rate <- define_fail_rate(
  duration = c(delay_effect_duration, 100),
  fail_rate = log(2) / median_ctrl,
  hr = median_ctrl / median_exp,
  dropout_rate = dropout_rate
)

# Other related parameters
alpha <- 0.025 # Type I error
beta <- 0.1 # Type II error
ratio <- 1 # Randomization ratio (experimental:control)

# Define cuttings of 2 IAs and 1 FA
# IA1
# The 1st interim analysis will occur at the later of the following 3 conditions:
# - At least 20 months have passed since the start of the study.
# - At least 100 events have occurred.
# - At least 20 months have elapsed after enrolling 200/400 subjects, with a
#   minimum of 20 months follow-up.
# However, if events accumulation is slow, we will wait for a maximum of 24 months.
ia1_cut <- create_cut(
  planned_calendar_time = 20,
  target_event_overall = 100,
  max_extension_for_target_event = 24,
  min_n_overall = 200,
  min_followup = 20
)

# IA2
# The 2nd interim analysis will occur at the later of the following 3 conditions:
# - At least 32 months have passed since the start of the study.
# - At least 200 events have occurred.
# - At least 10 months after IA1.
# However, if events accumulation is slow, we will wait for a maximum of 34 months.
ia2_cut <- create_cut(
  planned_calendar_time = 32,
  target_event_overall = 200,
  max_extension_for_target_event = 34,
  min_time_after_previous_analysis = 10
)

# FA
# The final analysis will occur at the later of the following 2 conditions:
# - At least 45 months have passed since the start of the study.
# - At least 350 events have occurred.
fa_cut <- create_cut(
  planned_calendar_time = 45,
  target_event_overall = 350
)

# Example 1: regular logrank test at all 3 analyses
sim_gs_n(
  n_sim = 3,
  sample_size = 400,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  test = wlr,
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
  weight = fh(rho = 0, gamma = 0)
)
#> Backend uses sequential processing.
#>   sim_id method          parameter analysis cut_date   n event   estimate
#> 1      1    WLR FH(rho=0, gamma=0)        1 24.00000 400   238 -14.341728
#> 2      1    WLR FH(rho=0, gamma=0)        2 32.00000 400   302 -25.890696
#> 3      1    WLR FH(rho=0, gamma=0)        3 47.15475 400   350 -31.002479
#> 4      2    WLR FH(rho=0, gamma=0)        1 24.00000 400   237 -28.176482
#> 5      2    WLR FH(rho=0, gamma=0)        2 32.00000 400   303 -27.607375
#> 6      2    WLR FH(rho=0, gamma=0)        3 45.20750 400   350 -31.998390
#> 7      3    WLR FH(rho=0, gamma=0)        1 24.00000 400   243  -8.083845
#> 8      3    WLR FH(rho=0, gamma=0)        2 32.00000 400   313 -18.518508
#> 9      3    WLR FH(rho=0, gamma=0)        3 45.00000 400   357 -15.358403
#>         se        z     info info0
#> 1 7.662756 1.871615 58.78481 59.25
#> 2 8.617322 3.004494 74.65232 75.50
#> 3 9.187449 3.374438 87.08857 87.50
#> 4 7.632234 3.691774 57.80591 59.25
#> 5 8.591683 3.213267 75.38614 75.75
#> 6 9.200596 3.477860 87.26857 87.50
#> 7 7.770050 1.040385 60.57613 60.75
#> 8 8.785069 2.107952 77.89776 78.25
#> 9 9.375340 1.638170 89.21569 89.25

# \donttest{
# Example 2: weighted logrank test by FH(0, 0.5) at all 3 analyses
sim_gs_n(
  n_sim = 3,
  sample_size = 400,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  test = wlr,
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
  weight = fh(rho = 0, gamma = 0.5)
)
#> Backend uses sequential processing.
#>   sim_id method            parameter analysis cut_date   n event   estimate
#> 1      1    WLR FH(rho=0, gamma=0.5)        1 24.00000 400   237 -14.458924
#> 2      1    WLR FH(rho=0, gamma=0.5)        2 32.00000 400   292 -16.365269
#> 3      1    WLR FH(rho=0, gamma=0.5)        3 47.13787 400   350 -17.603912
#> 4      2    WLR FH(rho=0, gamma=0.5)        1 24.00000 400   228  -6.217993
#> 5      2    WLR FH(rho=0, gamma=0.5)        2 32.00000 400   302 -14.879648
#> 6      2    WLR FH(rho=0, gamma=0.5)        3 45.00000 400   356 -19.941165
#> 7      3    WLR FH(rho=0, gamma=0.5)        1 24.00000 400   232 -18.771404
#> 8      3    WLR FH(rho=0, gamma=0.5)        2 32.00000 400   296 -22.955420
#> 9      3    WLR FH(rho=0, gamma=0.5)        3 45.00000 400   355 -31.838547
#>         se        z     info    info0
#> 1 4.195920 3.445948 16.89905 17.77081
#> 2 5.104649 3.205954 26.45812 26.70502
#> 3 6.091924 2.889713 38.31809 38.33033
#> 4 4.064090 1.529984 16.59433 16.63284
#> 5 5.300434 2.807251 28.37737 28.56080
#> 6 6.147001 3.244048 39.86389 39.86789
#> 7 4.094504 4.584536 16.45312 17.19726
#> 8 5.116139 4.486864 27.47104 27.62348
#> 9 5.976103 5.327643 39.70820 39.71523

# Example 3: weighted logrank test by MB(3) at all 3 analyses
sim_gs_n(
  n_sim = 3,
  sample_size = 400,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  test = wlr,
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
  weight = mb(delay = 3)
)
#> Backend uses sequential processing.
#>   sim_id method                       parameter analysis cut_date   n event
#> 1      1    WLR MB(delay = 3, max_weight = Inf)        1 24.00000 400   219
#> 2      1    WLR MB(delay = 3, max_weight = Inf)        2 32.00000 400   289
#> 3      1    WLR MB(delay = 3, max_weight = Inf)        3 48.05091 400   350
#> 4      2    WLR MB(delay = 3, max_weight = Inf)        1 24.00000 400   246
#> 5      2    WLR MB(delay = 3, max_weight = Inf)        2 32.00000 400   312
#> 6      2    WLR MB(delay = 3, max_weight = Inf)        3 45.00000 400   359
#> 7      3    WLR MB(delay = 3, max_weight = Inf)        1 24.00000 400   255
#> 8      3    WLR MB(delay = 3, max_weight = Inf)        2 32.00000 400   307
#> 9      3    WLR MB(delay = 3, max_weight = Inf)        3 45.73872 400   350
#>    estimate        se        z      info     info0
#> 1 -21.56149  9.127796 2.362180  83.37872  84.02876
#> 2 -19.06355 10.616555 1.795644 113.88678 113.93178
#> 3 -20.32039 11.727229 1.732753 138.69425 138.70858
#> 4 -35.30840  9.559822 3.693415  91.07733  93.04918
#> 5 -43.82029 10.800501 4.057246 118.99073 120.16949
#> 6 -47.37567 11.512562 4.115128 138.26765 138.66062
#> 7 -15.28769  9.794462 1.560850  96.65626  96.85809
#> 8 -16.06341 10.836071 1.482402 118.14499 118.25449
#> 9 -25.92916 11.576048 2.239897 135.18593 135.53621

# Example 4: weighted logrank test by early zero (6) at all 3 analyses
sim_gs_n(
  n_sim = 3,
  sample_size = 400,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  test = wlr,
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
  weight = early_zero(6)
)
#> Backend uses sequential processing.
#>   sim_id method                                parameter analysis cut_date   n
#> 1      1    WLR Xu 2017 with first 6 months of 0 weights        1 24.00000 400
#> 2      1    WLR Xu 2017 with first 6 months of 0 weights        2 32.00000 400
#> 3      1    WLR Xu 2017 with first 6 months of 0 weights        3 45.82652 400
#> 4      2    WLR Xu 2017 with first 6 months of 0 weights        1 24.00000 400
#> 5      2    WLR Xu 2017 with first 6 months of 0 weights        2 32.00000 400
#> 6      2    WLR Xu 2017 with first 6 months of 0 weights        3 49.28090 400
#> 7      3    WLR Xu 2017 with first 6 months of 0 weights        1 24.00000 400
#> 8      3    WLR Xu 2017 with first 6 months of 0 weights        2 32.00000 400
#> 9      3    WLR Xu 2017 with first 6 months of 0 weights        3 45.00000 400
#>   event  estimate       se        z     info info0
#> 1   246 -10.59987 5.155417 2.056065 26.22430 26.75
#> 2   303 -17.76063 6.355293 2.794621 40.39024 41.00
#> 3   350 -26.17475 7.141927 3.664942 52.12322 52.75
#> 4   224 -13.86612 5.032473 2.755329 24.31373 25.50
#> 5   294 -24.79518 6.505262 3.811558 41.51163 43.00
#> 6   350 -31.69640 7.360381 4.306354 56.36842 57.00
#> 7   236 -10.75188 5.038758 2.133835 25.91346 26.00
#> 8   302 -19.83549 6.346874 3.125239 42.47647 42.50
#> 9   353 -19.55176 7.165354 2.728653 54.99548 55.25

# Example 5: RMST at all 3 analyses
sim_gs_n(
  n_sim = 3,
  sample_size = 400,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  test = rmst,
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
  tau = 20
)
#> Backend uses sequential processing.
#>   sim_id method parameter analysis cut_date   n event    estimate        se
#> 1      1   RMST        20        1 24.00000 400   248  2.50301264 0.7518899
#> 2      1   RMST        20        2 32.00000 400   297  2.64873796 0.7354022
#> 3      1   RMST        20        3 49.48058 400   350  2.65073649 0.7360953
#> 4      2   RMST        20        1 24.00000 400   252  1.41183757 0.7595351
#> 5      2   RMST        20        2 32.00000 400   311  1.45891219 0.7383432
#> 6      2   RMST        20        3 45.00000 400   365  1.46871219 0.7379790
#> 7      3   RMST        20        1 24.00000 400   244 -0.17812793 0.7681930
#> 8      3   RMST        20        2 32.00000 400   309  0.03598688 0.7355750
#> 9      3   RMST        20        3 45.82981 400   350  0.04482531 0.7359993
#>             z
#> 1  3.32896146
#> 2  3.60175433
#> 3  3.60107786
#> 4  1.85881819
#> 5  1.97592691
#> 6  1.99018154
#> 7 -0.23187914
#> 8  0.04892347
#> 9  0.06090401

# Example 6: Milestone at all 3 analyses
sim_gs_n(
  n_sim = 3,
  sample_size = 400,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  test = milestone,
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
  ms_time = 10
)
#> Backend uses sequential processing.
#>   sim_id    method parameter analysis cut_date   n event  estimate        se
#> 1      1 milestone        10        1       24 400   243 0.2532512 0.1476521
#> 2      1 milestone        10        2       32 400   304 0.2474294 0.1465762
#> 3      1 milestone        10        3       45 400   354 0.2474294 0.1465762
#> 4      2 milestone        10        1       24 400   259 0.3482071 0.1423805
#> 5      2 milestone        10        2       32 400   315 0.3481048 0.1423502
#> 6      2 milestone        10        3       45 400   363 0.3481048 0.1423502
#> 7      3 milestone        10        1       24 400   242 0.1988469 0.1467060
#> 8      3 milestone        10        2       32 400   311 0.1704664 0.1430125
#> 9      3 milestone        10        3       45 400   368 0.1704664 0.1430125
#>          z
#> 1 1.715189
#> 2 1.688060
#> 3 1.688060
#> 4 2.445610
#> 5 2.445411
#> 6 2.445411
#> 7 1.355411
#> 8 1.191969
#> 9 1.191969
# }

# Warning: this example will be executable when we add info info0 to the milestone test
# Example 7: WLR with fh(0, 0.5) test at IA1,
# WLR with mb(6, Inf) at IA2, and milestone test at FA
ia1_test <- create_test(wlr, weight = fh(rho = 0, gamma = 0.5))
ia2_test <- create_test(wlr, weight = mb(delay = 6, w_max = Inf))
fa_test <- create_test(milestone, ms_time = 10)
if (FALSE) { # \dontrun{
sim_gs_n(
  n_sim = 3,
  sample_size = 400,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  test = list(ia1 = ia1_test, ia2 = ia2_test, fa = fa_test),
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut)
)
} # }

# WARNING: Multiple tests per cut will be enabled in a future version.
#          Currently does not work.
# Example 8: At IA1, we conduct 3 tests, LR, WLR with fh(0, 0.5), and RMST test.
# At IA2, we conduct 2 tests, LR and WLR with early zero (6).
# At FA, we conduct 2 tests, LR and milestone test.
ia1_test <- list(
  test1 = create_test(wlr, weight = fh(rho = 0, gamma = 0)),
  test2 = create_test(wlr, weight = fh(rho = 0, gamma = 0.5)),
  test3 = create_test(rmst, tau = 20)
)
ia2_test <- list(
  test1 = create_test(wlr, weight = fh(rho = 0, gamma = 0)),
  test2 = create_test(wlr, weight = early_zero(6))
)
fa_test <- list(
  test1 = create_test(wlr, weight = fh(rho = 0, gamma = 0)),
  test3 = create_test(milestone, ms_time = 20)
)
if (FALSE) { # \dontrun{
sim_gs_n(
  n_sim = 3,
  sample_size = 400,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  test = list(ia1 = ia1_test, ia2 = ia2_test, fa = fa_test),
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut)
)
} # }

# \donttest{
# Example 9: regular logrank test at all 3 analyses in parallel
plan("multisession", workers = 2)
sim_gs_n(
  n_sim = 3,
  sample_size = 400,
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  test = wlr,
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
  weight = fh(rho = 0, gamma = 0)
)
#> Using 2 cores with backend multisession
#>   sim_id method          parameter analysis cut_date   n event   estimate
#> 1      1    WLR FH(rho=0, gamma=0)        1       24 400   238 -12.713173
#> 2      1    WLR FH(rho=0, gamma=0)        2       32 400   298 -18.420232
#> 3      1    WLR FH(rho=0, gamma=0)        3       45 400   352 -23.852666
#> 4      2    WLR FH(rho=0, gamma=0)        1       24 400   252  -7.149188
#> 5      2    WLR FH(rho=0, gamma=0)        2       32 400   306 -14.941903
#> 6      2    WLR FH(rho=0, gamma=0)        3       45 400   356 -20.174083
#> 7      3    WLR FH(rho=0, gamma=0)        1       24 400   258 -17.450839
#> 8      3    WLR FH(rho=0, gamma=0)        2       32 400   319 -22.495646
#> 9      3    WLR FH(rho=0, gamma=0)        3       45 400   357 -27.876167
#>         se         z     info info0
#> 1 7.682310 1.6548633 59.23109 59.50
#> 2 8.593728 2.1434507 74.16443 74.50
#> 3 9.292240 2.5669447 87.81818 88.00
#> 4 7.926459 0.9019398 62.90079 63.00
#> 5 8.725428 1.7124550 76.17320 76.50
#> 6 9.346423 2.1584817 88.82022 89.00
#> 7 8.016516 2.1768607 63.94186 64.50
#> 8 8.887344 2.5312000 79.40439 79.75
#> 9 9.315263 2.9925261 88.82022 89.00
plan("sequential")

# Example 10: group sequential design with updated bounds -- efficacy only
x <- gs_design_ahr(analysis_time = 1:3*12) |> to_integer()
sim_gs_n(
  n_sim = 1,
  sample_size = max(x$analysis$n),
  enroll_rate = x$enroll_rate,
  fail_rate = x$fail_rate,
  test = wlr,
  cut = list(ia1 = create_cut(planned_calendar_time = x$analysis$time[1]),
             ia2 = create_cut(planned_calendar_time = x$analysis$time[2]),
             fa = create_cut(planned_calendar_time = x$analysis$time[3])),
  weight = fh(rho = 0, gamma = 0),
  original_design = x
)
#> Backend uses sequential processing.
#>   sim_id method          parameter analysis cut_date   n event   estimate
#> 1      1    WLR FH(rho=0, gamma=0)        1 12.00002 421    90  -2.835018
#> 2      1    WLR FH(rho=0, gamma=0)        2 23.99062 524   241  -9.823122
#> 3      1    WLR FH(rho=0, gamma=0)        3 35.93241 524   320 -17.474021
#>         se         z     info info0 planned_lower_bound planned_upper_bound
#> 1 4.742932 0.5977352 22.40000 22.50          -1.7052708            3.870248
#> 2 7.760276 1.2658212 59.95021 60.25           0.9601286            2.356655
#> 3 8.937393 1.9551587 79.38750 80.00           2.0047521            2.009758
#>   updated_lower_bound updated_upper_bound
#> 1          -2.0109114            4.074501
#> 2           0.9868377            2.356220
#> 3           1.9995136            2.006884

# Example 11: group sequential design with updated bounds -- efficacy & futility
x <- gs_design_ahr(
 alpha = 0.025, beta = 0.1, analysis_time = 1:3*12,
 upper = gs_spending_bound, upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025),
 lower = gs_spending_bound, lpar = list(sf = gsDesign::sfHSD, param = -4, total_spend = 0.01),
 test_upper = c(FALSE, TRUE, TRUE), test_lower = c(TRUE, FALSE, FALSE)) |> to_integer()
sim_gs_n(
  n_sim = 1,
  sample_size = max(x$analysis$n),
  enroll_rate = x$enroll_rate,
  fail_rate = x$fail_rate,
  test = wlr,
  cut = list(ia1 = create_cut(planned_calendar_time = x$analysis$time[1]),
             ia2 = create_cut(planned_calendar_time = x$analysis$time[2]),
             fa = create_cut(planned_calendar_time = x$analysis$time[3])),
  weight = fh(rho = 0, gamma = 0),
  original_design = x
)
#> Backend uses sequential processing.
#>   sim_id method          parameter analysis cut_date   n event   estimate
#> 1      1    WLR FH(rho=0, gamma=0)        1 11.95079 426    91  -1.846892
#> 2      1    WLR FH(rho=0, gamma=0)        2 23.95510 496   211 -12.647993
#> 3      1    WLR FH(rho=0, gamma=0)        3 35.96077 496   311 -19.050336
#>         se         z     info info0 planned_lower_bound planned_upper_bound
#> 1 4.767196 0.3874169 22.72527 22.75           -2.319759                  NA
#> 2 7.258625 1.7424778 52.22749 52.75                  NA            2.358356
#> 3 8.801474 2.1644485 77.24759 77.75                  NA            2.009328
#>   updated_lower_bound updated_upper_bound
#> 1           -2.361836                 Inf
#> 2                -Inf            2.450308
#> 3                -Inf            2.001302
# }
```
