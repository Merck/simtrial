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
#>   sim_id method          parameter analysis cut_date   n event  estimate
#> 1      1    WLR FH(rho=0, gamma=0)        1 24.00000 400   245 -26.39845
#> 2      1    WLR FH(rho=0, gamma=0)        2 32.00000 400   301 -34.97803
#> 3      1    WLR FH(rho=0, gamma=0)        3 48.63895 400   350 -41.42441
#> 4      2    WLR FH(rho=0, gamma=0)        1 24.00000 400   235 -20.37909
#> 5      2    WLR FH(rho=0, gamma=0)        2 32.00000 400   306 -30.55728
#> 6      2    WLR FH(rho=0, gamma=0)        3 45.00000 400   354 -32.66420
#> 7      3    WLR FH(rho=0, gamma=0)        1 24.00000 400   251 -16.50549
#> 8      3    WLR FH(rho=0, gamma=0)        2 32.00000 400   311 -22.97248
#> 9      3    WLR FH(rho=0, gamma=0)        3 45.00000 400   351 -32.13795
#>         se        z     info info0
#> 1 7.770184 3.397404 60.13878 61.25
#> 2 8.564565 4.084041 74.23256 75.25
#> 3 9.126837 4.538748 87.01714 87.50
#> 4 7.649520 2.664100 57.85532 58.75
#> 5 8.662696 3.527456 75.76471 76.50
#> 6 9.231999 3.538151 88.31921 88.50
#> 7 7.907949 2.087202 62.22311 62.75
#> 8 8.766078 2.620611 77.32476 77.75
#> 9 9.215294 3.487458 87.37322 87.75

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
#> 1      1    WLR FH(rho=0, gamma=0.5)        1 24.00000 400   238 -11.519460
#> 2      1    WLR FH(rho=0, gamma=0.5)        2 32.00000 400   302 -21.050307
#> 3      1    WLR FH(rho=0, gamma=0.5)        3 47.15475 400   350 -25.657580
#> 4      2    WLR FH(rho=0, gamma=0.5)        1 24.00000 400   237 -17.779288
#> 5      2    WLR FH(rho=0, gamma=0.5)        2 32.00000 400   303 -17.104511
#> 6      2    WLR FH(rho=0, gamma=0.5)        3 45.20750 400   350 -21.259844
#> 7      3    WLR FH(rho=0, gamma=0.5)        1 24.00000 400   243  -7.641473
#> 8      3    WLR FH(rho=0, gamma=0.5)        2 32.00000 400   313 -16.222460
#> 9      3    WLR FH(rho=0, gamma=0.5)        3 45.00000 400   357 -13.462925
#>         se        z     info    info0
#> 1 4.212840 2.734369 17.66159 18.06453
#> 2 5.287754 3.980954 28.10427 28.80492
#> 3 6.011098 4.268368 38.33550 38.53873
#> 4 4.159616 4.274261 17.21443 17.86992
#> 5 5.253254 3.255984 28.85187 28.85822
#> 6 6.029132 3.526187 38.42666 38.42690
#> 7 4.315878 1.770549 18.64245 18.87826
#> 8 5.485709 2.957222 30.52415 30.87758
#> 9 6.241763 2.156911 39.98650 39.99242

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
#> 1      1    WLR MB(delay = 3, max_weight = Inf)        1 24.00000 400   237
#> 2      1    WLR MB(delay = 3, max_weight = Inf)        2 32.00000 400   292
#> 3      1    WLR MB(delay = 3, max_weight = Inf)        3 47.13787 400   350
#> 4      2    WLR MB(delay = 3, max_weight = Inf)        1 24.00000 400   228
#> 5      2    WLR MB(delay = 3, max_weight = Inf)        2 32.00000 400   302
#> 6      2    WLR MB(delay = 3, max_weight = Inf)        3 45.00000 400   356
#> 7      3    WLR MB(delay = 3, max_weight = Inf)        1 24.00000 400   232
#> 8      3    WLR MB(delay = 3, max_weight = Inf)        2 32.00000 400   296
#> 9      3    WLR MB(delay = 3, max_weight = Inf)        3 45.00000 400   355
#>    estimate        se        z      info     info0
#> 1 -23.99092  9.119518 2.630723  81.75319  83.61133
#> 2 -27.23502 10.132880 2.687786 103.25527 104.18572
#> 3 -28.70345 11.105597 2.584593 125.58628 125.88236
#> 4 -15.71567  9.376773 1.676022  87.95612  88.29343
#> 5 -30.27144 10.886050 2.780755 118.83819 119.70573
#> 6 -37.58441 11.785381 3.189071 143.09412 143.37252
#> 7 -40.22868  9.270549 4.339406  84.05534  87.29638
#> 8 -47.10964 10.480221 4.495100 111.98258 113.59487
#> 9 -59.88361 11.354711 5.273900 137.01487 137.83879

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
#> 3      1    WLR Xu 2017 with first 6 months of 0 weights        3 48.05091 400
#> 4      2    WLR Xu 2017 with first 6 months of 0 weights        1 24.00000 400
#> 5      2    WLR Xu 2017 with first 6 months of 0 weights        2 32.00000 400
#> 6      2    WLR Xu 2017 with first 6 months of 0 weights        3 45.00000 400
#> 7      3    WLR Xu 2017 with first 6 months of 0 weights        1 24.00000 400
#> 8      3    WLR Xu 2017 with first 6 months of 0 weights        2 32.00000 400
#> 9      3    WLR Xu 2017 with first 6 months of 0 weights        3 45.73872 400
#>   event   estimate       se         z     info info0
#> 1   219  -7.233088 4.459505 1.6219486 20.17284 20.25
#> 2   289  -5.322163 6.090125 0.8739003 37.66887 37.75
#> 3   350  -6.283640 7.184216 0.8746452 52.15311 52.25
#> 4   246 -17.168891 5.269039 3.2584480 28.26087 28.75
#> 5   312 -23.808172 6.567445 3.6251801 45.08287 45.25
#> 6   359 -26.581367 7.266207 3.6582179 56.50000 56.50
#> 7   255  -7.792209 5.190782 1.5011629 27.46364 27.50
#> 8   307  -8.396865 6.324591 1.3276535 40.49383 40.50
#> 9   350 -16.086951 7.076497 2.2732929 50.87745 51.00

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
#>   sim_id method parameter analysis cut_date   n event  estimate        se
#> 1      1   RMST        20        1 24.00000 400   246 0.9620073 0.7703594
#> 2      1   RMST        20        2 32.00000 400   303 1.0704830 0.7434335
#> 3      1   RMST        20        3 45.82652 400   350 1.0754780 0.7446373
#> 4      2   RMST        20        1 24.00000 400   224 1.3565622 0.7654899
#> 5      2   RMST        20        2 32.00000 400   294 1.3234404 0.7136717
#> 6      2   RMST        20        3 49.28090 400   350 1.3267670 0.7112636
#> 7      3   RMST        20        1 24.00000 400   236 2.8786445 0.7466703
#> 8      3   RMST        20        2 32.00000 400   302 2.8336838 0.7139365
#> 9      3   RMST        20        3 45.00000 400   353 2.8416720 0.7127162
#>          z
#> 1 1.248777
#> 2 1.439918
#> 3 1.444298
#> 4 1.772149
#> 5 1.854411
#> 6 1.865366
#> 7 3.855308
#> 8 3.969098
#> 9 3.987102

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
#>   sim_id    method parameter analysis cut_date   n event     estimate        se
#> 1      1 milestone        10        1 24.00000 400   248  0.410514199 0.1481902
#> 2      1 milestone        10        2 32.00000 400   297  0.414543282 0.1478856
#> 3      1 milestone        10        3 49.48058 400   350  0.414543282 0.1478856
#> 4      2 milestone        10        1 24.00000 400   252  0.245726113 0.1434768
#> 5      2 milestone        10        2 32.00000 400   311  0.236741525 0.1427782
#> 6      2 milestone        10        3 45.00000 400   365  0.236741525 0.1427782
#> 7      3 milestone        10        1 24.00000 400   244 -0.008264004 0.1458867
#> 8      3 milestone        10        2 32.00000 400   309  0.016998485 0.1446776
#> 9      3 milestone        10        3 45.82981 400   350  0.016998485 0.1446776
#>             z
#> 1  2.77018499
#> 2  2.80313434
#> 3  2.80313434
#> 4  1.71265391
#> 5  1.65810723
#> 6  1.65810723
#> 7 -0.05664674
#> 8  0.11749215
#> 9  0.11749215
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
#> 1      1    WLR FH(rho=0, gamma=0)        1       24 400   243  -8.214034
#> 2      1    WLR FH(rho=0, gamma=0)        2       32 400   304 -17.213140
#> 3      1    WLR FH(rho=0, gamma=0)        3       45 400   354 -20.877501
#> 4      2    WLR FH(rho=0, gamma=0)        1       24 400   259 -23.598942
#> 5      2    WLR FH(rho=0, gamma=0)        2       32 400   315 -24.300718
#> 6      2    WLR FH(rho=0, gamma=0)        3       45 400   363 -29.136121
#> 7      3    WLR FH(rho=0, gamma=0)        1       24 400   242  -8.625253
#> 8      3    WLR FH(rho=0, gamma=0)        2       32 400   311 -16.179327
#> 9      3    WLR FH(rho=0, gamma=0)        3       45 400   368 -20.646537
#>         se        z     info info0
#> 1 7.785096 1.055097 60.66667 60.75
#> 2 8.690860 1.980603 75.67105 76.00
#> 3 9.327781 2.238207 88.39831 88.50
#> 4 8.018639 2.943011 63.93822 64.75
#> 5 8.818336 2.755703 78.46349 78.75
#> 6 9.410491 3.096132 90.63361 90.75
#> 7 7.773275 1.109603 60.29752 60.50
#> 8 8.798732 1.838825 77.45981 77.75
#> 9 9.487108 2.176273 91.69482 91.75
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
#> 1      1    WLR FH(rho=0, gamma=0)        1 12.00002 428   103  -5.027225
#> 2      1    WLR FH(rho=0, gamma=0)        2 23.99062 524   250 -18.698773
#> 3      1    WLR FH(rho=0, gamma=0)        3 35.93241 524   326 -33.341788
#>         se         z     info info0 planned_lower_bound planned_upper_bound
#> 1 5.073530 0.9908732 25.45631 25.75          -1.7052708            3.870248
#> 2 7.894615 2.3685479 61.47600 62.50           0.9601286            2.356655
#> 3 8.979590 3.7130636 79.73313 81.50           2.0047521            2.009758
#>   updated_lower_bound updated_upper_bound
#> 1          -1.7474033            3.870248
#> 2           0.9922834            2.356669
#> 3           1.9788284            2.003621

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
#> 1      1    WLR FH(rho=0, gamma=0)        1 11.95079 415    75  -8.370586
#> 2      1    WLR FH(rho=0, gamma=0)        2 23.95510 496   215 -24.615573
#> 3      1    WLR FH(rho=0, gamma=0)        3 35.96077 496   281 -24.634761
#>         se        z     info info0 planned_lower_bound planned_upper_bound
#> 1 4.328620 1.933777 18.00000 18.75           -2.319759                  NA
#> 2 7.306296 3.369091 52.15814 53.75                  NA            2.358356
#> 3 8.331563 2.956800 69.60142 70.25                  NA            2.009328
#>   updated_lower_bound updated_upper_bound
#> 1           -2.635304                 Inf
#> 2                -Inf            2.423161
#> 3                -Inf            1.990588
# }
```
