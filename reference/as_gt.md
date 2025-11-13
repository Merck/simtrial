# Convert summary table to a gt object

Convert summary table to a gt object

## Usage

``` r
as_gt(x, ...)

# S3 method for class 'simtrial_gs_wlr'
as_gt(
  x,
  title = "Summary of simulation results by WLR tests",
  subtitle = NULL,
  ...
)
```

## Arguments

- x:

  A object returned by
  [`summary()`](https://rdrr.io/r/base/summary.html).

- ...:

  Additional parameters (not used).

- title:

  Title of the gt table.

- subtitle:

  Subtitle of the gt table.

## Value

A gt table.

A gt table summarizing the simulation results.

## Examples

``` r
# Parameters for enrollment
enroll_rampup_duration <- 4 # Duration for enrollment ramp up
enroll_duration <- 16 # Total enrollment duration
enroll_rate <- gsDesign2::define_enroll_rate(
  duration = c(
    enroll_rampup_duration, enroll_duration - enroll_rampup_duration),
 rate = c(10, 30))

# Parameters for treatment effect
delay_effect_duration <- 3 # Delay treatment effect in months
median_ctrl <- 9 # Survival median of the control arm
median_exp <- c(9, 14) # Survival median of the experimental arm
dropout_rate <- 0.001
fail_rate <- gsDesign2::define_fail_rate(
  duration = c(delay_effect_duration, 100),
  fail_rate = log(2) / median_ctrl,
  hr = median_ctrl / median_exp,
  dropout_rate = dropout_rate)

# Other related parameters
alpha <- 0.025 # Type I error
beta <- 0.1 # Type II error
ratio <- 1 # Randomization ratio (experimental:control)

# Build a one-sided group sequential design
design <- gsDesign2::gs_design_ahr(
  enroll_rate = enroll_rate, fail_rate = fail_rate,
  ratio = ratio, alpha = alpha, beta = beta,
  analysis_time = c(12, 24, 36),
  upper = gsDesign2::gs_spending_bound,
  upar = list(sf = gsDesign::sfLDOF, total_spend = alpha),
  lower = gsDesign2::gs_b,
  lpar = rep(-Inf, 3))

# Define cuttings of 2 IAs and 1 FA
ia1_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[1]))
ia2_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[2]))
fa_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[3]))

# Run simulations
simulation <- sim_gs_n(
  n_sim = 3,
  sample_size = ceiling(design$analysis$n[3]),
  enroll_rate = design$enroll_rate,
  fail_rate = design$fail_rate,
  test = wlr,
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
  weight = fh(rho = 0, gamma = 0.5))
#> Backend uses sequential processing.

# Summarize simulations
simulation |>
 summary(bound = gsDesign::gsDesign(k = 3, test.type = 1, sfu = gsDesign::sfLDOF)$upper$bound) |>
 simtrial::as_gt()


  


Summary of simulation results by WLR tests
```

Weighted by FH(rho=0, gamma=0.5)

analysis

Time

N

Event

Crossing probability

1

12.14286

356

97

NA

2

24.71662

505

305

0.6666667

3

36.84514

505

405

1.0000000

\# Summarize simulations and compare with the planned design simulation
\|\> [summary](https://rdrr.io/r/base/summary.html)(design = design)
\|\> simtrial::as_gt()

[TABLE]
