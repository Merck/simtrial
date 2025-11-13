# Custom Fixed Design Simulations: A Tutorial on Writing Code from the Ground Up

``` r
library(gsDesign2)
library(simtrial)
library(dplyr)
library(gt)
library(doFuture)
library(tibble)
set.seed(2025)
```

The vignette [Simulate Fixed Designs with Ease via
sim_fixed_n](https://merck.github.io/simtrial/articles/sim_fixed_design_simple.html)
presents fixed design simulations using a single function call,
[`sim_fixed_n()`](https://merck.github.io/simtrial/reference/sim_fixed_n.md).
It offers a simple and straightforward process for running simulations
quickly.

If users are interested in the following aspects, we recommend
simulating from scratch rather than directly using
[`sim_fixed_n()`](https://merck.github.io/simtrial/reference/sim_fixed_n.md):

- Tests beyond the logrank test or Fleming-Harrington weighted logrank
  tests, such as modestly weighted logrank tests, RMST tests, and
  milestone tests.
- More complex cutoffs, such as analyzing data after 12 months of
  follow-up when at least 80% of the patient population is enrolled.
- Different dropout rates in the control and experimental groups.

The process for simulating from scratch is outlined in Steps 1 to 5
below.

## Step 1: Simulate time-to-event data

The
[`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md)
function allows for the simulation of a clinical trial with essentially
arbitrary patterns of enrollment, failure rates, and censoring. To
implement
[`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md),
you need to specify 5 design characteristics to simulate time-to-event
data:

1.  Sample Size (input as `n`).
2.  Stratified or Non-Stratified Designs (input as `stratum`).
3.  Randomization Ratio (input as `block`). The
    [`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md)
    function uses fixed block randomization.
4.  Enrollment Rate (input as `enroll_rate`). The
    [`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md)
    function supports piecewise enrollment, allowing the enrollment rate
    to be piecewise constant.
5.  Failure Rate (input as `fail_rate`) or time-to-event rate. The
    [`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md)
    function uses a piecewise exponential distribution for the failure
    rate, which makes it easy to define a distribution with changing
    failure rates over time. Specifically, in the \\j\\-th interval, the
    rate is denoted as \\\lambda_j \geq 0\\. We require that at least
    one interval has \\\lambda_j \> 0\\. There are two methods for
    defining the failure rate:

- Specify the failure rate by treatment group, stratum, and time period.
  An example can be found in Scenario b).
- Create a `fail_rate` using
  [`gsDesign2::define_fail_rate`](https://merck.github.io/gsDesign2/reference/define_fail_rate.html),
  and then convert it to the required format using
  [`to_sim_pw_surv()`](https://merck.github.io/simtrial/reference/to_sim_pw_surv.md).
  An example is provided in Scenario a).

1.  Dropout Rate (input as `dropout_Rate`). The
    [`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md)
    function accepts piecewise constant dropout rates, which may vary by
    treatment group. The configuration for dropout should be specified
    by treatment group, stratum, and time period, and setting up the
    dropout rate follows the same approach as the failure rate.

### Scenario a) The simplest scenario

We begin with the simplest implementation of
[`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md).
The following lines of code will generate 500 subjects using equal
randomization and an unstratified design.

``` r
n_sim <- 100
n <- 500
stratum <- data.frame(stratum = "All", p = 1)
block <- rep(c("experimental", "control"), 2)

enroll_rate <- define_enroll_rate(rate = 12, duration = n / 12)

fail_rate <- define_fail_rate(duration = c(6, Inf), fail_rate = log(2) / 10, 
                              hr = c(1, 0.7), dropout_rate = 0.0001)

uncut_data_a <- sim_pw_surv(n = n, stratum = stratum, block = block,
                            enroll_rate = enroll_rate,
                            fail_rate = to_sim_pw_surv(fail_rate)$fail_rate,
                            dropout_rate = to_sim_pw_surv(fail_rate)$dropout_rate)
```

The output of
[`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md)
is subject-level observations, including stratum, enrollment time for
the observation, treatment group the observation is randomized to,
failure time, dropout time, calendar time of enrollment plot the minimum
of failure time and dropout time ( `cte`), and an failure and dropout
indicator (`fail = 1` is a failure, `fail = 0` is a dropout).

``` r
uncut_data_a |> head() |> gt() |> tab_header("An Overview of Simulated TTE data")
```

| An Overview of Simulated TTE data |             |              |            |              |            |      |
|-----------------------------------|-------------|--------------|------------|--------------|------------|------|
| stratum                           | enroll_time | treatment    | fail_time  | dropout_time | cte        | fail |
| All                               | 0.04250966  | experimental | 1.1694497  | 6265.622     | 1.2119594  | 1    |
| All                               | 0.15597253  | control      | 73.5774306 | 21706.085    | 73.7334032 | 1    |
| All                               | 0.19363998  | experimental | 15.1356774 | 18096.246    | 15.3293174 | 1    |
| All                               | 0.23882703  | control      | 1.9020241  | 20844.870    | 2.1408512  | 1    |
| All                               | 0.27283752  | control      | 2.9475061  | 11058.946    | 3.2203437  | 1    |
| All                               | 0.29362231  | experimental | 0.5490203  | 7004.398     | 0.8426426  | 1    |

### Scenario b) Differential dropout rates

The dropout rate can differ between groups. For instance, in open-label
studies, the control group may experience a higher dropout rate. The
follow lines of code assumes the control group has a dropout rate of
0.002 for the first 10 months, which then decreases to 0.001 thereafter.
In contrast, the experimental group has a constant dropout rate of 0.001
throughout the study.

``` r
differential_dropout_rate <- data.frame(
  stratum = rep("All", 3), 
  period = c(1, 2, 1), 
  treatment = c("control", "control", "experimental"), 
  duration = c(10, Inf, Inf), 
  rate = c(0.002, 0.001, 0.001))

uncut_data_b <- sim_pw_surv(n = n, stratum = stratum, block = block,
                            enroll_rate = enroll_rate,
                            fail_rate = to_sim_pw_surv(fail_rate)$fail_rate,
                            dropout_rate = differential_dropout_rate)
```

### Scenario c) Stratified designs

The following code assumes there are two strata (biomarker-positive and
biomarker-negative) with equal prevalence of 0.5 for each. In the
control arm, the median survival time is 10 months for
biomarker-positive subjects and 8 months for biomarker-negative
subjects. For both strata, the hazard ratio is 1 for the first 3 months,
after which it decreases to 0.6 for biomarker-positive subjects and 0.8
for biomarker-negative subjects. The dropout rate is contently 0.001 for
both strata over time.

``` r
stratified_enroll_rate <- data.frame(
  stratum = c("Biomarker positive", "Biomarker negative"),
  rate = c(12, 12), 
  duration = c(1, 1))

stratified_fail_rate <- data.frame(
  stratum = c(rep("Biomarker positive", 3), rep("Biomarker negative", 3)), 
  period = c(1, 1, 2, 1, 1, 2), 
  treatment = rep(c("control", "experimental", "experimental"), 2),
  duration = c(Inf, 3, Inf, Inf, 3, Inf), 
  rate = c(# failure rate of biomarker positive subjects: control arm, exp arm period 1, exp arm period 2
           log(2) / 10, log(2) /10, log(2) / 10 * 0.6,
           # failure rate of biomarker negative subjects: control arm, exp arm period 1, exp arm period 2
           log(2) / 8, log(2) /8, log(2) / 8 * 0.8)
  )

stratified_dropout_rate <- data.frame(
  stratum = rep(c("Biomarker positive", "Biomarker negative"), each = 2),
  period = c(1, 1, 1, 1), 
  treatment = c("control", "experimental", "control", "experimental"),
  duration = rep(Inf, 4), 
  rate = rep(0.001, 4)
  )

uncut_data_c <- sim_pw_surv(n = n,
                            stratum = data.frame(stratum = c("Biomarker positive", "Biomarker negative"), 
                                                 p = c(0.5, 0.5)),
                            block = block,
                            enroll_rate = stratified_enroll_rate,
                            fail_rate = stratified_fail_rate,
                            dropout_rate = stratified_dropout_rate
                            )
```

### Scenario d) Multi-arm designs

Suppose you wish to have 3 arms: control, low-dose and high-dose. The
following code assumes the control arm has a median survival time of 10
months, the low-dose arm has a median survival time of 12 months, and
the high-dose arm has a median survival time of 14 months. The hazard
ratio for the low-dose arm is 0.8, and the hazard ratio for the
high-dose arm is 0.6. The dropout rate is 0.001 for all arms. Block size
is 7 with 3:2:2 randomization.

We begin by setting up enrollment, failure and dropout rates.

``` r
enroll_rate <- define_enroll_rate(rate = 12, duration = n / 12)

three_arm_fail_rate <- data.frame(
  stratum = "All",
  period = c(1, 1, 2, 1, 2), 
  treatment = c("control", "low-dose", "low-dose", "high-dose", "high-dose"),
  duration = c(Inf, 3, Inf, 3, Inf), 
  rate = c(# failure rate of control arm: period 1, period 2
           log(2) / 10,
           # failure rate of low-dose arm: period 1, period 2
           log(2) / c(10, 10 / .8),
           # failure rate of high-dose arm: period 1, period 2
           log(2) / c(10, 10 / .6)))

three_arm_dropout_rate <- data.frame(
  stratum = "All",
  period = c(1, 1, 1), 
  treatment = c("control", "low-dose", "high-dose"),
  duration = rep(Inf, 3), 
  rate = rep(0.001, 3))

uncut_data_d <- sim_pw_surv(n = n,
                            stratum = data.frame(stratum = "All"),
                            block = c(rep("control", 3), rep("low-dose", 2), rep("high-dose", 2)),
                            enroll_rate = enroll_rate,
                            fail_rate = three_arm_fail_rate,
                            dropout_rate = three_arm_dropout_rate)
```

For illustration purposes, we will focus on scenario b) for the
following discussion.

``` r
uncut_data <- uncut_data_b
```

## Step 2: Cut data

The
[`get_analysis_date()`](https://merck.github.io/simtrial/reference/get_analysis_date.md)
derives analysis date for interim/final analysis given multiple
conditions, see [the help page of `get_analysis_date()` at the pkgdown
website](https://merck.github.io/simtrial/reference/get_analysis_date.html).

Users can cut for analysis at the 24th month and there are 300 events,
whichever arrives later. This is equivalent to `timing_type = 4` in
[`sim_fixed_n()`](https://merck.github.io/simtrial/reference/sim_fixed_n.md).

``` r
cut_date_a <- get_analysis_date(data = uncut_data,
                                planned_calendar_time = 24,
                                target_event_overall = 300)
```

Users can also cut by the maximum of targeted 300 event and minimum
follow-up 12 months. This is equivalent to `timing_type = 5` in
[`sim_fixed_n()`](https://merck.github.io/simtrial/reference/sim_fixed_n.md).

``` r
cut_date_b <- get_analysis_date(data = uncut_data,
                                min_followup = 12,
                                target_event_overall = 300)
```

Users can cut data when there are 300 events, with maximum time
extension to reach targeted events of 24 months. This is not enabled in
`timing_type` of
[`sim_fixed_n()`](https://merck.github.io/simtrial/reference/sim_fixed_n.md).

``` r
cut_date_c <- get_analysis_date(data = uncut_data,
                                max_extension_for_target_event = 12,
                                target_event_overall = 300)
```

Users can cut data after 12 months followup when 80% of the patients are
enrolled in the overall population as below. This is not enabled in
`timing_type` of
[`sim_fixed_n()`](https://merck.github.io/simtrial/reference/sim_fixed_n.md).

``` r
cut_date_d <- get_analysis_date(data = uncut_data,
                                min_n_overall = 100 * 0.8,
                                min_followup = 12)
```

More examples are available in the [reference page of
`get_analysis_date()`](https://merck.github.io/simtrial/reference/get_analysis_date.html)
For illustration purposes, we will focus on scenario d) for the
following discussion.

``` r
cut_date <- cut_date_d
cat("The cutoff date is ", round(cut_date, 2))
```

    ## The cutoff date is  20.19

``` r
cut_data <- uncut_data |> cut_data_by_date(cut_date)
cut_data |> head() |> gt() |> tab_header(paste0("An Overview of TTE data Cut at ", round(cut_date, 2), "Months"))
```

| An Overview of TTE data Cut at 20.19Months |       |         |              |
|--------------------------------------------|-------|---------|--------------|
| tte                                        | event | stratum | treatment    |
| 19.288154                                  | 1     | All     | control      |
| 20.029362                                  | 0     | All     | experimental |
| 5.400825                                   | 1     | All     | control      |
| 8.723331                                   | 1     | All     | experimental |
| 19.860289                                  | 0     | All     | control      |
| 4.589687                                   | 1     | All     | experimental |

## Step 3: Run tests

The simtrial package provides many options for testing methods,
including (weighted) logrank tests, RMST test, milestone test, and
MaxComboi test, see the \[Section “Compute p-values/test statistics” at
the pkgdown reference page\]
(<https://merck.github.io/simtrial/reference/index.html#compute-p-values-test-statistics>).

The following code lists all possible tests available in simtrial. Users
can select one of the tests listed above or combine the testing results
to make comparisons across tests. For demonstration purposes, we will
aggregate all tests together.

``` r
# Logrank test
sim_res_lr <- cut_data |> wlr(weight = fh(rho = 0, gamma = 0))

# weighted logrank test by Fleming-Harrington weights
sim_res_fh <- cut_data |> wlr(weight = fh(rho = 0, gamma = 0.5))

# Modestly weighted logrank test
sim_res_mb <- cut_data |> wlr(weight = mb(delay = Inf, w_max = 2))

# Weighted logrank test by Xu 2017's early zero weights
sim_res_xu <- cut_data |> wlr(weight = early_zero(early_period = 3))

# RMST test
sim_res_rmst <- cut_data |> rmst(tau = 10)

# Milestone test
sim_res_ms <- cut_data |> milestone(ms_time = 10)

# Maxcombo tests comboing multiple weighted logrank test with Fleming-Harrington weights
sim_res_mc <- cut_data |> maxcombo(rho = c(0, 0), gamma = c(0, 0.5))
```

The output of the tests mentioned above are lists including:

- The testing method employed (WLR, RMST, milestone, or MaxCombo), which
  can be accessed using `sim_res_rmst$method`.
- The parameters associated with the testing method. For instance, the
  RMST test parameter is 10, indicating that the RMST is evaluated at
  month 10. You can find this information using
  `sim_res_rmst$parameter`.
- The point estimate and standard error for the testing method used. For
  example, the point estimate for RMST represents the survival
  difference between the experimental group and the control group. This
  estimate can be retrieved with `sim_res_rmst$estimate` and
  `sim_res_rmst$se`.
- The Z-score for the testing method, accessible via `sim_res_rmst$z`.
  Please note that the Z-score is not provided for the MaxCombo test;
  instead, the p-value is reported (`sim_res_mc$p_value`).

``` r
sim_res <- tribble(
  ~Method, ~Parameter, ~Z, ~Estimate, ~SE, ~`P value`,
  sim_res_lr$method, sim_res_lr$parameter, sim_res_lr$z, sim_res_lr$estimate, sim_res_lr$se, pnorm(-sim_res_lr$z),
  sim_res_fh$method, sim_res_fh$parameter, sim_res_fh$z, sim_res_fh$estimate, sim_res_fh$se, pnorm(-sim_res_fh$z),
  sim_res_mb$method, sim_res_mb$parameter, sim_res_mb$z, sim_res_mb$estimate, sim_res_mb$se, pnorm(-sim_res_mb$z),
  sim_res_xu$method, sim_res_xu$parameter, sim_res_xu$z, sim_res_xu$estimate, sim_res_xu$se, pnorm(-sim_res_xu$z),
  sim_res_rmst$method, sim_res_rmst$parameter|> as.character(), sim_res_rmst$z, sim_res_rmst$estimate, sim_res_rmst$se, pnorm(-sim_res_rmst$z),
  sim_res_ms$method, sim_res_ms$parameter |> as.character(), sim_res_ms$z, sim_res_ms$estimate, sim_res_ms$se, pnorm(-sim_res_ms$z),
  sim_res_mc$method, sim_res_mc$parameter, NA, NA, NA, sim_res_mc$p_value
  ) 

sim_res |> gt() |> tab_header("One Simulation Results")
```

| One Simulation Results |                                          |          |             |           |             |
|------------------------|------------------------------------------|----------|-------------|-----------|-------------|
| Method                 | Parameter                                | Z        | Estimate    | SE        | P value     |
| WLR                    | FH(rho=0, gamma=0)                       | 1.788225 | -9.3460253  | 5.2264273 | 0.036869897 |
| WLR                    | FH(rho=0, gamma=0.5)                     | 2.360073 | -6.8225888  | 2.8908376 | 0.009135661 |
| WLR                    | MB(delay = Inf, max_weight = 2)          | 2.130368 | -16.9085058 | 7.9368946 | 0.016570624 |
| WLR                    | Xu 2017 with first 3 months of 0 weights | 2.600908 | -11.0553977 | 4.2505921 | 0.004648873 |
| RMST                   | 10                                       | 1.128125 | 0.5589405   | 0.4954596 | 0.129633491 |
| milestone              | 10                                       | 1.940708 | 0.4501018   | 0.2319266 | 0.026146861 |
| MaxCombo               | FH(0, 0) + FH(0, 0.5)                    | NA       | NA          | NA        | 0.012855297 |

## Step 4: Perform the above single simulation repeatedly

We will now merge Steps 1 to 3 into a single function named `one_sim()`,
which facilitates a single simulation run. The construction of
`one_sim()` involves copying all the lines of code from Steps 1 to 3.

``` r
one_sim <- function(sim_id = 1, 
                    # arguments from Step 1: design characteristic
                    n, stratum, enroll_rate, fail_rate, dropout_rate, block, 
                    # arguments from Step 2； cutting method
                    min_n_overall, min_followup,
                    # arguments from Step 3； testing method
                    fh, mb, xu, rmst, ms, mc
                    ) {
    # Step 1: simulate a time-to-event data
    uncut_data <- sim_pw_surv(
      n = n,
      stratum = stratum,
      block = block,
      enroll_rate = enroll_rate,
      fail_rate = fail_rate,
      dropout_rate = dropout_rate) 
    
    ## Step 2: Cut data
    cut_date <- get_analysis_date(min_n_overall = min_n_overall, min_followup = min_followup, data = uncut_data)
    cut_data <- uncut_data |> cut_data_by_date(cut_date)
    
    # Step 3: Run tests
    sim_res_lr <- cut_data |> wlr(weight = fh(rho = 0, gamma = 0))
    sim_res_fh <- cut_data |> wlr(weight = fh(rho = fh$rho, gamma = fh$gamma))
    sim_res_mb <- cut_data |> wlr(weight = mb(delay = mb$delay, w_max = mb$w_max))
    sim_res_xu <- cut_data |> wlr(weight = early_zero(early_period = xu$early_period))
    sim_res_rmst <- cut_data |> rmst(tau = rmst$tau)
    sim_res_ms <- cut_data |> milestone(ms_time = ms$ms_time)
    sim_res_mc <- cut_data |> maxcombo(rho = mc$rho, gamma = mc$gamma)
    
    sim_res <- tribble(
      ~`Sim ID`, ~Method, ~Parameter, ~Z, ~Estimate, ~SE, ~`P value`,
      sim_id, sim_res_lr$method, sim_res_lr$parameter, sim_res_lr$z, sim_res_lr$estimate, sim_res_lr$se, pnorm(-sim_res_lr$z),
      sim_id, sim_res_fh$method, sim_res_fh$parameter, sim_res_fh$z, sim_res_fh$estimate, sim_res_fh$se, pnorm(-sim_res_fh$z),
      sim_id, sim_res_mb$method, sim_res_mb$parameter, sim_res_mb$z, sim_res_mb$estimate, sim_res_mb$se, pnorm(-sim_res_mb$z),
      sim_id, sim_res_xu$method, sim_res_xu$parameter, sim_res_xu$z, sim_res_xu$estimate, sim_res_xu$se, pnorm(-sim_res_xu$z),
      sim_id, sim_res_rmst$method, sim_res_rmst$parameter|> as.character(), sim_res_rmst$z, sim_res_rmst$estimate, sim_res_rmst$se, pnorm(-sim_res_rmst$z),
      sim_id, sim_res_ms$method, sim_res_ms$parameter |> as.character(), sim_res_ms$z, sim_res_ms$estimate, sim_res_ms$se, pnorm(-sim_res_ms$z),
      sim_id, sim_res_mc$method, sim_res_mc$parameter, NA, NA, NA, sim_res_mc$p_value
  ) 
      
    return(sim_res)
}
```

After that, we will execute `one_sim()` multiple times using parallel
computation. The following lines of code uses 2 workers to run 100
simulations.

``` r
set.seed(2025)

plan("multisession", workers = 2)
ans <- foreach(
  sim_id = seq_len(n_sim),
  .errorhandling = "stop",
  .options.future = list(seed = TRUE)
  ) %dofuture% {
    ans_new <- one_sim(
      sim_id = sim_id, 
      # arguments from Step 1: design characteristic
      n = n, 
      stratum = stratum, 
      enroll_rate = enroll_rate, 
      fail_rate = to_sim_pw_surv(fail_rate)$fail_rate, 
      dropout_rate = differential_dropout_rate, 
      block = block, 
      # arguments from Step 2； cutting method
      min_n_overall = 500 * 0.8,
      min_followup = 12,
      # arguments from Step 3； testing method
      fh = list(rho = 0, gamma = 0.5), 
      mb = list(delay = Inf, w_max = 2), 
      xu = list(early_period = 3), 
      rmst = list(tau = 10), 
      ms = list(ms_time = 10), 
      mc = list(rho = c(0, 0), gamma = c(0, 0.5))
      )
                              
    ans_new
  }

ans <- data.table::rbindlist(ans)

plan("sequential")
```

The output from the parallel computation resembles the output of
`sim_fix_n()` described in the vignette [Simulate Fixed Designs with
Ease via
sim_fixed_n](https://merck.github.io/simtrial/articles/sim_fixed_design_simple.html).
Each row in the output corresponds to the simulation results for each
testing method per each repeation.

``` r
ans |> head() |> gt() |> tab_header("Overview Each Simulation results")
```

| Overview Each Simulation results |           |                                          |           |             |            |            |
|----------------------------------|-----------|------------------------------------------|-----------|-------------|------------|------------|
| Sim ID                           | Method    | Parameter                                | Z         | Estimate    | SE         | P value    |
| 1                                | WLR       | FH(rho=0, gamma=0)                       | 1.9005840 | -18.1176359 | 9.5326679  | 0.02867826 |
| 1                                | WLR       | FH(rho=0, gamma=0.5)                     | 2.1478743 | -12.7153128 | 5.9199519  | 0.01586187 |
| 1                                | WLR       | MB(delay = Inf, max_weight = 2)          | 2.0945611 | -32.4645715 | 15.4994626 | 0.01810501 |
| 1                                | WLR       | Xu 2017 with first 3 months of 0 weights | 1.9926533 | -16.4165941 | 8.2385601  | 0.02314971 |
| 1                                | RMST      | 10                                       | 0.7072865 | 0.2184929   | 0.3089170  | 0.23969421 |
| 1                                | milestone | 10                                       | 0.9610773 | 0.1280086   | 0.1331928  | 0.16825665 |

## Step 5: Summarize simulations

Using the 100 parallel simulations provided above, users can summarize
the simulated power and compare it across different testing methods with
some data manipulation using `dplyr`. Please note that the power
calculation for the MaxCombo test differs from the other tests, as it
does not report a Z-score.

``` r
ans_non_mc <- ans |>
  filter(Method != "MaxCombo") |>
  group_by(Method, Parameter) %>% 
  summarise(Power = mean(Z > -qnorm(0.025))) |>
  ungroup()

ans_mc <- ans |>
  filter(Method == "MaxCombo") |>
  summarize(Power = mean(`P value` < 0.025), Method = "MaxCombo", Parameter = "FH(0, 0) + FH(0, 0.5)") 

ans_non_mc |>
  union(ans_mc) |>
  gt() |>
  tab_header("Summary from 100 simulations")
```

| Summary from 100 simulations |                                          |       |
|------------------------------|------------------------------------------|-------|
| Method                       | Parameter                                | Power |
| RMST                         | 10                                       | 0.06  |
| WLR                          | FH(rho=0, gamma=0)                       | 0.39  |
| WLR                          | FH(rho=0, gamma=0.5)                     | 0.53  |
| WLR                          | MB(delay = Inf, max_weight = 2)          | 0.54  |
| WLR                          | Xu 2017 with first 3 months of 0 weights | 0.48  |
| milestone                    | 10                                       | 0.13  |
| MaxCombo                     | FH(0, 0) + FH(0, 0.5)                    | 0.52  |

### References
