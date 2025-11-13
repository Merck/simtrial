# Simulate Group Sequential Designs with Ease via sim_gs_n

``` r
library(gsDesign2)
library(simtrial)
library(dplyr)
library(gt)

set.seed(2025)
```

The
[`sim_gs_n()`](https://merck.github.io/simtrial/reference/sim_gs_n.md)
function simulates group sequential designs with fixed sample size and
multiple analyses. There are many advantages of calling
[`sim_gs_n()`](https://merck.github.io/simtrial/reference/sim_gs_n.md)
directly.

- It is simple, which allows for a single function call to perform
  thousands of simulations.
- It allows for a variety of testing methods, such as logrank, weighted
  logrank test, RMST, and milestone test, via the `test = ...` argument.
- It automatically implements a parallel computation backend, allowing
  users to achieve shorter running times.
- It enables flexible data cutting methods for analyses,
  specifically: (1) planned calendar time, (2) targeted events, (3)
  maximum time extension to reach targeted events, (4) planned minimum
  time after the previous analysis, (5) minimal follow-up time after
  specified enrollment fraction; various combinations of these cut
  methods can also be specified.

The process for simulating via
[`sim_gs_n()`](https://merck.github.io/simtrial/reference/sim_gs_n.md)
is outlined in Steps 1 to 3 below.

## Step 1: Define design paramaters

To run simulations for a group sequential design, several design
characteristics are required. The following code creates a design for an
unstratified 2-arm trial with equal randomization. Enrollment is
targeted to last for 12 months at a constant enrollment rate. The
control arm is specified as exponential with a median of 10 months. The
experimental arm distribution has a piecewise exponential distribution
with a delay of 3 months with no benefit relative to control (HR = 1)
followed by a hazard ratio of 0.6 thereafter. Additionally, there is an
exponential dropout rate of 0.001 per month (or unit of time). The set
up of these parameters is similar to the vignette [Simulate Fixed
Designs with Ease via
sim_fixed_n](https://merck.github.io/simtrial/articles/sim_fixed_design_simple.html).
The total sample size is derived for 90% power.

``` r
stratum <- data.frame(stratum = "All", p = 1)
block <- rep(c("experimental", "control"), 2)
# enrollment rate will be updated later, 
# multiplied by a constant to get targeted power
enroll_rate <- data.frame(stratum = "All", rate = 1, duration = 12)
fail_rate <- data.frame(stratum = "All",
                        duration = c(3, Inf), fail_rate = log(2) / 10, 
                        hr = c(1, 0.6), dropout_rate = 0.001)
# Derive design using the average hazard ratio method
x <- gs_design_ahr(enroll_rate = enroll_rate, fail_rate = fail_rate,
                   analysis_time = c(12, 24, 36), alpha = 0.025, beta = 0.1,
                   # spending function for upper bound
                   upper = gs_spending_bound, 
                   upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025),
                   # Fixed lower bound
                   lower = gs_b,
                   lpar = rep(-Inf, 3)) |> to_integer()

sample_size <- x$analysis$n |> max()
event <- x$analysis$event
eff_bound <- x$bound$z[x$bound$bound == "upper"]
```

``` r
cat(paste("The total sample size is ", sample_size, "\n", sep = ''))
```

    ## The total sample size is 362

``` r
cat("The number of events at IA1, IA2 and FA are:", event, "\n")
```

    ## The number of events at IA1, IA2 and FA are: 106 227 287

``` r
cat("The efficacy bounds at IA1, IA2 and FA are:", round(eff_bound, 3), "\n")
```

    ## The efficacy bounds at IA1, IA2 and FA are: 3.508 2.269 2.023

``` r
cat("Targeted analysis times:", round(x$analysis$time, 1), "\n")
```

    ## Targeted analysis times: 12 24 35.9

Now we get the updated planned enrollment rate from the design to
achieve the above targeted sample size.

``` r
enroll_rate <- x$enroll_rate
enroll_rate
```

    ##   stratum     rate duration
    ## 1     All 30.16667       12

There are additional parameters required for a group sequential design
simulation as demonstrated below. One is the testing method. We focus on
logrank here. Users can change these to other tests of interest; in
comments below we demonstrate logrank and modestly weighted logrank test
(Magirr and Burman (2019)) and Fleming-Harrington (Harrington and
Fleming (1982)) tests; how to set up group sequential designs for these
alternate tests is beyond the scope of this article. More testing
methods are available at [reference page of
`simtrial`](https://merck.github.io/simtrial/reference/index.html#compute-p-values-test-statistics).

``` r
# Example for logrank
weight <- fh(rho = 0, gamma = 0)
test <- wlr
# Example for Modestly Weighted Logrank Test (Magirr-Burman)
# weight <- mb(delay = Inf, w_max = 2)
# Example for Fleming-Harrington(0, 0.5)
# weight <- fh(rho = 0, gamma = 0.5)
```

The final step is to specify how when data is cut off for each analysis
in the group sequential design. The
[`create_cut()`](https://merck.github.io/simtrial/reference/create_cut.md)
function includes 5 options for the analysis cutoff:

1.  planned calendar time,
2.  targeted events,
3.  maximum time extension to reach targeted events,
4.  planned minimum time after the previous analysis, and
5.  minimal follow-up time after specified enrollment fraction.

More details and examples are available at the [help
page](https://merck.github.io/simtrial/reference/get_analysis_date.html).

A straightforward method for cutting analyses is based on events. For
instance, the following code specifies 2 interim analyses and 1 final
analysis cut when 106, 227, 287 events occur. **In this event-driven
approach, there is no need to update the efficacy boundary.**

``` r
ia1_cut <- create_cut(target_event_overall = event[1])
ia2_cut <- create_cut(target_event_overall = event[2])
fa_cut <- create_cut(target_event_overall = event[3])

cut <- list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut)
```

Users can set more complex cutting. For example,

- The first interim analysis occurs at the targeted IA 1 analyusis time
  or when at least the IA 1 targeted events are observed, which is
  later. However, if the target number of events is not reached, we will
  wait a maximum of 16 months from start of enrollment.
- The second interim analysis is targeted to take place at the targeted
  time or targeted events for IA 2, whichever is later. Additionally,
  this interim analysis will be scheduled at least 10 months after the
  first interim analysis, but no later than 28 months from start of
  enrollment.
- The final analysis is set for the targeted final analysis time or
  targeted events at the final analysis, whichever is later. It must be
  at least 6 months after IA 2.

**Please keep in mind that if the cut is not event-driven, boundary
updates are necessary; this is not covered here.** You can find more
information about [the boundary updates in the boundary update
vignette](https://merck.github.io/simtrial/articles/sim_fixed_design_simple.html).
In this vignette, we use the event-driven cut from above for
illustrative purposes.

``` r
ia1_cut <- create_cut(
  planned_calendar_time = round(x$analysis$time[1]), 
  target_event_overall = x$analysis$event[1],
  max_extension_for_target_event = 16)

ia2_cut <- create_cut(
  planned_calendar_time = round(x$analysis$time[2]),
  target_event_overall = x$analysis$event[2],
  min_time_after_previous_analysis = 10, 
  max_extension_for_target_event = 28)

fa_cut <- create_cut(
  planned_calendar_time = round(x$analysis$time[3]),
  min_time_after_previous_analysis = 6,
  target_event_overall = x$analysis$event[3])

cut <- list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut)
```

## Step 2: Run `sim_gs_n()`

Now that we have set up the design characteristics in Step 1, we can
proceed to run
[`sim_gs_n()`](https://merck.github.io/simtrial/reference/sim_gs_n.md)
for a specified number of simulations. This function automatically
utilizes a parallel computing backend, which helps reduce the running
time.

``` r
n_sim <- 100 # Number of simulated trials
sim_res <- sim_gs_n(
  n_sim = n_sim,
  sample_size = sample_size, stratum = stratum, block = block,
  enroll_rate = enroll_rate, fail_rate = fail_rate,
  test = test, weight = weight, cut = cut)
```

The output of `sim_gs_n` is a data frame with one row per simulation per
analysis. We show results for the first 2 simulated trials here. The
estimate column is the *sum(0 - E)* for the logrank test; the se column
is its standard error estimated under the null hypothesis. The `z`
column is the test statistic for the logrank test (estimate / se). The
`info` and `info0` columns are the information at the current analysis
under the alternate and null hypotheses, respectively.

``` r
sim_res |> head(n = 6) |> gt() |> tab_header("Overview Each Simulation results") |>
  fmt_number(columns = c(5, 8:12), decimals = 2)
```

| Overview Each Simulation results |        |                    |          |          |     |       |          |      |      |       |       |
|----------------------------------|--------|--------------------|----------|----------|-----|-------|----------|------|------|-------|-------|
| sim_id                           | method | parameter          | analysis | cut_date | n   | event | estimate | se   | z    | info  | info0 |
| 1                                | WLR    | FH(rho=0, gamma=0) | 1        | 12.53    | 362 | 106   | −4.00    | 5.13 | 0.78 | 26.46 | 26.50 |
| 1                                | WLR    | FH(rho=0, gamma=0) | 2        | 24.68    | 362 | 227   | −23.11   | 7.46 | 3.10 | 55.82 | 56.75 |
| 1                                | WLR    | FH(rho=0, gamma=0) | 3        | 35.99    | 362 | 287   | −36.82   | 8.27 | 4.45 | 70.68 | 71.75 |
| 2                                | WLR    | FH(rho=0, gamma=0) | 1        | 12.19    | 348 | 106   | −5.20    | 5.15 | 1.01 | 26.26 | 26.50 |
| 2                                | WLR    | FH(rho=0, gamma=0) | 2        | 24.66    | 362 | 227   | −8.33    | 7.53 | 1.11 | 56.50 | 56.75 |
| 2                                | WLR    | FH(rho=0, gamma=0) | 3        | 37.57    | 362 | 287   | −18.96   | 8.44 | 2.24 | 71.11 | 71.75 |

## Step 3: Summarize simulations

With the 100 simulations provided, users can summarize the simulated
power and compare it to the target power of 90% as follows:

``` r
sim_res |>
  left_join(data.frame(analysis = 1:3, eff_bound = eff_bound)) |>
  group_by(analysis) |>
  summarize(`Mean time` = mean(cut_date), `sd(time)` = sd(cut_date), `Simulated power` = mean(z >= eff_bound)) |>
  ungroup() |>
  mutate(`Asymptotic power` = x$bound$probability[x$bound$bound == "upper"]) |>
  gt() |>
  tab_header("Summary of 100 simulations") |> 
  fmt_number(columns = 2, decimals = 1) |>
  fmt_number(columns = 3:5, decimals = 2)
```

| Summary of 100 simulations |           |          |                 |                  |
|----------------------------|-----------|----------|-----------------|------------------|
| analysis                   | Mean time | sd(time) | Simulated power | Asymptotic power |
| 1                          | 12.0      | 0.79     | 0.03            | 0.01             |
| 2                          | 23.8      | 1.53     | 0.66            | 0.66             |
| 3                          | 35.7      | 2.04     | 0.87            | 0.90             |

### References

Harrington, David P, and Thomas R Fleming. 1982. “A Class of Rank Test
Procedures for Censored Survival Data.” *Biometrika* 69 (3): 553–66.

Magirr, Dominic, and Carl-Fredrik Burman. 2019. “Modestly Weighted
Logrank Tests.” *Statistics in Medicine* 38 (20): 3782–90.
