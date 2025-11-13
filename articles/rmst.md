# Restricted mean survival time (RMST)

## Introduction

Restricted mean survival time (RMST) is defined as the area under the
survival curve up to a specific time point. It can be interpreted as the
average survival time during a defined time period ranging from time 0
to a specific follow-up time point, which is a straightforward and
clinically meaningful way to interpret the contrast in survival between
groups.

The RMST may provide valuable information for comparing two survival
curves when the proportional hazards assumption is not met, such as in
cases of *crossing* or *delayed separation* of survival curves.

## RMST vs. logrank

The log-rank test calculates the test statistics using the survival rate
at each time point, and then summarizes them to test the equality of the
survival curves as a whole for the entire follow-up period.

A comparison of the RMST between two survival curves provides an
estimate of the duration of time gained or lost associated with the
exposure.

Although RMST has an advantage over the hazard ratio, a previous study
showed that the difference in RMST often has operating characteristics
similar to the log-rank test under the proportional hazards assumption
(Royston and Parmar 2013). However, in the case of crossing survival
curves, the efficacy of an intervention may be demonstrated by showing a
difference in RMST between the two curves, although the log-rank test
may fail to detect a difference because of the occurrence of
non-proportional hazards.

## Estimation of RMST in a single arm at a single time point

Assume the event time as \\T\\, and the survival function as \\S(t) =
Pr(T\>t)\\. The restricted mean survival time at a pre-specified cutoff
time point \\\tau\\ is \\ \text{RMST}(\tau) = E\[\min (T, \tau)\] =
\int\_{0}^{\tau} S(u) d u. \\ Suppose there are \\D\\ events, with
distinct observed event times at \\t_1 \< t_2 \< \ldots \<t_D\\. For any
\\i = 1, \ldots, D\\, let \\Y_i\\ be the number of risk just prior to
\\t_i\\, and let \\d_i\\ be the number of subjects that fail at \\t_i\\.
The Kaplan-Meier (product-limit) estimate of the survival function at
\\t_i\\ is

\\ \hat{S}(t_i) = \prod\_{j=1}^{i} \left( 1-\frac{d\_{j}}{Y\_{j}}
\right) \\ Based on the definition and formula above,
\\\text{RMST}(\tau)\\ can be estimated by \\ \widehat{\text{RMST}}(\tau)
= \int\_{0}^{\tau} \hat{S}(t) d t = \sum\_{i=1}^{L\_{\tau}}
\hat{S}\left(t\_{i-1}\right)\left(t\_{i}-t\_{i-1}\right) +
\hat{S}\left(t\_{L\_{\tau}}\right)\left(\tau-t\_{L\_{\tau}}\right), \\
where \\L\_{\tau}\\ is the number of \\t_i\\ values that are less than
\\\tau\\.

The standard error of \\\widehat{\text{RMST}}(\tau)\\ can be estimated
by \\ \hat{\sigma} = \widehat{\text{Var}}(\widehat{\text{RMST}}(\tau)) =
\sqrt{\sum\_{i=1}^{L\_\tau} \frac{d\_{i}
A\_{i}^{2}}{Y\_{i}\left(Y\_{i}-d\_{i}\right)}} \\ where \\ A\_{i} =
\int\_{t_i}^{\tau} \hat{S}(t) d t = \sum\_{j=i}^{L\_\tau}
\hat{S}(t\_{j}) (t\_{j+1}-t\_{j}) +
\hat{S}(t\_{L\_\tau})(\tau-t\_{L\_\tau}) \\ and
\\m=\sum\_{j=1}^{L\_\tau} d\_{j}\\. The \\(1-\alpha)\\ confidence
interval for \\\text{RMST}\\ can be calculated as \\ \left\[
\widehat{\operatorname{RMST}}(\tau) - z\_{\alpha/2}\hat{\sigma}, \\\\
\widehat{\operatorname{RMST}}(\tau) + z\_{\alpha/2}\hat{\sigma} \right\]
\\ where \\\alpha\\ is a predefined significant level, and
\\z\_{\alpha/2}\\ is the upper \\1-\alpha/2\\ critical value for the
standard normal distribution.

``` r
# Simulate NPH data from the piecewise model
library(simtrial)
# Table display
library(gt)
```

``` r
data(ex1_delayed_effect)
data_single_arm <- ex1_delayed_effect[ex1_delayed_effect$trt == 1, ]
simtrial:::rmst_single_arm(
  time_var = data_single_arm$month,
  event_var = data_single_arm$evntd,
  tau = 10
) |> gt()
```

| cutoff_time | group        | rmst     | variance   | std       | lcl      | ucl      | event |
|-------------|--------------|----------|------------|-----------|----------|----------|-------|
| 10          | Single Group | 6.495175 | 0.05711322 | 0.2389837 | 6.026776 | 6.963575 | 127   |

### Estimation of RMST differences in 2 arms at a single time point

Let \\\text{RMST}\_{1}(\tau)\\ and \\\text{RMST}\_{2}(\tau)\\ be the
RMST of treatment group 1 and 2 at predefined time \\\tau\\, the RMST
difference between the 2 treatment groups (\\\theta\\) can be defined as
\\ \theta = \text{RMST}\_1(\tau) - \text{RMST}\_2(\tau). \\ The expected
value of \\\theta\\ is \\E(\theta) = E\[\text{RMST}\_{1}(\tau)\] -
E\[\text{RMST}\_{2}(\tau)\]\\. If two treatment groups are independent,
the variance of \\\theta\\ is \\ \text{Var}(\theta) = \sigma\_{1}^{2} +
\sigma\_{2}^{2} \\

Similarly, the \\(1-\alpha)\\ confidence interval for RMST difference
between 2 groups can be calculated as: \\ \left\[ \hat{\theta} -
z\_{\alpha/2}\sqrt{\hat{\sigma}\_1^2 + \hat{\sigma}\_2^2}, \\\\
\hat{\theta} + z\_{\alpha/2}\sqrt{\hat{\sigma}\_1^2 + \hat{\sigma}\_2^2}
\right\]. \\

``` r
tau <- 10

data(ex1_delayed_effect)

ex1_delayed_effect |>
  rmst(
    var_label_tte = "month",
    var_label_event = "evntd",
    var_label_group = "trt",
    tau = 10,
    reference = "0"
  )
#> $method
#> [1] "RMST"
#> 
#> $parameter
#> [1] 10
#> 
#> $estimate
#> [1] 0.8650493
#> 
#> $se
#> [1] 0.3900344
#> 
#> $z
#> [1] 2.21788
```

The R package survRM2 (Uno et al. 2022) performs two-sample comparisons
using RMST as a summary measure of the survival time distribution. Three
kinds of between-group contrast metrics (i.e., the difference in RMST,
the ratio of RMST and the ratio of the restricted mean time lost (RMTL))
are computed. It performs an ANCOVA-type covariate adjustment as well as
unadjusted analyses for those measures. We use this R package as
validation of
[`simtrial::rmst()`](https://merck.github.io/simtrial/reference/rmst.md).

``` r
verify <- survRM2::rmst2(
  time = ex1_delayed_effect$month,
  status = ex1_delayed_effect$evntd,
  arm = ex1_delayed_effect$trt,
  tau = tau,
  alpha = 0.05
)

verify$RMST.arm1$rmst[1] - verify$RMST.arm0$rmst[1]
#>      Est. 
#> 0.8650493
```

### References

Royston, Patrick, and Mahesh KB Parmar. 2013. “Restricted Mean Survival
Time: An Alternative to the Hazard Ratio for the Design and Analysis of
Randomized Trials with a Time-to-Event Outcome.” *BMC Medical Research
Methodology* 13 (1): 1–15.

Uno, Hajime, Lu Tian, Miki Horiguchi, Angel Cronin, Chakib Battioui, and
James Bell. 2022. “survRM2: Comparing Restricted Mean Survival Time.”
<https://CRAN.R-project.org/package=survRM2>.
