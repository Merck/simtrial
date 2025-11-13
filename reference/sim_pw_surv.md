# Simulate a stratified time-to-event outcome randomized trial

`sim_pw_surv()` enables simulation of a clinical trial with essentially
arbitrary patterns of enrollment, failure rates and censoring. The
piecewise exponential distribution allows a simple method to specify a
distribution and enrollment pattern where the enrollment, failure, and
dropout rate changes over time. While the main purpose may be to
generate a trial that can be analyzed at a single point in time or using
group sequential methods, the routine can also be used to simulate an
adaptive trial design. Enrollment, failure, and dropout rates are
specified by treatment group, stratum and time period. Fixed block
randomization is used; blocks must include treatments provided in
failure and dropout specification. Default arguments are set up to allow
very simple implementation of a non-proportional hazards assumption for
an unstratified design.

## Usage

``` r
sim_pw_surv(
  n = 100,
  stratum = data.frame(stratum = "All", p = 1),
  block = c(rep("control", 2), rep("experimental", 2)),
  enroll_rate = data.frame(rate = 9, duration = 1),
  fail_rate = data.frame(stratum = rep("All", 4), period = rep(1:2, 2), treatment =
    c(rep("control", 2), rep("experimental", 2)), duration = rep(c(3, 1), 2), rate =
    log(2)/c(9, 9, 9, 18)),
  dropout_rate = data.frame(stratum = rep("All", 2), period = rep(1, 2), treatment =
    c("control", "experimental"), duration = rep(100, 2), rate = rep(0.001, 2))
)
```

## Arguments

- n:

  Number of observations. If length(n) \> 1, the length is taken to be
  the number required.

- stratum:

  A data frame with stratum specified in `stratum`, probability
  (incidence) of each stratum in `p`.

- block:

  Vector of treatments to be included in each block. Also used to
  calculate the attribute "ratio" (for more details see the section
  Value below).

- enroll_rate:

  Enrollment rates; see details and examples.

- fail_rate:

  Failure rates; see details and examples; note that treatments need to
  be the same as input in block.

- dropout_rate:

  Dropout rates; see details and examples; note that treatments need to
  be the same as input in block.

## Value

A data frame with the following variables for each observation:

- `stratum`: Stratum for the observation.

- `enroll_time`: Enrollment time for the observation.

- `treatment`: Treatment group; this will be one of the values in the
  input `block`.

- `fail_time`: Failure time generated using
  [`rpwexp()`](https://merck.github.io/simtrial/reference/rpwexp.md).

- `dropout_time`: Dropout time generated using
  [`rpwexp()`](https://merck.github.io/simtrial/reference/rpwexp.md).

- `cte`: Calendar time of enrollment plus the minimum of failure time
  and dropout time.

- `fail`: Indicator that `cte` was set using failure time; i.e., 1 is a
  failure, 0 is a dropout.

The data frame also has the attribute "ratio", which is calculated as
the number of "experimental" treatments divided by the number of
"control" treatments from the input argument `block`.

## Examples

``` r
library(dplyr)

# Example 1
sim_pw_surv(n = 20)
#>    stratum enroll_time    treatment  fail_time dropout_time       cte fail
#> 1      All  0.02142839 experimental  0.3408276     98.58106  0.362256    1
#> 2      All  0.02391720      control 10.6974702    126.02694 10.721387    1
#> 3      All  0.12132222      control 12.2867886    222.00638 12.408111    1
#> 4      All  0.33721875 experimental 48.4943098   3940.42836 48.831529    1
#> 5      All  0.34688264      control  5.5165464    747.72368  5.863429    1
#> 6      All  0.40649620      control 46.7612306   1035.36115 47.167727    1
#> 7      All  0.47819083 experimental  4.6297940    591.70021  5.107985    1
#> 8      All  0.61818734 experimental 25.3045387   1074.12731 25.922726    1
#> 9      All  0.75531457      control 12.3506497    280.66035 13.105964    1
#> 10     All  0.77138700 experimental  8.2785049    460.90846  9.049892    1
#> 11     All  0.78564859      control 27.5026277    159.48258 28.288276    1
#> 12     All  0.86767983 experimental  1.6351292    480.78425  2.502809    1
#> 13     All  0.98040092      control 24.6939029   2312.83176 25.674304    1
#> 14     All  1.09322355      control  2.0823743    195.25467  3.175598    1
#> 15     All  1.37734513 experimental 19.0080405     80.18976 20.385386    1
#> 16     All  1.42045481 experimental  0.7551894    278.42173  2.175644    1
#> 17     All  1.66141361      control  2.5583299    204.01239  4.219743    1
#> 18     All  1.73594205 experimental 32.7148660   1376.63546 34.450808    1
#> 19     All  1.92466435 experimental  0.1472751    176.17803  2.071939    1
#> 20     All  1.96150812      control 10.4512454   2664.66156 12.412754    1

# Example 2
# 3:1 randomization
sim_pw_surv(
  n = 20,
  block = c(rep("experimental", 3), "control")
)
#>    stratum enroll_time    treatment   fail_time dropout_time        cte fail
#> 1      All   0.3217242 experimental   9.4928294   742.229574   9.814554    1
#> 2      All   0.3565670 experimental  12.4943694   685.656819  12.850936    1
#> 3      All   0.3830329      control   7.9711034  3357.369819   8.354136    1
#> 4      All   0.4244951 experimental  41.0686377  1134.605316  41.493133    1
#> 5      All   0.4762787 experimental 116.8281381  1848.415820 117.304417    1
#> 6      All   0.4887699      control  28.9180215   146.013236  29.406791    1
#> 7      All   0.6136093 experimental   1.1196077   162.507928   1.733217    1
#> 8      All   0.6891445 experimental  44.2434604     4.963087   5.652231    0
#> 9      All   0.7894425 experimental   4.5214569   456.380132   5.310899    1
#> 10     All   1.0188476      control   5.8872020  1544.620162   6.906050    1
#> 11     All   1.0209447 experimental   6.9775576   276.484136   7.998502    1
#> 12     All   1.2873632 experimental  14.8730426   426.240944  16.160406    1
#> 13     All   1.4099907 experimental   0.1715598    30.568471   1.581550    1
#> 14     All   1.4468614 experimental  10.0687941  3025.328402  11.515656    1
#> 15     All   1.4532458      control   2.4995260  1411.869284   3.952772    1
#> 16     All   1.5900233 experimental  17.4354083  1147.682266  19.025432    1
#> 17     All   1.6577903 experimental  99.0519123    76.763221  78.421012    0
#> 18     All   1.8990358 experimental  19.0166778   440.688450  20.915714    1
#> 19     All   2.1388782      control  19.5929947   286.458360  21.731873    1
#> 20     All   2.2620772 experimental   3.6437051  1873.164324   5.905782    1

# Example 3
# Simulate 2 stratum; will use defaults for blocking and enrollRates
sim_pw_surv(
  n = 20,
  # 2 stratum,30% and 70% prevalence
  stratum = data.frame(stratum = c("Low", "High"), p = c(.3, .7)),
  fail_rate = data.frame(
    stratum = c(rep("Low", 4), rep("High", 4)),
    period = rep(1:2, 4),
    treatment = rep(c(
      rep("control", 2),
      rep("experimental", 2)
    ), 2),
    duration = rep(c(3, 1), 4),
    rate = c(.03, .05, .03, .03, .05, .08, .07, .04)
  ),
  dropout_rate = data.frame(
    stratum = c(rep("Low", 2), rep("High", 2)),
    period = rep(1, 4),
    treatment = rep(c("control", "experimental"), 2),
    duration = rep(1, 4),
    rate = rep(.001, 4)
  )
)
#>    stratum enroll_time    treatment fail_time dropout_time       cte fail
#> 1     High  0.05272381 experimental  9.066633    494.34540  9.119357    1
#> 2      Low  0.13566515 experimental 78.546249    128.42704 78.681914    1
#> 3     High  0.26260540 experimental 54.966515     21.84553 22.108138    0
#> 4     High  0.33482645      control  7.189508   2008.07518  7.524335    1
#> 5      Low  0.48250535      control 14.974977    517.85966 15.457482    1
#> 6     High  0.68305515      control  7.064524    513.65581  7.747579    1
#> 7     High  0.81710026      control 21.438898   1164.23533 22.255998    1
#> 8     High  1.07033538      control 19.713971   2212.00385 20.784306    1
#> 9      Low  1.41332672      control 19.401860   2560.64452 20.815187    1
#> 10    High  1.61703621 experimental 17.517759    699.12197 19.134795    1
#> 11    High  1.61704018 experimental  7.688537    788.59361  9.305577    1
#> 12     Low  1.68979094 experimental  8.225572     89.70154  9.915363    1
#> 13     Low  1.73695233      control  4.916503   1871.11105  6.653455    1
#> 14     Low  1.78646504      control  6.955286    276.59924  8.741751    1
#> 15    High  1.79986228 experimental 16.726290    382.11254 18.526153    1
#> 16     Low  1.80771972 experimental 40.504842    653.13538 42.312562    1
#> 17     Low  1.90859178 experimental  4.846215    170.22468  6.754807    1
#> 18    High  1.99798037      control 43.722842     23.72062 25.718601    0
#> 19    High  2.00729138      control 15.319417    661.56404 17.326709    1
#> 20    High  2.02787066 experimental 22.863487    134.13149 24.891357    1
# Example 4
# If you want a more rectangular entry for a data.frame
fail_rate <- bind_rows(
  data.frame(stratum = "Low", period = 1, treatment = "control", duration = 3, rate = .03),
  data.frame(stratum = "Low", period = 1, treatment = "experimental", duration = 3, rate = .03),
  data.frame(stratum = "Low", period = 2, treatment = "experimental", duration = 3, rate = .02),
  data.frame(stratum = "High", period = 1, treatment = "control", duration = 3, rate = .05),
  data.frame(stratum = "High", period = 1, treatment = "experimental", duration = 3, rate = .06),
  data.frame(stratum = "High", period = 2, treatment = "experimental", duration = 3, rate = .03)
)

dropout_rate <- bind_rows(
  data.frame(stratum = "Low", period = 1, treatment = "control", duration = 3, rate = .001),
  data.frame(stratum = "Low", period = 1, treatment = "experimental", duration = 3, rate = .001),
  data.frame(stratum = "High", period = 1, treatment = "control", duration = 3, rate = .001),
  data.frame(stratum = "High", period = 1, treatment = "experimental", duration = 3, rate = .001)
)

sim_pw_surv(
  n = 12,
  stratum = data.frame(stratum = c("Low", "High"), p = c(.3, .7)),
  fail_rate = fail_rate,
  dropout_rate = dropout_rate
)
#>    stratum enroll_time    treatment fail_time dropout_time        cte fail
#> 1      Low   0.1321033 experimental 49.865443   1562.73537  49.997547    1
#> 2     High   0.3007465 experimental 18.103313    117.95707  18.404059    1
#> 3     High   0.4960628      control 13.752074    908.18929  14.248136    1
#> 4     High   0.5495319 experimental 36.791951    599.65573  37.341483    1
#> 5      Low   0.7205039 experimental 22.158442   1352.21930  22.878946    1
#> 6     High   0.8137968      control 34.253180    930.71607  35.066977    1
#> 7     High   0.8824684      control  8.837317     13.13534   9.719785    1
#> 8     High   0.9442470 experimental 53.139230   1493.85426  54.083477    1
#> 9     High   0.9865067 experimental 57.256334   1265.50109  58.242841    1
#> 10     Low   1.2318692      control 41.671776    403.53936  42.903645    1
#> 11    High   1.2658735      control 56.696444     46.61400  47.879875    0
#> 12     Low   1.4219050      control 98.788315    553.18286 100.210220    1
```
