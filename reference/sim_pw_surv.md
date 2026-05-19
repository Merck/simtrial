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
#> 1      All   0.1366848      control  3.8388328    847.48951  3.975518    1
#> 2      All   0.3524424 experimental  1.0412043   1901.83435  1.393647    1
#> 3      All   0.3738708      control  3.3195804     58.16193  3.693451    1
#> 4      All   0.3763596 experimental  4.2301976    160.37688  4.606557    1
#> 5      All   0.4737646      control  5.1065099   1375.31437  5.580275    1
#> 6      All   0.6896611 experimental 32.7491726    197.03324 33.438834    1
#> 7      All   0.6993250      control  3.7912304     11.34259  4.490555    1
#> 8      All   0.7589386 experimental  2.2875405    804.91681  3.046479    1
#> 9      All   0.8306332 experimental 18.3949403    126.02694 19.225574    1
#> 10     All   0.9706297      control  0.3408276     98.58106  1.311457    1
#> 11     All   1.1077569      control 25.7471549   3940.42836 26.854912    1
#> 12     All   1.1238294 experimental 21.5735772    222.00638 22.697407    1
#> 13     All   1.1380910      control  3.8148970    591.70021  4.952988    1
#> 14     All   1.2201222 experimental  8.0330927    747.72368  9.253215    1
#> 15     All   1.3328433      control 14.1522694   1074.12731 15.485113    1
#> 16     All   1.4456659 experimental 90.5224612   1035.36115 91.968127    1
#> 17     All   1.7297875 experimental 21.7012994    280.66035 23.431087    1
#> 18     All   1.7728972      control  5.6392525    460.90846  7.412150    1
#> 19     All   2.0138560      control  1.6351292    480.78425  3.648985    1
#> 20     All   2.0883844 experimental 52.0052554    159.48258 54.093640    1

# Example 2
# 3:1 randomization
sim_pw_surv(
  n = 20,
  block = c(rep("experimental", 3), "control")
)
#>    stratum enroll_time    treatment   fail_time dropout_time        cte fail
#> 1      All  0.08448291 experimental  34.5768873   847.822120  34.661370    1
#> 2      All  0.10368791 experimental   3.5338628   255.836968   3.637551    1
#> 3      All  0.42541209      control   5.7220114   192.504381   6.147423    1
#> 4      All  0.46025494 experimental   9.4928294   742.229574   9.953084    1
#> 5      All  0.48672076 experimental  12.4943694   685.656819  12.981090    1
#> 6      All  0.52818303 experimental  41.0686377  1134.605316  41.596821    1
#> 7      All  0.57996659 experimental 116.8281381  1848.415820 117.408105    1
#> 8      All  0.59245783      control  24.3216439  1508.981005  24.914102    1
#> 9      All  0.71729722 experimental   1.1196077   162.507928   1.836905    1
#> 10     All  0.79283241 experimental  44.2434604     4.963087   5.755919    0
#> 11     All  0.89313041      control   7.9711034  3357.369819   8.864234    1
#> 12     All  1.12253551 experimental   4.5214569   456.380132   5.643992    1
#> 13     All  1.12463257 experimental   6.9775576   276.484136   8.102190    1
#> 14     All  1.39105109      control  28.9180215   146.013236  30.309073    1
#> 15     All  1.51367859 experimental  14.8730426   426.240944  16.386721    1
#> 16     All  1.55054933 experimental   0.1715598    30.568471   1.722109    1
#> 17     All  1.55693371 experimental  10.0687941  3025.328402  11.625728    1
#> 18     All  1.69371116      control   5.8872020  1544.620162   7.580913    1
#> 19     All  1.76147824 experimental  17.4354083  1147.682266  19.196887    1
#> 20     All  2.00272370 experimental  99.0519123    76.763221  78.765945    0

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
#>    stratum enroll_time    treatment  fail_time dropout_time       cte fail
#> 1     High   0.3141897 experimental  3.7017709    2288.6606  4.015961    1
#> 2      Low   0.9300871 experimental 66.9358395     661.5640 67.865927    1
#> 3     High   0.9828109 experimental 30.6210009     790.7104 31.603812    1
#> 4      Low   1.0657523      control  3.7685408     688.7488  4.834293    1
#> 5     High   1.1926925      control  7.3043175     485.1607  8.497010    1
#> 6     High   1.2649136      control  0.4369107     475.1619  1.701824    1
#> 7      Low   1.4125925      control  2.9900514     910.0930  4.402644    1
#> 8     High   1.6131423      control  9.8640246    1625.1118 11.477167    1
#> 9     High   1.7471874 experimental 41.5564604     397.5415 43.303648    1
#> 10    High   2.0004225 experimental 18.7568570     759.0516 20.757279    1
#> 11     Low   2.3434139 experimental 17.1218605    2356.3875 19.465274    1
#> 12    High   2.5471233      control 10.9824201    1487.1176 13.529543    1
#> 13    High   2.5471273      control  5.9014068    3407.8273  8.448534    1
#> 14     Low   2.6198781 experimental 38.8078443     246.7672 41.427722    1
#> 15     Low   2.6670395      control 14.2627076     185.8252 16.929747    1
#> 16     Low   2.7165522 experimental 73.7334616    1215.1453 76.450014    1
#> 17    High   2.7299494 experimental  9.0666334    1004.5395 11.796583    1
#> 18     Low   2.7378069      control  4.6044936     287.7643  7.342300    1
#> 19     Low   2.8386789 experimental  0.7906874     145.3865  3.629366    1
#> 20    High   2.9280675      control  2.6826298    1135.5534  5.610697    1
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
#>    stratum enroll_time    treatment fail_time dropout_time      cte fail
#> 1     High   0.2006600 experimental  18.52521   1352.21930 18.72587    1
#> 2     High   0.2630211      control  11.06366    117.95707 11.32668    1
#> 3     High   0.3951244 experimental  31.24363   1250.15328 31.63875    1
#> 4      Low   0.5637675 experimental  84.13295   2834.82220 84.69672    1
#> 5      Low   0.7590838 experimental  20.59329    908.18929 21.35238    1
#> 6     High   0.8125530      control  12.66199    599.65573 13.47454    1
#> 7     High   0.9835250      control  23.87517   1493.85426 24.85870    1
#> 8     High   1.0768179 experimental  12.77229   2963.64946 13.84911    1
#> 9      Low   1.1454895      control  31.02387     13.13534 14.28083    0
#> 10    High   1.2072681 experimental  49.09118    403.53936 50.29845    1
#> 11    High   1.2495278      control  33.68354   1265.50109 34.93307    1
#> 12    High   1.4948903      control  36.15380    687.60368 37.64869    1
```
