# Zero early weighting function

Zero early weighting function

## Usage

``` r
early_zero(early_period, fail_rate = NULL)
```

## Arguments

- early_period:

  The initial delay period where weights increase; after this, weights
  are constant at the final weight in the delay period.

- fail_rate:

  Failure rate

## Value

A list of parameters of the zero early weighting function

## References

Xu, Z., Zhen, B., Park, Y., & Zhu, B. (2017). "Designing therapeutic
cancer vaccine trials with delayed treatment effect."

## Examples

``` r
library(gsDesign2)
#> 
#> Attaching package: ‘gsDesign2’
#> The following object is masked from ‘package:simtrial’:
#> 
#>     as_gt

# Example 1: Unstratified ----
sim_pw_surv(n = 200) |>
  cut_data_by_event(125) |>
  wlr(weight = early_zero(early_period = 2))
#> $method
#> [1] "WLR"
#> 
#> $parameter
#> [1] "Xu 2017 with first 2 months of 0 weights"
#> 
#> $estimate
#> [1] -19.03831
#> 
#> $se
#> [1] 4.794997
#> 
#> $z
#> [1] 3.970454
#> 
#> $info
#> [1] 22.73958
#> 
#> $info0
#> [1] 24
#> 

# Example 2: Stratified ----
n <- 500
# Two strata
stratum <- c("Biomarker-positive", "Biomarker-negative")
prevalence_ratio <- c(0.6, 0.4)

# Enrollment rate
enroll_rate <- define_enroll_rate(
  stratum = rep(stratum, each = 2),
  duration = c(2, 10, 2, 10),
  rate = c(c(1, 4) * prevalence_ratio[1], c(1, 4) * prevalence_ratio[2])
)
enroll_rate$rate <- enroll_rate$rate * n / sum(enroll_rate$duration * enroll_rate$rate)

# Failure rate
med_pos <- 10 # Median of the biomarker positive population
med_neg <- 8 # Median of the biomarker negative population
hr_pos <- c(1, 0.7) # Hazard ratio of the biomarker positive population
hr_neg <- c(1, 0.8) # Hazard ratio of the biomarker negative population
fail_rate <- define_fail_rate(
  stratum = rep(stratum, each = 2),
  duration = c(3, 1000, 4, 1000),
  fail_rate = c(log(2) / c(med_pos, med_pos, med_neg, med_neg)),
  hr = c(hr_pos, hr_neg),
  dropout_rate = 0.01
)

# Simulate data
temp <- to_sim_pw_surv(fail_rate) # Convert the failure rate
set.seed(2023)

sim_pw_surv(
  n = n, # Sample size
  # Stratified design with prevalence ratio of 6:4
  stratum = data.frame(stratum = stratum, p = prevalence_ratio),
  # Randomization ratio
  block = c("control", "control", "experimental", "experimental"),
  enroll_rate = enroll_rate, # Enrollment rate
  fail_rate = temp$fail_rate, # Failure rate
  dropout_rate = temp$dropout_rate # Dropout rate
) |>
  cut_data_by_event(125) |>
  wlr(weight = early_zero(early_period = 2, fail_rate = fail_rate))
#> $method
#> [1] "WLR"
#> 
#> $parameter
#> [1] "Xu 2017 with first 2 months of 0 weights"
#> 
#> $estimate
#> [1] 1.207753
#> 
#> $se
#> [1] 1.133941
#> 
#> $z
#> [1] -1.065093
#> 
#> $info
#> [1] 1.285211
#> 
#> $info0
#> [1] 1.298506
#> 
```
