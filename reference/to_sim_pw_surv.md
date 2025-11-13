# Convert enrollment and failure rates from `sim_fixed_n()` to `sim_pw_surv()` format

`to_sim_pw_surv()` converts failure rates and dropout rates entered in
the simpler format for
[`sim_fixed_n()`](https://merck.github.io/simtrial/reference/sim_fixed_n.md)
to that used for
[`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md).
The `fail_rate` argument for
[`sim_fixed_n()`](https://merck.github.io/simtrial/reference/sim_fixed_n.md)
requires enrollment rates, failure rates hazard ratios and dropout rates
by stratum for a 2-arm trial,
[`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md)
is in a more flexible but less obvious but more flexible format. Since
[`sim_fixed_n()`](https://merck.github.io/simtrial/reference/sim_fixed_n.md)
automatically analyzes data and
[`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md)
just produces a simulation dataset, the latter provides additional
options to analyze or otherwise evaluate individual simulations in ways
that
[`sim_fixed_n()`](https://merck.github.io/simtrial/reference/sim_fixed_n.md)
does not.

## Usage

``` r
to_sim_pw_surv(
  fail_rate = data.frame(stratum = "All", duration = c(3, 100), fail_rate = log(2)/c(9,
    18), hr = c(0.9, 0.6), dropout_rate = rep(0.001, 2))
)
```

## Arguments

- fail_rate:

  Piecewise constant control group failure rates, hazard ratio for
  experimental vs. control, and dropout rates by stratum and time
  period.

## Value

A list of two data frame components formatted for
[`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md):
`fail_rate` and `dropout_rate`.

## Examples

``` r
# Example 1
# Convert standard input
to_sim_pw_surv()
#> $fail_rate
#>   stratum period    treatment duration       rate
#> 1     All      1      control        3 0.07701635
#> 2     All      2      control      100 0.03850818
#> 3     All      1 experimental        3 0.06931472
#> 4     All      2 experimental      100 0.02310491
#> 
#> $dropout_rate
#>   stratum period    treatment duration  rate
#> 1     All      1      control        3 0.001
#> 2     All      2      control      100 0.001
#> 3     All      1 experimental        3 0.001
#> 4     All      2 experimental      100 0.001
#> 

# Stratified example
fail_rate <- data.frame(
  stratum = c(rep("Low", 3), rep("High", 3)),
  duration = rep(c(4, 10, 100), 2),
  fail_rate = c(
    .04, .1, .06,
    .08, .16, .12
  ),
  hr = c(
    1.5, .5, 2 / 3,
    2, 10 / 16, 10 / 12
  ),
  dropout_rate = .01
)

x <- to_sim_pw_surv(fail_rate)

# Do a single simulation with the above rates
# Enroll 300 patients in ~12 months at constant rate
sim <- sim_pw_surv(
  n = 300,
  stratum = data.frame(stratum = c("Low", "High"), p = c(.6, .4)),
  enroll_rate = data.frame(duration = 12, rate = 300 / 12),
  fail_rate = x$fail_rate,
  dropout_rate = x$dropout_rate
)

# Cut after 200 events and do a stratified logrank test
sim |>
  cut_data_by_event(200) |> # Cut data
  wlr(weight = fh(rho = 0, gamma = 0)) # Stratified logrank
#> $method
#> [1] "WLR"
#> 
#> $parameter
#> [1] "FH(rho=0, gamma=0)"
#> 
#> $estimate
#> [1] 2.340246
#> 
#> $se
#> [1] 6.972734
#> 
#> $z
#> [1] -0.3356281
#> 
#> $info
#> [1] 49.92
#> 
#> $info0
#> [1] 50
#> 
```
