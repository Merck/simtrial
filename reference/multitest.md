# Perform multiple tests on trial data cutting

WARNING: This experimental function is a work-in-progress. The function
arguments and/or returned output format may change as we add additional
features.

## Usage

``` r
multitest(data, ...)
```

## Arguments

- data:

  Trial data cut by
  [`cut_data_by_event()`](https://merck.github.io/simtrial/reference/cut_data_by_event.md)
  or
  [`cut_data_by_date()`](https://merck.github.io/simtrial/reference/cut_data_by_date.md)

- ...:

  One or more test functions. Use
  [`create_test()`](https://merck.github.io/simtrial/reference/create_test.md)
  to change the default arguments of each test function.

## Value

A list of test results, one per test. If the test functions are named in
the call to `multitest()`, the returned list uses the same names.

## See also

[`create_test()`](https://merck.github.io/simtrial/reference/create_test.md)

## Examples

``` r
trial_data <- sim_pw_surv(n = 200)
trial_data_cut <- cut_data_by_event(trial_data, 150)

# create cutting test functions
wlr_partial <- create_test(wlr, weight = fh(rho = 0, gamma = 0))
rmst_partial <- create_test(rmst, tau = 20)
maxcombo_partial <- create_test(maxcombo, rho = c(0, 0), gamma = c(0, 0.5))

multitest(
  data = trial_data_cut,
  wlr = wlr_partial,
  rmst = rmst_partial,
  maxcombo = maxcombo_partial
)
#> $wlr
#> $wlr$method
#> [1] "WLR"
#> 
#> $wlr$parameter
#> [1] "FH(rho=0, gamma=0)"
#> 
#> $wlr$estimate
#> [1] -21.51717
#> 
#> $wlr$se
#> [1] 5.971903
#> 
#> $wlr$z
#> [1] 3.603067
#> 
#> $wlr$info
#> [1] 36.96
#> 
#> $wlr$info0
#> [1] 37.5
#> 
#> 
#> $rmst
#> $rmst$method
#> [1] "RMST"
#> 
#> $rmst$parameter
#> [1] 20
#> 
#> $rmst$estimate
#> [1] 3.504274
#> 
#> $rmst$se
#> [1] 1.064673
#> 
#> $rmst$z
#> [1] 3.291409
#> 
#> 
#> $maxcombo
#> $maxcombo$method
#> [1] "MaxCombo"
#> 
#> $maxcombo$parameter
#> [1] "FH(0, 0) + FH(0, 0.5)"
#> 
#> $maxcombo$z
#> [1] -3.603067 -3.567134
#> 
#> $maxcombo$p_value
#> [1] 0.0002371133
#> 
#> 
```
