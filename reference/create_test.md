# Create a cutting test function

Create a cutting test function for use with
[`sim_gs_n()`](https://merck.github.io/simtrial/reference/sim_gs_n.md)

## Usage

``` r
create_test(test, ...)
```

## Arguments

- test:

  A test function such as
  [`wlr()`](https://merck.github.io/simtrial/reference/wlr.md),
  [`maxcombo()`](https://merck.github.io/simtrial/reference/maxcombo.md),
  or [`rmst()`](https://merck.github.io/simtrial/reference/rmst.md)

- ...:

  Arguments passed to the cutting test function

## Value

A function that accepts a data frame of simulated trial data and returns
a test result

## See also

[`sim_gs_n()`](https://merck.github.io/simtrial/reference/sim_gs_n.md),
[`create_cut()`](https://merck.github.io/simtrial/reference/create_cut.md)

## Examples

``` r
# Simulate trial data
trial_data <- sim_pw_surv()

# Cut after 150 events
trial_data_cut <- cut_data_by_event(trial_data, 150)

# Create a cutting test function that can be used by sim_gs_n()
regular_logrank_test <- create_test(wlr, weight = fh(rho = 0, gamma = 0))

# Test the cutting
regular_logrank_test(trial_data_cut)
#> $method
#> [1] "WLR"
#> 
#> $parameter
#> [1] "FH(rho=0, gamma=0)"
#> 
#> $estimate
#> [1] -16.60282
#> 
#> $se
#> [1] 4.370912
#> 
#> $z
#> [1] 3.798481
#> 
#> $info
#> [1] 23.1828
#> 
#> $info0
#> [1] 23.25
#> 

# The results are the same as directly calling the function
stopifnot(all.equal(
  regular_logrank_test(trial_data_cut),
  wlr(trial_data_cut, weight = fh(rho = 0, gamma = 0))
))
```
