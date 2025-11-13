# Fleming-Harrington weighting function

Fleming-Harrington weighting function

## Usage

``` r
fh(rho = 0, gamma = 0)
```

## Arguments

- rho:

  Non-negative number. `rho = 0, gamma = 0` is equivalent to regular
  logrank test.

- gamma:

  Non-negative number. `rho = 0, gamma = 0` is equivalent to regular
  logrank test.

## Value

A list of parameters of the Fleming-Harrington weighting function

## Examples

``` r
sim_pw_surv(n = 200) |>
  cut_data_by_event(100) |>
  wlr(weight = fh(rho = 0, gamma = 1))
#> $method
#> [1] "WLR"
#> 
#> $parameter
#> [1] "FH(rho=0, gamma=1)"
#> 
#> $estimate
#> [1] -2.354546
#> 
#> $se
#> [1] 1.596507
#> 
#> $z
#> [1] 1.474812
#> 
#> $info
#> [1] 2.598126
#> 
#> $info0
#> [1] 2.625333
#> 
```
