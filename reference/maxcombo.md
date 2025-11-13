# MaxCombo test

WARNING: This experimental function is a work-in-progress. The function
arguments will change as we add additional features.

## Usage

``` r
maxcombo(
  data = cut_data_by_event(sim_pw_surv(n = 200), 150),
  rho = c(0, 0, 1),
  gamma = c(0, 1, 1),
  return_variance = FALSE,
  return_corr = FALSE
)
```

## Arguments

- data:

  A TTE dataset.

- rho:

  Numeric vector. Must be greater than or equal to zero. Must be the
  same length as `gamma`.

- gamma:

  Numeric vector. Must be greater than or equal to zero. Must be the
  same length as `rho`.

- return_variance:

  A logical flag that, if `TRUE`, adds columns estimated variance for
  weighted sum of observed minus expected; see details; Default:
  `FALSE`.

- return_corr:

  A logical flag that, if `TRUE`, adds columns estimated correlation for
  weighted sum of observed minus expected; see details; Default:
  `FALSE`.

## Value

A list containing the test method (`method`), parameters of this test
method (`parameter`), point estimate of the treatment effect
(`estimate`), standardized error of the treatment effect (`se`), Z-score
of each test of the MaxCombo (`z`), p-values (`p_value`) and the
correlation matrix of each tests in MaxCombo (begin with `v`)

## See also

[`wlr()`](https://merck.github.io/simtrial/reference/wlr.md),
[`rmst()`](https://merck.github.io/simtrial/reference/rmst.md),
[`milestone()`](https://merck.github.io/simtrial/reference/milestone.md)

## Examples

``` r
sim_pw_surv(n = 200) |>
  cut_data_by_event(150) |>
  maxcombo(rho = c(0, 0), gamma = c(0, 1), return_corr = TRUE)
#> $method
#> [1] "MaxCombo"
#> 
#> $parameter
#> [1] "FH(0, 0) + FH(0, 1)"
#> 
#> $z
#> [1] -1.429270 -1.925974
#> 
#> $corr
#>          v1        v2
#> 1 1.0000000 0.8572258
#> 2 0.8572258 1.0000000
#> 
#> $p_value
#> [1] 0.03992553
#> 
```
