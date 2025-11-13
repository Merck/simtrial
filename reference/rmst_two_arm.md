# Calculate RMST difference

Calculate RMST difference

## Usage

``` r
rmst_two_arm(
  time_var,
  event_var,
  group_var,
  trunc_time,
  reference = sort(unique(group_var))[1],
  alpha = 0.05
)
```

## Arguments

- time_var:

  A numeric vector of follow up time.

- event_var:

  A numeric or integer vector of the status indicator; 0=alive 1=event.

- group_var:

  A vector of treatment groups.

- trunc_time:

  A numeric vector of pre-defined cut-off time point(s).

- reference:

  Group name of reference group for RMST comparison. Default is the
  first group name by alphabetical order.

- alpha:

  A numeric value of the significant level for RMST confidence interval.
  Default is 0.05.

## Value

A list of 2 data frames of RMST calculations:

- `rmst_per_arm`: the calculation results per group.

- `rmst_diff`: the calculation results of RMST differences.

## Examples

``` r
data(ex1_delayed_effect)
with(
  ex1_delayed_effect,
  simtrial:::rmst_two_arm(
    time_var = month,
    event_var = evntd,
    group_var = trt,
    trunc_time = 6,
    reference = "0",
    alpha = 0.05
  )
)
#> $rmst_per_arm
#>   cutoff_time group     rmst   variance       std      lcl      ucl event
#> 1           6     0 4.340067 0.02902105 0.1703557 4.006176 4.673958    68
#> 2           6     1 4.552177 0.01607455 0.1267854 4.303682 4.800672   100
#> 
#> $rmst_diff
#>   cutoff_time group rmst_diff  variance       std        lcl       ucl
#> 1           6 1 - 0 0.2121097 0.0450956 0.2123572 -0.2041029 0.6283222
#> 
```
