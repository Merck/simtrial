# Milestone test for two survival curves

Milestone test for two survival curves

## Usage

``` r
milestone(data, ms_time, test_type = c("log-log", "naive"))
```

## Arguments

- data:

  Data frame containing at least 3 columns:

  - `tte` - Time to event.

  - `event` - Event indicator.

  - `treatment` - Grouping variable.

- ms_time:

  Milestone analysis time.

- test_type:

  Method to build the test statistics. There are 2 options:

  - `"naive"`: a naive approach by dividing the KM survival difference
    by its standard derivations, see equation (1) of Klein, J. P.,
    Logan, B., Harhoff, M., & Andersen, P. K. (2007).

  - `"log-log"`: a log-log transformation of the survival, see
    equation (3) of Klein, J. P., Logan, B., Harhoff, M., &
    Andersen, P. K. (2007).

## Value

A list frame containing:

- `method` - The method, always `"milestone"`.

- `parameter` - Milestone time point.

- `estimate` - Survival difference between the experimental and control
  arm.

- `se` - Standard error of the control and experimental arm.

- `z` - Test statistics.

## References

Klein, J. P., Logan, B., Harhoff, M., & Andersen, P. K. (2007).
"Analyzing survival curves at a fixed point in time." *Statistics in
Medicine*, 26(24), 4505–4519.

## Examples

``` r
cut_data <- sim_pw_surv(n = 200) |>
  cut_data_by_event(150)

cut_data |>
  milestone(10, test_type = "log-log")
#> $method
#> [1] "milestone"
#> 
#> $parameter
#> [1] 10
#> 
#> $estimate
#> [1] 0.3090112
#> 
#> $se
#> [1] 0.2084854
#> 
#> $z
#> [1] 1.482171
#> 

cut_data |>
  milestone(10, test_type = "naive")
#> $method
#> [1] "milestone"
#> 
#> $parameter
#> [1] 10
#> 
#> $estimate
#> [1] 0.105612
#> 
#> $se
#> [1] 0.07066159
#> 
#> $z
#> [1] 1.494617
#> 
```
