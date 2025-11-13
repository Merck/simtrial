# Create a cutting function

Create a cutting function for use with
[`sim_gs_n()`](https://merck.github.io/simtrial/reference/sim_gs_n.md)

## Usage

``` r
create_cut(...)
```

## Arguments

- ...:

  Arguments passed to
  [`get_analysis_date()`](https://merck.github.io/simtrial/reference/get_analysis_date.md)

## Value

A function that accepts a data frame of simulated trial data and returns
a cut date

## See also

[`get_analysis_date()`](https://merck.github.io/simtrial/reference/get_analysis_date.md),
[`sim_gs_n()`](https://merck.github.io/simtrial/reference/sim_gs_n.md)

## Examples

``` r
# Simulate trial data
trial_data <- sim_pw_surv()

# Create a cutting function that applies the following 2 conditions:
# - At least 45 months have passed since the start of the study
# - At least 300 events have occurred
cutting <- create_cut(
  planned_calendar_time = 45,
  target_event_overall = 350
)

# Cut the trial data
cutting(trial_data)
#> [1] 77.87317
```
