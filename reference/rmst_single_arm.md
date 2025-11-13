# Calculate RMST for a single cut-off time point

Calculate RMST for a single cut-off time point

## Usage

``` r
rmst_single_arm(
  time_var,
  event_var,
  tau,
  group_label = "Single Group",
  alpha = 0.05
)
```

## Arguments

- time_var:

  A numeric vector of follow up time.

- event_var:

  A numeric or integer vector of the status indicator; 0=alive 1=event.

- tau:

  A value of pre-defined cut-off time point.

- group_label:

  A character of customized treatment group name.

- alpha:

  A numeric value of the significant level for RMST confidence interval.
  Default is 0.05.

## Value

A data frame of

- Cutoff time: same as `tau`;

- Group label: same as `group_label`;

- Estimated RMST;

- Variance, standard error, and CIs of the estimated RMST;

- Number of events.

## Examples

``` r
data(ex1_delayed_effect)
data_single_arm <- ex1_delayed_effect[ex1_delayed_effect$trt == 1, ]
simtrial:::rmst_single_arm(
  time_var = data_single_arm$month,
  event_var = data_single_arm$evntd,
  tau = 10,
  group_label = "Treatment 1",
  alpha = 0.05
)
#>   cutoff_time       group     rmst   variance       std      lcl      ucl event
#> 1          10 Treatment 1 6.495175 0.05711322 0.2389837 6.026776 6.963575   127
```
