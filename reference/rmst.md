# RMST difference of 2 arms

RMST difference of 2 arms

## Usage

``` r
rmst(
  data,
  tau = 10,
  var_label_tte = "tte",
  var_label_event = "event",
  var_label_group = "treatment",
  formula = NULL,
  reference = "control",
  alpha = 0.05
)
```

## Arguments

- data:

  A time-to-event dataset with a column `tte` indicating the survival
  time and a column of `event` indicating whether it is event or censor.

- tau:

  RMST analysis time.

- var_label_tte:

  Column name of the TTE variable.

- var_label_event:

  Column name of the event variable.

- var_label_group:

  Column name of the grouping variable.

- formula:

  (default: `NULL`) A formula that indicates the TTE, event, and group
  variables using the syntax `Surv(tte, event) ~ group)` (see Details
  below for more information). This is an alternative to specifying the
  variables as strings. If a formula is provided, the values passed to
  `var_label_tte`, `var_label_event`, and `var_label_group` are ignored.

- reference:

  A group label indicating the reference group.

- alpha:

  Type I error.

## Value

The z statistics.

## Details

The argument `formula` is provided as a convenience to easily specify
the TTE, event, and grouping variables using the syntax
`Surv(tte, event) ~ group)`.
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) is from the
{survival} package
([`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)). You
can also explicitly name the arguments passed to
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html), for example the
following is equivalent `Surv(event = event, time = tte) ~ group)`. Note
however that the function
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) is never actually
executed. Similarly, any other functions applied in the formula are also
ignored, thus you shouldn't apply any transformation functions such as
[`log()`](https://rdrr.io/r/base/Log.html) since these will have no
effect.

## Examples

``` r
data(ex1_delayed_effect)
rmst(
  data = ex1_delayed_effect,
  var_label_tte = "month",
  var_label_event = "evntd",
  var_label_group = "trt",
  tau = 10,
  reference = "0"
)
#> $method
#> [1] "RMST"
#> 
#> $parameter
#> [1] 10
#> 
#> $estimate
#> [1] 0.8650493
#> 
#> $se
#> [1] 0.3900344
#> 
#> $z
#> [1] 2.21788
#> 

# Formula interface
rmst(
  data = ex1_delayed_effect,
  formula = Surv(month, evntd) ~ trt,
  tau = 10,
  reference = "0"
)
#> $method
#> [1] "RMST"
#> 
#> $parameter
#> [1] 10
#> 
#> $estimate
#> [1] 0.8650493
#> 
#> $se
#> [1] 0.3900344
#> 
#> $z
#> [1] 2.21788
#> 
```
