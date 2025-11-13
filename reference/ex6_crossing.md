# Time-to-event data example 6 for non-proportional hazards working group

Survival objects reverse-engineered datasets from published Kaplan-Meier
curves. Individual trials are de-identified since the data are only
approximations of the actual data. Data are intended to evaluate methods
and designs for trials where non-proportional hazards may be anticipated
for outcome data.

## Usage

``` r
data(ex6_crossing)
```

## Format

Data frame with 4 variables:

- `id`: Sequential numbering of unique identifiers.

- `month`: Time-to-event.

- `event`: 1 for event, 0 for censored.

- `trt`: 1 for experimental, 0 for control.

## References

Lin, Ray S., Ji Lin, Satrajit Roychoudhury, Keaven M. Anderson, Tianle
Hu, Bo Huang, Larry F Leon, Jason J.Z. Liao, Rong Liu, Xiaodong Luo,
Pralay Mukhopadhyay, Rui Qin, Kay Tatsuoka, Xuejing Wang, Yang Wang,
Jian Zhu, Tai-Tsang Chen, Renee Iacona & Cross-Pharma Non-proportional
Hazards Working Group. 2020. Alternative analysis methods for time to
event endpoints under nonproportional hazards: A comparative analysis.
*Statistics in Biopharmaceutical Research* 12(2): 187–198.

## See also

[ex1_delayed_effect](https://merck.github.io/simtrial/reference/ex1_delayed_effect.md),
[ex2_delayed_effect](https://merck.github.io/simtrial/reference/ex2_delayed_effect.md),
[ex3_cure_with_ph](https://merck.github.io/simtrial/reference/ex3_cure_with_ph.md),
[ex4_belly](https://merck.github.io/simtrial/reference/ex4_belly.md),
[ex5_widening](https://merck.github.io/simtrial/reference/ex5_widening.md)

## Examples

``` r
library(survival)

data(ex6_crossing)
km1 <- with(ex6_crossing, survfit(Surv(month, evntd) ~ trt))
km1
#> Call: survfit(formula = Surv(month, evntd) ~ trt)
#> 
#>         n events median 0.95LCL 0.95UCL
#> trt=0 145    111  10.66    8.83    12.5
#> trt=1 145    113   9.92    7.38    14.3
plot(km1)
```
