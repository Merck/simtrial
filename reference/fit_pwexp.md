# Piecewise exponential survival estimation

Computes survival function, density function, -2 \* log-likelihood based
on input dataset and intervals for piecewise constant failure rates.
Initial version assumes observations are right censored or events only.

## Usage

``` r
fit_pwexp(
  srv = Surv(time = ex1_delayed_effect$month, event = ex1_delayed_effect$evntd),
  intervals = array(3, 3)
)
```

## Arguments

- srv:

  Input survival object (see
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html));
  note that only 0 = censored, 1 = event for
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html).

- intervals:

  Vector containing positive values indicating interval lengths where
  the exponential rates are assumed. Note that a final infinite interval
  is added if any events occur after the final interval specified.

## Value

A matrix with rows containing interval length, estimated rate, -2 \*
log-likelihood for each interval.

## Examples

``` r
# Use default arguments for delayed effect example dataset (ex1_delayed_effect)
library(survival)

# Example 1
rateall <- fit_pwexp()
rateall
#>   intervals     ttot event       rate     m2ll
#> 1         3 937.1785    97 0.10350216 634.0236
#> 2         3 605.3572    71 0.11728612 446.3257
#> 3         3 346.8482    30 0.08649317 206.8614
#> 4       Inf 254.1148    20 0.07870458 141.6822

# Example 2
# Estimate by treatment effect
rate1 <- with(subset(ex1_delayed_effect, trt == 1), fit_pwexp(Surv(month, evntd)))
rate0 <- with(subset(ex1_delayed_effect, trt == 0), fit_pwexp(Surv(month, evntd)))

rate1
#>   intervals     ttot event       rate      m2ll
#> 1         3 620.4375    64 0.10315302 418.75734
#> 2         3 415.8482    36 0.08657005 248.16970
#> 3         3 256.2053    19 0.07415927 136.85853
#> 4       Inf 205.4186    13 0.06328542  97.76261
rate0
#>   intervals      ttot event      rate      m2ll
#> 1         3 316.74106    33 0.1041861 215.26408
#> 2         3 189.50899    35 0.1846878 188.23619
#> 3         3  90.64288    11 0.1213554  68.39871
#> 4       Inf  48.69624     7 0.1437483  41.15568
rate1$rate / rate0$rate
#> [1] 0.9900847 0.4687372 0.6110917 0.4402517

# Chi-square test for (any) treatment effect (8 - 4 parameters = 4 df)
pchisq(sum(rateall$m2ll) - sum(rate1$m2ll + rate0$m2ll),
  df = 4,
  lower.tail = FALSE
)
#> [1] 0.006424744

# Compare with logrank
survdiff(formula = Surv(month, evntd) ~ trt, data = ex1_delayed_effect)
#> Call:
#> survdiff(formula = Surv(month, evntd) ~ trt, data = ex1_delayed_effect)
#> 
#>         N Observed Expected (O-E)^2/E (O-E)^2/V
#> trt=0 121       86     67.7      4.97      7.35
#> trt=1 240      132    150.3      2.24      7.35
#> 
#>  Chisq= 7.3  on 1 degrees of freedom, p= 0.007 

# Example 3
# Simple model with 3 rates same for each for 3 months,
# different for each treatment after months
rate1a <- with(subset(ex1_delayed_effect, trt == 1), fit_pwexp(Surv(month, evntd), 3))
rate0a <- with(subset(ex1_delayed_effect, trt == 0), fit_pwexp(Surv(month, evntd), 3))
rate1a$rate / rate0a$rate
#> [1] 0.9900847 0.4808339

m2ll0 <- rateall$m2ll[1] + rate1a$m2ll[2] + rate0a$m2ll[2]
m2ll1 <- sum(rate0$m2ll) + sum(rate1$m2ll)

# As a measure of strength, chi-square examines improvement in likelihood
pchisq(m2ll0 - m2ll1, df = 5, lower.tail = FALSE)
#> [1] 0.741822
```
