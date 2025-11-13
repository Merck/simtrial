# Magirr and Burman weighting function

Magirr and Burman weighting function

## Usage

``` r
mb(delay = 4, w_max = Inf)
```

## Arguments

- delay:

  The initial delay period where weights increase; after this, weights
  are constant at the final weight in the delay period.

- w_max:

  Maximum weight to be returned. Set `delay = Inf`, `w_max = 2` to be
  consistent with recommendation of Magirr (2021).

## Value

A list of parameters of the Magirr and Burman weighting function

## Details

Magirr and Burman (2019) proposed a weighted logrank test to have better
power than the logrank test when the treatment effect is delayed, but to
still maintain good power under a proportional hazards assumption. In
Magirr (2021), (the equivalent of) a maximum weight was proposed as
opposed to a fixed time duration over which weights would increase. The
weights for some early interval specified by the user are the inverse of
the combined treatment group empirical survival distribution; see
details. After this initial period, weights are constant at the maximum
of the previous weights. Another advantage of the test is that under
strong null hypothesis that the underlying survival in the control group
is greater than or equal to underlying survival in the experimental
group, Type I error is controlled as the specified level.

We define \\t^\*\\ to be the input variable `delay`. This specifies an
initial period during which weights increase. We also set a maximum
weight \\w\_{\max}\\. To define specific weights, we let \\S(t)\\ denote
the Kaplan-Meier survival estimate at time \\t\\ for the combined data
(control plus experimental treatment groups). The weight at time \\t\\
is then defined as \$\$w(t)=\min(w\_{\max}, S(\min(t, t^\*))^{-1}).\$\$

## References

Magirr, Dominic, and Carl‐Fredrik Burman. 2019. "Modestly weighted
logrank tests." *Statistics in Medicine* 38 (20): 3782–3790.

Magirr, Dominic. 2021. "Non‐proportional hazards in immuno‐oncology: Is
an old perspective needed?" *Pharmaceutical Statistics* 20 (3): 512–527.

## Examples

``` r
sim_pw_surv(n = 200) |>
  cut_data_by_event(100) |>
  wlr(weight = mb(delay = 8, w_max = Inf))
#> $method
#> [1] "WLR"
#> 
#> $parameter
#> [1] "MB(delay = 8, max_weight = Inf)"
#> 
#> $estimate
#> [1] -12.94058
#> 
#> $se
#> [1] 6.689765
#> 
#> $z
#> [1] 1.934385
#> 
#> $info
#> [1] 44.70922
#> 
#> $info0
#> [1] 45.36046
#> 
```
