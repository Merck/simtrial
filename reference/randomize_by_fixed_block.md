# Permuted fixed block randomization

Fixed block randomization. The `block` input should repeat each
treatment code the number of times it is to be included within each
block. The final block will be a partial block if `n` is not an exact
multiple of the block length.

## Usage

``` r
randomize_by_fixed_block(n = 10, block = c(0, 0, 1, 1))
```

## Arguments

- n:

  Sample size to be randomized.

- block:

  Vector of treatments to be included in each block.

## Value

A treatment group sequence (vector) of length `n` with treatments from
`block` permuted within each block having block size equal to the length
of `block`.

## Examples

``` r
library(dplyr)

# Example 1
# 2:1 randomization with block size 3, treatments "A" and "B"
data.frame(x = 1:10) |> mutate(Treatment = randomize_by_fixed_block(block = c("A", "B", "B")))
#>     x Treatment
#> 1   1         B
#> 2   2         B
#> 3   3         A
#> 4   4         A
#> 5   5         B
#> 6   6         B
#> 7   7         B
#> 8   8         A
#> 9   9         B
#> 10 10         A

# Example 2
# Stratified randomization
data.frame(stratum = c(rep("A", 10), rep("B", 10))) |>
  group_by(stratum) |>
  mutate(Treatment = randomize_by_fixed_block())
#> # A tibble: 20 × 2
#> # Groups:   stratum [2]
#>    stratum Treatment
#>    <chr>       <dbl>
#>  1 A               1
#>  2 A               1
#>  3 A               0
#>  4 A               0
#>  5 A               1
#>  6 A               0
#>  7 A               0
#>  8 A               1
#>  9 A               1
#> 10 A               0
#> 11 B               1
#> 12 B               0
#> 13 B               1
#> 14 B               0
#> 15 B               0
#> 16 B               1
#> 17 B               0
#> 18 B               1
#> 19 B               1
#> 20 B               0
```
