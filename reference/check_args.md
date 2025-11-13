# Check argument types, length, or dimension

Check argument types, length, or dimension

## Usage

``` r
check_args(arg, type, length = NULL, dim = NULL)
```

## Arguments

- arg:

  An argument to be checked.

- type:

  A character vector of candidate argument type.

- length:

  A numeric value of argument length or `NULL`.

- dim:

  A numeric vector of argument dimension or `NULL`.

## Value

Check failure detailed error message.

## Details

If `type`, `length` or `dim` is `NULL`, the corresponding check will not
be executed.

## Specification

The contents of this section are shown in PDF user manual only.

## Examples

``` r
if (FALSE) { # \dontrun{
tbl <- as.data.frame(matrix(1:9, nrow = 3))
simtrial:::check_args(arg = tbl, type = c("data.frame"))

vec <- c("a", "b", "c")
simtrial:::check_args(arg = vec, type = c("character"), length = 3)
} # }
```
