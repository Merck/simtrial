# simtrial

<!-- badges: start -->
[![R build status](https://github.com/Merck/simtrial/workflows/R-CMD-check/badge.svg)](https://github.com/Merck/simtrial/actions)
[![Codecov test coverage](https://codecov.io/gh/Merck/simtrial/branch/main/graph/badge.svg)](https://codecov.io/gh/Merck/simtrial?branch=main)
<!-- badges: end -->

## Installation

You can install from GitHub:

```r
remotes::install_github("Merck/simtrial")
```

## Overview

`simtrial` is a small package built initially to focus on evaluating weighted logrank tests and combination tests based on such tests. The intent is to use `tidyverse` (data wrangling) programming procedures and to have a package that is easy to qualify for use in a regulated environment.

Initial areas of focus are:

- Generating time-to-event data for stratified trials using piecewise constant enrollment and piecewise exponential failure rates. Both proportional and non-proportional hazards are supported.
Under proportional hazards, the assumptions are along the lines of those used by Lachin and Foulkes as implemented in the `gsDesign` for deriving group sequential designs.
- Setting up data cutoffs for (interim and final) analyses.
- Support for weighted logrank tests with arbitrary weighting schemes, specifically supporting the Fleming-Harrington set of tests, including the logrank test.

## Future developments

Expectations for future development include:

- Provide a test suite to document that the package is fit for use in a regulatory environment.
- Further examples.
