# simtrial

<!-- badges: start -->
[![R-CMD-check](https://github.com/Merck/simtrial/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Merck/simtrial/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/Merck/simtrial/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Merck/simtrial?branch=main)
![CRAN status](https://www.r-pkg.org/badges/version/simtrial)
<!-- badges: end -->

simtrial is a fast and extensible clinical trial simulation framework
for time-to-event endpoints.

## Installation

You can install from GitHub:

```r
remotes::install_github("Merck/simtrial")
```

## Overview

simtrial is an R package built initially to focus on evaluating weighted
logrank tests and combination tests based on such tests.
simtrial is designed with a core philosophy of basing most computations on
efficient table transformations and to have a package that is easy to qualify
for use in regulated environments.
It utilizes the blazingly fast data.table for tabular data processing,
enhanced by C++ implementations to ensure optimal performance.

Initial areas of focus are:

- Generating time-to-event data for stratified trials using piecewise constant
  enrollment and piecewise exponential failure rates.
  Both proportional and non-proportional hazards are supported.
  Under proportional hazards, the assumptions are along the lines of those
  used by Lachin and Foulkes as implemented in
  [gsDesign](https://keaven.github.io/gsDesign/) for deriving
  group sequential designs.
- Setting up data cutoffs for (interim and final) analyses.
- Support for weighted logrank tests with arbitrary weighting schemes,
  specifically supporting the Fleming-Harrington set of tests,
  including the logrank test.

## Future developments

Expectations for future development include:

- Provide a test suite to document that the package is fit for use in a
  regulatory environment.
- Further examples.
