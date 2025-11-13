# simtrial

simtrial is a fast and extensible clinical trial simulation framework
for time-to-event endpoints.

## Installation

The easiest way to get simtrial is to install from CRAN:

``` r
install.packages("simtrial")
```

Alternatively, to use a new feature or get a bug fix, you can install
the development version of simtrial from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("Merck/simtrial")
```

## Overview

simtrial is intended to be a general purpose tool for simulating fixed,
group sequential or adaptive clinical trials. It allows stratified
populations and flexible parameters for generating enrollment, event
times, dropout times. It takes care of bookkeeping to enable easily
going from data generation to creating analysis datasets for evaluation
of standard or innovative designs and testing procedures. For a single
endpoint, it will easily generate trials with multiple arms (e.g., a
single or multiple experimental arms versus a common control) and
multiple study populations (e.g., overall population and biomarker
positive). While tools are built into the package for logrank and
weighted logrank tests, arbitrary testing and estimation procedures are
easily applied. In addition to weighted logrank tests, we support
combinations of weighted logrank tests (e.g., the MaxCombo test). The
package used piecewise constant enrollment, failure and dropout rates as
a simple model able to approximate arbitrary distributions easily. This
model also enables simulating non-proportional hazards assumptions that
are transparent for users to explain to non-statistical collaborators.

simtrial is designed with a core philosophy of basing most computations
on efficient table transformations and to have a package that is easy to
qualify for use in regulated environments. It utilizes the blazingly
fast data.table for tabular data processing, enhanced by C++
implementations to ensure optimal performance. However, it does not
require the user to be a data.table or C++ user.

Initial areas of focus are:

- Generating time-to-event data for stratified trials using piecewise
  constant enrollment and piecewise exponential failure rates. Both
  proportional and non-proportional hazards are supported. Under
  proportional hazards, the assumptions are along the lines of those
  used by Lachin and Foulkes as implemented in
  [gsDesign](https://keaven.github.io/gsDesign/) for deriving group
  sequential designs.
- Setting up data cutoffs for (interim and final) analyses.
- Support for weighted logrank tests with arbitrary weighting schemes,
  specifically supporting the Fleming-Harrington set of tests, including
  the logrank test.

## Future developments

Expectations for future development include:

- Provide a test suite to document that the package is fit for use in a
  regulatory environment.
- Further examples.
