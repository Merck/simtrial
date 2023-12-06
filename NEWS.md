# simtrial 0.3.1

## Significant user-visible changes

- API: function names, argument names changed.
  We recommend to checking out the
  [function reference page](https://merck.github.io/simtrial/reference/)
  to see the latest naming scheme and search the merged
  [pull requests](https://github.com/Merck/simtrial/pulls?q=is%3Apr+is%3Aclosed)
  to see detailed change history
  (thanks, @LittleBeannie, @lili-ling-msd, @XintongLi2023).
- Dataset names are updated.

## Improvements

- The table backend has been rewritten to use data.table for
  optimal computational performance. On average, 3x to 5x speedup can be
  observed compared to the previous implementation
  (thanks, @jdblischak, #111).
- Use `%dofuture%` to follow modern standard (thanks, @cmansch, #77).
- `rpwexp()` inverse CDF method, naive method implementations unexported
  as internal functions (#174).
  The C++ implementation is from #15 (thanks, @jianxiaoyang).

## New features

- New function `early_zero_weight()` is added as a weighting function
  for early data removal (thanks, @LittleBeannie, #123).
- New function `get_analysis_date()` is added to enable the calculation
  of interim/final analysis dates based on various conditions
  (thanks, @LittleBeannie, #122).

## Documentation

- Add a vignette to demonstrate the parallelization workflow and
  coding best practices (thanks, @cmansch, #113 and #134).

## Bug fixes

# simtrial 0.2.2

GitHub release in February 2023.

This is the version that enables parallel computation in `simfix()`.

# simtrial 0.2.1

GitHub release in May 2022.

This version supports the _Biometrical Journal_ paper
"A unified framework for weighted parametric group sequential design (WPGSD)"
by Keaven M. Anderson, Zifang Guo, Jing Zhao, and Linda Z. Sun.

# simtrial 0.2.0

Internal development release in August 2020.

- Updated vignettes and website.
- Prepared for Regulatory/Industry training session in September.

# simtrial 0.1.7.9004

Internal development release in February 2020.

- Added `wMB()` to compute Magirr-Burman weights.
- Added vignette to demonstrate working with different weighting schemes.
- Replaced `Depends` with `Imports` in `DESCRIPTION`.

# simtrial 0.1.7.9003

Internal development release in November 2019.

- Incorporated new functions to simplify use (`simfix()`, `simfix2simPWSurv()`, `pMaxCombo()`).
- Removed `hgraph()` with intent to put it into a release of gsDesign.
- Limited to 2 essential vignettes.
- Added continuous integration/continuous deployment (YAML) and pkgdown for website generation.
- Limited dependencies to those that are essential; this removed some convenience functions not related to core package functionality.
