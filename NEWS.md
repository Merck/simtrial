# simtrial 0.3.1

This release introduces significant changes to the API, improves simulation
performance substantially, and adds new features and documentation.

## Significant user-visible changes

- Complete overhaul of the API. Function and argument names now use
  snake case for consistency and readability. See the
  [function reference](https://merck.github.io/simtrial/reference/)
  for the updated naming scheme. Detailed change history is available in the
  [merged pull requests](https://github.com/Merck/simtrial/pulls?q=is%3Apr+is%3Aclosed)
  (thanks, @LittleBeannie, @lili-ling-msd, and @XintongLi2023).
- Dataset names updated to snake case (thanks, @nanxstats, #164).
- The base pipe operator is now used throughout the package.
  The magrittr pipe is no longer re-exported (thanks, @nanxstats, #146).

## Improvements

- Rewritten table backend for simtrial functions using data.table,
  achieving a 3x to 5x speedup compared to the previous implementation
  (thanks, @jdblischak, #111).
- `sim_fixed_n()` now utilizes the `%dofuture%` operator for parallelization,
  enhancing flexibility and reproducibility (thanks, @cmansch, #110).
- `rpwexp()` adopts the inverse CDF method for random number generation,
  with the naive methods now as internal functions
  (thanks, @jianxiaoyang, #15 and #174).
- `sim_fixed_n()` is optimized to skip Breslow's method in the absence of ties
  (thanks, @jdblischak, #130).
- The internal function for computing Z statistics in Fleming-Harrington
  weighted logrank tests is now named `wlr_z_stat()` (thanks, @elong0527, #105).

## New features

- `early_zero_weight()` is added as a weighting function for early data removal
  (thanks, @LittleBeannie, #123).
- `get_analysis_date()` is added to calculate interim/final analysis dates
  under various conditions (thanks, @LittleBeannie, #122).

## Documentation

- New `vignette("workflow")` providing an overview of data manipulations
  involved in TTE simulations (thanks, @keaven, #99).
- New `vignette("parallel")` demonstrating the parallelization workflow and
  coding best practices (thanks, @cmansch, #113 and #134).

## Miscellaneous

- Added a hex sticker logo with a generative art design for the package
  (thanks, @keaven, #158).

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
