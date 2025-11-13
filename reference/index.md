# Package index

## Distributions

- [`rpwexp()`](https://merck.github.io/simtrial/reference/rpwexp.md) :
  The piecewise exponential distribution
- [`rpwexp_enroll()`](https://merck.github.io/simtrial/reference/rpwexp_enroll.md)
  : Generate piecewise exponential enrollment
- [`fit_pwexp()`](https://merck.github.io/simtrial/reference/fit_pwexp.md)
  : Piecewise exponential survival estimation

## Simulate data under the piecewise model

- [`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md)
  : Simulate a stratified time-to-event outcome randomized trial

- [`sim_fixed_n()`](https://merck.github.io/simtrial/reference/sim_fixed_n.md)
  : Simulation of fixed sample size design for time-to-event endpoint

- [`sim_gs_n()`](https://merck.github.io/simtrial/reference/sim_gs_n.md)
  : Simulate group sequential designs with fixed sample size

- [`to_sim_pw_surv()`](https://merck.github.io/simtrial/reference/to_sim_pw_surv.md)
  :

  Convert enrollment and failure rates from
  [`sim_fixed_n()`](https://merck.github.io/simtrial/reference/sim_fixed_n.md)
  to
  [`sim_pw_surv()`](https://merck.github.io/simtrial/reference/sim_pw_surv.md)
  format

## Cut data

- [`cut_data_by_date()`](https://merck.github.io/simtrial/reference/cut_data_by_date.md)
  : Cut a dataset for analysis at a specified date
- [`cut_data_by_event()`](https://merck.github.io/simtrial/reference/cut_data_by_event.md)
  : Cut a dataset for analysis at a specified event count
- [`get_cut_date_by_event()`](https://merck.github.io/simtrial/reference/get_cut_date_by_event.md)
  : Get date at which an event count is reached
- [`get_analysis_date()`](https://merck.github.io/simtrial/reference/get_analysis_date.md)
  : Derive analysis date for interim/final analysis given multiple
  conditions
- [`create_cut()`](https://merck.github.io/simtrial/reference/create_cut.md)
  : Create a cutting function

## Compute p-values/test statistics

- [`counting_process()`](https://merck.github.io/simtrial/reference/counting_process.md)
  : Process survival data into counting process format
- [`rmst()`](https://merck.github.io/simtrial/reference/rmst.md) : RMST
  difference of 2 arms
- [`milestone()`](https://merck.github.io/simtrial/reference/milestone.md)
  : Milestone test for two survival curves
- [`wlr()`](https://merck.github.io/simtrial/reference/wlr.md) :
  Weighted logrank test
- [`maxcombo()`](https://merck.github.io/simtrial/reference/maxcombo.md)
  : MaxCombo test
- [`create_test()`](https://merck.github.io/simtrial/reference/create_test.md)
  : Create a cutting test function
- [`multitest()`](https://merck.github.io/simtrial/reference/multitest.md)
  : Perform multiple tests on trial data cutting

## Summarize simulations

- [`summary(`*`<simtrial_gs_wlr>`*`)`](https://merck.github.io/simtrial/reference/summary.md)
  : Summary of group sequential simulations.
- [`as_gt()`](https://merck.github.io/simtrial/reference/as_gt.md) :
  Convert summary table to a gt object

## Randomization algorithms

- [`randomize_by_fixed_block()`](https://merck.github.io/simtrial/reference/randomize_by_fixed_block.md)
  : Permuted fixed block randomization

## Weight functions for WLR test

- [`mb()`](https://merck.github.io/simtrial/reference/mb.md) : Magirr
  and Burman weighting function
- [`fh()`](https://merck.github.io/simtrial/reference/fh.md) :
  Fleming-Harrington weighting function
- [`early_zero()`](https://merck.github.io/simtrial/reference/early_zero.md)
  : Zero early weighting function

## Example datasets

- [`ex1_delayed_effect`](https://merck.github.io/simtrial/reference/ex1_delayed_effect.md)
  : Time-to-event data example 1 for non-proportional hazards working
  group
- [`ex2_delayed_effect`](https://merck.github.io/simtrial/reference/ex2_delayed_effect.md)
  : Time-to-event data example 2 for non-proportional hazards working
  group
- [`ex3_cure_with_ph`](https://merck.github.io/simtrial/reference/ex3_cure_with_ph.md)
  : Time-to-event data example 3 for non-proportional hazards working
  group
- [`ex4_belly`](https://merck.github.io/simtrial/reference/ex4_belly.md)
  : Time-to-event data example 4 for non-proportional hazards working
  group
- [`ex5_widening`](https://merck.github.io/simtrial/reference/ex5_widening.md)
  : Time-to-event data example 5 for non-proportional hazards working
  group
- [`ex6_crossing`](https://merck.github.io/simtrial/reference/ex6_crossing.md)
  : Time-to-event data example 6 for non-proportional hazards working
  group
- [`mb_delayed_effect`](https://merck.github.io/simtrial/reference/mb_delayed_effect.md)
  : Simulated survival dataset with delayed treatment effect
