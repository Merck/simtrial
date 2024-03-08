#  Copyright (c) 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
#  All rights reserved.
#
#  This file is part of the simtrial program.
#
#  simtrial is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Simulate group sequential designs with fixed sample size
#'
#' WARNING: This experimental function is a work-in-progress. The function
#' arguments will change as we add additional features.
#'
#' @inheritParams sim_fixed_n
#' @param test A test function such as [wlr()],
#'   [maxcombo()], or [rmst()]. The simulated data set is
#'   passed as the first positional argument to the test function provided.
#' @param cutting A list of cutting functions created by [create_cutting()],
#'   see examples.
#' @param seed Random seed.
#' @param ... Arguments passed to the test function provided by the argument
#'   `test`.
#'
#' @return A data frame summarizing the simulation ID, analysis date,
#'   z statistics or p-values.
#'
#' @export
#'
#' @examples
#' library(gsDesign2)
#'
#' # Parameters for enrollment
#' enroll_rampup_duration <- 4 # Duration for enrollment ramp up
#' enroll_duration <- 16 # Total enrollment duration
#' enroll_rate <- define_enroll_rate(
#'   duration = c(
#'     enroll_rampup_duration,
#'     enroll_duration - enroll_rampup_duration
#'   ),
#'   rate = c(10, 30)
#' )
#'
#' # Parameters for treatment effect
#' delay_effect_duration <- 3 # Delay treatment effect in months
#' median_ctrl <- 9 # Survival median of the control arm
#' median_exp <- c(9, 14) # Survival median of the experimental arm
#' dropout_rate <- 0.001
#' fail_rate <- define_fail_rate(
#'   duration = c(delay_effect_duration, 100),
#'   fail_rate = log(2) / median_ctrl,
#'   hr = median_ctrl / median_exp,
#'   dropout_rate = dropout_rate
#' )
#'
#' # Other related parameters
#' alpha <- 0.025 # Type I error
#' beta <- 0.1 # Type II error
#' ratio <- 1 # Randomization ratio (experimental:control)
#'
#' # Define cuttings of 2 IAs and 1 FA
#' # IA1
#' # The 1st interim analysis will occur at the later of the following 3 conditions:
#' # - At least 20 months have passed since the start of the study.
#' # - At least 100 events have occurred.
#' # - At least 20 months have elapsed after enrolling 200/400 subjects, with a
#' #   minimum of 20 months follow-up.
#' # However, if events accumulation is slow, we will wait for a maximum of 24 months.
#' ia1 <- create_cutting(
#'   planned_calendar_time = 20,
#'   target_event_overall = 100,
#'   max_extension_for_target_event = 24,
#'   min_n_overall = 200,
#'   min_followup = 20
#' )
#'
#' # IA2
#' # The 2nd interim analysis will occur at the later of the following 3 conditions:
#' # - At least 32 months have passed since the start of the study.
#' # - At least 250 events have occurred.
#' # - At least 10 months after IA1.
#' # However, if events accumulation is slow, we will wait for a maximum of 34 months.
#' ia2 <- create_cutting(
#'   planned_calendar_time = 32,
#'   target_event_overall = 200,
#'   max_extension_for_target_event = 34,
#'   min_time_after_previous_analysis = 10
#' )
#'
#' # FA
#' # The final analysis will occur at the later of the following 2 conditions:
#' # - At least 45 months have passed since the start of the study.
#' # - At least 300 events have occurred.
#' fa <- create_cutting(
#'   planned_calendar_time = 45,
#'   target_event_overall = 350
#' )
#'
#' # Test 1: regular logrank test
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr,
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024,
#'   weight = fh(rho = 0, gamma = 0)
#' )
#'
#' # Test 2: weighted logrank test by FH(0, 0.5)
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr,
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024,
#'   weight = fh(rho = 0, gamma = 0.5)
#' )
#'
#' # Test 3: weighted logrank test by MB(3)
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr,
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024,
#'   weight = mb(delay = 3)
#' )
#'
#' # Test 4: weighted logrank test by early zero (6)
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr,
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024,
#'   weight = early_zero(6)
#' )
#'
#' # Test 5: RMST
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = rmst,
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024,
#'   tau = 20
#' )
#'
#' # Test 6: Milestone
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = milestone,
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024,
#'   ms_time = 10
#' )
#'
#' # Test 7: MaxCombo (WLR-FH(0,0) + WLR-FH(0, 0.5))
#' # for all analyses
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = maxcombo,
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024,
#'   rho = c(0, 0),
#'   gamma = c(0, 0.5)
#' )
#'
#' # Test 8: MaxCombo (WLR-FH(0,0.5) + milestone(10))
#' # for all analyses
#' \dontrun{
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = maxcombo(test1 = wlr, test2 = milestone),
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024,
#'   test1_par = list(weight = fh(rho = 0, gamma = 0.5)),
#'   test2_par = list(ms_time = 10)
#' )
#' }
#'
#' # Test 9: MaxCombo (WLR-FH(0,0) at IAs
#' # and WLR-FH(0,0) + milestone(10) + WLR-MB(4,2) at FA)
#' \dontrun{
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = list(ia1 = wlr, ia2 = wlr, fa = maxcombo),
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024,
#'   test_par = list(
#'     ia1 = list(weight = fh(rho = 0, gamma = 0)),
#'     ia2 = list(weight = fh(rho = 0, gamma = 0)),
#'     ia3 = list(
#'       test1_par = list(weight = fh(rho = 0, gamma = 0)),
#'       test2_par = list(ms_time = 10),
#'       test3_par = list(delay = 4, w_max = 2)
#'     )
#'   )
#' )
#' }
sim_gs_n <- function(
    n_sim = 1000,
    sample_size = 500,
    stratum = data.frame(stratum = "All", p = 1),
    enroll_rate = data.frame(duration = c(2, 2, 10), rate = c(3, 6, 9)),
    fail_rate = data.frame(
      stratum = "All",
      duration = c(3, 100),
      fail_rate = log(2) / c(9, 18),
      hr = c(.9, .6),
      dropout_rate = rep(.001, 2)
    ),
    block = rep(c("experimental", "control"), 2),
    test = wlr,
    cutting = NULL,
    seed = 2024,
    ...) {
  # Input checking
  # TODO

  # Simulate for `n_sim` times
  ans <- NULL
  for (sim_id in seq_len(n_sim)) {
    set.seed(seed + sim_id)
    # Generate data
    simu_data <- sim_pw_surv(
      n = sample_size,
      stratum = stratum,
      block = block,
      enroll_rate = enroll_rate,
      fail_rate = to_sim_pw_surv(fail_rate)$fail_rate,
      dropout_rate = to_sim_pw_surv(fail_rate)$dropout_rate
    )

    # Initialize the cut date of IA(s) and FA
    n_analysis <- length(cutting)
    cut_date <- rep(-100, n_analysis)
    ans_1sim <- NULL

    for (i_analysis in seq_len(n_analysis)) {
      # Get cut date
      cut_date[i_analysis] <- cutting[[i_analysis]](data = simu_data)

      # Cut the data
      simu_data_cut <- simu_data |> cut_data_by_date(cut_date[i_analysis])

      # Test
      ans_1sim_new <- test(simu_data_cut, ...)
      ans_1sim_new$analysis <- i_analysis
      ans_1sim_new$cut_date <- cut_date[i_analysis]
      ans_1sim_new$sim_id <- sim_id
      ans_1sim_new$n <- nrow(simu_data_cut)
      ans_1sim_new$event <- sum(simu_data_cut$event)

      # rbind simulation results for all IA(s) and FA in 1 simulation
      ans_1sim <- rbind(ans_1sim, ans_1sim_new)
    }

    ans <- rbind(ans, ans_1sim)
  }
  return(ans)
}

#' Create a cutting function
#'
#' Create a cutting function for use with \code{\link{sim_gs_n}}
#'
#' @param ... Arguments passed to \code{\link{get_analysis_date}}
#'
#' @return A function that accepts a data frame of simulated trial data and
#'   returns a cut date
#'
#' @export
#'
#' @seealso \code{\link{get_analysis_date}}, \code{\link{sim_gs_n}}
#'
#' @examples
#' # Simulate trial data
#' trial_data <- sim_pw_surv()
#'
#' # Create a cutting function that applies the following 2 conditions:
#' # - At least 45 months have passed since the start of the study
#' # - At least 300 events have occurred
#' cutting <- create_cutting(
#'   planned_calendar_time = 45,
#'   target_event_overall = 350
#' )
#'
#' # Cut the trial data
#' cutting(trial_data)
create_cutting <- function(...) {
  function(data) {
    get_analysis_date(data, ...)
  }
}
