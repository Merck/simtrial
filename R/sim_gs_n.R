#  Copyright (c) 2023 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
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

#' Simulate group sequantial designs with fixed sample size
#' @inheritParams sim_fixed_n
#' @param test a functional call of the test such as \code{wlr()} or \code{maxcombo()}
#' @param cutting a functional call of the cutting for IA(s) and FA, see examples
#' @param seed random seed
#'
#' @return a data frame summaring the simulation ID, analysis date, z statistics or p-values
#' @export
#'
#' @examples
#' library(gsDesign2)
#'
#' # parameters for enrollment
#' enroll_rampup_duration <- 4 # duration for enrollment ramp up
#' enroll_duration <- 16       # total enrollment duration
#' enroll_rate <- define_enroll_rate(duration = c(enroll_rampup_duration,
#'                                                enroll_duration - enroll_rampup_duration),
#'                                   rate = c(10, 30))
#'
#' # parameters for treatment effect
#' delay_effect_duration <- 3  # delay treatment effect in months
#' median_col <- 9             # survival median of the control arm
#' median_exp <- c(9, 14)      # survival median of the experimental arm
#' dropout_rate <- 0.001
#' fail_rate <- define_fail_rate(duration = c(delay_effect_duration, 100),
#'                               hr = median_col / median_exp,
#'                               dropout_rate = dropout_rate)
#'
#' # other related parameters
#' alpha <- 0.025              # type I error
#' beta <- 0.1                 # type II error
#' ratio <- 1                  # randomization ratio (exp:col)
#'
#' # Define cuttings of 2 IAs and 1 FA
#' # IA1
#' # The 1st interim analysis will occur at the later of the following 3 conditions:
#' # - At least 20 months have passed since the start of the study
#' # - At least 100 events have occurred
#' # - At least 20 months have elapsed after enrolling 200/400 subjects, with a minimum of 20 months follow-up
#' # However, if events accumulation is slow, we will wait for a maximum of 24 months.
#' ia1 <- get_analysis_date(data,
#'                          planned_calendar_time = 20,
#'                          target_event_overall = 100,
#'                          max_extension_for_target_event = 24,
#'                          min_n_overall = 200,
#'                          min_followup = 20) |> quote()
#'
#' # IA2
#' # The 2nd interim analysis will occur at the later of the following 3 conditions:
#' # - At least 32 months have passed since the start of the study
#' # - At least 250 events have occurred
#' # - At least 10 months after IA1
#' # However, if events accumulation is slow, we will wait for a maximum of 34 months.
#' ia2 <- get_analysis_date(data,
#'                          planned_calendar_time = 32,
#'                          target_event_overall = 200,
#'                          max_extension_for_target_event = 34,
#'                          min_time_after_previous_analysis = 10) |> quote()
#'
#' # FA
#' # The final analysis will occur at the later of the following 2 conditions:
#' # - At least 45 months have passed since the start of the study
#' # - At least 300 events have occurred
#' fa <- get_analysis_date(data,
#'                        planned_calendar_time = 45,
#'                        target_event_overall = 350) |> quote()
#'
#' # Test 1: regular logrank test
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr(data, weight = fh(rho = 0, gamma = 0)) |> quote(),
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024)
#'
#' # Test 2: weighted logrank test by FH(0, 0.5)
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr(data, weight = fh(rho = 0, gamma = 0.5)) |> quote(),
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024)
#'
#'
#' # Test 3: weighted logrank test by MB(6)
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr(data, weight = mb(delay = 3)) |> quote(),
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024)
#'
#' # Test 4: weighted logrank test by early zero (6)
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr(data, weight = early_zero(6)) |> quote(),
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024)
#'
#' # Test 5: weighted logrank test by RMST
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr(data, weight = rmst(10)) |> quote(),
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024)
#'
#' # Test 6: weighted logrank test by RMST
#' sim_gs_n(
#'   n_sim = 3,
#'   sample_size = 400,
#'   enroll_rate = enroll_rate,
#'   fail_rate = fail_rate,
#'   test = wlr(data, weight = rmst(10)) |> quote(),
#'   cutting = list(ia1 = ia1, ia2 = ia2, fa = fa),
#'   seed = 2024)
sim_gs_n <- function(
    # number of simulations
  n_sim = 1000,
  # sample size
  sample_size = 500,
  # multinomial probability distribution for stratum enrollment
  stratum = data.frame(stratum = "All", p = 1),
  # enrollment rates
  enroll_rate = data.frame(duration = c(2, 2, 10), rate = c(3, 6, 9)),
  # failure rates
  fail_rate = data.frame(
    stratum = "All",
    duration = c(3, 100),
    fail_rate = log(2) / c(9, 18),
    hr = c(.9, .6),
    dropout_rate = rep(.001, 2)
  ),
  # fixed block randomization specification
  block = rep(c("experimental", "control"), 2),
  # default is to to logrank testing
  # but alternative tests (such as rmst, maxcombo) can be specified
  test = wlr(weight = fh(rho = 0, gamma = 0)) |> quote(),
  # cutting for IA(s) and FA
  cutting = NULL,
  # random seed
  seed = 2024
){
  # input checking
  # TODO

  # simulate for n_sim times
  ans <- NULL
  for (sim_id in 1:n_sim) {
    set.seed(seed + sim_id)
    # generate data
    simu_data <- sim_pw_surv(
      n = sample_size,
      stratum = stratum,
      block = block,
      enroll_rate = enroll_rate,
      fail_rate = to_sim_pw_surv(fail_rate)$fail_rate,
      dropout_rate = to_sim_pw_surv(fail_rate)$dropout_rate)

    # initialize the cut date of IA(s) and FA
    n_analysis <- length(cutting)
    cut_date <- rep(-100, n_analysis)
    ans_1sim <- NULL

    for (i_analysis in 1:n_analysis) {

      # get cut date
      if (i_analysis < n_analysis) {
        cut_date[i_analysis] <- cutting[[paste0("ia", i_analysis)]] |> eval(envir = env(data = simu_data))
      } else {
        cut_date[i_analysis] <- cutting[["fa"]] |> eval(envir = env(data = simu_data))
      }

      # cut the data
      simu_data_cut <- simu_data |> cut_data_by_date(cut_date[i_analysis])

      # test
      ans_1sim_new <- eval(test, envir = env(data = simu_data_cut)) |>
        mutate(analysis = i_analysis,
               cut_date = cut_date[i_analysis],
               sim_id = sim_id)

      # rbind simulation results for all IA(s) and FA in 1 simulation
      ans_1sim <- rbind(ans_1sim, ans_1sim_new)
    }

    ans <- rbind(ans, ans_1sim)
  }
  return(ans)
}
