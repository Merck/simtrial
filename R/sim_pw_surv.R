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

#' Simulate a stratified time-to-event outcome randomized trial
#'
#' `sim_pw_surv()` enables simulation of a clinical trial with
#' essentially arbitrary patterns of enrollment, failure rates and censoring.
#' The piecewise exponential distribution allows a simple method to specify
#' a distribution and enrollment pattern where the enrollment, failure,
#' and dropout rate changes over time.
#' While the main purpose may be to generate a trial that can be analyzed
#' at a single point in time or using group sequential methods,
#' the routine can also be used to simulate an adaptive trial design.
#' Enrollment, failure, and dropout rates are specified by treatment group,
#' stratum and time period.
#' Fixed block randomization is used; blocks must include treatments provided
#' in failure and dropout specification.
#' Default arguments are set up to allow very simple implementation of
#' a non-proportional hazards assumption for an unstratified design.
#'
#' @param n Number of observations.
#'   If length(n) > 1, the length is taken to be the number required.
#' @param stratum A tibble with stratum specified in `stratum`,
#'   probability (incidence) of each stratum in `p`.
#' @param block Vector of treatments to be included in each block.
#' @param enroll_rate Enrollment rates; see details and examples.
#' @param fail_rate Failure rates; see details and examples;
#'   note that treatments need to be the same as input in block.
#' @param dropout_rate Dropout rates; see details and examples;
#'   note that treatments need to be the same as input in block.
#'
#' @return A data frame with the following variables for each observation:
#' - `stratum`.
#' - `enroll_time`: Enrollment time for the observation.
#' - `Treatment`: Treatment group; this will be one of the values
#'   in the input `block`.
#' - `fail_time`: Failure time generated using [rpwexp()].
#' - `dropout_time`: Dropout time generated using [rpwexp()].
#' - `cte`: Calendar time of enrollment plot the minimum of
#'   failure time and dropout time.
#' - `fail`: Indicator that `cte` was set using failure time;
#'   i.e., 1 is a failure, 0 is a dropout.
#'
#' @importFrom data.table ":=" .N data.table setDF setDT setorderv
#' @importFrom dplyr group_by mutate
#' @importFrom tibble tibble
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr)
#'
#' # Example 1
#' sim_pw_surv(n = 20)
#'
#' # Example 2
#' # 3:1 randomization
#' sim_pw_surv(
#'   n = 20,
#'   block = c(rep("experimental", 3), "control")
#' )
#'
#' # Example 3
#' # Simulate 2 stratum; will use defaults for blocking and enrollRates
#' sim_pw_surv(
#'   n = 20,
#'   # 2 stratum,30% and 70% prevalence
#'   stratum = tibble(stratum = c("Low", "High"), p = c(.3, .7)),
#'   fail_rate = tibble(
#'     stratum = c(rep("Low", 4), rep("High", 4)),
#'     period = rep(1:2, 4),
#'     treatment = rep(c(
#'       rep("control", 2),
#'       rep("experimental", 2)
#'     ), 2),
#'     duration = rep(c(3, 1), 4),
#'     rate = c(.03, .05, .03, .03, .05, .08, .07, .04)
#'   ),
#'   dropout_rate = tibble(
#'     stratum = c(rep("Low", 2), rep("High", 2)),
#'     period = rep(1, 4),
#'     treatment = rep(c("control", "experimental"), 2),
#'     duration = rep(1, 4),
#'     rate = rep(.001, 4)
#'   )
#' )
#' # Example 4
#' # If you want a more rectangular entry for a tibble
#' fail_rate <- bind_rows(
#'   tibble(stratum = "Low", period = 1, treatment = "control", duration = 3, rate = .03),
#'   tibble(stratum = "Low", period = 1, treatment = "experimental", duration = 3, rate = .03),
#'   tibble(stratum = "Low", period = 2, treatment = "experimental", duration = 3, rate = .02),
#'   tibble(stratum = "High", period = 1, treatment = "control", duration = 3, rate = .05),
#'   tibble(stratum = "High", period = 1, treatment = "experimental", duration = 3, rate = .06),
#'   tibble(stratum = "High", period = 2, treatment = "experimental", duration = 3, rate = .03)
#' )
#'
#' dropout_rate <- bind_rows(
#'   tibble(stratum = "Low", period = 1, treatment = "control", duration = 3, rate = .001),
#'   tibble(stratum = "Low", period = 1, treatment = "experimental", duration = 3, rate = .001),
#'   tibble(stratum = "High", period = 1, treatment = "control", duration = 3, rate = .001),
#'   tibble(stratum = "High", period = 1, treatment = "experimental", duration = 3, rate = .001)
#' )
#'
#' sim_pw_surv(
#'   n = 12,
#'   stratum = tibble(stratum = c("Low", "High"), p = c(.3, .7)),
#'   fail_rate = fail_rate,
#'   dropout_rate = dropout_rate
#' )
sim_pw_surv <- function(
    n = 100,
    stratum = tibble(stratum = "All", p = 1),
    block = c(rep("control", 2), rep("experimental", 2)),
    enroll_rate = tibble(rate = 9, duration = 1),
    fail_rate = tibble(
      stratum = rep("All", 4),
      period = rep(1:2, 2),
      treatment = c(rep("control", 2), rep("experimental", 2)),
      duration = rep(c(3, 1), 2),
      rate = log(2) / c(9, 9, 9, 18)
    ),
    dropout_rate = tibble(
      stratum = rep("All", 2),
      period = rep(1, 2),
      treatment = c("control", "experimental"),
      duration = rep(100, 2),
      rate = rep(.001, 2)
    )) {
  # Start table by generating stratum and enrollment times
  x <- data.table(stratum = sample(
    x = stratum$stratum,
    size = n,
    replace = TRUE,
    prob = stratum$p
  ))
  x[, enroll_time := rpw_enroll(n, enroll_rate)]
  # The awkward back and forth ordering is to maintain 1:1 parity with
  # dplyr::group_by() for backwards compatibility. group_by() sorts by the
  # grouping variable and then returns the rows to their original positions.
  # This is mainly for testing for backwards compatibility. Since the
  # treatments are assigned randomly by group, it would still be statistically
  # valid without this ordering
  setorderv(x, "stratum")
  x[, treatment := randomize_by_fixed_block(n = .N, block = block), by = "stratum"]
  setorderv(x, "enroll_time")

  # Generate time to failure and time to dropout
  unique_stratum <- unique(x$stratum)
  unique_treatment <- unique(x$treatment)
  x[, fail_time := 0]
  x[, dropout_time := 0]

  for (sr in unique_stratum) {
    for (tr in unique_treatment) {
      indx <- x$stratum == sr & x$treatment == tr
      x$fail_time[indx] <- rpwexpinvRcpp(
        n = sum(indx),
        fail_rate = fail_rate[fail_rate$stratum == sr & fail_rate$treatment == tr, , drop = FALSE]
      )
      x$dropout_time[indx] <- rpwexpinvRcpp(
        n = sum(indx),
        fail_rate = dropout_rate[dropout_rate$stratum == sr & dropout_rate$treatment == tr, , drop = FALSE]
      )
    }
  }

  # Set calendar time-to-event and failure indicator
  ans <- setDT(x)
  ans[, cte := pmin(dropout_time, fail_time) + enroll_time]
  ans[, fail := (fail_time <= dropout_time) * 1]

  return(setDF(ans))
}
