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

#' Conversion of enrollment and failure rates from `sim_fixed_n()` to
#' `sim_pw_surv()` format
#'
#' `simfix2simpwsurv()` converts failure rates and dropout rates entered in
#' the simpler format for [sim_fixed_n()] to that used for [sim_pw_surv()].
#' The `fail_rate` argument for [sim_fixed_n()] requires enrollment rates,
#' failure rates hazard ratios and dropout rates by stratum for a 2-arm trial,
#' [sim_pw_surv()] is in a more flexible but less obvious but more flexible
#' format. Since [sim_fixed_n()] automatically analyzes data and [sim_pw_surv()]
#' just produces a simulation dataset, the latter provides additional options
#' to analyze or otherwise evaluate individual simulations in ways that
#' [sim_fixed_n()] does not.
#'
#' @param fail_rate Piecewise constant control group failure rates,
#'   hazard ratio for experimental vs. control,
#'   and dropout rates by stratum and time period.
#'
#' @return A list of two data frame components formatted for
#'   [sim_pw_surv()]: `fail_rate` and `dropout_rate`.
#'
#' @importFrom tibble tibble
#' @importFrom data.table ":=" .N as.data.table setDF
#'
#' @export
#'
#' @examples
#' library(tidyr)
#' library(dplyr)
#' library(tibble)
#'
#' # Example 1
#' # Convert standard input
#' simfix2simpwsurv()
#'
#' # Stratified example
#' fail_rate <- tibble(
#'   stratum = c(rep("Low", 3), rep("High", 3)),
#'   duration = rep(c(4, 10, 100), 2),
#'   fail_rate = c(
#'     .04, .1, .06,
#'     .08, .16, .12
#'   ),
#'   hr = c(
#'     1.5, .5, 2 / 3,
#'     2, 10 / 16, 10 / 12
#'   ),
#'   dropout_rate = .01
#' )
#'
#' x <- simfix2simpwsurv(fail_rate)
#'
#' # Do a single simulation with the above rates
#' # Enroll 300 patients in ~12 months at constant rate
#' sim <- sim_pw_surv(
#'   n = 300,
#'   stratum = tibble(stratum = c("Low", "High"), p = c(.6, .4)),
#'   enroll_rate = tibble(duration = 12, rate = 300 / 12),
#'   fail_rate = x$fail_rate,
#'   dropout_rate = x$dropout_rate
#' )
#'
#' # Cut after 200 events and do a stratified logrank test
#' dat <- sim %>%
#'   cut_data_by_event(200) %>% # cut data
#'   counting_process(arm = "experimental") %>% # convert format for tenFH
#'   wlr(rho_gamma = tibble(rho = 0, gamma = 0)) # stratified logrank
simfix2simpwsurv <- function(
    # Failure rates as in sim_fixed_n()
    fail_rate = tibble(
      stratum = "All",
      duration = c(3, 100),
      fail_rate = log(2) / c(9, 18),
      hr = c(.9, .6),
      dropout_rate = rep(.001, 2)
    )) {
  # Put failure rates into sim_pw_surv format
  fr_control <- as.data.table(fail_rate)
  fr_control[, `:=`(
    treatment = "control",
    rate = fail_rate,
    period = seq_len(.N)
  ), by = "stratum"]

  fr_experimental <- as.data.table(fail_rate)
  fr_experimental[, `:=`(
    treatment = "experimental",
    rate = fail_rate * hr,
    period = seq_len(.N)
  ), by = "stratum"]

  fr <- rbind(fr_control, fr_experimental)
  fr <- fr[, c("stratum", "period", "treatment", "duration", "rate")]

  # Put dropout rates into sim_pw_surv format
  dr_control <- as.data.table(fail_rate)
  dr_control[, `:=`(
    treatment = "control",
    rate = dropout_rate,
    period = seq_len(.N)
  ), by = "stratum"]

  dr_experimental <- dr_control
  dr_experimental$treatment <- "experimental"

  dr <- rbind(dr_control, dr_experimental)
  dr <- dr[, c("stratum", "period", "treatment", "duration", "rate")]

  list(fail_rate = setDF(fr), dropout_rate = setDF(dr))
}
