#  Copyright (c) 2022 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
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

#' @import dplyr
#' @importFrom tibble tibble
NULL

#' Conversion of enrollment and failure rates from sim_fixed_n() to sim_pw_surv() format
#'
#' `simfix2simPWSurv()` converts failure rates and dropout rates entered in the simpler
#' format for `sim_fixed_n()` to that used for `simtrial::sim_pw_surv()`.
#' The `failRates` argument for `sim_fixed_n()` requires enrollment rates, failure rates
#' hazard ratios and dropout rates by strata for a 2-arm trial, `simtrial::sim_pw_surv()`
#' is in a more flexible but less obvious but more flexible format.
#' Since `sim_fixed_n()` automatically analyzes data and `simtrial::sim_pw_surv()` just produces
#' a simulation dataset, the latter provides additional options to analyze or otherwise evaluate
#' individual simulations in ways that `sim_fixed_n()` does not.
#' @param failRates Piecewise constant control group failure rates, hazard ratio for experimental vs control,
#'  and dropout rates by stratum and time period.
#' @return A \code{list} of two `tibble` components formatted for `simtrial::sim_pw_surv()`:
#' `failRates` and `dropoutRates`.
#'
#' @examples
#' library(tidyr)
#' library(dplyr)
#' library(tibble)
#'
#' # example 1
#' # Convert standard input
#' simfix2simPWSurv()
#'
#' # Stratified example
#' failRates <- tibble(Stratum = c(rep("Low", 3),rep("High", 3)),
#'                     duration = rep(c(4, 10, 100), 2),
#'                     failRate = c(.04, .1, .06,
#'                                  .08,.16,.12),
#'                     hr = c(1.5, .5, 2/3,
#'                            2, 10/16, 10/12),
#'                     dropoutRate =.01)
#'
#' x <- simfix2simPWSurv(failRates)
#'
#' # Do a single simulation with the above rates
#' # Enroll 300 patients in ~12 months at constant rate
#' sim <- sim_pw_surv(n = 300,
#'                  strata = tibble(Stratum = c("Low","High"), p = c(.6, .4)),
#'                  enroll_rate = tibble(duration = 12, rate = 300 / 12),
#'                  failRates = x$failRates,
#'                  dropoutRates = x$dropoutRates)
#'
#' # Cut after 200 events and do a stratified logrank test
#' dat <- sim %>%
#'   cut_data_by_event(200) %>%              # cut data
#'   counting_process(arm = "Experimental") %>%  # convert format for wlr
#'   wlr(rg = tibble(rho=0,gamma=0))    # stratified logrank
#'
#' @export
#'
simfix2simPWSurv <- function(
  # failure rates as in sim_fixed_n()
  failRates = tibble(Stratum = "All",
                     duration = c(3, 100),
                     failRate = log(2) / c(9, 18),
                     hr = c(.9, .6),
                     dropoutRate = rep(.001, 2))
  ){
  # put failure rates into sim_pw_surv format
  fr <- rbind(failRates %>%
                group_by(Stratum) %>%
                mutate(Treatment = "Control",
                       rate = failRate, period = 1:n()) %>%
                ungroup(),
              failRates %>%
                group_by(Stratum) %>%
                mutate(Treatment = "Experimental",
                       rate = failRate * hr,
                       period = 1:n()) %>%
                ungroup()
              ) %>%
    select("Stratum", "period", "Treatment", "duration", "rate")

  # put dropout rates into sim_pw_surv format
  dr <- failRates %>%
    group_by(Stratum) %>%
    mutate(Treatment = "Control",
           rate = dropoutRate,
           period = 1:n()) %>%
    select("Stratum", "period", "Treatment", "duration", "rate") %>%
    ungroup()

  dr <- rbind(dr,
              dr %>% mutate(Treatment = "Experimental"))

  return(list(failRates = fr, dropoutRates = dr))
}
