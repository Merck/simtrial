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

#' The piecewise exponential distribution
#'
#' The piecewise exponential distribution allows a simple method to specify
#' a distribution where the hazard rate changes over time.
#' It is likely to be useful for conditions where failure rates change,
#' but also for simulations where there may be a delayed treatment effect
#' or a treatment effect that that is otherwise changing
#' (for example, decreasing) over time.
#' `rpwexp()` is to support simulation of both the Lachin and Foulkes (1986)
#' sample size method for (fixed trial duration) as well as the
#' Kim and Tsiatis (1990) method (fixed enrollment rates and either
#' fixed enrollment duration or fixed minimum follow-up);
#' see [gsDesign::nSurv()].
#'
#' Using the `cumulative = TRUE` option, enrollment times that piecewise
#' constant over time can be generated.
#'
#' @param n Number of observations to be generated.
#' @param fail_rate A data frame containing `duration` and `rate` variables.
#'   `rate` specifies failure rates during the corresponding interval duration
#'   specified in `duration`. The final interval is extended to be infinite
#'   to ensure all observations are generated.
#'
#' @return The generated random numbers.
#'
#' @export
#'
#' @examples
#' # Example 1
#' # Exponential failure times
#' x <- rpwexp(
#'   n = 10000,
#'   fail_rate = data.frame(rate = 5, duration = 1)
#' )
#' plot(sort(x), (10000:1) / 10001,
#'   log = "y", main = "Exponential simulated survival curve",
#'   xlab = "Time", ylab = "P{Survival}"
#' )
#'
#' # Example 2
#'
#' # Get 10k piecewise exponential failure times.
#' # Failure rates are 1 for time 0 to 0.5, 3 for time 0.5 to 1, and 10 for > 1.
#' # Intervals specifies duration of each failure rate interval
#' # with the final interval running to infinity.
#' x <- rpwexp(
#'   n = 1e4,
#'   fail_rate = data.frame(rate = c(1, 3, 10), duration = c(.5, .5, 1))
#' )
#' plot(sort(x), (1e4:1) / 10001,
#'   log = "y", main = "PW Exponential simulated survival curve",
#'   xlab = "Time", ylab = "P{Survival}"
#' )
rpwexp <- function(
    n = 100,
    fail_rate = data.frame(duration = c(1, 1), rate = c(10, 20))) {
  rpwexp_inverse_cdf_cpp(n = n, fail_rate = fail_rate)
}
