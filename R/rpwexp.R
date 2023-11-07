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
#' Kim and Tsiatis(1990) method (fixed enrollment rates and either
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
  n_rate <- nrow(fail_rate)

  if (n_rate == 1) {
    # Set failure time to Inf if 0 failure rate
    if (fail_rate$rate == 0) {
      ans <- rep(Inf, n)
    } else {
      # Generate exponential failure time if non-0 failure rate
      ans <- stats::rexp(n, fail_rate$rate)
    }
  } else {
    # Start of first failure rate interval
    start_time <- 0
    # Ends of failure rate interval
    end_time <- cumsum(fail_rate$duration)
    # Initiate vector for failure times
    ans <- rep(0, n)
    # Index for event times not yet reached
    indx <- rep(TRUE, n)

    for (i in 1:n_rate) {
      # Number of event times left to generate
      n_event_left <- sum(indx)

      # Stop if event is arrived
      if (n_event_left == 0) {
        break
      }

      # Set failure time to Inf for interval i if 0 fail rate
      if (fail_rate$rate[i] == 0) {
        ans[indx] <- start_time + rep(Inf, n_event_left)
      } else {
        # Generate exponential failure time for interval i if non-0 failure rate
        ans[indx] <- start_time + stats::rexp(n_event_left, fail_rate$rate[i])
      }

      # Skip this for last interval as all remaining times are generated there
      if (i < n_rate) {
        start_time <- end_time[i]
        # Update index of event times not yet reached
        indx <- (ans > end_time[i])
      }
    }
  }

  ans
}
