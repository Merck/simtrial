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

#' Generate piecewise exponential enrollment
#'
#' With piecewise exponential enrollment rate generation any enrollment rate
#' distribution can be easily approximated.
#' `rpw_enroll()` is to support simulation of both the Lachin and Foulkes (1986)
#' sample size method for (fixed trial duration) as well as the
#' Kim and Tsiatis(1990) method (fixed enrollment rates and either
#' fixed enrollment duration or fixed minimum follow-up);
#' see [gsDesign::nSurv()].
#'
#' @param n Number of observations.
#'   Default of `NULL` yields random enrollment size.
#' @param enroll_rate A tibble containing period duration (`duration`)
#'   and enrollment rate (`rate`). for specified enrollment periods.
#'   If necessary, last period will be extended to ensure enrollment
#'   of specified `n`.
#'
#' @return A vector of random enrollment times.
#'
#' @importFrom data.table ":=" .N as.data.table last
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' # Example 1
#' # Piecewise uniform (piecewise exponential inter-arrival times) for 10k patients enrollment
#' # Enrollment rates of 5 for time 0-100, 15 for 100-300, and 30 thereafter
#' x <- rpw_enroll(
#'   n = 1e5,
#'   enroll_rate = tibble(
#'     rate = c(5, 15, 30),
#'     duration = c(100, 200, 100)
#'   )
#' )
#' plot(x, 1:1e5,
#'   main = "Piecewise uniform enrollment simulation",
#'   xlab = "Time",
#'   ylab = "Enrollment"
#' )
#'
#' # Example 2
#' # Exponential enrollment
#' x <- rpw_enroll(
#'   n = 1e5,
#'   enroll_rate = tibble(rate = .03, duration = 1)
#' )
#' plot(x, 1:1e5,
#'   main = "Simulated exponential inter-arrival times",
#'   xlab = "Time",
#'   ylab = "Enrollment"
#' )
rpw_enroll <- function(
    n = NULL,
    enroll_rate = tibble(duration = c(1, 2), rate = c(2, 5))) {
  # Take care of the simple case first if it is exponential enrollment
  if (nrow(enroll_rate) == 1) {
    # Stop with error message if only 1 enrollment period and the
    # enrollment rate is less or equal with 0
    if (enroll_rate$rate <= 0) {
      stop("rpw_enroll: please specify > 0 enrollment rate, otherwise enrollment cannot finish.")
    }
    # Otherwise, return inter-arrival exponential times
    else {
      ans <- cumsum(stats::rexp(n = n, rate = enroll_rate$rate))
      return(ans)
    }
  }

  # Build `y` summarizes the start/end time, period order, etc.
  y <- as.data.table(enroll_rate)
  y[, period := seq_len(.N)]
  y[, finish := cumsum(duration)]
  y[, lambda := duration * rate]
  y[, origin := c(0, finish[-length(finish)])]
  y[, N := stats::rpois(n = .N, lambda = lambda)]

  # Deal with extreme cases where none randomized in fixed intervals
  if (sum(y$N) == 0) {
    if (is.null(n)) {
      ans <- NULL
      return(ans)
    }

    if (last(enroll_rate$rate) <= 0) {
      # Stop with error message if enrollment has not finished but enrollment rate for last period is less or equal with 0
      stop("rpw_enroll: please specify > 0 enrollment rate for the last period; otherwise enrollment cannot finish.")
    } else {
      # Otherwise, return inter-arrival exponential times
      ans <- cumsum(stats::rexp(n = n, rate = last(enroll_rate$rate))) + last(y$finish)
      return(ans)
    }
  }

  # Generate sorted uniform observations for Poisson count for each interval
  z <- y[, .(enroll_time = sort(stats::runif(n = N, min = origin, max = finish))), by = "period"]

  # If n not specified, return generated times
  if (is.null(n)) {
    ans <- z$enroll_time
    return(ans)
  }

  # If n already achieved, return first n observations
  if (nrow(z) >= n) {
    ans <- z$enroll_time[1:n]
    return(ans)
  }

  # After specified finite intervals, add required additional observations with
  # exponential inter-arrival times
  n_add <- n - nrow(z)
  # Stop with error message if enrollment has not finished but
  # enrollment rate for last period is less or equal with 0
  if (last(enroll_rate$rate) <= 0) {
    stop("rpw_enroll: please specify > 0 enrollment rate for the last period; otherwise enrollment cannot finish.")
  }
  # Otherwise, return inter-arrival exponential times
  else {
    ans <- c(
      z$enroll_time,
      cumsum(stats::rexp(n_add, rate = last(enroll_rate$rate))) + last(y$finish)
    )
    return(ans)
  }
}
