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

#' @import tibble
#' @import dplyr
NULL

#' The Piecewise Exponential Distribution
#'
#' The piecewise exponential distribution allows a simple method to specify a distribtuion
#' where the hazard rate changes over time. It is likely to be useful for conditions where
#' failure rates change, but also for simulations where there may be a delayed treatment
#' effect or a treatment effect that that is otherwise changing (e.g., decreasing) over time.
#' \code{rpwexp()} is to support simulation of both the Lachin and Foulkes (1986) sample size
#' method for (fixed trial duration) as well as the Kim and Tsiatis(1990) method
#' (fixed enrollment rates and either fixed enrollment duration or fixed minimum follow-up);
#' see \code{\link[gsDesign:nSurv]{gsDesign}}.
#'
#' Using the \code{cumulative=TRUE} option, enrollment times that piecewise constant over
#' time can be generated.
#'
#' @param n Number of observations to be generated.
#' @param failRates A tibble containing \code{duration} and \code{rate} variables.
#' \code{rate} specifies failure rates during the corresponding interval duration
#' specified in \code{duration}. The final interval is extended to be infinite
#' to ensure all observations are generated.
#' @examples
#' # get 10k piecewise exponential failure times
#' # failure rates are 1 for time 0-.5, 3 for time .5 - 1 and 10 for >1.
#' # intervals specifies duration of each failure rate interval
#' # with the final interval running to infinity
#' x <- rpwexp(10000, failRates=tibble::tibble(rate = c(1, 3, 10), duration = c(.5,.5,1)))
#' plot(sort(x),(10000:1)/10001,log="y", main="PW Exponential simulated survival curve",
#' xlab="Time",ylab="P{Survival}")
#' # exponential failure times
#' x <- rpwexp(10000, failRates=tibble::tibble(rate = 5, duration=1))
#'
#' plot(sort(x),(10000:1)/10001,log="y", main="Exponential simulated survival curve",
#'      xlab="Time",ylab="P{Survival}")
#'
#' @export
rpwexp <- function(n=100,
                   failRates=tibble(duration=c(1,1),rate=c(10,20))
                  ){
  n_rates <- nrow(failRates)
  if (n_rates == 1){
    # set failure time to Inf if 0 failure rate
    if(failRates$rate == 0) times = rep(Inf, n)
    # generate exponential failure time if non-0 failure rate
    else times = stats::rexp(n,failRates$rate)
  }else{
    starttime <- 0 # start of first failure rate interval
    finish <- cumsum(failRates$duration) # ends of failure rate interval
    times <- rep(0,n) # initiate vector for failure times
    indx <- rep(TRUE,n) # index for event times not yet reached
    for(i in 1:n_rates){
      nindx <- sum(indx) # number of event times left to generate
      if (nindx==0) break # stop if finished
      # set failure time to Inf for inveral i if 0 fail rate
      if (failRates$rate[i] == 0) times[indx] = starttime + rep(Inf, nindx)
      # generate exponential failure time for interval i if non-0 faiurel rate
      else times[indx] <- starttime + stats::rexp(nindx,failRates$rate[i])
      if (i < n_rates){ # skip this for last interval as all remaining times are generated there
        starttime <- finish[i]
        indx <- (times > finish[i]) # update index of event times not yet reached
      }
    }
  }
  times
}
