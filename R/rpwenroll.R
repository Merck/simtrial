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
#' @import tidyr
NULL

#' Generate Piecewise Exponential Enrollment
#'
#' With piecewise exponential enrollment rate generation any enrollment rate distribution can be easily approximated.
#' \code{rpwenroll()} is to support simulation of both the Lachin and Foulkes (1986) sample size
#' method for (fixed trial duration) as well as the Kim and Tsiatis(1990) method
#' (fixed enrollment rates and either fixed enrollment duration or fixed minimum follow-up);
#' see \code{\link[gsDesign:nSurv]{gsDesign}}.
#'
#' @param n Number of observations.
#' Default of \code{NULL} yields random enrollment size.
#' @param enrollRates A tibble containing period duration (\code{duration}) and enrollment rate (\code{rate})
#' for specified enrollment periods.
#' If necessary, last period will be extended to ensure enrollment of specified \code{n}.
#' @return A vector of random enrollment times.
#' @examples
#' # piecewise uniform (piecewise exponential inter-arrival times) for 10k patients enrollment
#' # enrollment rates of 5 for time 0-100, 15 for 100-300, and 30 thereafter
#' x <- rpwenroll(n=10000, enrollRates=tibble::tibble(rate = c(5, 15, 30), duration = c(100,200,100)))
#' plot(x,1:10000,
#'      main="Piecewise uniform enrollment simulation",xlab="Time",
#'      ylab="Enrollment")
#' # exponential enrollment
#' x <- rpwenroll(10000, enrollRates=tibble::tibble(rate = .03, duration = 1))
#' plot(x,1:10000,main="Simulated exponential inter-arrival times",
#'      xlab="Time",ylab="Enrollment")
#'
#' @export
rpwenroll <- function(n = NULL,
                      enrollRates = tibble::tibble(duration=c(1,2),rate=c(2,5))
){# take care of the simple case first
  if(nrow(enrollRates)==1) {
    # stop with error message if only 1 enrollment period and the enrollment rate is less or equal with 0
    if (enrollRates$rate <=0) stop("Please specify > 0 enrollment rate, otherwise enrollment cannot finish.")
    # otherwise, return inter-arrival exponential times
    else return(cumsum(stats::rexp(n=n,rate=enrollRates$rate)))
  }
  y <- enrollRates %>%
          dplyr::mutate(period=dplyr::row_number(),
                        finish=cumsum(duration),
                        lambda=duration*rate,
                        origin=dplyr::lag(finish,default=0)) %>%
          group_by(period) %>%
          dplyr::mutate(N=stats::rpois(n=1,lambda=lambda))
  # deal with extreme cases where none randomized in fixed intervals
  if (sum(y$N)==0){
    if (is.null(n)) return(NULL)
    # stop with error message if enrollment has not finished but enrollment rate for last period is less or equal with 0
    if (dplyr::last(enrollRates$rate) <=0) stop("Please specify > 0 enrollment rate for the last period; otherwise enrollment cannot finish.")
    # otherwise, return inter-arrival exponential times
    else return(cumsum(stats::rexp(n,rate=dplyr::last(enrollRates$rate)))+dplyr::last(y$finish))
  }
  # generate sorted uniform observations for Poisson count for each interval
  z <- tidyr::expand(y,enrollTime=sort(stats::runif(n=N,min=origin,max=finish)))
  if (is.null(n)) return(z$enrollTime) # if n not specified, return generated times
  if (nrow(z) >= n) return(z$enrollTime[1:n]) # if n already achieved, return first n observations
  # after specified finite intervals, add required additional observations with
  # exponential inter-arrival times
  nadd <- n-nrow(z)
  # stop with error message if enrollment has not finished but enrollment rate for last period is less or equal with 0
  if (dplyr::last(enrollRates$rate) <=0) stop("Please specify > 0 enrollment rate for the last period; otherwise enrollment cannot finish.")
  # Otherwise, return inter-arrival exponential times
  else return(c(z$enrollTime, cumsum(stats::rexp(nadd,rate=dplyr::last(enrollRates$rate)))+dplyr::last(y$finish)))
}
