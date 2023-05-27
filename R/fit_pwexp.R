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

#' @import survival
NULL

#' Piecewise exponential survival estimation
#'
#' Computes survival function, density function, -2*log-likelihood based
#' on input dataset and intervals for piecewise constant failure rates.
#' Initial version assumes observations are right censored or events only.
#'
#' @param srv input survival object (see \code{Surv});
#' note that only 0 = censored, 1 = event for \code{Surv}
#' @param intervals Vector containing positive values indicating interval lengths where the
#' exponential rates are assumed.
#' Note that a final infinite interval is added if any events occur after the final interval
#' specified.
#'
#' @return A matrix with rows containing interval length, estimated rate,
#' -2*log-likelihood for each interval.
#'
#' @examples
#' # use default arguments for delayed effect example dataset (Ex1delayedEffect)
#' library(survival)
#'
#' # example 1
#' rateall <- fit_pwexp()
#' rateall
#'
#' # example 2
#' # Estimate by treatment effect
#' rate1 <- with(subset(Ex1delayedEffect, trt == 1), fit_pwexp(Surv(month, evntd)))
#' rate0 <- with(subset(Ex1delayedEffect, trt == 0), fit_pwexp(Surv(month, evntd)))
#'
#' rate1
#' rate0
#' rate1$rate / rate0$rate
#'
#' # chi-square test for (any) treatment effect (8 - 4 parameters = 4 df)
#' pchisq(sum(rateall$m2ll) - sum(rate1$m2ll + rate0$m2ll),
#'   df = 4,
#'   lower.tail = FALSE
#' )
#'
#' # compare with logrank
#' survdiff(formula = Surv(month, evntd) ~ trt, data = Ex1delayedEffect)
#'
#' # example 3
#' # simple model with 3 rates same for each for 3 months,
#' # different for each treatment after months
#' rate1a <- with(subset(Ex1delayedEffect, trt == 1), fit_pwexp(Surv(month, evntd), 3))
#' rate0a <- with(subset(Ex1delayedEffect, trt == 0), fit_pwexp(Surv(month, evntd), 3))
#' rate1a$rate / rate0a$rate
#'
#' m2ll0 <- rateall$m2ll[1] + rate1a$m2ll[2] + rate0a$m2ll[2]
#' m2ll1 <- sum(rate0$m2ll) + sum(rate1$m2ll)
#'
#' # as a measure of strength, chi-square examines improvement in likelihood
#' pchisq(m2ll0 - m2ll1, df = 5, lower.tail = FALSE)
#'
#' @export
#'
fit_pwexp <- function(
    srv = Surv(time = Ex1delayedEffect$month, event = Ex1delayedEffect$evntd),
    intervals = array(3, 3)) {
  if (!is.Surv(srv)) {
    stop("fit_pwexp: srv must be a survival object!")
  }

  # only allow status 0,1
  xx <- data.frame(time = srv[, "time"], status = srv[, "status"])
  if (nrow(subset(xx, status != 0 & status != 1))) {
    stop("fit_pwexp: srv may only have status values of 0 or 1!")
  }

  # check for late observation after sum(intervals)
  if (nrow(subset(xx, time > sum(intervals) & status > 0)) > 0) {
    intervals <- c(intervals, Inf)
  }

  times <- c(0, cumsum(intervals))

  ans <- NULL
  for (i in 1:length(intervals)) {
    dat <- subset(xx, time > times[i])
    dat$status[dat$time > times[i + 1]] <- 0
    dat$time[dat$time > times[i + 1]] <- times[i + 1]
    dat$time <- dat$time - times[i]

    event <- sum(dat$status)
    ttot <- sum(dat$time)
    rate <- event / ttot

    if (ttot > 0) {
      ans_new <- data.frame(
        intervals = intervals[i],
        ttot = ttot,
        event = event,
        rate = rate,
        m2ll = 2 * (rate * ttot - event * log(rate))
      )
      ans <- rbind(ans, ans_new)
    }
  }

  return(ans)
}
