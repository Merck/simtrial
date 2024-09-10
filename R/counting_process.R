#  Copyright (c) 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
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

#' Process survival data into counting process format
#'
#' Produces a data frame that is sorted by stratum and time.
#' Included in this is only the times at which one or more event occurs.
#' The output dataset contains stratum, TTE (time-to-event),
#' at risk count, and count of events at the specified TTE
#' sorted by stratum and TTE.
#'
#' @details
#' The function only considered two group situation.
#'
#' The tie is handled by the Breslow's Method.
#'
#' @param x A data frame with no missing values and contain variables:
#'   - `stratum`: Stratum.
#'   - `treatment`: Treatment group.
#'   - `tte`: Observed time.
#'   - `event`: Binary event indicator, `1` represents event,
#'     `0` represents censoring.
#' @param arm Value in the input `treatment` column that indicates
#'   treatment group value.
#'
#' @return
#' A data frame grouped by `stratum` and sorted within stratum by `tte`.
#' Remain rows with at least one event in the population, at least one subject
#' is at risk in both treatment group and control group.
#' Other variables in this represent the following within each stratum at
#' each time at which one or more events are observed:
#' - `event_total`: Total number of events
#' - `event_trt`: Total number of events at treatment group
#' - `n_risk_total`: Number of subjects at risk
#' - `n_risk_trt`: Number of subjects at risk in treatment group
#' - `s`: Left-continuous Kaplan-Meier survival estimate
#' - `o_minus_e`: In treatment group, observed number of events minus expected
#'   number of events. The expected number of events is estimated by assuming
#'   no treatment effect with hypergeometric distribution with parameters total
#'   number of events, total number of events at treatment group and number of
#'   events at a time. (Same assumption of log-rank test under the null
#'   hypothesis)
#' - `var_o_minus_e`: Variance of `o_minus_e` under the same assumption.
#'
#' @importFrom data.table ":=" as.data.table setDF uniqueN
#'
#' @export
#'
#' @details
#' The output produced by [counting_process()] produces a
#' counting process dataset grouped by stratum and sorted within stratum
#' by increasing times where events occur. The object is assigned the class
#' "counting_process". It also has the attributes "n_ctrl" and "n_exp",
#' which are the totals of the control and experimental treatments,
#' respectively, from the input time-to-event data.
#'
#' @examples
#' # Example 1
#' x <- data.frame(
#'   stratum = c(rep(1, 10), rep(2, 6)),
#'   treatment = rep(c(1, 1, 0, 0), 4),
#'   tte = 1:16,
#'   event = rep(c(0, 1), 8)
#' )
#' counting_process(x, arm = 1)
#'
#' # Example 2
#' x <- sim_pw_surv(n = 400)
#' y <- cut_data_by_event(x, 150) |> counting_process(arm = "experimental")
#' # Weighted logrank test (Z-value and 1-sided p-value)
#' z <- sum(y$o_minus_e) / sqrt(sum(y$var_o_minus_e))
#' c(z, pnorm(z))
counting_process <- function(x, arm) {
  unique_treatment <- unique(x$treatment)

  if (length(unique_treatment) > 2) {
    stop("counting_process: expected two groups.")
  }

  if (!arm %in% unique_treatment) {
    stop("counting_process: arm is not a valid treatment group value.")
  }

  if (!all(unique(x$event) %in% c(0, 1))) {
    stop("counting_process: event indicator must be 0 (censoring) or 1 (event).")
  }

  # initilize the number of subjects at risk
  ans <- as.data.table(x)
  ans <- ans[order(tte, decreasing = TRUE), ]
  ans[, one := 1]
  ans[, `:=`(
    n_risk_total = cumsum(one),
    n_risk_trt = cumsum(treatment == arm)
  ), by = "stratum"]

  # Handling ties using Breslow's method
  if (uniqueN(ans[, .(stratum, tte)]) < nrow(ans)) { # ties
    ans[, mtte := -tte]
    ans <- ans[, .(
      event_total = sum(event),
      event_trt = sum((treatment == arm) * event),
      tte = tte[1],
      n_risk_total = max(n_risk_total),
      n_risk_trt = max(n_risk_trt)
    ), by = c("stratum", "mtte")]
    ans[, mtte := NULL]
  } else { # no ties
    ans <- ans[, .(
      stratum,
      event_total = event,
      event_trt = (treatment == arm) * event,
      tte,
      n_risk_total,
      n_risk_trt
    )]
  }

  # Keep calculation for observed time with at least one event,
  # at least one subject is at risk in both treatment group and control group.
  ans <- ans[event_total > 0 & n_risk_total - n_risk_trt > 0 & n_risk_trt > 0, ]
  ans[, s := 1 - event_total / n_risk_total]
  ans <- ans[order(stratum, tte), ]
  # Left continuous Kaplan-Meier Estimator
  ans[, s := c(1, cumprod(s)[-length(s)]), by = "stratum"]
  # Observed events minus Expected events in treatment group
  ans[, o_minus_e := event_trt - n_risk_trt / n_risk_total * event_total]
  # Variance of o_minus_e
  ans[, var_o_minus_e := (n_risk_total - n_risk_trt) *
    n_risk_trt * event_total * (n_risk_total - event_total) /
    n_risk_total^2 / (n_risk_total - 1)]

  setDF(ans)
  class(ans) <- c("counting_process", class(ans))

  attr(ans, "ratio") <- attributes(x)$ratio

  return(ans)
}
