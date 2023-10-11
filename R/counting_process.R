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

#' Process survival data into counting process format
#'
#' Produces a tibble that is sorted by stratum and time.
#' Included in this is only the times at which one or more event occurs.
#' The output dataset contains stratum, tte (time-to-event),
#' at risk count and count of events at the specified tte
#' sorted by stratum and tte.
#'
#' The function only considered two group situation.
#'
#' The tie is handled by the Breslow's Method.
#'
#' @param x A tibble with no missing values and contain variables:
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
#'
#' - `events`: Total number of events
#' - `n_event_tol`: Total number of events at treatment group
#' - `n_risk_tol`: Number of subjects at risk
#' - `n_risk_trt`: Number of subjects at risk in treatment group
#' - `S`: Left-continuous Kaplan-Meier survival estimate
#' - `o_minus_e`: In treatment group, observed number of events minus expected
#'   number of events. The expected number of events is estimated by assuming
#'   no treatment effect with hypergeometric distribution with parameters total
#'   number of events, total number of events at treatment group and number of
#'   events at a time. (Same assumption of log-rank test under the null
#'   hypothesis)
#' - `var_o_minus_e`: Variance of `o_minus_e` under the same assumption.
#'
#' @importFrom data.table ":=" as.data.table setDF
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' # Example 1
#' x <- tibble(
#'   stratum = c(rep(1, 10), rep(2, 6)),
#'   treatment = rep(c(1, 1, 0, 0), 4),
#'   tte = 1:16,
#'   event = rep(c(0, 1), 8)
#' )
#' counting_process(x, arm = 1)
#'
#' # Example 2
#' x <- sim_pw_surv(n = 400)
#' y <- cut_data_by_event(x, 150) %>% counting_process(arm = "experimental")
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

  ans <- as.data.table(x)
  ans <- ans[order(tte, decreasing = TRUE), ]
  ans[, one := 1]
  ans[, `:=`(
    n_risk_tol = cumsum(one),
    n_risk_trt = cumsum(treatment == arm)
  ), by = "stratum"]

  # Handling ties using Breslow's method
  ans[, mtte := -tte]
  ans <- ans[, .(
    events = sum(event),
    n_event_tol = sum((treatment == arm) * event),
    tte = tte[1],
    n_risk_tol = max(n_risk_tol),
    n_risk_trt = max(n_risk_trt)
  ), by = c("stratum", "mtte")]

  # Keep calculation for observed time with at least one event,
  # at least one subject is at risk in both treatment group and control group.
  ans <- ans[events > 0 & n_risk_tol - n_risk_trt > 0 & n_risk_trt > 0, ]
  ans[, mtte := NULL]
  ans[, s := 1 - events / n_risk_tol]
  ans <- ans[order(stratum, tte), ]
  # Left continuous Kaplan-Meier Estimator
  ans[, s := c(1, cumprod(s)[-length(s)]), by = "stratum"]
  # Observed events minus Expected events in treatment group
  ans[, o_minus_e := n_event_tol - n_risk_trt / n_risk_tol * events]
  # Variance of o_minus_e
  ans[, var_o_minus_e := (n_risk_tol - n_risk_trt) *
    n_risk_trt * events * (n_risk_tol - events) /
    n_risk_tol^2 / (n_risk_tol - 1)]

  return(setDF(ans))
}
