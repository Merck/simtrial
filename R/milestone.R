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

#' Milestone test for two survival curves
#'
#' @param data Data frame containing at least 3 columns:
#'   - `tte` - Time to event.
#'   - `event` - Event indicator.
#'   - `treatment` - Grouping variable.
#' @param ms_time Milestone analysis time.
#'
#' @return A data frame containing:
#' - `method` - The method, always `"milestone"`.
#' - `z` - Test statistics.
#' - `ms_time` - Milestone time point.
#' - `surv0` - Survival rate of the control arm.
#' - `surv1` - Survival rate of the experimental arm.
#' - `surv_diff` - Survival difference between the experimental and control arm.
#' - `std_err0` - Standard error of the control arm.
#' - `std_err1` - Standard error of the experimental arm.
#' @references Klein, John P., et al. "Analyzing survival curves at a fixed point in time."
#' Statistics in medicine 26.24 (2007): 4505-4519.
#'
#' @export
#'
#' @examples
#' sim_pw_surv(n = 200) |>
#'   cut_data_by_event(150) |>
#'   milestone(10)
milestone <- function(data, ms_time) {
  # Fit into KM curves
  fit <- survfit(Surv(tte, event) ~ treatment, data = data)
  fit_res <- summary(fit, time = ms_time, extend = TRUE)

  # Survival difference
  surv_col <- fit_res$surv[1]
  surv_exp <- fit_res$surv[2]
  diff_survival <- surv_exp - surv_col

  # Indicator whether the std is NA or not
  var_col <- fit_res$std.err[1]^2
  var_exp <- fit_res$std.err[2]^2

  na_col <- is.na(var_col)
  na_exp <- is.na(var_exp)

  sigma2_col <- var_col / (surv_col^2)
  sigma2_exp <- var_exp / (surv_exp^2)

  # Calculate the test statistics
  if (na_col + na_exp == 2) {
    z <- -Inf
  } else {
    z <- (log(-log(surv_exp)) - log(-log(surv_col)))^2 / (sigma2_exp / (log(surv_exp))^2 + sigma2_col / (log(surv_col))^2)
  }

  ans <- data.frame(
    method = "milestone", z = z, ms_time = ms_time,
    surv_col = surv_col, surv_exp = surv_exp,
    surv_diff = diff_survival,
    std_err_col = fit_res$std.err[1], std_err_exp = fit_res$std.err[2]
  )
  return(ans)
}
