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
#' @return A data frame containing the method (\code{method}), always "milestone",
#' test statistics (\code{z}), milestone time point (\code{ms_time}),
#' survival rate of the control arm (\code{surv0}),
#' survival rate of the experimental arm (\code{surv1}),
#' survival difference between the experimental and control arm (\code{surv_diff}),
#' standard error of the control arm (\code{std_err0}),
#' standard error of the experimental arm (\code{std_err1}).
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
  diff_survival <- fit_res$surv[2] - fit_res$surv[1]

  # Indicator whether the std is NA or not
  na_col <- is.na(fit_res$std.err[1])
  na_exp <- is.na(fit_res$std.err[2])

  # Calculate the test statistics
  if (na_col + na_exp == 2) {
    z <- -Inf
  } else {
    term1 <- if (na_col) 0 else fit_res$std.err[1]
    term2 <- if (na_exp) 0 else fit_res$std.err[2]
    var_survival <- term1^2 + term2^2
    z <- diff_survival / sqrt(var_survival)
  }

  ans <- data.frame(method = "milestone", z = z, ms_time = ms_time,
                    surv0 = fit_res$surv[1], surv1 = fit_res$surv[2],
                    surv_diff = diff_survival,
                    std_err0 = fit_res$std.err[1], std_err1 = fit_res$std.err[2])
  return(ans)
}
