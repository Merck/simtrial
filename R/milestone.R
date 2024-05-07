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
#' @param test_type Method to build the test statistics.
#' There are 2 options:
#'   - `"native"`: a native approach by dividing the KM survival difference by its standard derivations,
#'   see equation (1) of Klein, J. P., Logan, B., Harhoff, M., & Andersen, P. K. (2007).
#'   - `"log-log"`: a log-log transformation of the survival, see equation (3) of
#'   Klein, J. P., Logan, B., Harhoff, M., & Andersen, P. K. (2007).
#' @return A list frame containing:
#'   - `method` - The method, always `"milestone"`.
#'   - `parameter` - Milestone time point.
#'   - `estimate` - Survival difference between the experimental and control arm.
#'   - `se` - Standard error of the control and experimental arm.
#'   - `z` - Test statistics.
#'
#' @references
#' Klein, J. P., Logan, B., Harhoff, M., & Andersen, P. K. (2007).
#' "Analyzing survival curves at a fixed point in time."
#' _Statistics in Medicine_, 26(24), 4505--4519.
#'
#' @export
#'
#' @examples
#' cut_data <- sim_pw_surv(n = 200) |>
#'   cut_data_by_event(150)
#'
#' cut_data |>
#'   milestone(10, test_type = "log-log")
#'
#' cut_data |>
#'   milestone(10, test_type = "naive")
milestone <- function(data, ms_time, test_type = c("log-log", "naive")) {
  test_type <- match.arg(test_type)

  # Fit into KM curves
  fit <- survfit(Surv(tte, event) ~ treatment, data = data)
  fit_res <- summary(fit, time = ms_time, extend = TRUE)

  # Survival difference
  surv_ctrl <- fit_res$surv[1]
  surv_exp <- fit_res$surv[2]
  surv_diff <- surv_exp - surv_ctrl

  # Indicator whether the std is NA or not
  var_ctrl <- fit_res$std.err[1]^2
  var_exp <- fit_res$std.err[2]^2

  na_ctrl <- is.na(var_ctrl)
  na_exp <- is.na(var_exp)

  sigma2_ctrl <- var_ctrl / (surv_ctrl^2)
  sigma2_exp <- var_exp / (surv_exp^2)

  # Calculate the test statistics
  ans <- list()
  ans$method <- "milestone"
  ans$parameter <- ms_time

  if (na_ctrl + na_exp == 2) {
    z <- -Inf
  } else {
    if (test_type == "naive") {
      z_numerator <- surv_diff
      z_denominator <- surv_exp * sqrt(sigma2_exp) + surv_ctrl * sqrt(sigma2_ctrl)
    } else if (test_type == "log-log") {
      z_numerator <- log(-log(surv_exp)) - log(-log(surv_ctrl))
      z_denominator <- sqrt(sigma2_exp) / log(surv_exp) + sqrt(sigma2_ctrl) / log(surv_ctrl)
    }
  }

  ans$estimate <- z_numerator
  ans$se <- z_denominator
  ans$z <- z_numerator / z_denominator

  return(ans)
}
