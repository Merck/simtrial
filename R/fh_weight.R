#  Copyright (c) 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
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

#' Fleming-Harrington weighted logrank tests
#'
#' With output from the function [counting_process()].
#'
#' @param x A [counting_process()]-class data frame with a counting process
#'   dataset.
#' @param rho A numerical value greater than equal to zero,
#' to specify one Fleming-Harrington weighted logrank test
#' @param gamma A numerical value greater than equal to zero,
#' to specify one Fleming-Harrington weighted logrank test
#'
#' @return A data frame with `rho` and `gamma` as input and the FH weights
#' for the data in `x`.
#'
#' @importFrom data.table data.table merge.data.table setDF
#' @noRd
fh_weight <- function(
    x = sim_pw_surv(n = 200) |>
      cut_data_by_event(150) |>
      counting_process(arm = "experimental"),
    rho = 0,
    gamma = 1) {
  # Input checking ----
  if (length(rho) != 1 || length(gamma) != 1) {
    stop("fh_weight: please input single numerical values of rho and gamma.")
  }

  if (!is.data.frame(x)) {
    stop("fh_weight: x in `fh_weight()` must be a data frame.")
  }

  if (!("s" %in% names(x))) {
    stop("fh_weight: x column names in `fh_weight()` must contain s.")
  }

  if (!("o_minus_e" %in% names(x))) {
    stop("fh_weight: x column names in `fh_weight()` must contain o_minus_e.")
  }

  if (!("var_o_minus_e" %in% names(x))) {
    stop("fh_weight: x column names in `fh_weight()` must contain var_o_minus_e.")
  }

  x$weight <- x$s^rho * (1 - x$s)^gamma

  return(x)
}
