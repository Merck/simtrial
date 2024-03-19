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

#' Maxcombo test
#'
#' WARNING: This experimental function is a work-in-progress. The function
#' arguments will change as we add additional features.
#'
#' @param data a tte dataset
#' @param rho Numeric vector passed to [fh_weight()]. Must be greater
#'   than or equal to zero. Must be the same length as `gamma`.
#' @param gamma Numeric vector passed to [fh_weight()]. Must be
#'   greater than or equal to zero. Must be the same length as `rho`.
#'
#' @return pvalues
#' @export
#'
#' @seealso [fh_weight()]
#'
#' @examples
#' sim_pw_surv(n = 200) |>
#'   cut_data_by_event(150) |>
#'   maxcombo(rho = c(0, 0), gamma = c(0, 0.5))
maxcombo <- function(data, rho, gamma) {
  stopifnot(
    is.numeric(rho), is.numeric(gamma),
    rho >= 0, gamma >= 0,
    length(rho) == length(gamma)
  )

  ans <- data |>
    counting_process(arm = "experimental") |>
    fh_weight(
      rho_gamma = data.frame(rho = rho, gamma = gamma),
      return_corr = TRUE
    )

  ans <- data.frame(p_value = pvalue_maxcombo(ans))

  return(ans)
}
