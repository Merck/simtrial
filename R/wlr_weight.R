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

#' Fleming-Harrington weighting function
#'
#' @param rho Non-negative number. \code{rho = 0, gamma = 0} is equivalent to regular logrank test.
#' @param gamma Non-negative number. \code{rho = 0, gamma = 0} is equivalent to regular logrank test.
#'
#' @export
#' @return A list of parameters of the Fleming-Harrington weighting function
#' @examples
#' fh(rho = 0, gamma = 0.5)
fh <- function(rho = 0, gamma = 0){
  structure(list(rho = rho, gamma = gamma), class = c("list", "fh", "wlr"))
}

#' Magirr and Burman weighting function
#'
#' @param delay The initial delay period where weights increase;
#'   after this, weights are constant at the final weight in the delay period.
#' @param w_max Maximum weight to be returned.
#'   Set `delay = Inf`, `w_max = 2` to be consistent with recommendation of
#'   Magirr (2021).
#'
#' @return A list of parameters of the Magirr and Burman weighting function
#' @export
#'
#' @examples
#' mb(delay = 6, w_max = 2)
mb <- function(delay = 4, w_max = Inf){
  structure(list(delay = delay, w_max = w_max), class = c("list", "mb", "wlr"))
}

#' Zero early weighting function
#'
#' @param early_period The initial delay period where weights increase;
#'   after this, weights are constant at the final weight in the delay period.
#'
#' @return A list of parameters of the zero early weighting function
#' @references
#' Xu, Z., Zhen, B., Park, Y., & Zhu, B. (2017).
#' "Designing therapeutic cancer vaccine trials with delayed treatment effect."
#' @export
#'
#' @examples
#' early_zero(6)
early_zero <- function(early_period){
  structure(list(early_period = early_period), class = c("list", "early_period", "wlr"))
}
