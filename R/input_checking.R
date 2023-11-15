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

#' Check if the input `x` is a single non-negative number (or missing)
#'
#' @param x a scalar value
#' @param label the label of `x`
#'
#' @return an error message or nothing
#' @noRd
#' @examples
#' input_check_scalar(x = 100, label = "my_x")
input_check_scalar <- function(x = NA, label = "x") {
  if (length(x) == 1 && is.numeric(x) && x >= 0) {
    return(invisible())
  }
  if (length(x) == 1 && is.na(x)) {
    return(invisible())
  }

  stop(paste0(label, " must be a single non-negative number (or NA)"))
}

#' Check if the input `x` is a vector of positive numbers (missing values allowed)
#'
#' @param x a vector
#' @param label the label of `x`
#'
#' @return an error message or nothing
#' @noRd
#' @examples
#' input_check_vector(x = 1:3, label = "my_x")
#' input_check_vector(x = c(1, 2, NA), label = "my_x")
input_check_vector <- function(x = NA, label = "x") {
  missing <- is.na(x)
  positive <- is.numeric(x) & x > 0
  if (all(missing | positive)) {
    return(invisible())
  }

  stop(paste0(label, " must be a vector with only positive numbers and missing values"))
}
