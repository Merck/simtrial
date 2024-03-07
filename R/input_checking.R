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

#' Check if the input `x` is a single non-negative number (or missing)
#'
#' @param x a scalar value
#' @param require_whole_number Must the numbers be whole numbers (default: FALSE)
#'
#' @return an error message or nothing
#' @noRd
#' @examples
#' a <- 100
#' input_check_scalar(a)
input_check_scalar <- function(x = NA, require_whole_number = FALSE) {
  label <- deparse(substitute(x))
  if (length(x) == 1 && is.na(x)) {
    return(invisible())
  }

  passed <- length(x) == 1 && is.numeric(x) && x >= 0
  if (!passed) {
    stop(paste0(label, " must be a single non-negative number (or NA)"))
  }

  if (require_whole_number) {
    if (!is_whole_number(x)) {
      stop(paste0(label, " must be a single non-negative whole number (or NA)"))
    }
  }

  return(invisible())
}

#' Check if the input `x` is a vector of positive numbers (missing values allowed)
#'
#' @param x a vector
#' @param require_whole_number Must the numbers be whole numbers (default: FALSE)
#'
#' @return an error message or nothing
#' @noRd
#' @examples
#' a <- 1:3
#' input_check_vector(a)
#' b <- c(1, 2, NA)
#' input_check_vector(b)
input_check_vector <- function(x = NA, require_whole_number = FALSE) {
  label <- deparse(substitute(x))
  missing <- is.na(x)
  positive <- is.numeric(x) & x > 0

  if (require_whole_number) {
    whole_number <- is_whole_number(x)
    passed <- all(missing | (positive & whole_number))
  } else {
    passed <- all(missing | positive)
  }

  if (passed) {
    return(invisible())
  }

  if (require_whole_number) {
    stop(paste0(label, " must be a vector with only positive whole numbers and missing values"))
  } else {
    stop(paste0(label, " must be a vector with only positive numbers and missing values"))
  }
}

#' Test if numbers are whole numbers
#'
#' @param x a numeric vector
#' @param tol tolerance
#'
#' @return TRUE, FALSE, or NA
#' @seealso [base::is.integer()]
#' @noRd
#' @examples
#' x <- c(1.1, -1.1, 0, 2, NA)
#' is_whole_number(x)
#' ## [1] FALSE FALSE  TRUE  TRUE    NA
is_whole_number <- function(x, tol = .Machine$double.eps^0.5) {
  return(abs(x - round(x)) < tol)
}
