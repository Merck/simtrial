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

#' Time-to-event data example 1 for non-proportional hazards working group
#'
#' Survival objects reverse-engineered datasets from published Kaplan-Meier
#' curves.
#' Individual trials are de-identified since the data are only
#' approximations of the actual data.
#' Data are intended to evaluate methods and designs for trials where
#' non-proportional hazards may be anticipated for outcome data.
#'
#' @docType data
#'
#' @usage data(Ex1delayedEffect)
#'
#' @format
#' Data frame with 4 variables:
#' - `id`: Sequential numbering of unique identifiers.
#' - `month`: Time-to-event.
#' - `event`: 1 for event, 0 for censored.
#' - `trt`: 1 for experimental, 0 for control.
#'
#' @keywords datasets
#'
#' @references TBD
#'
#' @seealso
#' [Ex2delayedEffect],
#' [Ex3curewithph],
#' [Ex4belly],
#' [Ex5widening],
#' [Ex6crossing]
#'
#' @examples
#' library(survival)
#'
#' data(Ex1delayedEffect)
#' km1 <- with(Ex1delayedEffect, survfit(Surv(month, evntd) ~ trt))
#' km1
#' plot(km1)
#' with(subset(Ex1delayedEffect, trt == 1), survfit(Surv(month, evntd) ~ trt))
#' with(subset(Ex1delayedEffect, trt == 0), survfit(Surv(month, evntd) ~ trt))
"Ex1delayedEffect"
