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

#' Magirr and Burman modestly weighted logrank tests
#'
#' Computes Magirr-Burman weights and adds them to a dataset created by
#' [counting_process()].
#' These weights can then be used to compute a z-statistic for the
#' modestly weighted logrank test proposed.
#'
#' @param x A [counting_process()]-class data frame with a counting process
#'   dataset.
#' @param delay The initial delay period where weights increase;
#'   after this, weights are constant at the final weight in the delay period.
#' @param w_max Maximum weight to be returned.
#'   Set `delay = Inf`, `w_max = 2` to be consistent with recommendation of
#'   Magirr (2021).
#'
#' @return A list.
#' @importFrom data.table ":=" as.data.table data.table fifelse merge.data.table setDF
mb_weight <- function(x, delay = 4, w_max = Inf) {
  # check input failure rate assumptions
  if (!is.data.frame(x)) {
    stop("x in `mb_weight()` must be a data frame")
  }

  # check input delay
  if (!is.numeric(delay)) {
    stop("delay in `mb_weight()` must be a non-negative number")
  }

  if (!delay >= 0) {
    stop("delay in `mb_weight()` must be a non-negative number")
  }

  if (!is.numeric(w_max)) {
    stop("w_max (maximum weight) in `mb_weight()` must be a positive number")
  }

  if (!delay > 0) {
    stop("w_max (maximum weight) in `mb_weight()` must be a positive number")
  }

  if (max(names(x) == "stratum") != 1) {
    stop("x column names in `mb_weight()` must contain stratum")
  }

  if (max(names(x) == "tte") != 1) {
    stop("x column names in `mb_weight()` must contain tte")
  }

  if (max(names(x) == "s") != 1) {
    stop("x column names in `mb_weight()` must contain s")
  }

  ans <- list()
  ans$method <- "WLR"
  ans$parameter <- paste0("MB(delay = ", delay, ", max_weight = ", w_max, ")")
  # Compute max weight by stratum
  x2 <- as.data.table(x)
  # Make sure you don't lose any stratum!
  tbl_all_stratum <- data.table(stratum = unique(x2$stratum))

  # Look only up to delay time
  res <- x2[tte <= delay, ]
  # Weight before delay specified as 1/S
  res <- res[, .(max_weight = max(1 / s)), by = "stratum"]
  # Get back stratum with no records before delay ends
  res <- res[tbl_all_stratum, on = "stratum"]
  # `max_weight` is 1 when there are no records before delay ends
  res[, max_weight := fifelse(is.na(max_weight), 1, max_weight)]
  # Cut off weights at w_max
  res[, max_weight := pmin(w_max, max_weight)]
  # Now merge max_weight back to stratified dataset
  res <- merge.data.table(res, x2, by = "stratum", all = TRUE)
  # Weight is min of max_weight and 1/S which will increase up to delay
  res[, mb_weight := pmin(max_weight, 1 / s)]
  res[, max_weight := NULL]

  setDF(res)

  ans$estimate <- sum(res$o_minus_e * res$mb_weight)
  ans$se <- sqrt(sum(res$var_o_minus_e * res$mb_weight^2))
  ans$z <- ans$estimate / ans$se

  return(ans)
}
