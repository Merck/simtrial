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

#' Zero early weight for weighted logrank tests
#'
#' @param x A [counting_process()]-class data frame with a counting process dataset.
#' @param early_period The initial delay period where weights increase;
#'   after this, weights are constant at the final weight in the delay period.
#' @param fail_rate A data frame record the failure rate.
#'
#' @return A data frame. The column `weight` contains the weights for the
#'   early zero weighted logrank test for the data in `x`.
#'
#' @importFrom data.table ":=" as.data.table fifelse merge.data.table setDF
#' @noRd
early_zero_weight <- function(x, early_period = 4, fail_rate = NULL) {
  ans <- as.data.table(x)
  n_stratum <- length(unique(ans$stratum))

  # If it is unstratified design
  if (n_stratum == 1) {
    ans[, weight := fifelse(tte < early_period, 0, 1)]
  } else {
    if (is.null(fail_rate)) {
      stop("For stratified design to use `early_zero_weight()`, `fail_rate` can't be `NULL`.")
    }
    fail_rate <- as.data.table(fail_rate)
    # require 2 piece failure rate per stratum
    two_piece_fr <- fail_rate[, .(check = .N == 2), by = "stratum"]
    if (!all(two_piece_fr$check)) {
      stop("`early_zero_weight()` only allows delayed treatment effect, that is, 2 piece failure rate with HR = 1 at the first period.")
    }

    late_hr <- fail_rate[hr != 1, .(stratum, hr)]
    delay_change_time <- fail_rate[hr == 1, .(stratum, duration)]

    ans <- merge.data.table(ans, late_hr, by = "stratum", all.x = TRUE)
    ans <- merge.data.table(ans, delay_change_time, by = "stratum", all.x = TRUE)
    ans[, weight := fifelse(tte < duration, 0, log(hr))]
  }

  setDF(ans)

  return(ans)
}
