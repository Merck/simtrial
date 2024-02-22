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
#' @return A data frame. The column `mb_weight` contains the weights for the
#'   Magirr-Burman weighted logrank test for the data in `x`.
#'
#' @details
#' Magirr and Burman (2019) proposed a weighted logrank test to have better
#' power than the logrank test when the treatment effect is delayed,
#' but to still maintain good power under a proportional hazards assumption.
#' In Magirr (2021), (the equivalent of) a maximum weight was proposed
#' as opposed to a fixed time duration over which weights would increase.
#' The weights for some early interval specified by the user are the inverse
#' of the combined treatment group empirical survival distribution; see details.
#' After this initial period, weights are constant at the maximum of the
#' previous weights. Another advantage of the test is that under strong
#' null hypothesis that the underlying survival in the control group is
#' greater than or equal to underlying survival in the experimental group,
#' Type I error is controlled as the specified level.
#'
#' We define \eqn{t^*} to be the input variable `delay`.
#' This specifies an initial period during which weights increase.
#' We also set a maximum weight \eqn{w_{\max}}.
#' To define specific weights, we let \eqn{S(t)} denote the Kaplan-Meier
#' survival estimate at time \eqn{t} for the combined data
#' (control plus experimental treatment groups).
#' The weight at time \eqn{t} is then defined as
#' \deqn{w(t)=\min(w_{\max}, S(\min(t, t^*))^{-1}).}
#'
#' @references
#' Magirr, Dominic, and Carl‐Fredrik Burman. 2019.
#' "Modestly weighted logrank tests."
#' _Statistics in Medicine_ 38 (20): 3782--3790.
#'
#' Magirr, Dominic. 2021.
#' "Non‐proportional hazards in immuno‐oncology: Is an old perspective needed?"
#' _Pharmaceutical Statistics_ 20 (3): 512--527.
#'
#' @importFrom data.table ":=" as.data.table data.table fifelse merge.data.table setDF
#'
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' # Use default enrollment and event rates at cut at 100 events
#' # For transparency, may be good to set either `delay` or `w_max` to `Inf`
#' x <- sim_pw_surv(n = 200) |>
#'   cut_data_by_event(125) |>
#'   counting_process(arm = "experimental")
#'
#' # Example 1
#' # Compute Magirr-Burman weights with `delay = 6`
#' ZMB <- x |>
#'   mb_weight(delay = 6, w_max = Inf) |>
#'   summarize(
#'     S = sum(o_minus_e * mb_weight),
#'     V = sum(var_o_minus_e * mb_weight^2),
#'     z = S / sqrt(V)
#'   )
#'
#' # Compute p-value of modestly weighted logrank of Magirr-Burman
#' pnorm(ZMB$z)
#'
#' # Example 2
#' # Now compute with maximum weight of 2 as recommended in Magirr, 2021
#' ZMB2 <- x |>
#'   mb_weight(delay = Inf, w_max = 2) |>
#'   summarize(
#'     S = sum(o_minus_e * mb_weight),
#'     V = sum(var_o_minus_e * mb_weight^2),
#'     z = S / sqrt(V)
#'   )
#'
#' # Compute p-value of modestly weighted logrank of Magirr-Burman
#' pnorm(ZMB2$z)
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

  # Compute max weight by stratum
  x2 <- as.data.table(x)
  # Make sure you don't lose any stratum!
  tbl_all_stratum <- data.table(stratum = unique(x2$stratum))

  # Look only up to delay time
  ans <- x2[tte <= delay, ]
  # Weight before delay specified as 1/S
  ans <- ans[, .(max_weight = max(1 / s)), by = "stratum"]
  # Get back stratum with no records before delay ends
  ans <- ans[tbl_all_stratum, on = "stratum"]
  # `max_weight` is 1 when there are no records before delay ends
  ans[, max_weight := fifelse(is.na(max_weight), 1, max_weight)]
  # Cut off weights at w_max
  ans[, max_weight := pmin(w_max, max_weight)]
  # Now merge max_weight back to stratified dataset
  ans <- merge.data.table(ans, x2, by = "stratum", all = TRUE)
  # Weight is min of max_weight and 1/S which will increase up to delay
  ans[, mb_weight := pmin(max_weight, 1 / s)]
  ans[, max_weight := NULL]

  setDF(ans)
  return(ans)
}
