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

#' Fleming-Harrington weighted logrank tests
#'
#' With output from the function [counting_process()].
#'
#' @param x A [counting_process()]-class `tibble` with a counting process dataset.
#' @param rho_gamma A `tibble` with variables `rho` and `gamma`, both greater
#'   than equal to zero, to specify one Fleming-Harrington weighted logrank test
#'   per row; Default: `tibble(rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))`.
#' @param return_variance A logical flag that, if `TRUE`, adds columns
#'   estimated variance for weighted sum of observed minus expected;
#'   see details; Default: `FALSE`.
#'
#' @return A `tibble` with `rho_gamma` as input and the FH test statistic
#'   for the data in `x`. (`z`, a directional square root of the usual
#'   weighted logrank test); if variance calculations are specified
#'   (for example, to be used for covariances in a combination test),
#'   then this will be returned in the column `Var`.
#'
#' @details
#' The input value `x` produced by [counting_process()] produces a
#' counting process dataset grouped by stratum and sorted within stratum
#' by increasing times where events occur.
#' - \eqn{z} - Standardized normal Fleming-Harrington weighted logrank test.
#' - \eqn{i} - Stratum index.
#' - \eqn{d_i} - Number of distinct times at which events occurred in
#'   stratum \eqn{i}.
#' - \eqn{t_{ij}} - Ordered times at which events in stratum
#'   \eqn{i}, \eqn{j = 1, 2, \ldots, d_i} were observed;
#'   for each observation, \eqn{t_{ij}} represents the time post study entry.
#' - \eqn{O_{ij.}} - Total number of events in stratum \eqn{i} that occurred
#'   at time \eqn{t_{ij}}.
#' - \eqn{O_{ije}} - Total number of events in stratum \eqn{i} in the
#'   experimental treatment group that occurred at time \eqn{t_{ij}}.
#' - \eqn{N_{ij.}} - Total number of study subjects in stratum \eqn{i}
#'   who were followed for at least duration.
#' - \eqn{E_{ije}} - Expected observations in experimental treatment group
#'   given random selection of \eqn{O_{ij.}} from those in
#'   stratum \eqn{i} at risk at time \eqn{t_{ij}}.
#' - \eqn{V_{ije}} - Hypergeometric variance for \eqn{E_{ije}} as
#'   produced in `Var` from [counting_process()].
#' - \eqn{N_{ije}} - Total number of study subjects in
#'   stratum \eqn{i} in the experimental treatment group
#'   who were followed for at least duration \eqn{t_{ij}}.
#' - \eqn{E_{ije}} - Expected observations in experimental group in
#'   stratum \eqn{i} at time \eqn{t_{ij}} conditioning on the overall number
#'   of events and at risk populations at that time and sampling at risk
#'   observations without replacement:
#'   \deqn{E_{ije} = O_{ij.} N_{ije}/N_{ij.}}
#' - \eqn{S_{ij}} - Kaplan-Meier estimate of survival in combined
#'   treatment groups immediately prior to time \eqn{t_{ij}}.
#' - \eqn{\rho, \gamma} - Real parameters for Fleming-Harrington test.
#' - \eqn{X_i} - Numerator for signed logrank test in stratum \eqn{i}
#'   \deqn{X_i = \sum_{j=1}^{d_{i}} S_{ij}^\rho(1-S_{ij}^\gamma)(O_{ije}-E_{ije})}
#' - \eqn{V_{ij}} - Variance used in denominator for Fleming-Harrington
#'   weighted logrank tests
#'   \deqn{V_i = \sum_{j=1}^{d_{i}} (S_{ij}^\rho(1-S_{ij}^\gamma))^2V_{ij})}
#'   The stratified Fleming-Harrington weighted logrank test is then computed as:
#'   \deqn{z = \sum_i X_i/\sqrt{\sum_i V_i}.}
#'
#' @importFrom tibble tibble
#' @import dplyr
#'
#' @export
#'
#' @examples
#' library(tidyr)
#'
#' # Use default enrollment and event rates at cut at 100 events
#' x <- sim_pw_surv(n = 200) %>%
#'   cut_data_by_event(100) %>%
#'   counting_process(arm = "experimental")
#'
#' # Compute logrank (FH(0,0)) and FH(0,1)
#' wlr(x, rho_gamma = tibble(rho = c(0, 0), gamma = c(0, 1)))
wlr <- function(
    x = sim_pw_surv(n = 200) %>%
      cut_data_by_event(150) %>%
      counting_process(arm = "experimental"),
    rho_gamma = tibble(
      rho = c(0, 0, 1, 1),
      gamma = c(0, 1, 0, 1)
    ),
    return_variance = FALSE) {
  # Check input failure rate assumptions
  if (!is.data.frame(x)) {
    stop("wlr: x in `wlr()` must be a data frame.")
  }

  if (!("s" %in% names(x))) {
    stop("wlr: x column names in `wlr()` must contain s.")
  }

  if (!("o_minus_e" %in% names(x))) {
    stop("wlr: x column names in `wlr()` must contain o_minus_e.")
  }

  if (!("var_o_minus_e" %in% names(x))) {
    stop("wlr: x column names in `wlr()` must contain var_o_minus_e.")
  }

  # Get minimal columns from counting_process item
  xx <- x %>%
    ungroup() %>%
    select(s, o_minus_e, var_o_minus_e)

  rho_gamma$z <- rep(0, nrow(rho_gamma))

  if (return_variance) {
    rho_gamma$Var <- rep(0, nrow(rho_gamma))
  }

  for (i in 1:nrow(rho_gamma)) {
    y <- xx %>%
      mutate(
        weight = s^rho_gamma$rho[i] * (1 - s)^rho_gamma$gamma[i],
        weighted_o_minus_e = weight * o_minus_e,
        weighted_var = weight^2 * var_o_minus_e
      ) %>%
      summarize(
        weighted_var = sum(weighted_var),
        weighted_o_minus_e = sum(weighted_o_minus_e)
      )

    rho_gamma$z[i] <- y$weighted_o_minus_e / sqrt(y$weighted_var)

    if (return_variance) {
      rho_gamma$Var[i] <- y$weighted_var
    }
  }

  rho_gamma
}
