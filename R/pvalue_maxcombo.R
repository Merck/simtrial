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

#' MaxCombo p-value
#'
#' Computes p-values for the MaxCombo test based on output from [tenFHcorr()].
#' This is still in an experimental stage and is intended for use with
#' the [sim_fixed_n()] trial simulation routine.
#' However, it can also be used to analyze clinical trial data such as
#' that provided in the ADaM ADTTE format.
#'
#' @param z A dataset output from [tenFHcorr()]; see examples.
#' @param dummy_var A dummy input that allows [dplyr::group_map()] to be used to
#'   compute p-values for multiple simulations.
#' @param algorithm This is passed directly to the `algorithm` argument
#'   in [mvtnorm::pmvnorm()].
#'
#' @return A numeric p-value.
#'
#' @importFrom mvtnorm pmvnorm GenzBretz
#' @importFrom dplyr select starts_with
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr)
#'
#' # Example 1
#' x <- sim_fixed_n(
#'   n_sim = 1,
#'   timing_type = 5,
#'   rho_gamma = tibble(
#'     rho = c(0, 0, 1),
#'     gamma = c(0, 1, 1)
#'   )
#' )
#' head(x)
#' pvalue_maxcombo(x)
#'
#' # Example 2
#' # Only use cuts for events, events + min follow-up
#' xx <- sim_fixed_n(
#'   n_sim = 100,
#'   timing_type = 5,
#'   rho_gamma = tibble(
#'     rho = c(0, 0, 1),
#'     gamma = c(0, 1, 1)
#'   )
#' )
#' head(xx)
#'
#' # MaxCombo power estimate for cutoff at max of targeted events, minimum follow-up
#' p <- xx %>%
#'   group_by(sim) %>%
#'   group_map(pvalue_maxcombo) %>%
#'   unlist()
#' mean(p < .025)
pvalue_maxcombo <- function(
    z,
    dummy_var,
    algorithm = mvtnorm::GenzBretz(maxpts = 50000, abseps = 0.00001)) {
  ans <- (1 - mvtnorm::pmvnorm(
    lower = rep(
      z$z %>% min() %>% as.numeric(),
      nrow(z)
    ),
    corr = z %>%
      select(starts_with("V")) %>%
      data.matrix(),
    algorithm = algorithm
  )[1]
  ) %>% as.numeric()

  ans
}
