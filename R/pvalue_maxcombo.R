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
#' Computes p-values for the MaxCombo test based on output from [fh_weight()].
#' This is still in an experimental stage and is intended for use with
#' the [sim_fixed_n()] trial simulation routine.
#' However, it can also be used to analyze clinical trial data such as
#' that provided in the ADaM ADTTE format.
#'
#' @param z A dataset output from [fh_weight()]; see examples.
#' @param algorithm This is passed directly to the `algorithm` argument
#'   in [mvtnorm::pmvnorm()].
#'
#' @return A numeric p-value.
#'
#' @importFrom mvtnorm pmvnorm GenzBretz
#'
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' # Example 1
#' x <- sim_fixed_n(
#'   n_sim = 1,
#'   timing_type = 5,
#'   rho_gamma = data.frame(
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
#'   rho_gamma = data.frame(
#'     rho = c(0, 0, 1),
#'     gamma = c(0, 1, 1)
#'   )
#' )
#' head(xx)
#'
#' # MaxCombo power estimate for cutoff at max of targeted events, minimum follow-up
#' p <- xx %>%
#'   group_by(sim) %>%
#'   group_map(~ pvalue_maxcombo(.x)) %>%
#'   unlist()
#' mean(p < .025)
pvalue_maxcombo <- function(
    z,
    algorithm = mvtnorm::GenzBretz(maxpts = 50000, abseps = 0.00001)) {
  zmin <- as.numeric(min(z$z))
  lower_limits <- rep.int(zmin, times = nrow(z))

  correlation_matrix <- z[, grep("^v", colnames(z))]
  correlation_matrix <- data.matrix(correlation_matrix)

  distribution <- mvtnorm::pmvnorm(
    lower = lower_limits,
    corr = correlation_matrix,
    algorithm = algorithm
  )

  ans <- 1 - distribution[1]

  as.numeric(ans)
}
