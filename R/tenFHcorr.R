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

#' Fleming-Harrington weighted logrank tests plus correlations
#'
#' Correlations can be used with [mvtnorm::pmvnorm()] to compute
#' p-value for MaxCombo, the maximum of the specified
#' Fleming-Harrington tests.
#'
#' @param x A [counting_process()]-class `tibble` with a counting process dataset.
#' @param rho_gamma A `tibble` with variables `rho` and `gamma`, both greater
#'   than equal to zero, to specify one Fleming-Harrington weighted logrank test
#'   per row.
#' @param corr Logical; if `TRUE` (default), return correlation matrix;
#'   otherwise, return covariance matrix.
#'
#' @return A `tibble` with `rho_gamma` as input, the FH test statistics
#'   specified for the data in `z`, and the correlation or
#'   covariance matrix for these tests in variables starting with `v`.
#'
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr left_join
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr)
#'
#' # Use default enrollment and event rates at cut of 100 events
#' x <- sim_pw_surv(n = 200) %>%
#'   cut_data_by_event(100) %>%
#'   counting_process(arm = "experimental")
#'
#' # Compute logrank (FH(0,0)) and FH(0,1)
#' x <- x %>% tenFHcorr(rho_gamma = tibble(
#'   rho = c(0, 0),
#'   gamma = c(0, 1)
#' ))
#'
#' # Compute p-value for MaxCombo
#' library(mvtnorm)
#' 1 - pmvnorm(
#'   lower = rep(min(x$z), nrow(x)),
#'   corr = data.matrix(select(x, -c(rho, gamma, z))),
#'   algorithm = GenzBretz(maxpts = 50000, abseps = 0.00001)
#' )[1]
#'
#' # Check that covariance is as expected
#' x <- sim_pw_surv(n = 200) %>%
#'   cut_data_by_event(100) %>%
#'   counting_process(arm = "experimental")
#'
#' x %>% tenFHcorr(
#'   rho_gamma = tibble(
#'     rho = c(0, 0),
#'     gamma = c(0, 1)
#'   ),
#'   corr = FALSE
#' )
#'
#' # Off-diagonal element should be variance in following
#' x %>% tenFHcorr(
#'   rho_gamma = tibble(
#'     rho = 0,
#'     gamma = .5
#'   ),
#'   corr = FALSE
#' )
#'
#' # Compare off diagonal result with wlr()
#' x %>% wlr(rho_gamma = tibble(rho = 0, gamma = .5))
tenFHcorr <- function(
    x = sim_pw_surv(n = 200) %>%
      cut_data_by_event(100) %>%
      counting_process(arm = "experimental"),
    rho_gamma = tibble(
      rho = c(0, 0, 1, 1),
      gamma = c(0, 1, 0, 1)
    ),
    corr = TRUE) {
  n_weight <- nrow(rho_gamma)

  # Get average rho and gamma for FH covariance matrix
  # We want ave_rho[i,j] = (rho[i] + rho[j])/2
  # and     ave_gamma[i,j] = (gamma[i] + gamma[j])/2
  ave_rho <- (matrix(rho_gamma$rho, nrow = n_weight, ncol = n_weight, byrow = FALSE) +
    matrix(rho_gamma$rho, nrow = n_weight, ncol = n_weight, byrow = TRUE)
  ) / 2
  ave_gamma <- (matrix(rho_gamma$gamma, nrow = n_weight, ncol = n_weight) +
    matrix(rho_gamma$gamma, nrow = n_weight, ncol = n_weight, byrow = TRUE)
  ) / 2

  # Convert back to tibble
  rg_new <- tibble(rho = as.numeric(ave_rho), gamma = as.numeric(ave_gamma))
  # Get unique values of rho, gamma
  rg_unique <- rg_new %>% unique()

  # Compute FH statistic for unique values
  # and merge back to full set of pairs
  rg_fh <- rg_new %>% left_join(
    wlr(x, rg_unique, return_variance = TRUE),
    by = c("rho" = "rho", "gamma" = "gamma")
  )

  # Get z statistics for input rho, gamma combinations
  z <- rg_fh$z[(0:(n_weight - 1)) * n_weight + 1:n_weight]

  # Get correlation matrix
  cov_mat <- matrix(rg_fh$Var, nrow = n_weight, byrow = TRUE)

  if (corr) {
    corr_mat <- stats::cov2cor(cov_mat)
  } else {
    corr_mat <- cov_mat
  }

  colnames(corr_mat) <- paste("v", 1:ncol(corr_mat), sep = "")

  # Return combined values
  ans <- cbind(rho_gamma, z, as_tibble(corr_mat))
  ans
}
