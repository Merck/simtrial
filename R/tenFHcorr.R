#  Copyright (c) 2022 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
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

#' @importFrom tibble tibble as_tibble
#' @import dplyr
NULL

#' @title Fleming-Harrington Weighted Logrank Tests plus Correlations
#'
#' @description
#' Correlations can be used with \code{mvtnorm::pmvnorm} to compute
#' p-value for MaxCombo, the maximum of the specifed
#' Fleming-Harrington tests
#'
#' @param x a \code{counting_process}-class \code{tibble} with a counting process dataset
#' @param rho_gamma a \code{tibble} with variables \code{rho} and \code{gamma}, both greater than equal
#' to zero, to specify one Fleming-Harrington weighted logrank test per row
#' @param corr a logical; if TRUE (default), return correlation matrix; otherwise, return covariance matrix
#'
#' @return a `tibble` with \code{rho_gamma} as input, the FH test statistics specified
#' for the data in \code{Z}, and the correlation or covariance matrix for these tests in variables starting
#' with \code{V}
#'
#' @examples
#' library(tidyr)
#' library(dplyr)
#'
#' # Use default enrollment and event rates at cut of 100 events
#' x <- sim_pw_surv(n = 200) %>%
#'   cut_data_by_event(100) %>%
#'   counting_process(arm = "Experimental")
#'
#' # compute logrank (FH(0,0)) and FH(0,1)
#' x <- x %>% tenFHcorr(rho_gamma = tibble(rho = c(0, 0),
#'                                  gamma = c(0, 1)))
#'
#' # compute p-value for MaxCombo
#' library(mvtnorm)
#' 1 - pmvnorm(lower = rep(min(x$Z), nrow(x)),
#'             corr = data.matrix(select(x, -c(rho, gamma, Z))),
#'             algorithm = GenzBretz(maxpts = 50000, abseps = 0.00001))[1]
#'
#' # check that covariance is as expected
#' x <- sim_pw_surv(n = 200) %>%
#'   cut_data_by_event(100) %>%
#'   counting_process(arm = "Experimental")
#'
#' x %>% tenFHcorr(rho_gamma = tibble(rho = c(0, 0),
#'                             gamma = c(0, 1)),
#'                 corr = FALSE)
#'
#' # Off-diagonal element should be variance in following
#' x %>% tenFHcorr(rho_gamma = tibble(rho = 0,
#'                             gamma =.5),
#'                 corr = FALSE)
#'
#' # compare off diagonal result with wlr()
#' x %>% wlr(rho_gamma = tibble(rho = 0, gamma =.5))
#'
#' @export
#' @rdname tenFHcorr
tenFHcorr <- function(x = sim_pw_surv(n = 200) %>%
                            cut_data_by_event(100) %>%
                            counting_process(arm = "Experimental"),
                      rho_gamma = tibble(rho = c(0, 0, 1, 1),
                                  gamma = c(0, 1, 0, 1)),
                      corr = TRUE
){

  n_weight <- nrow(rho_gamma)

  # Get average rho and gamma for FH covariance matrix
  # We want ave_rho[i,j] = (rho[i] + rho[j])/2
  # and     ave_gamma[i,j] = (gamma[i] + gamma[j])/2
  ave_rho <- (matrix(rho_gamma$rho, nrow = n_weight, ncol = n_weight, byrow = FALSE) +
                matrix(rho_gamma$rho, nrow = n_weight, ncol = n_weight, byrow = TRUE)
              )/2
  ave_gamma <- (matrix(rho_gamma$gamma, nrow = n_weight, ncol = n_weight) +
                  matrix(rho_gamma$gamma,nrow = n_weight,ncol = n_weight, byrow = TRUE)
               )/2

  # Convert back to tibble
  rg_new <- tibble(rho = as.numeric(ave_rho), gamma = as.numeric(ave_gamma))
  # get unique values of rho, gamma
  rg_unique <- rg_new %>% unique()

  # compute FH statistic for unique values
  # and merge back to full set of pairs
  rg_fh <- rg_new %>% left_join(wlr(x, rg_unique, returnVariance = TRUE),
                                by = c("rho" = "rho","gamma" = "gamma"))

  # get Z statistics for input rho, gamma combinations
  Z <- rg_fh$Z[(0:(n_weight - 1)) * n_weight + 1:n_weight]

  # get correlation matrix
  cov_mat <- matrix(rg_fh$Var, nrow = n_weight, byrow = TRUE)

  if (corr){
    corr_mat <- stats::cov2cor(cov_mat)
  } else{
    corr_mat <- cov_mat
  }

  names(corr_mat) <- paste("V", 1:ncol(corr_mat), sep = "")

  # return combined values
  ans <- cbind(rho_gamma, Z, as_tibble(corr_mat))
  return(ans)
}
