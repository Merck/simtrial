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
#' @param x A [counting_process()]-class data frame with a counting process
#'   dataset.
#' @param rho_gamma A data frame with variables `rho` and `gamma`, both greater
#'   than equal to zero, to specify one Fleming-Harrington weighted logrank test
#'   per row; Default: `data.frame(rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))`.
#' @param return_variance A logical flag that, if `TRUE`, adds columns
#'   estimated variance for weighted sum of observed minus expected;
#'   see details; Default: `FALSE`.
#' @param return_corr A logical flag that, if `TRUE`, adds columns
#'   estimated correlation for weighted sum of observed minus expected;
#'   see details; Default: `FALSE`.
#'
#' @return A data frame with `rho_gamma` as input and the FH test statistic
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
#' @importFrom data.table data.table merge.data.table setDF
#'
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' # Example 1
#' # Use default enrollment and event rates at cut at 100 events
#' x <- sim_pw_surv(n = 200) |>
#'   cut_data_by_event(100) |>
#'   counting_process(arm = "experimental")
#'
#' # Compute logrank FH(0, 1)
#' fh_weight(x, rho_gamma = data.frame(rho = 0, gamma = 1))
#' fh_weight(x, rho_gamma = data.frame(rho = 0, gamma = 1), return_variance = TRUE)
#'
#' # Compute the corvariance between FH(0, 0), FH(0, 1) and FH(1, 0)
#' fh_weight(x, rho_gamma = data.frame(rho = c(0, 0, 1), gamma = c(0, 1, 0)))
#' fh_weight(x, rho_gamma = data.frame(rho = c(0, 0, 1), gamma = c(0, 1, 0)), return_variance = TRUE)
#' fh_weight(x, rho_gamma = data.frame(rho = c(0, 0, 1), gamma = c(0, 1, 0)), return_corr = TRUE)
#'
#' # Example 2
#' # Use default enrollment and event rates at cut of 100 events
#' set.seed(123)
#' x <- sim_pw_surv(n = 200) |>
#'   cut_data_by_event(100) |>
#'   counting_process(arm = "experimental") |>
#'   fh_weight(rho_gamma = data.frame(rho = c(0, 0), gamma = c(0, 1)), return_corr = TRUE)
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
#' x <- sim_pw_surv(n = 200) |>
#'   cut_data_by_event(100) |>
#'   counting_process(arm = "experimental")
#'
#' x |> fh_weight(
#'   rho_gamma = data.frame(
#'     rho = c(0, 0),
#'     gamma = c(0, 1)
#'   ),
#'   return_variance = TRUE
#' )
#'
#' # Off-diagonal element should be variance in following
#' x |> fh_weight(
#'   rho_gamma = data.frame(
#'     rho = 0,
#'     gamma = .5
#'   ),
#'   return_variance = TRUE
#' )
#'
#' # Compare off diagonal result with fh_weight()
#' x |> fh_weight(rho_gamma = data.frame(rho = 0, gamma = .5))
fh_weight <- function(
    x = sim_pw_surv(n = 200) |>
      cut_data_by_event(150) |>
      counting_process(arm = "experimental"),
    rho_gamma = data.frame(
      rho = c(0, 0, 1, 1),
      gamma = c(0, 1, 0, 1)
    ),
    return_variance = FALSE,
    return_corr = FALSE) {
  n_weight <- nrow(rho_gamma)

  # Check input failure rate assumptions
  if (!is.data.frame(x)) {
    stop("fh_weight: x in `fh_weight()` must be a data frame.")
  }

  if (!("s" %in% names(x))) {
    stop("fh_weight: x column names in `fh_weight()` must contain s.")
  }

  if (!("o_minus_e" %in% names(x))) {
    stop("fh_weight: x column names in `fh_weight()` must contain o_minus_e.")
  }

  if (!("var_o_minus_e" %in% names(x))) {
    stop("fh_weight: x column names in `fh_weight()` must contain var_o_minus_e.")
  }

  if (return_variance && return_corr) {
    stop("fh_weight: can't report both covariance and correlation for MaxCombo test.")
  }

  if (return_corr && n_weight == 1) {
    stop("fh_weight: can't report the correlation for a single WLR test.")
  }

  if (n_weight == 1) {
    ans <- wlr_z_stat(x, rho_gamma = rho_gamma, return_variance = return_variance)
  } else {
    # Get average rho and gamma for FH covariance matrix
    # We want ave_rho[i,j]   = (rho[i] + rho[j])/2
    # and     ave_gamma[i,j] = (gamma[i] + gamma[j])/2
    ave_rho <- (matrix(rho_gamma$rho, nrow = n_weight, ncol = n_weight, byrow = FALSE) +
      matrix(rho_gamma$rho, nrow = n_weight, ncol = n_weight, byrow = TRUE)
    ) / 2
    ave_gamma <- (matrix(rho_gamma$gamma, nrow = n_weight, ncol = n_weight) +
      matrix(rho_gamma$gamma, nrow = n_weight, ncol = n_weight, byrow = TRUE)
    ) / 2

    # Convert to data.table
    rg_new <- data.table(rho = as.numeric(ave_rho), gamma = as.numeric(ave_gamma))
    # Get unique values of rho, gamma
    rg_unique <- unique(rg_new)

    # Compute FH statistic for unique values
    # and merge back to full set of pairs
    rg_fh <- merge.data.table(
      x = rg_new,
      y = wlr_z_stat(x, rho_gamma = rg_unique, return_variance = TRUE),
      by = c("rho", "gamma"),
      all.x = TRUE,
      sort = FALSE
    )

    # Get z statistics for input rho, gamma combinations
    z <- rg_fh$z[(0:(n_weight - 1)) * n_weight + 1:n_weight]

    # Get correlation matrix
    cov_mat <- matrix(rg_fh$var, nrow = n_weight, byrow = TRUE)

    if (return_corr) {
      corr_mat <- stats::cov2cor(cov_mat)
      colnames(corr_mat) <- paste("v", 1:ncol(corr_mat), sep = "")
      ans <- cbind(rho_gamma, z, as.data.frame(corr_mat))
    } else if (return_variance) {
      corr_mat <- cov_mat
      colnames(corr_mat) <- paste("v", 1:ncol(corr_mat), sep = "")
      ans <- cbind(rho_gamma, z, as.data.frame(corr_mat))
    } else if (return_corr + return_corr == 0) {
      corr_mat <- NULL
      ans <- cbind(rho_gamma, z)
    }
  }

  setDF(ans)
  return(ans)
}

# Build an internal function to compute the Z statistics
# under a sequence of rho and gamma of WLR.
wlr_z_stat <- function(x, rho_gamma, return_variance) {
  ans <- rho_gamma

  xx <- x[, c("s", "o_minus_e", "var_o_minus_e")]

  ans$z <- rep(0, nrow(rho_gamma))

  if (return_variance) {
    ans$var <- rep(0, nrow(rho_gamma))
  }

  for (i in 1:nrow(rho_gamma)) {
    weight <- xx$s^rho_gamma$rho[i] * (1 - xx$s)^rho_gamma$gamma[i]
    weighted_o_minus_e <- weight * xx$o_minus_e
    weighted_var <- weight^2 * xx$var_o_minus_e

    weighted_var_total <- sum(weighted_var)
    weighted_o_minus_e_total <- sum(weighted_o_minus_e)


    ans$z[i] <- weighted_o_minus_e_total / sqrt(weighted_var_total)

    if (return_variance) {
      ans$var[i] <- weighted_var_total
    }
  }

  return(ans)
}
