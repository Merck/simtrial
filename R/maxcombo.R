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

#' Maxcombo test
#'
#' WARNING: This experimental function is a work-in-progress. The function
#' arguments will change as we add additional features.
#'
#' @param data a tte dataset
#' @param rho Numeric vector passed to [fh_weight()]. Must be greater
#'   than or equal to zero. Must be the same length as `gamma`.
#' @param gamma Numeric vector passed to [fh_weight()]. Must be
#'   greater than or equal to zero. Must be the same length as `rho`.
#'
#' @return pvalues
#' @export
#'
#' @seealso [fh_weight()]
#'
#' @examples
#' sim_pw_surv(n = 200) |>
#'   cut_data_by_event(150) |>
#'   maxcombo(rho = c(0, 0), gamma = c(0, 1))
maxcombo <- function(data = sim_pw_surv(n = 200) |> cut_data_by_event(150),
                     rho = c(0, 0, 1),
                     gamma = c(0, 1, 1),
                     return_variance = FALSE,
                     return_corr = FALSE) {
  # Input checking ----
  n_weight <- length(rho)

  if(any(!is.numeric(rho)) || any(!is.numeric(gamma)) || any(rho < 0) || any(gamma < 0)) {
    stop("maxcombo: please input positive values of rho and gamma.")
  }

  if (length(rho) != length(gamma)) {
    stop("maxcombo: please input rho and gamma of equal length.")
  }

  if (length(rho) == 1) {
    stop("maxcombo: please input multiple rho and gamma to make it a MaxCombo test.")
  }

  if (return_variance && return_corr) {
    stop("maxcombo: can't report both covariance and correlation for MaxCombo test.")
  }

  if (return_corr && n_weight == 1) {
    stop("maxcombo: can't report the correlation for a single WLR test.")
  }

  # Get average rho and gamma for FH covariance matrix ----
  # We want ave_rho[i,j]   = (rho[i] + rho[j])/2
  # and     ave_gamma[i,j] = (gamma[i] + gamma[j])/2
  ave_rho <- (matrix(rho, nrow = n_weight, ncol = n_weight, byrow = FALSE) +
                matrix(rho, nrow = n_weight, ncol = n_weight, byrow = TRUE)
  ) / 2
  ave_gamma <- (matrix(gamma, nrow = n_weight, ncol = n_weight) +
                  matrix(gamma, nrow = n_weight, ncol = n_weight, byrow = TRUE)
  ) / 2

  # Convert to all original and average rho/gamma into a data.table ----
  rg_new <- data.table(rho = as.numeric(ave_rho), gamma = as.numeric(ave_gamma))
  rg_unique <- unique(rg_new)

  # Compute WLR FH Z-score and variance for all original and average rho/gamma ----
  # Compute FH statistic for unique values
  res <- list()
  x <- data |> counting_process(arm = "experimental")
  # Loop for each WLR-FH test with 1 rho and 1 gamma
  for (i in seq_len(nrow(rg_unique))) {
    weight <- x$s^rg_unique$rho[i] * (1 - x$s)^rg_unique$gamma[i]
    weighted_o_minus_e <- weight * x$o_minus_e
    weighted_var <- weight^2 * x$var_o_minus_e

    weighted_var_total <- sum(weighted_var)
    weighted_o_minus_e_total <- sum(weighted_o_minus_e)

    res$rho[i] <- rg_unique$rho[i]
    res$gamma[i] <- rg_unique$gamma[i]
    res$estimation[i] <- weighted_o_minus_e_total
    res$se[i] <- sqrt(weighted_var_total)
    res$var[i] <- weighted_var_total
    res$z[i] <- res$estimation[i] / res$se[i]
  }

  # Merge back to full set of pairs ----
  rg_fh <- merge.data.table(
    x = rg_new,
    y = res |> as.data.frame(),
    by = c("rho", "gamma"),
    all.x = TRUE,
    sort = FALSE)

  # Tidy outputs ----
  ans <- list()
  ans$method <- "Maxcombo"
  temp <- data.frame(rho = rho, gamma = gamma) %>% mutate(x= paste0("FH(", rho, ", ", gamma, ")"))
  ans$parameter <- paste(temp$x, collapse = " + ")
  ans$estimation <- NULL
  ans$se <- NULL

  # Get z statistics for input rho, gamma combinations
  ans$z <- rg_fh$z[(0:(n_weight - 1)) * n_weight + 1:n_weight]

  # Get correlation matrix
  cov_mat <- matrix(rg_fh$var, nrow = n_weight, byrow = TRUE)
  corr_mat <- stats::cov2cor(cov_mat)

  cov_mat <- as.data.frame(cov_mat)
  corr_mat <- as.data.frame(corr_mat)
  colnames(cov_mat) <- paste("v", seq_len(ncol(cov_mat)), sep = "")
  colnames(corr_mat) <- paste("v", seq_len(ncol(corr_mat)), sep = "")

  if (return_corr) {
    ans$corr <- corr_mat
  } else if (return_variance) {
    ans$var <- cov_mat
  }

  # Get p-values
  res_df <- as.data.frame(cbind(ans$z, corr_mat))
  colnames(res_df)[1] <- "z"
  ans$p_value <- pvalue_maxcombo(res_df)

  return(ans)
}
