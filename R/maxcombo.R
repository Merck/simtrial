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
#' @param data a tte dataset
#' @param test1 maxcombo test1
#' @param test2 maxcombo test2
#' @param ... additional tests
#'
#' @return pvalues
#' @export
#'
#' @examples
#' sim_pw_surv(n = 200) |>
#'   cut_data_by_event(150) |>
#'   maxcombo(test1 = wlr(data, rho = 0, gamma = 0) |> quote(),
#'            test2 = wlr(data, rho = 0, gamma = 0.5) |> quote())
maxcombo <- function(data, test1, test2, ...){
  all_args <- match.call(expand.dots = FALSE)
  args <- all_args[-1]  # Exclude the first element (function name)

  n_test <- length(args) - 1
  rho_vector <- NULL
  gamma_vector <- NULL

  for (i in seq_len(n_test)) {
    test_i <- get(paste0("test", i))
    rho_vector <- c(rho_vector, test_i$rho)
    gamma_vector <- c(gamma_vector, test_i$gamma)
  }

  ans <- data |>
    counting_process(arm = "experimental") |>
    fh_weight(
      rho_gamma = data.frame(rho = rho_vector, gamma = gamma_vector),
      return_corr = TRUE)

  ans <-  data.frame(p_value = pvalue_maxcombo(ans))

  return(ans)
}
