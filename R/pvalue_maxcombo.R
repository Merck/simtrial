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

#' @import dplyr
#' @import tibble
#' @import mvtnorm
NULL

#' MaxCombo p-value
#'
#' \code{pMaxCombo()} computes p-values for the MaxCombo test
#' based on output from \code{simtrial::tenFHcorr()}.
#' This is still in an experimental stage and is intended for use with
#' the \code{simtrial::simfix()} trial simulation routine.
#' However, it can also be used to analyze clinical trial data such as that provided in the
#' ADaM ADTTE format.
#' @param Z a dataset output from \code{tenFHcorr()}; see examples.
#' @param dummyvar a dummy input that allows \code{group_map()} to be used to
#' compute p-values for multiple simulations.
#' @param algorithm This is passed directly to the \code{algorithm} argument in the \code{mvtnorm::pmvnorm()}
#'
#' @return A numeric p-value
#'
#' @examples
#' library(tidyr)
#' library(tibble)
#' library(dplyr)
#'
#' # example 1
#' x <- simfix(nsim = 1,
#'             timingType = 5,
#'             rg = tibble(rho = c(0, 0, 1),
#'                         gamma = c(0, 1, 1)))
#' head(x)
#' pMaxCombo(x)
#'
#' # example 2
#' # Only use cuts for events, events + min follow-up
#' xx <- simfix(nsim = 100,
#'              timingType = 5,
#'              rg = tibble(rho = c(0, 0, 1),
#'                          gamma = c(0, 1, 1)))
#' head(xx)
#' # MaxCombo power estimate for cutoff at max of targeted events, minimum follow-up
#' p <- xx %>% group_by(Sim) %>% group_map(pMaxCombo) %>% unlist()
#' mean(p < .025)
#'
#' @export
#'
pvalue_maxcombo<- function(Z,
                      dummyvar,
                      algorithm = GenzBretz(maxpts = 50000, abseps = 0.00001)){

  ans <- (1 - mvtnorm::pmvnorm(lower = rep(Z$Z %>% min() %>% as.numeric(),
                                           nrow(Z)),
                               corr = Z %>%
                                        select(starts_with("V")) %>%
                                        data.matrix(),
                               algorithm = algorithm)[1]
          ) %>% as.numeric()

  return(ans)
}
