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

#' @import tibble
#' @import dplyr
NULL

#' Magirr and Burman Modestly Weighted Logrank Tests
#'
#' Magirr and Burman (2019) have proposed a weighted logrank test to have better power than
#' the logrank test when the treatment effect is delayed, but to still maintain good power under
#' a proportional hazards assumption.
#' In Magirr (2021), (the equivalent of) a maximum weight was proposed as opposed to a fixed
#' time duration over which weights would increase.
#' The weights for some early interval specified by the user are the inverse of the combined treatment group
#' empirical survival distribution; see details.
#' After this initial period, weights are constant at the maximum of the previous weights.
#' Another advantage of the test is that under strong null hypothesis that the underlying survival in the control group
#' is greater than or equal to underlying survival in the experimental group,
#' Type I error is controlled as the specified level.
#'
#' This function computes Magirr-Burman weights and adds them to a dataset created by the \code{counting_process()} function.
#' These weights can then be used to compute a Z-statistic for the modestly weighted logrank test proposed.
#'
#' @param x a \code{counting_process}-class \code{tibble} with a counting process dataset
#' @param delay The initial delay period where weights increase;
#' after this, weights are constant at the final weight in the delay period
#' @param wmax Maximum weight to be returned.
#' set \code{delay = Inf, wmax = 2} to be consistent with recommendation of
#' Magirr (2021).
#'
#' @return a vector with weights for the Magirr-Burman weighted logrank test
#' for the data in \code{x}
#'
#' @details
#' We define \eqn{t^*} to be the input variable \code{delay}.
#' This specifies an initial period during which weights increase.
#' We also set a maximum weight \eqn{w_max}.
#' To define specific weights, we let \eqn{S(t)} denote the Kaplan-Meier survival estimate at time \eqn{t}
#' for the combined data (control plus experimental treatment groups).
#' The weight at time \eqn{t} is then defined as
#' \deqn{w(t)=\min(w_max, S(\min(t,t^*))^{-1}).}
#'
#' @references
#' Magirr, Dominic, and Carl‐Fredrik Burman.
#' "Modestly weighted logrank tests."
#' \emph{Statistics in Medicine} 38.20 (2019): 3782--3790.
#'
#' Magirr, Dominic.
#' "Non‐proportional hazards in immuno‐oncology: Is an old perspective needed?"
#' \emph{Pharmaceutical Statistics} 20.3 (2021): 512--527.
#'
#' @examples
#' library(tidyr)
#' library(dplyr)
#'
#' # Use default enrollment and event rates at cut at 100 events
#' # For transparency, may be good to set either `delay` or `wmax` to Inf`
#' x <- sim_pw_surv(n = 200) %>%
#'   cut_data_by_event(125) %>%
#'   counting_process(arm = "Experimental")
#'
#' # example 1
#' # compute Magirr-Burman weights with `delay = 6`
#' ZMB <- x %>%
#'   mb_weight(delay = 6, wmax = Inf) %>%
#'   summarize(S = sum(o_minus_e * mb_weight),
#'             V = sum(var_o_minus_e * mb_weight^2),
#'             Z = S / sqrt(V))
#'
#' # Compute p-value of modestly weighted logrank of Magirr-Burman
#' pnorm(ZMB$Z)
#'
#' # example 2
#' # Now compute with maximum weight of 2 as recommended in Magirr, 2021
#' ZMB2 <- x %>%
#'   mb_weight(delay = Inf, wmax = 2) %>%
#'   summarize(S = sum(o_minus_e * mb_weight),
#'             V = sum(var_o_minus_e * mb_weight^2),
#'             Z = S / sqrt(V))
#'
#' # Compute p-value of modestly weighted logrank of Magirr-Burman
#' pnorm(ZMB2$Z)
#'
#' @export
mb_weight <- function(x, delay = 4, wmax = Inf)
{
  # check input failure rate assumptions
  if(!is.data.frame(x)){
    stop("x in `mb_weight()` must be a data frame")
  }

  # check input delay
  if(!is.numeric(delay)){
    stop("delay in `mb_weight()` must be a non-negative number")
  }

  if(!delay >= 0){
    stop("delay in `mb_weight()` must be a non-negative number")
  }

  if(!is.numeric(wmax)){
    stop("wmax (maximum weight) in `mb_weight()` must be a positive number")
  }

  if(!delay > 0){
    stop("wmax (maximum weight) in `mb_weight()` must be a positive number")
  }

  if(max(names(x)=="stratum") != 1){
    stop("x column names in `mb_weight()` must contain stratum")
  }

  if(max(names(x)=="tte") != 1){
    stop("x column names in `mb_weight()` must contain tte")
  }

  if(max(names(x)=="S") != 1){
    stop("x column names in `mb_weight()` must contain S")
  }

  # Compute max weight by stratum
  x2 <- x %>% group_by(stratum)
  # Make sure you don't lose any stratum!
  tbl_all_stratum <- x2 %>% summarize()

  ans <- x2 %>%
    # look only up to delay time
    filter(tte <= delay) %>%
    # weight before delay specified as 1/S
    summarize(max_weight = max(1/S)) %>%
    # get back stratum with no records before delay ends
    right_join(tbl_all_stratum, by = "stratum") %>%
    # max_weight is 1 when there are no records before delay ends
    mutate(max_weight = tidyr::replace_na(max_weight, 1)) %>%
    # Cut off weights at wmax
    mutate(max_weight = pmin(wmax, max_weight)) %>%
    # Now merge max_weight back to stratified dataset
    full_join(x2, by = "stratum") %>%
    # Weight is min of max_weight and 1/S which will increase up to delay
    mutate(mb_weight = pmin(max_weight, 1/S)) %>%
    select(-max_weight)

  return(ans)
}
