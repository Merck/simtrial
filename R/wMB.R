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
#' The weights for some early interval specified by the user are the inverse of the combined treatment group
#' empirical survival distribution; see details.
#' After this initial period, weights are constant at the maximum of the previous weights.
#' Another advantage of the test is that under strong null hypothesis that the underlying survival in the control group
#' is greater than or equal to underlying survival in the experimental group,
#' Type I error is controlled as the specified level.
#'
#' This function computes Magirr-Burman weights and adds them to a dataset created by the \code{tensurv()} function.
#' These weights can then be used to compute a Z-statistic for the modestly weighted logrank test proposed.
#'
#' @param x a \code{tensurv}-class \code{tibble} with a counting process dataset
#' @param delay The initial delay period where weights increase;
#' after this, weights are constant at the final weigh in the delay period
#' @return a vector with weights for the Magirr-Burman weighted logrank test
#' for the data in \code{x}
#' @details
#' We define \eqn{t^*} to be the input variable \code{delay}.
#' This specifies an initial period during which weights increase.
#' To define specific weights, we let \eqn{S(t)} denote the Kaplan-Meier survival estimate at time \eqn{t}
#' for the combined data (control plus experimental treatment groups).
#' The weight at time \eqn{t} is then defined as
#' \deqn{w(t)=S(\min(t,t^*))^{-1}.}
#'
#' @references
#' Magirr, Dominic, and Carl‐Fredrik Burman.
#' "Modestly weighted logrank tests."
#' \emph{Statistics in Medicine} 38.20 (2019): 3782-3790.
#'
#' @examples
#' library(tidyr)
#' library(dplyr)
#' # Use default enrollment and event rates at cut at 100 events
#' x <- simPWSurv(n=200) %>% cutDataAtCount(125) %>% tensurv(txval="Experimental")
#' # compute Magirr-Burman weights with
#' ZMB <- x %>% wMB(6) %>%
#'              summarize(S=sum(OminusE*wMB),V=sum(Var*wMB^2),Z=S/sqrt(V))
#' # Compute p-value of modestly weighted logrank of Magirr-Burman
#' pnorm(ZMB$Z)
#' @export
wMB <- function(x, delay = 4)
{
  # check input failure rate assumptions
  if(!is.data.frame(x)){stop("x in `wMB()` must be a data frame")}

  # check input delay
  if(!is.numeric(delay)){stop("delay in `wMB()` must be a non-negative number")}
  if(!delay >= 0){stop("delay in `wMB()` must be a non-negative number")}

  if(max(names(x)=="Stratum") != 1){stop("x column names in `wMB()` must contain Stratum")}
  if(max(names(x)=="tte") != 1){stop("x column names in `wMB()` must contain tte")}
  if(max(names(x)=="S") != 1){stop("x column names in `wMB()` must contain S")}


  # Compute max weight by stratum
  x2 <- x %>% group_by(Stratum)
  allstrat <- x2 %>% summarize()                        # Make sure you don't lose any strata!
  maxwgt <- x2 %>%
            filter(tte <= delay) %>%                    # look only up to delay time
            summarize(maxwgt = max(1/S)) %>%            # weight before delay specified as 1/S
            right_join(allstrat, by = "Stratum") %>%    # get back strata with no records before delay ends
            mutate(maxwgt = replace_na(maxwgt, 1)) %>%  # maxwgt is 1 when there are no records before delay ends
            full_join(x2,by="Stratum") %>%              # Now merge maxwgt back to stratified dataset
            mutate(wMB=pmin(maxwgt,1/S)) %>%            # Weight is min of maxwgt and 1/S which will increase up to delay
            select(-maxwgt)
  return(maxwgt)
}
