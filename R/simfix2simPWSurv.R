#  Copyright (c) 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
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
NULL

#' Conversion of enrollment and failure rates from simfix() to simPWSurv() format
#'
#' `simfix2simPWSurv()` converts failure rates and dropout rates entered in the simpler
#' format for `simfix()` to that used for `simtrial::simPWSurv()`.
#' The `failRates` argument for `simfix()` requires enrollment rates, failure rates
#' hazard ratios and dropout rates by strata for a 2-arm trial, `simtrial::simPWSurv()`
#' is in a more flexible but less obvious but more flexible format.
#' Since `simfix()` automatically analyzes data and `simtrial::simPWSurv()` just produces
#' a simulation dataset, the latter provides additional options to analyze or otherwise evaluate
#' individual simulations in ways that `simfix()` does not.
#' @param failRates Piecewise constant control group failure rates, hazard ratio for experimental vs control,
#'  and dropout rates by stratum and time period.
#' @return A \code{list} of two `tibble` components formatted for `simtrial::simPWSurv()`:
#' `failRates` and `dropoutRates`.
#' @examples
#' library(tidyr)
#' library(dplyr)
#' # Convert standard input
#' simfix2simPWSurv()
#' # Stratified example
#' failRates <- tibble::tibble(Stratum=c(rep("Low",3),rep("High",3)),
#'                             duration=rep(c(4,10,100),2),
#'                             failRate=c(.04,.1,.06,
#'                                        .08,.16,.12),
#'                             hr=c(1.5,.5,2/3,
#'                                  2, 10/16, 10/12),
#'                             dropoutRate=.01
#' )
#' x <- simfix2simPWSurv(failRates)
#' # Do a single simulation with the above rates
#' # Enroll 300 patients in ~12 months at constant rate
#' sim <-
#'     simPWSurv(n=300,
#'           enrollStrata=tibble::tibble(Stratum=c("Low","High"),p=c(.6,.4)),
#'           enrollRates=tibble::tibble(duration=12,rate=300/12),
#'           failRates=x$failRates,
#'           dropoutRates=x$dropoutRates)
#' # Cut after 200 events and do a stratified logrank test
#' dat <- sim %>%
#'        cutDataAtCount(200) %>%            # cut data
#'        tensurv(txval="Experimental") %>%  # convert format for tenFH
#'        tenFH(rg=tibble(rho=0,gamma=0))    # stratified logrank
#' @export
simfix2simPWSurv <-
  function(# failure rates as in simfix()
  failRates=tibble::tibble(Stratum="All",
                           duration=c(3,100),
                           failRate=log(2)/c(9,18),
                           hr=c(.9,.6),
                           dropoutRate=rep(.001,2))
){# put failure rates into simPWSurv format
  fr <- rbind(failRates %>% group_by(Stratum) %>% mutate(Treatment="Control",rate=failRate,period=1:n()) %>%
                ungroup(),
              failRates %>% group_by(Stratum) %>% mutate(Treatment="Experimental",rate=failRate*hr,period=1:n()) %>%
                ungroup()
  ) %>% select("Stratum","period","Treatment","duration","rate")
  # put dropout rates into simPWSurv format
  dr <-  failRates %>%  group_by(Stratum) %>% mutate(Treatment="Control",rate=dropoutRate,period=1:n()) %>%
    select("Stratum","period","Treatment","duration","rate") %>%
    ungroup()
  dr <- rbind(dr, dr %>% mutate(Treatment="Experimental"))
  return(list(failRates=fr,dropoutRates=dr))
}
