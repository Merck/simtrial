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

#' Simulate a stratified time-to-event outcome randomized trial
#'
#' \code{simPWSurv} enables simulation of a clinical trial with essentially arbitrary
#' patterns of enrollment, failure rates and censoring.
#' The piecewise exponential distribution allows a simple method to specify a distribtuion
#' and enrollment pattern
#' where the enrollment, failure and dropout rate changes over time.
#' While the main purpose may be to generate a trial that can be analyzed at a single point in time or
#' using group sequential methods, the routine can also be used to simulate an adaptive trial design.
#' Enrollment, failure and dropout rates are specified by treatment group, stratum and time period.
#' Fixed block randomization is used; blocks must include treatments provided in failure and dropout
#' specification.
#' Default arguments are set up to allow very simple implementation of a non-proportional hazards assumption
#' for an unstratified design.
#'
#' @param n Number of observations.
#' If length(n) > 1, the length is taken to be the number required.
#' @param enrollStrata A tibble with strata specified in `Stratum`, probability (incidence) of each stratum
#' in `p`
#' @param block Vector of treatments to be included in each  block
#' @param enrollRates Enrollment rates; see details and examples
#' @param failRates Failure rates; see details and examples; note that treatments need
#' to be the same as input in block
#' @param dropoutRates Dropout rates; see details and examples; note that treatments need
#' to be the same as input in block
#'
#' @return a \code{tibble} with the following variables for each observation
#' \code{Stratum},
#' \code{enrollTime} (enrollment time for the observation),
#' \code{Treatment} (treatment group; this will be one of the values in the input \code{block}),
#' \code{failTime} (failure time generated using \code{rpwexp()}),
#' \code{dropoutTime} (dropout time generated using \code{rpwexp()}),
#' \code{cte} (calendar time of enrollment plot the minimum of failure time and dropout time),
#' \code{fail} (indicator that \code{cte} was set using failure time; i.e., 1 is a failure, 0 is a dropout).
#' @examples
#' library(dplyr)
#' # Tests
#'  simPWSurv(n=20)
#'  # 3:1 randomization
#'  simPWSurv(n=20,block=c(rep("Experimental",3),"Control"))
#'
#' # Simulate 2 strata; will use defaults for blocking and enrollRates
#' simPWSurv(n=20,
#'           # 2 strata,30% and 70% prevalence
#'           enrollStrata=tibble::tibble(Stratum=c("Low","High"),p=c(.3,.7)),
#'           failRates=tibble::tibble(Stratum=c(rep("Low",4),rep("High",4)),
#'                                    period=rep(1:2,4),
#'                                    Treatment=rep(c(rep("Control",2),rep("Experimental",2)),2),
#'                                    duration=rep(c(3,1),4),
#'                                    rate=c(.03,.05,.03,.03,.05,.08,.07,.04)),
#'           dropoutRates=tibble::tibble(Stratum=c(rep("Low",2),rep("High",2)),
#'                                       period=rep(1,4),
#'                                       Treatment=rep(c("Control","Experimental"),2),
#'                                       duration=rep(1,4),
#'                                       rate=rep(.001,4))
#')
#'
#'# If you want a more rectangular entry for a tibble
#'failRates <- bind_rows(
#'    tibble(Stratum="Low" ,period=1,Treatment="Control"     ,duration=3,rate=.03),
#'    tibble(Stratum="Low" ,period=1,Treatment="Experimental",duration=3,rate=.03),
#'    tibble(Stratum="Low" ,period=2,Treatment="Experimental",duration=3,rate=.02),
#'    tibble(Stratum="High",period=1,Treatment="Control"     ,duration=3,rate=.05),
#'    tibble(Stratum="High",period=1,Treatment="Experimental",duration=3,rate=.06),
#'    tibble(Stratum="High",period=2,Treatment="Experimental",duration=3,rate=.03)
#')
#'dropoutRates <- bind_rows(
#'    tibble(Stratum="Low" ,period=1,Treatment="Control"     ,duration=3,rate=.001),
#'    tibble(Stratum="Low" ,period=1,Treatment="Experimental",duration=3,rate=.001),
#'    tibble(Stratum="High",period=1,Treatment="Control"     ,duration=3,rate=.001),
#'    tibble(Stratum="High",period=1,Treatment="Experimental",duration=3,rate=.001)
#')
#'simPWSurv(n=12,enrollStrata=tibble(Stratum=c("Low","High"),p=c(.3,.7)),
#'          failRates=failRates,dropoutRates=dropoutRates)
#' @export
simPWSurv <- function(n=100,
                      enrollStrata=tibble::tibble(Stratum="All",p=1),
                      block=c(rep("Control",2),rep("Experimental",2)),
                      enrollRates=tibble::tibble(rate=9,
                                                 duration=1),
                      failRates=tibble::tibble(Stratum=rep("All",4),
                                               period=rep(1:2,2),
                                               Treatment=c(rep("Control",2), rep("Experimental",2)),
                                               duration=rep(c(3,1),2),
                                               rate=log(2)/c(9,9,9,18)),
                      dropoutRates=tibble::tibble(Stratum=rep("All",2),
                                              period=rep(1,2),
                                              Treatment=c("Control","Experimental"),
                                              duration=rep(100,2),
                                              rate=rep(.001,2))
                      ){
# start tibble by generating strata and enrollment times
    #return(
    x<-  tibble::tibble(Stratum=sample(x=enrollStrata$Stratum,size=n,replace=TRUE,prob=enrollStrata$p)) %>%
         mutate(enrollTime=rpwenroll(n, enrollRates)) %>%
         group_by(Stratum) %>% mutate(Treatment=fixedBlockRand(n=n(),block=block))  %>% # assign treatment
         # generate time to failure and time to dropout
         dplyr::group_by(Stratum,Treatment)
    utr <- unique(x$Treatment)
    usr <- unique(x$Stratum)
    x$failTime <- 0
    x$dropoutTime <- 0
    for(sr in usr){for(tr in utr){
      indx <- x$Stratum==sr & x$Treatment==tr
      x$failTime[indx] <- rpwexp(n=sum(indx),failRates=filter(failRates,Stratum==sr&Treatment==tr))
      x$dropoutTime[indx] <- rpwexp(n=sum(indx),failRates=filter(dropoutRates,Stratum==sr&Treatment==tr))
    }}
    # set calendar time-to-event and failure indicator
    return(x %>% mutate(cte=pmin(dropoutTime,failTime)+enrollTime,
                        fail=(failTime <= dropoutTime)*1))
}
