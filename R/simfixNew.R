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
#' @import survival
NULL

#' Simulation of fixed sample size design for time-to-event endpoint
#'
#' `simfix()` provide simulations of a single endpoint two-arm trial
#' where the enrollment, hazard ratio, and failure and dropout rates change over time.
#' @param nsim Number of simulations to perform.
#' @param sampleSize Total sample size per simulation.
#' @param targetEvents Targeted event count for analysis.
#' @param strata A tibble with strata specified in `Stratum`, probability (incidence) of each stratum in `p`.
#' @param enrollRates Piecewise constant enrollment rates by time period.
#' Note that these are overall population enrollment rates and the `strata` argument controls the
#' random distribution between strata.
#' @param failRates Piecewise constant control group failure rates, hazard ratio for experimental vs control,
#'  and dropout rates by stratum and time period.
#' @param totalDuration Total follow-up from start of enrollment to data cutoff.
#' @param block As in `simtrial::simPWSurv()`. Vector of treatments to be included in each block.
#' @param timingType A numeric vector determining data cutoffs used; see details.
#' Default is to include all available cutoff methods.
#' @param rg As in `simtrial::tenFHCorr()`.
#' A \code{tibble} with variables \code{rho} and \code{gamma}, both greater than equal
#' to zero, to specify one Fleming-Harrington weighted logrank test per row.
#'
#' @details \code{timingType} has up to 5 elements indicating different options for data cutoff.
#' 1 uses the planned study duration, 2 the time the targeted event count is achieved,
#' 3 the planned minimum follow-up after enrollment is complete,
#' 4 the maximum of planned study duration and targeted event count cuts (1 and 2),
#' 5 the maximum of targeted event count and minimum follow-up cuts (2 and 3).
#' @return A \code{tibble} including columns \code{Events} (event count), \code{lnhr} (log-hazard ratio),
#' \code{Z} (normal test statistic; < 0 favors experimental) cut (text describing cutoff used),
#' \code{Duration} (duration of trial at cutoff for analysis) and \code{sim} (sequential simulation id).
#' One row per simulated dataset per cutoff specified in \code{timingType}, per test statistic specified.
#' If multiple Fleming-Harrington tests are specified in \code{rg}, then columns {rho,gamma}
#' are also included.
#' @examples
#' library(tidyr)
#' library(dplyr)
#' # Show output structure
#' simfix(nsim=3)
#' # Example with 2 tests: logrank and FH(0,1)
#' simfix(nsim=1,rg=tibble::tibble(rho=0,gamma=c(0,1)))
#' # Power by test
#' # Only use cuts for events, events + min follow-up
#' xx <- simfix(nsim=100,timingType=c(2,5),rg=tibble::tibble(rho=0,gamma=c(0,1)))
#' # Get power approximation for FH, data cutoff combination
#' xx %>% group_by(cut,rho,gamma) %>% summarise(mean(Z<=qnorm(.025)))
#' # MaxCombo power estimate for cutoff at max of targeted events, minimum follow-up
#' p <- xx %>%  filter(cut != "Targeted events") %>% group_by(Sim) %>% group_map(pMaxCombo)
#' p <- unlist(p)
#' mean(p<.025)
#' # MaxCombo estimate for targeted events cutoff
#' p <- unlist(xx %>%  filter(cut == "Targeted events") %>% group_by(Sim) %>% group_map(pMaxCombo))
#' mean(p<.025)
#' @export
simfixNew <- function(nsim=1000,
                   sampleSize=500, # sample size
                   targetEvents=350,  # targeted total event count
                   # multinomial probability distribution for strata enrollment
                   strata = tibble::tibble(Stratum = "All", p = 1),
                   # enrollment rates as in AHR()
                   enrollRates=tibble::tibble(duration=c(2,2,10),
                                              rate=c(3,6,9)),
                   # failure rates as in AHR()
                   failRates=tibble::tibble(Stratum="All",
                                            duration=c(3,100),
                                            failRate=log(2)/c(9,18),
                                            hr=c(.9,.6),
                                            dropoutRate=rep(.001,2)),
                   totalDuration=30, # total planned trial duration; single value
                   block=rep(c("Experimental","Control"),2), # Fixed block randomization specification
                   timingType=1:5, # select desired cutoffs for analysis (default is all types)
                   # default is to to logrank testing, but one or more Fleming-Harrington tests
                   # can be specified
                   rg=tibble::tibble(rho=0,gamma=0)
){# check input values
  # check input enrollment rate assumptions
  if(max(names(enrollRates)=="duration") != 1){stop("enrollRates column names in `simfix()` must contain duration")}
  if(max(names(enrollRates)=="rate") != 1){stop("enrollRates column names in `simfix()` must contain  rate")}

  # check input failure rate assumptions
  if(max(names(failRates)=="Stratum") != 1){stop("failRates column names in `simfix()` must contain Stratum")}
  if(max(names(failRates)=="duration") != 1){stop("failRates column names in `simfix()` must contain duration")}
  if(max(names(failRates)=="failRate") != 1){stop("failRates column names in `simfix()` must contain failRate")}
  if(max(names(failRates)=="hr") != 1){stop("failRates column names in `simfix()` must contain hr")}
  if(max(names(failRates)=="dropoutRate") != 1){stop("failRates column names in `simfix()` must contain dropoutRate")}

  # check input trial durations
  if(!is.numeric(totalDuration)){stop("totalDuration in `simfix()` must be a single positive number")}
  if(!is.vector(totalDuration)){stop("totalDuration in `simfix()` must be a single positive number")}
  if(length(totalDuration) != 1){stop("totalDuration in `simfix()` must be a single positive number")}
  if(!min(totalDuration) > 0){stop("totalDuration in `simfix()` must be a single positive number")}

  strata2 <- names(table(failRates$Stratum))
  if(nrow(strata)!= length(strata2)){stop("Stratum in `simfix()` must be the same in strata and failRates")}
  for(s in strata$Stratum){
    if(max(strata2==s) != 1){stop("Stratum in `simfix()` must be the same in strata and failRates")}
  }

  if(!nsim > 0){stop("nsim in `simfix()` must be positive integer")}
  if(length(nsim) != 1){stop("nsim in `simfix()` must be positive integer")}
  if(nsim != ceiling(nsim)){stop("nsim in `simfix()` must be positive integer")}

  if(!targetEvents > 0){stop("targetEvents in `simfix()` must be positive")}
  if(length(targetEvents) != 1){stop(("targetEvents in `simfix()` must be positive"))}

  if(!sampleSize > 0){stop("sampleSize in `simfix()` must be positive")}
  if(length(sampleSize) != 1){stop("sampleSize in `simfix()` must be positive")}

  nstrata <- nrow(strata)
  doAnalysis <- function(d,rg,nstrata){
    if (nrow(rg)==1){Z = tibble::tibble(Z=(d %>%
                                           simtrial::tensurv(txval="Experimental") %>%
                                           simtrial::tenFH(rg=rg)
    )$Z
    )
    }else Z = d %>%
              simtrial::tensurv(txval="Experimental") %>%
              simtrial::tenFHcorr(rg=rg,corr=TRUE)
    r <- cbind(tibble::tibble(Events = sum(d$event),
                              lnhr = as.numeric(ifelse(nstrata>1,
                                     survival::coxph(survival::Surv(tte,event)~
                                                       (Treatment=="Experimental")+
                                                       survival::strata(Stratum),data=d)$coefficients,
                                     survival::coxph(survival::Surv(tte,event)~
                                                       (Treatment=="Experimental"),data=d)$coefficients))
                             ),
                Z
    )
    r
  }
  # compute minimum planned follow-up time
  minFollow <- max(0,totalDuration - sum(enrollRates$duration))
  # put failure rates into simPWSurv format
  xx <- simfix2simPWSurv(failRates)
  fr <- xx$failRates
  dr <- xx$dropoutRates
  results <- NULL
  for(i in 1:nsim){
    sim <- simtrial::simPWSurvNew(n = sampleSize,
                               strata = strata,
                               enrollRates = enrollRates,
                               failRates = fr,
                               dropoutRates = dr,
                               block = block)
    # study date that targeted event rate achieved
    tedate <- sim %>% getCutDateForCount(targetEvents)
    # study data that targeted minimum follow-up achieved
    tmfdate <- max(sim$enrollTime) + minFollow
    # Compute tests for all specified cutoff options
    r1 <- NULL ; r2 <- NULL; r3 <- NULL;
    tests <- rep(FALSE,3)
    ## Total planned trial duration or max of this and targeted events
    if (max(1 == timingType)) tests[1] <- TRUE # planned duration cutoff
    if (max(2 == timingType)) tests[2] <- TRUE # targeted events cutoff
    if (max(3 == timingType)) tests[3] <- TRUE # minimum follow-up duration target
    if (max(4 == timingType)){ # max of planned duration, targeted events
      if (tedate > totalDuration){
        tests[2] <- TRUE
      }else tests[1] <- TRUE
    }
    if (max(5 == timingType)){ # max of minimum follow-up, targeted events
      if (tedate > tmfdate){
        tests[2] <- TRUE
      }else tests[3] <- TRUE
    }
    if (tests[1]){ # Total duration cutoff
      d <- sim %>% cutData(totalDuration)
      r1 <- d %>% doAnalysis(rg,nstrata)
    }
    if (tests[2]){ # targeted events cutoff
      d <- sim %>% cutData(tedate)
      r2 <- d %>% doAnalysis(rg,nstrata)
    }
    if (tests[3]){ # minimum follow-up cutoff
      d <- sim %>% cutData(tmfdate)
      r3 <- d %>% doAnalysis(rg,nstrata)
    }
    addit <- NULL
    # planned duration cutoff
    if (max(1 == timingType)) addit <- rbind(addit,
                                             r1 %>% mutate(cut="Planned duration",Duration=totalDuration))
    # targeted events cutoff
    if (max(2 == timingType)) addit <- rbind(addit, r2 %>% mutate(cut="Targeted events",Duration=tedate))
    # minimum follow-up duration target
    if (max(3 == timingType)) addit <- rbind(addit, r3 %>% mutate(cut="Minimum follow-up",Duration=tmfdate))
    if (max(4 == timingType)){ # max of planned duration, targeted events
      if (tedate > totalDuration){
        addit <- rbind(addit, r2 %>% mutate(cut="Max(planned duration, event cut)",Duration=tedate))
      }else addit <- rbind(addit, r1 %>% mutate(cut="Max(planned duration, event cut)",Duration=totalDuration))
    }
    if (max(5 == timingType)){ # max of minimum follow-up, targeted events
      if (tedate > tmfdate){
        addit <- rbind(addit, r2 %>% mutate(cut="Max(min follow-up, event cut)",Duration=tedate))
      }else addit <- rbind(addit, r3 %>% mutate(cut="Max(min follow-up, event cut)",Duration=tmfdate))
    }
    results <- rbind(results, addit %>% mutate(Sim=i))
  }
  results
}
