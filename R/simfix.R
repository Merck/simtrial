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
#' @importFrom tibble tibble
#' @import survival
#' @import foreach
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
#' @param seed Optional. Initial seed for simulations
#'
#' @details \code{timingType} has up to 5 elements indicating different options for data cutoff.
#' 1 uses the planned study duration, 2 the time the targeted event count is achieved,
#' 3 the planned minimum follow-up after enrollment is complete,
#' 4 the maximum of planned study duration and targeted event count cuts (1 and 2),
#' 5 the maximum of targeted event count and minimum follow-up cuts (2 and 3).
#'
#' @return A \code{tibble} including columns \code{Events} (event count), \code{lnhr} (log-hazard ratio),
#' \code{Z} (normal test statistic; < 0 favors experimental) cut (text describing cutoff used),
#' \code{Duration} (duration of trial at cutoff for analysis) and \code{sim} (sequential simulation id).
#' One row per simulated dataset per cutoff specified in \code{timingType}, per test statistic specified.
#' If multiple Fleming-Harrington tests are specified in \code{rg}, then columns {rho,gamma}
#' are also included.
#'
#' @examples
#' library(tidyr)
#' library(dplyr)
#' library(doParallel)
#' library(tibble)
#'
#' # example 1
#' # Show output structure
#' simfix(nsim = 3)
#'
#' # example 2
#' # Example with 2 tests: logrank and FH(0,1)
#' simfix(nsim = 1,rg = tibble(rho = 0, gamma = c(0, 1)))
#'
#' # example 3
#' # Power by test
#' # Only use cuts for events, events + min follow-up
#' xx <- simfix(nsim = 100,
#'              timingType = c(2, 5),
#'              rg = tibble(rho = 0, gamma = c(0, 1)))
#' # Get power approximation for FH, data cutoff combination
#' xx %>%
#'   group_by(cut, rho, gamma) %>%
#'   summarise(mean(Z <= qnorm(.025)))
#'
#' # MaxCombo power estimate for cutoff at max of targeted events, minimum follow-up
#' p <- xx %>%
#'   filter(cut != "Targeted events") %>%
#'   group_by(Sim) %>%
#'   group_map(value_maxcombo) %>%
#'   unlist()
#'
#' mean(p < .025)
#'
#' # MaxCombo estimate for targeted events cutoff
#' p <- xx %>%
#'   filter(cut == "Targeted events") %>%
#'   group_by(Sim) %>%
#'   group_map(value_maxcombo) %>%
#'   unlist()
#'
#' mean(p < .025)
#'
#' # example 3
#' # Use two cores
#' registerDoParallel(2)
#' simfix(nsim = 10, seed = 2022)
#' stopImplicitCluster()
#' registerDoSEQ()
#'
#' @export
#'
simfix <- function(nsim = 1000,
                   sampleSize = 500, # sample size
                   targetEvents = 350,  # targeted total event count
                   # multinomial probability distribution for strata enrollment
                   strata = tibble(Stratum = "All", p = 1),
                   # enrollment rates as in AHR()
                   enrollRates = tibble(duration = c(2, 2, 10), rate = c(3, 6, 9)),
                   # failure rates as in AHR()
                   failRates = tibble(Stratum = "All",
                                      duration = c(3, 100),
                                      failRate = log(2) / c(9, 18),
                                      hr = c(.9, .6),
                                      dropoutRate = rep(.001, 2)),
                   # total planned trial duration; single value
                   totalDuration = 30,
                   # Fixed block randomization specification
                   block = rep(c("Experimental", "Control"), 2),
                   # select desired cutoffs for analysis (default is all types)
                   timingType = 1:5,
                   # default is to to logrank testing, but one or more Fleming-Harrington tests
                   # can be specified
                   rg = tibble(rho = 0, gamma = 0),
                   seed = NULL
                   ){
  # check input values
  # check input enrollment rate assumptions
  if(!("duration" %in% names(enrollRates)) ){
    stop("simfix: enrollRates column names in `simfix()` must contain duration!")
  }

  if(!("rate" %in% names(enrollRates))){
    stop("simfix: enrollRates column names in `simfix()` must contain  rate!")
  }


  # check input failure rate assumptions
  if(!("Stratum" %in% names(failRates))){
    stop("simfix: failRates column names in `simfix()` must contain Stratum!")
  }

  if(!("duration" %in% names(failRates))){
    stop("simfix: failRates column names in `simfix()` must contain duration!")
  }

  if(!("failRate" %in% names(failRates))){
    stop("simfix: failRates column names in `simfix()` must contain failRate!")
  }

  if(!("hr" %in% names(failRates))){
    stop("simfix: failRates column names in `simfix()` must contain hr")
  }

  if(!("dropoutRate" %in% names(failRates))){
    stop("simfix: failRates column names in `simfix()` must contain dropoutRate")
  }

  # check input trial duration
  if(!is.numeric(totalDuration)){
    stop("simfix: totalDuration in `simfix()` must be a single positive number!")
  }

  if(!is.vector(totalDuration)){
    stop("simfix: totalDuration in `simfix()` must be a single positive number!")
  }

  if(length(totalDuration) != 1){
    stop("simfix: totalDuration in `simfix()` must be a single positive number!")
  }

  if(!min(totalDuration) > 0){
    stop("simfix: totalDuration in `simfix()` must be a single positive number!")
  }

  # check stratum
  strata2 <- names(table(failRates$Stratum))
  if(nrow(strata) != length(strata2)){
    stop("simfix: Stratum in `simfix()` must be the same in strata and failRates!")
  }

  if(any(is.na(match(strata$Stratum, strata2))) | any(is.na(match(strata2, strata$Stratum)))){
    stop("simfix: Stratum in `simfix()` must be the same in strata and failRates!")
  }

  # check nsim
  if(nsim <= 0){
    stop("simfix: nsim in `simfix()` must be positive integer!")
  }

  if(length(nsim) != 1){
    stop("simfix: nsim in `simfix()` must be positive integer!")
  }

  if(nsim != ceiling(nsim)){
    stop("simfix: nsim in `simfix()` must be positive integer!")
  }

  # check targetEvents
  if(targetEvents <= 0){
    stop("simfix: targetEvents in `simfix()` must be positive!")
  }

  if(length(targetEvents) != 1){
    stop(("simfix: targetEvents in `simfix()` must be positive!"))
  }

  # check sampleSize
  if(sampleSize <= 0){
    stop("simfix: sampleSize in `simfix()` must be positive!")
  }

  if(length(sampleSize) != 1){
    stop("simfix: sampleSize in `simfix()` must be positive")
  }

  # check seed
  if(is.null(seed)) {
    setSeed <- FALSE
  } else {
    if (!is.numeric(seed)){stop("simfix: seed in `simfix()` must be a number")}
    setSeed <- TRUE
  }

  n_stratum <- nrow(strata)

  # build a function to calculate Z and log-hr
  doAnalysis <- function(d, rg, n_stratum){
    if (nrow(rg) == 1){
      Z <- tibble(Z = (d %>% tensurv(txval = "Experimental") %>% tenFH(rg = rg))$Z)
    } else{
      Z <- d %>% tensurv(txval = "Experimental") %>% tenFHcorr(rg = rg, corr = TRUE)
    }

    ans <- tibble(
      Events = sum(d$event),
      lnhr = ifelse(n_stratum > 1,
                    survival::coxph(survival::Surv(tte, event) ~ (Treatment == "Experimental") + survival::strata(Stratum), data = d)$coefficients,
                    survival::coxph(survival::Surv(tte, event) ~ (Treatment == "Experimental"), data = d)$coefficients
                    ) %>% as.numeric()
      )

    ans <- cbind(ans, Z)
    return(ans)
  }

  # compute minimum planned follow-up time
  minFollow <- max(0, totalDuration - sum(enrollRates$duration))

  # put failure rates into simPWSurv format
  temp <- simfix2simPWSurv(failRates)
  fr <- temp$failRates
  dr <- temp$dropoutRates
  results <- NULL

  # parallel
  `%op%` <- get_operator()
  results <- foreach::foreach(
    i = seq_len(nsim),
    .combine = "rbind",
    .errorhandling = "pass"
    ) %op% {

    # set random seed
    if (setSeed){
      set.seed(seed + i - 1)
    }

    # generate piecewise data
    sim <- simPWSurv(n = sampleSize,
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
    r1 <- NULL
    r2 <- NULL
    r3 <- NULL

    tests <- rep(FALSE, 3)

    ## Total planned trial duration or max of this and targeted events
    if (1 %in% timingType){
      tests[1] <- TRUE       # planned duration cutoff
    }

    if (2 %in% timingType){
      tests[2] <- TRUE       # targeted events cutoff
    }

    if (3 %in% timingType){
      tests[3] <- TRUE       # minimum follow-up duration target
    }

    if (4 %in% timingType){  # max of planned duration, targeted events
      if (tedate > totalDuration){
        tests[2] <- TRUE
      }else{
        tests[1] <- TRUE
      }
    }

    if (5 %in% timingType){  # max of minimum follow-up, targeted events
      if (tedate > tmfdate){
        tests[2] <- TRUE
      }else{
        tests[3] <- TRUE
      }
    }

    # Total duration cutoff
    if (tests[1]){
      d <- sim %>% cutData(totalDuration)
      r1 <- d %>% doAnalysis(rg, n_stratum)
    }

    # targeted events cutoff
    if (tests[2]){
      d <- sim %>% cutData(tedate)
      r2 <- d %>% doAnalysis(rg, n_stratum)
    }

    # minimum follow-up cutoff
    if (tests[3]){
      d <- sim %>% cutData(tmfdate)
      r3 <- d %>% doAnalysis(rg, n_stratum)
    }

    addit <- NULL
    # planned duration cutoff
    if (1 %in% timingType){
      addit <- rbind(addit,
                     r1 %>% mutate(cut = "Planned duration",
                                   Duration = totalDuration))
    }

    # targeted events cutoff
    if (2 %in% timingType){
      addit <- rbind(addit,
                     r2 %>% mutate(cut = "Targeted events",
                                   Duration = tedate))
    }

    # minimum follow-up duration target
    if (3 %in% timingType){
      addit <- rbind(addit,
                     r3 %>% mutate(cut = "Minimum follow-up",
                                   Duration = tmfdate))
    }

    # max of planned duration, targeted events
    if (4 %in% timingType){
      if (tedate > totalDuration){
        addit <- rbind(addit,
                       r2 %>% mutate(cut = "Max(planned duration, event cut)",
                                     Duration = tedate))
      }else{
        addit <- rbind(addit,
                       r1 %>% mutate(cut = "Max(planned duration, event cut)",
                                     Duration = totalDuration))
      }
    }

    # max of minimum follow-up, targeted events
    if (5 %in% timingType){
      if (tedate > tmfdate){
        addit <- rbind(addit,
                       r2 %>% mutate(cut = "Max(min follow-up, event cut)",
                                     Duration = tedate))
      }else{
        addit <- rbind(addit,
                       r3 %>% mutate(cut = "Max(min follow-up, event cut)",
                                     Duration = tmfdate))
      }
    }

    addit %>% mutate(Sim = i)
    }

  return(results)
}

# Get operator (parallel or serial)
get_operator <- function() {
  is_par <- foreach::getDoParWorkers() > 1
  if (is_par) {
    message("Using ", foreach::getDoParWorkers(), " cores with backend ", foreach::getDoParName())
    res <- foreach::`%dopar%`
  } else {
    res <- foreach::`%do%`
  }
  res
}
