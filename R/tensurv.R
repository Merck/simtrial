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

#' Process Survival Data into Counting Process Format
#'
#' Produces a tibble that is sorted by stratum and time.
#' Included in this is only the times at which one or more event occurs.
#' The output dataset contains Stratum, tte (time-to-event), at risk count and count of events at the specified tte
#' sorted by Stratum and tte.
#'
#' The function only considered two group situation.
#'
#' The tie is handled by the Breslow's Method.
#'
#' @param x a tibble with no missing values and contain variables
#'
#'     \code{Stratum} (Stratum)
#'
#'     \code{Treatment} (Treatment group)
#'
#'     \code{tte} (Observed time)
#'
#'     \code{event} (Binary event indicator, 1 represents event, 0 represents censoring)
#'
#' @param txval value in the input \code{Treatment} column that indicates treatment group value.
#'
#'
#' @return A \code{tibble} grouped by
#' \code{Stratum} and sorted within strata by \code{tte}.
#' Remain rows with at least one event in the population, at least one subject is at risk in both treatment group and control group.
#' Other variables in this represent the following within each stratum at each time at which one or more
#' events are observed:
#'
#' \code{events} (Total number of events)
#'
#' \code{txevents} (Total number of events at treatment group)
#'
#' \code{atrisk} (Number of subjects at risk)
#'
#' \code{txatrisk} (Number of subjects at risk in treatment group)
#'
#' \code{S} (Left-continuous Kaplan-Meier survival estimate)
#'
#' \code{OminusE} (In treatment group, observed number of events minus expected number of events.
#'           The expected number of events is estimated by assuming no treatment effect with hypergeometric distribution with
#'           parameters total number of events, total number of events at treatment group and number of events at a time.
#'           (Same assumption of log-rank test under the null hypothesis)
#'
#' \code{Var} (variance of OminusE under the same assumption).
#'
#' @examples
#' library(dplyr)
#'
#' # Example 1
#' x=tibble(Stratum = c(rep(1,10),rep(2,6)),
#' Treatment = rep(c(1,1,0,0),4),
#' tte = 1:16,
#' event= rep(c(0,1),8))
#'
#' tensurv(x, txval=1)
#'
#' # Example 2
#' x <- simPWSurv(n=400)
#' y <- cutDataAtCount(x,150) %>% tensurv(txval = "Experimental")
#' # weighted logrank test (Z-value and 1-sided p-value)
#' z <- sum(y$OminusE)/sqrt(sum(y$Var))
#' c(z,pnorm(z))
#'
#' @export
tensurv <- function(x, txval){

    u.trt = unique(x$Treatment)
    if(length(u.trt) > 2){
      stop("Expected two groups")
    }

    if(! txval %in% u.trt){
      stop("txval is not a valid treatment group value")
    }

    if(! all(unique(x$event) %in% c(0,1) ) ){
      stop("Event indicator must be 0 (censoring) or 1 (event)")
    }
    x %>% group_by(Stratum) %>% arrange(desc(tte)) %>%
          mutate(one=1,
                 atrisk=cumsum(one),
                 txatrisk=cumsum(Treatment==txval)) %>%
          # Handling ties using Breslow's method
          group_by(Stratum, mtte=desc(tte)) %>%
          dplyr::summarise(events=sum(event),
                           txevents=sum((Treatment==txval)*event),
                           tte=first(tte),
                           atrisk=max(atrisk),
                           txatrisk=max(txatrisk)) %>%
          # Keep calculation for observed time with at least one event, at least one subject is
          # at risk in both treatment group and control group.
          filter(events>0, atrisk-txatrisk>0, txatrisk>0) %>%
          select(-mtte) %>%
          mutate(s=1-events/atrisk) %>%
          arrange(Stratum, tte) %>%
          group_by(Stratum) %>%
          mutate(S=lag(cumprod(s), default=1),               # left continuous Kaplan-Meier Estimator
                 OminusE=txevents-txatrisk/atrisk*events,    #  Observed events minus Expected events in treatment group
                 Var=(atrisk-txatrisk)*txatrisk*events*(atrisk-events)/atrisk^2/(atrisk-1)) %>% #Variance of OminusE
                 select(-s)
}

