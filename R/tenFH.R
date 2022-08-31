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

#' Fleming-Harrington Weighted Logrank Tests
#'
#' With output from the function \code{tensurv}
#' @param x a \code{tensurv}-class \code{tibble} with a counting process dataset
#' @param rg a \code{tibble} with variables \code{rho} and \code{gamma}, both greater than equal
#' to zero, to specify one Fleming-Harrington weighted logrank test per row;
#' Default: tibble(rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1))
#' @param returnVariance a logical flag that, if true, adds columns
#' estimated variance for weighted sum of observed minus expected; see details; Default: FALSE
#' @return a `tibble` with \code{rg} as input and the FH test statistic
#' for the data in \code{x}
#' (\code{Z}, a directional square root of the usual weighted logrank test);
#' if variance calculations are specified (e.g., to be used for covariances in a combination test),
#' the this will be returned in the column \code{Var}
#' @details
#' The input value \code{x} produced by \code{tensurv()} produces a counting process dataset
#' grouped by strata and sorted within strata by increasing times where events occur.
#' \itemize{
#' \item \eqn{Z} - standardized normal Fleming-Harrington weighted logrank test
#' \item \eqn{i}  - stratum index
#' \item \eqn{d_i} - number of distinct times at which events occurred in stratum \eqn{i}
#' \item \eqn{t_{ij}} - ordered times at which events in stratum \eqn{i}, \eqn{j=1,2,\ldots d_i} were observed;
#' for each observation, \eqn{t_{ij}} represents the time post study entry
#' \item \eqn{O_{ij.}} - total number of events in stratum \eqn{i} that occurred at time \eqn{t_{ij}}
#' \item \eqn{O_{ije}} - total number of events in stratum \eqn{i} in the experimental treatment group that occurred
#' at time \eqn{t_{ij}}
#' \item \eqn{N_{ij.}} - total number of study subjects in stratum \eqn{i} who were followed for at least duration
#' \item \eqn{E_{ije}} - expected observations in experimental treatment group given random selection of \eqn{O_{ij.}}
#' from those in stratum \eqn{i} at risk at time \eqn{t_{ij}}
#' \item \eqn{V_{ije}} - hypergeometric variance for \eqn{E_{ije}} as produced in \code{Var}
#' from the \code{tensurv()} routine
#' \item \eqn{N_{ije}} - total number of study subjects in stratum \eqn{i} in the experimental treatment group
#' who were followed for at least duration \eqn{t_{ij}}
#' \item \eqn{E_{ije}} - expected observations in experimental group in stratum \eqn{i} at time \eqn{t_{ij}}
#' conditioning on the overall number of events and at risk populations at that time and sampling at risk
#' observations without replacement:
#' \deqn{E_{ije} = O_{ij.} N_{ije}/N_{ij.}}
#' \item \eqn{S_{ij}} - Kaplan-Meier estimate of survival in combined treatment groups immediately prior
#' to time \eqn{t_{ij}}
#' \item \eqn{\rho, \gamma} - real parameters for Fleming-Harrington test
#' \item \eqn{X_i} - Numerator for signed logrank test in stratum \eqn{i}
#' \deqn{X_i = \sum_{j=1}^{d_{i}} S_{ij}^\rho(1-S_{ij}^\gamma)(O_{ije}-E_{ije})}
#' \item \eqn{V_{ij}} - variance used in denominator for Fleming-Harrington weighted logrank tests
#' \deqn{V_i = \sum_{j=1}^{d_{i}} (S_{ij}^\rho(1-S_{ij}^\gamma))^2V_{ij})}
#'
#' The stratified Fleming-Harrington weighted logrank test is then computed as:
#' \deqn{Z = \sum_i X_i/\sqrt{\sum_i V_i}}
#' }
#' @examples
#' library(tidyr)
#' # Use default enrollment and event rates at cut at 100 events
#' x <- simPWSurv(n=200) %>% cutDataAtCount(100) %>% tensurv(txval="Experimental")
#' # compute logrank (FH(0,0)) and FH(0,1)
#' tenFH(x,rg=tibble(rho=c(0,0),gamma=c(0,1)))
#' @export

tenFH <- function(x=simPWSurv(n=200) %>% cutDataAtCount(150) %>% tensurv(txval = "Experimental"),
                  rg=tibble(rho=c(0,0,1,1),gamma=c(0,1,0,1)),
                  returnVariance=FALSE){
  # check input failure rate assumptions
  if(!is.data.frame(x)){stop("x in `tenFH()` must be a data frame")}
  if(max(names(x)=="S") != 1){stop("x column names in `tenFH()` must contain S")}
  if(max(names(x)=="OminusE") != 1){stop("x column names in `tenFH()` must contain OminusE ")}
  if(max(names(x)=="Var") != 1){stop("x column names in `tenFH()` must contain Var")}

  # get minimal columns from tensurv item
  xx <- x %>%
        ungroup() %>%
        select(S,OminusE,Var)
  rg$Z <- rep(0,nrow(rg))
  if (returnVariance) rg$Var <- rep(0,nrow(rg))
  for(i in 1:nrow(rg)){
    y <- xx %>% mutate(w=S^rg$rho[i]*(1-S)^rg$gamma[i],
                       wOminusE=w*OminusE,
                       wVar=w^2*Var) %>%
                summarize(wVar=sum(wVar),wOminusE=sum(wOminusE))
    rg$Z[i] <- y$wOminusE/sqrt(y$wVar)
    if (returnVariance) rg$Var[i] <- y$wVar
  }
  rg
}
#' @rdname tenFH
#' @export
