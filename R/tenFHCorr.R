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

#' @import tibble
#' @import dplyr
NULL

#' @title Fleming-Harrington Weighted Logrank Tests plus Correlations
#'
#' @description
#' Correlations can be used with \code{mvtnorm::pmvnorm} to compute
#' p-value for MaxCombo, the maximum of the specifed
#' Fleming-Harrington tests
#'
#' @param x a \code{tensurv}-class \code{tibble} with a counting process dataset
#' @param rg a \code{tibble} with variables \code{rho} and \code{gamma}, both greater than equal
#' to zero, to specify one Fleming-Harrington weighted logrank test per row
#' @param corr a logical; if TRUE (default), return correlation matrix; otherwise, return covariance matrix
#' @return a `tibble` with \code{rg} as input, the FH test statistics specified
#' for the data in \code{Z}, and the correlation or covariance matrix for these tests in variables starting
#' with \code{V}
#' @examples
#' library(tidyr)
#' library(dplyr)
#' # Use default enrollment and event rates at cut of 100 events
#' x <- simPWSurv(n=200) %>% cutDataAtCount(100) %>% tensurv(txval="Experimental")
#' # compute logrank (FH(0,0)) and FH(0,1)
#' x <- tenFHcorr(rg=tibble(rho=c(0,0),gamma=c(0,1)),x=x)
#' # compute p-value for MaxCombo
#' library(mvtnorm)
#' 1-pmvnorm(lower=rep(min(x$Z),nrow(x)),corr=data.matrix(select(x,-c(rho,gamma,Z))),
#' algorithm=GenzBretz(maxpts=50000,abseps=0.00001))[1]
#' # check that covariance is as expected
#' x <- simPWSurv(n=200) %>%
#'          cutDataAtCount(100) %>%
#'          tensurv(txval="Experimental")
#' x %>% tenFHcorr(rg=tibble(rho=c(0,0),gamma=c(0,1)),corr=FALSE)
#' # Off-diagonal element should be variance in following
#' x %>% tenFHcorr(rg=tibble(rho=0,gamma=.5),corr=FALSE)
#' # compare off diagonal result with tenFH()
#' x %>% tenFH(rg=tibble(rho=0,gamma=.5))
#' @export

tenFHcorr <- function(x=simPWSurv(n=200) %>% cutDataAtCount(100) %>%
                        tensurv(txval = "Experimental"),
                      rg=tibble(rho=c(0,0,1,1),gamma=c(0,1,0,1)),
                      corr=TRUE
                     ){
  # Get average rho and gamma for FH covariance matrix
  # We want rhoave[i,j] = (rho[i]+rho[j])/2
  # and     gamave[i,j] = (gamma[i]+gamma[j])/2
  nr <- nrow(rg)
  rhoave <- (matrix(rg$rho,nrow=nr,ncol=nr)+matrix(rg$rho,nrow=nr,ncol=nr,byrow=TRUE))/2
  gamave <- (matrix(rg$gamma,nrow=nr,ncol=nr)+matrix(rg$gamma,nrow=nr,ncol=nr,byrow=TRUE))/2
  # Convert back to tibble
  rg2 <- tibble(rho=as.numeric(rhoave), gamma=as.numeric(gamave))
  # get unique values of rho, gamma
  rgu <- rg2 %>% unique()
  # compute FH statistic for unique values
  # and merge back to full set of pairs
  # rgFH <- tenFH(x,rgu,returnVariance=TRUE) %>% right_join(rg2,by=c("rho"="rho","gamma"="gamma"))
  # FIXED by KA TO GET SORT CORRECT 8/21/2020
  rgFH <- rg2 %>% left_join(tenFH(x,rgu,returnVariance=TRUE),by=c("rho"="rho","gamma"="gamma"))
  # get Z statistics for input rho, gamma combinations
  Z <- rgFH$Z[(0:(nrow(rg)-1))*nrow(rg)+1:nrow(rg)]
  # get correlation matrix
  c <- matrix(rgFH$Var,nrow=nrow(rg),byrow=TRUE)
  if (corr) c <- stats::cov2cor(c)
  names(c) <- paste("V",1:ncol(c),sep="")
  # return combined values
  mat = cbind(rg, Z, tibble(c))
  return(mat) 
}
#' @rdname tenFHcorr
#' @export
