## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- warning=FALSE, message=FALSE---------------------------------------
library(simtrial)
library(knitr)
library(tibble)
library(dplyr)

## ------------------------------------------------------------------------
fixedBlockRand(n=10,block=c("A","Dog","Cat","Cat"))

## ------------------------------------------------------------------------
fixedBlockRand(n=20)

## ------------------------------------------------------------------------
rpwenroll(n = 20, enrollRates = tibble(duration = c(1, 2), 
                                       rate = c(2,5)))

## ----fig.width=6---------------------------------------------------------
x <- rpwexp(10000, 
            failRates=tibble(rate = c(1, 3, 10), 
                             duration = c(.5,.5,1)))
plot(sort(x),(10000:1)/10001,log="y", 
     main="PW Exponential simulated survival curve",
     xlab="Time", ylab="P{Survival}")

## ------------------------------------------------------------------------
strata <- tibble(Stratum=c("Negative","Positive"), p=c(.5,.5))

block <- c(rep("Control",2),rep("Experimental",2))

enrollRates <- tibble(rate=c(3, 6, 9), duration=c(3,2,1))

failRates <- tibble(Stratum=c(rep("Negative",4),rep("Positive",4)),
                    period=rep(1:2,4),
                    Treatment=rep(c(rep("Control",2), rep("Experimental",2)),2),
                    duration=rep(c(3,1),4),
                    rate=log(2)/c(4,9,4.5,10,4,9,8,18))
dropoutRates <- tibble(Stratum=c(rep("Negative",4),rep("Positive",4)),
                       period=rep(1:2,4),
                       Treatment=rep(c(rep("Control",2), rep("Experimental",2)),2),
                       duration=rep(c(3,1),4),
                       rate=rep(c(.001,.001),4))

## ------------------------------------------------------------------------
x <- simPWSurv(n=400,
              strata = strata,
              block = block,
              enrollRates = enrollRates,
              failRates=failRates,
              dropoutRates=dropoutRates)
head(x) %>% kable(digits=2)

## ------------------------------------------------------------------------
y <- cutData(x,cutDate=5)
head(y) %>% kable(digits=2)

## ------------------------------------------------------------------------
cut50Positive <- getCutDateForCount(filter(x,Stratum=="Positive"),50)
y50Positive <- cutData(x,cut50Positive)
with(y50Positive,table(Stratum,event))

## ------------------------------------------------------------------------
y150 <- cutDataAtCount(x,150)
table(y150$event,y150$Treatment)

## ------------------------------------------------------------------------
ten150 <- tensurv(y150,txval="Experimental")
head(ten150) %>% kable(digits=2)

## ------------------------------------------------------------------------
z <- with(ten150,sum(OminusE)/sqrt(sum(Var)))
c(z,pnorm(z))

## ------------------------------------------------------------------------
xx <- mutate(ten150,w=S*(1-S)^2)
z <- with(xx,sum(OminusE*w)/sum(sqrt(Var*w^2)))
c(z,pnorm(z))

## ------------------------------------------------------------------------
tenFH(x=ten150,rg=tibble(rho=c(0,0,1,1),gamma=c(0,1,0,1))) %>% kable(digits=2)

## ----message=FALSE-------------------------------------------------------
x <- ten150 %>% tenFHcorr(rg=tibble(rho=c(0,0,1,1),gamma=c(0,1,0,1)))
x %>% kable(digits=2)

## ------------------------------------------------------------------------
# compute p-value for MaxCombo
pMaxCombo(x)

## ------------------------------------------------------------------------
strata <- tibble(Stratum = "All", p = 1)
enrollRates <- tibble(duration = c(2, 2, 10), 
                      rate = c(3, 6, 9)
                     ) 
failRates <- tibble(Stratum = "All", 
                    duration = c(3, 100),
                    failRate = log(2)/c(9, 18), 
                    hr = c(0.9, 0.6), 
                    dropoutRate = rep(0.001, 2)
                   )
block <- rep(c("Experimental", "Control"), 2)
rg <- tibble(rho = 0, gamma = 0)

## ------------------------------------------------------------------------
simfix(nsim = 2,                  # Number of simulations
       sampleSize = 500,          # Trial sample size
       targetEvents = 350,        # Targeted events at analysis
       strata = strata,           # Study strata
       enrollRates = enrollRates, # Enrollment rates
       failRates = failRates,     # Failure rates
       totalDuration = 30,        # Planned trial duration 
       block = block,             # Block for treatment
       timingType = 1:5,          # Use all possible data cutoff methods
       rg = rg                    # FH test(s) to use; in this case, logrank
) %>% kable(digits=2)

