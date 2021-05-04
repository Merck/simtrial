## ----setup, include = FALSE,message=FALSE,warning=FALSE----
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

options( width = 58 )


## ---- message=FALSE, warning=FALSE----------------------
library(simtrial)
library(dplyr)
library(survival)

## ---- message=FALSE, warning=FALSE----------------------
studyDuration = 36
sampleSize <- 200
enrollRates <- tibble::tibble(duration = 12, rate = 200/12)
failRates <- tibble::tribble(
  ~Stratum,  ~duration, ~failRate, ~hr, ~dropoutRate,
  "All",     6,         log(2)/15,  1,  0,
  "All",     36,        log(2)/15,  .7, 0
)

## ---- message=FALSE,warning=FALSE,fig.width=7.5, fig.height=4----
set.seed(7783)
xpar <- simfix2simPWSurv(failRates)
MBdelay <- simPWSurv(n = sampleSize, 
                     strata = tibble::tibble(Stratum = "All", p = 1),
                     block = c(rep("Control", 2), rep("Experimental", 2)),
                     enrollRates = enrollRates,
                     failRates = xpar$failRates, 
                     dropoutRates = xpar$dropoutRates) %>% 
           cutData(studyDuration)
fit <- survfit(Surv(tte,event)~Treatment,data=MBdelay)
plot(fit,col=1:2,mark="|", xaxt="n")
axis(1, xaxp=c(0, 36, 6))

## -------------------------------------------------------
xx <- MBdelay %>% 
    tensurv(txval="Experimental") %>%
    tenFHcorr(rg=tibble(rho=c(0,0,1),gamma=c(0,1,1))) %>%
    mutate(p=pnorm(Z))
xx

## -------------------------------------------------------
xx %>% pMaxCombo()

## -------------------------------------------------------
ZMB <-  MBdelay %>% 
        tensurv(txval="Experimental") %>% 
        wMB(6) %>% 
        summarize(S=sum(OminusE*wMB),V=sum(Var*wMB^2),Z=S/sqrt(V))
# Compute p-value of modestly weighted logrank of Magirr-Burman
pnorm(ZMB$Z)

## -------------------------------------------------------
xx <- MBdelay %>% 
    tensurv(txval="Experimental") %>%
    tenFHcorr(rg=tibble(rho=c(0,0,.5),gamma=c(0,.5,.5))) %>%
    mutate(p=pnorm(Z))
xx

## -------------------------------------------------------
xx %>% pMaxCombo()

## ---- message=FALSE, warning=FALSE----------------------
studyDuration = 5
sampleSize <- 2000
enrollDuration <- .0001
enrollRates <- tibble::tibble(duration = enrollDuration, rate = sampleSize/enrollDuration)
failRates <- tibble::tibble(Stratum="All",
                 failRate=0.25,
                 dropoutRate=0,
                 hr=c(4/.25,.19/.25),
                 duration=c(.1,4.9)
)

## ---- message=FALSE,warning=FALSE,fig.width=7.5, fig.height=4----
set.seed(7783)
xpar <- simfix2simPWSurv(failRates)
FHwn <- simPWSurv(n = sampleSize, 
                     strata = tibble::tibble(Stratum = "All", p = 1),
                     block = c(rep("Control", 2), rep("Experimental", 2)),
                     enrollRates = enrollRates,
                     failRates = xpar$failRates, 
                     dropoutRates = xpar$dropoutRates) %>% 
           cutData(studyDuration)
fit <- survfit(Surv(tte,event)~Treatment,data=FHwn)
plot(fit,col=1:2,mark="|", xaxt="n")
axis(1, xaxp=c(0, 36, 6))

## -------------------------------------------------------
xx <- FHwn %>% 
    tensurv(txval="Experimental") %>%
    tenFHcorr(rg=tibble(rho=c(0,0,1),gamma=c(0,1,1))) %>%
    mutate(p=pnorm(Z))
xx

## -------------------------------------------------------
xx %>% pMaxCombo()

## -------------------------------------------------------
ZMB <-  FHwn %>% 
        tensurv(txval="Experimental") %>% 
        wMB(6) %>% 
        summarize(S=sum(OminusE*wMB),V=sum(Var*wMB^2),Z=S/sqrt(V))
# Compute p-value of modestly weighted logrank of Magirr-Burman
pnorm(ZMB$Z)

## -------------------------------------------------------
xx <- FHwn %>% 
    tensurv(txval="Experimental") %>%
    tenFHcorr(rg=tibble(rho=c(0,0,.5),gamma=c(0,.5,.5))) %>%
    mutate(p=pnorm(Z))
xx

## -------------------------------------------------------
xx %>% pMaxCombo()

