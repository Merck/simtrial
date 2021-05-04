## ----setup, include = FALSE,message=FALSE,warning=FALSE----
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

options( width = 58 )

## ----message=FALSE,warning=FALSE------------------------
library(simtrial)
library(knitr)
library(dplyr)
x <- simfix(nsim=1,timingType=5,rg=tibble::tibble(rho=c(0,0,1),gamma=c(0,1,1)))
x %>% kable(digits=2)

## ----warning=FALSE,message=FALSE------------------------
pMaxCombo(x)

## ----message=FALSE,warning=FALSE,cache=FALSE------------
s <- simPWSurv(n=100)
head(s) %>% kable(digits=2)

## ----warning=FALSE,message=FALSE------------------------
x <- s %>% cutDataAtCount(75)
head(x) %>% kable(digits=2)

## ----warning=FALSE,message=FALSE------------------------
Z <- s %>% cutDataAtCount(75) %>% 
           tensurv(txval="Experimental") %>% 
           tenFHcorr(rg=tibble(rho=c(0,0,1,1),gamma=c(0,1,0,1)))
Z %>% kable(digits=2)

## ----warning=FALSE,message=FALSE------------------------
pMaxCombo(Z)

## ----warning=FALSE,message=FALSE------------------------
pMaxCombo(Z %>% select(-c(V1,V4)) %>% filter((rho==0 & gamma==1) | (rho==1 & gamma==0)))

## ----warning=FALSE,message=FALSE------------------------
library(survival)
head(aml) %>% kable()

## ----warning=FALSE,message=FALSE------------------------
x <- aml %>% transmute(tte=time,event=status,Stratum="All",Treatment=as.character(x))
head(x) %>% kable()

## ----warning=FALSE,message=FALSE------------------------
x %>% 
  tensurv(txval="Maintained") %>% 
  tenFHcorr(rg=tibble(rho=0,gamma=c(0,1))) %>% 
  pMaxCombo()

## ----cache=FALSE,warning=FALSE,message=FALSE------------
# Only use cut events + min follow-up
xx <- simfix(nsim=100,timingType=5,rg=tibble(rho=c(0,0,1),gamma=c(0,1,1)))
# MaxCombo power estimate for cutoff at max of targeted events, minimum follow-up
p <- unlist(xx %>%  group_by(Sim) %>% group_map(pMaxCombo))
mean(p<.001)

## ----cache=FALSE,warning=FALSE,message=FALSE------------
# Only use cuts for events and events + min follow-up
xx <- simfix(nsim=100,timingType=c(2,5),rg=tibble(rho=0,gamma=c(0,1)))
head(xx) %>% kable(digits=2)

## ----warning=FALSE,message=FALSE------------------------
# subset to targeted events cutoff tests
p <- unlist(xx %>% filter(cut=="Targeted events") %>% group_by(Sim) %>% group_map(pMaxCombo))
mean(p<.025)

## ----warning=FALSE,message=FALSE------------------------
# subset to targeted events cutoff tests
p <- unlist(xx %>% filter(cut!="Targeted events") %>% group_by(Sim) %>% group_map(pMaxCombo))
mean(p<.025)

