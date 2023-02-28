
strata <- tibble::tibble(Stratum=c("Low","High"),p=c(.4,.6))

block <- c(rep("Control",2),rep("Experimental",2))

enrollRates = tibble::tibble(duration = c(5,195), rate = c(100,3000))

failRates <- bind_rows(
  tibble::tibble(Stratum="Low" ,period=1,Treatment="Control"     ,duration=3,rate=.03),
  tibble::tibble(Stratum="Low" ,period=2,Treatment="Control"     ,duration=297,rate=.03),
  tibble::tibble(Stratum="Low" ,period=1,Treatment="Experimental",duration=3,rate=.03),
  tibble::tibble(Stratum="Low" ,period=2,Treatment="Experimental",duration=297,rate=.02),
  tibble::tibble(Stratum="High",period=1,Treatment="Control"     ,duration=3,rate=.05),
  tibble::tibble(Stratum="High",period=2,Treatment="Control"     ,duration=297,rate=.05),
  tibble::tibble(Stratum="High",period=1,Treatment="Experimental",duration=3,rate=.06),
  tibble::tibble(Stratum="High",period=2,Treatment="Experimental",duration=297,rate=.03)
)
dropoutRates <- bind_rows(
  tibble::tibble(Stratum="Low" ,period=1,Treatment="Control"     ,duration=300,rate=.001),
  tibble::tibble(Stratum="Low" ,period=1,Treatment="Experimental",duration=300,rate=.001),
  tibble::tibble(Stratum="High",period=1,Treatment="Control"     ,duration=300,rate=.001),
  tibble::tibble(Stratum="High",period=1,Treatment="Experimental",duration=300,rate=.001)
)
set.seed(1)
x <- simPWSurv(n=400000,
               strata = strata,
               block = block,
               enrollRates = enrollRates,
               failRates=failRates,
               dropoutRates=dropoutRates)

#prepare to test block
block1<-x%>%filter(Stratum=='Low')
block2<-x%>%filter(Stratum=='High')
bktest1 <- c()
j=1
for (i in seq(1,floor(nrow(block1)/4))){
  j=4*i-3
  bktest1[i]<-sum(stringr::str_count(block1$Treatment[j:(j+3)], "Control"))
}
j=1
bktest2<-0
for (i in seq(1,floor(nrow(block2)/4))){
  j=4*i-3
  bktest2[i]<-sum(stringr::str_count(block2$Treatment[j:(j+3)], "Control"))
}

#prepare to test failRates
y <- cut_data_by_date(x,cutDate=300)
intervals<-c(3)
rate00 <- with(subset(y,Treatment=='Control'|Stratum=='Low'), pwexpfit(Surv(tte,event),intervals))
rate01 <- with(subset(y,Treatment=='Control'|Stratum=='High'), pwexpfit(Surv(tte,event),intervals))
rate10 <- with(subset(y,Treatment=='Experimental'|Stratum=='Low'), pwexpfit(Surv(tte,event),intervals))
rate11 <- with(subset(y,Treatment=='Experimental'|Stratum=='High'), pwexpfit(Surv(tte,event),intervals))
ratetest<- c(rate00$rate,rate10$rate,rate01$rate, rate11$rate)
xevent<-bind_rows(rate00, rate01,rate10,rate11)

testthat::test_that("Strata percentage calculated from simulated dataset must be within
                    the tolerance=0.002 of strata in setup (0.4,0.6)",{
  expect_equal(object=c(sum(stringr::str_count(x$Stratum, "Low"))/400000,
                        sum(stringr::str_count(x$Stratum, "High"))/400000),
               expected=c(0.4, 0.6), tolerance=0.002)
})

testthat::test_that("block calculated from simulated dataset equals size of 4 with 1:1
                    randomization, which is 2 for each arm",{
  expect_equal(object=bktest1, expected=rep(2,length(bktest1)))
  expect_equal(object=bktest2, expected=rep(2,length(bktest2)))
})

testthat::test_that("failRates calculated from simulated dataset must be within the
                    tolerance=0.1 of failRates in setting",{
  expect_equal(object=ratetest, expected=failRates$rate, tolerance=0.1)})

testthat::test_that("DropoutRates calculated from simulated dataset must be within
                    the tolerance=0.0005 of DropoutRates=0.001 in setup",{
  duration=300
  drtest<-0
  for (i in 1:duration){
    drtest[i]=sum(x$dropoutTime<=i&x$dropoutTime>(i-1))/400000
  }
  expect_equal(object=drtest, expected=rep(0.001,300), tolerance=0.001)})

testthat::test_that("enrollRates calculated from simulated dataset must be within
                    the relative tolerance=0.05 of enrollRates in setup",{
  duration=300
  entest<-0
  for (i in 1:duration){
    entest[i]=sum(x$enrollTime<=i&x$enrollTime>(i-1))
  }
  entest1<-entest[entest!=0]
  entestexp<-c(rep(100,5), rep(3000,length(entest1)-5))
  entest2<-(entest1-entestexp)/entestexp
  expect_equal(object=entest2, expected=rep(0,length(entest1)), tolerance=0.05)})

#check the arguments, by changing n, the actual number of events changes
set.seed(2468)
z <- simPWSurv(n=300000,
               strata = strata,
               block = block,
               enrollRates = enrollRates,
               failRates=failRates,
               dropoutRates=dropoutRates)

y1 <- cut_data_by_date(z,cutDate=300)
intervals<-c(3)
rate00 <- with(subset(y1,Treatment=='Control'|Stratum=='Low'), pwexpfit(Surv(tte,event),intervals))
rate01 <- with(subset(y1,Treatment=='Control'|Stratum=='High'), pwexpfit(Surv(tte,event),intervals))
rate10 <- with(subset(y1,Treatment=='Experimental'|Stratum=='Low'), pwexpfit(Surv(tte,event),intervals))
rate11 <- with(subset(y1,Treatment=='Experimental'|Stratum=='High'), pwexpfit(Surv(tte,event),intervals))
zevent<-bind_rows(rate00, rate01,rate10,rate11)

testthat::test_that("The actual number of events changes by changing total sample size",{
  expect_false(unique(xevent$event==zevent$event))})

