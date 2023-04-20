#study design using gsDesign
alpha<-0.025
gamma<- c(5,5,47)
R <- c(1,1,9)
median.c <- 7
hr <- 0.65
dropout <- 0.05
#power=0.9 with target events=227
PE <- gsDesign::nSurv(alpha = alpha, beta = c(1-0.9), sided = 1, lambdaC = log(2)/median.c, hr = hr,
                      eta = -log(1-dropout)/12, gamma = gamma, R = R, T=18)
#power=0.93 with duration=18
PE <- gsDesign::nSurv(alpha = alpha, beta = c(1-0.93), sided = 1, lambdaC = log(2)/median.c, hr = hr,
                      eta = -log(1-dropout)/12, gamma = gamma, R = R, T=18)

#test for power comparing sim_fixed_n results with simple study design
set.seed(1234)
test2<-sim_fixed_n(nsim=100,
              sample_size=434,
              target_event=227,
              stratum = tibble::tibble(stratum = "All", p = 1),
              enroll_rate=tibble::tibble(duration=c(1,1,9),
                                         rate= c(5,5,47)),
              fail_rate=tibble::tibble(stratum="All",
                                       duration=c(100),
                                       fail_rate=log(2)/7,
                                       hr=0.65,
                                       dropoutRate=-log(1-0.05)/12),
              totalDuration=18,
              block=rep(c("Experimental","Control"),2),
              timingType=1:5,
              rg=tibble::tibble(rho=0,gamma=0)
)
#load("./fixtures/test_data_simfix.Rdata")
testthat::test_that("test for sim_fixed_n power comparing to gsDesign results with fixed duration in timingType=1",{
  tt1test<-subset(test2,test2$cut=='Planned duration',select=c(Events, lnhr, Z,Duration,Sim))
  expect_equal(object=sum(as.integer(tt1test$Z<(-1.96)))/100, expected=0.93, tolerance=0.02)
})

testthat::test_that("test for sim_fixed_n power comparing to gsDesign results with target events in timingType=2",{
  tt2test<-subset(test2,test2$cut=='Targeted events',select=c(Events, lnhr, Z,Duration,Sim))
  expect_equal(object=sum(as.integer(tt2test$Z<(-1.96)))/100, expected=0.90, tolerance=0.02)
})

testthat::test_that("test for events in the correct directions in timingType=3 comparing to timingType=2",{
  tt2test<-subset(test2,test2$cut=='Targeted events',select=c(Events, lnhr, Z,Duration,Sim))
  tt3test<-subset(test2,test2$cut=='Minimum follow-up',select=c(Events, lnhr, Z,Duration,Sim))
  ttvalue<-0
  for (i in 1:nrow(tt3test)){
    if ((tt3test$Duration[i]>tt2test$Duration[i])&(tt3test$Events[i]>=tt2test$Events[i])) {ttvalue[i]=1 }
    else if ((tt3test$Duration[i]<=tt2test$Duration[i])&(tt3test$Events[i]<=tt2test$Events[i])) {ttvalue[i]=1 }
    else {ttvalue[i]=0}
  }
  expect_equal(object=unique(ttvalue),expected=1)}
)

testthat::test_that("test for timingType=4 outputs using timingType 1 and 2 output",{
  tt1test<-subset(test2,test2$cut=='Planned duration',select=c(Events, lnhr, Z,Duration,Sim))
  tt2test<-subset(test2,test2$cut=='Targeted events',select=c(Events, lnhr, Z,Duration,Sim))
  tt4test<-subset(test2,test2$cut=='Max(planned duration, event cut)',select=c(Events, lnhr, Z,Duration,Sim))
  tt4event=0
  for (i in 1:nrow(tt4test)){
    if (tt1test$Duration[i]<tt2test$Duration[i]) {tt4event[i]=tt2test$Events[i]}
    else {tt4event[i]=tt1test$Events[i]}
  }
  expect_equal(object=tt4event,expected=tt4test$Events)
})


testthat::test_that("test for timingType=5 outputs using timingType 2 and 3 output",{
  tt2test<-subset(test2,test2$cut=='Targeted events',select=c(Events, lnhr, Z,Duration,Sim))
  tt3test<-subset(test2,test2$cut=='Minimum follow-up',select=c(Events, lnhr, Z,Duration,Sim))
  tt5test<-subset(test2,test2$cut=='Max(min follow-up, event cut)',select=c(Events, lnhr, Z,Duration,Sim))
  tt5event=0
  for (i in 1:nrow(tt5test)){
    if (tt2test$Duration[i]<tt3test$Duration[i]) {tt5event[i]=tt3test$Events[i]}
    else {tt5event[i]=tt2test$Events[i]}
  }
  expect_equal(object=tt5event,expected=tt5test$Events)
})

