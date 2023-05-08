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
test2<-sim_fixed_n(n_sim=100,
              sample_size=434,
              target_event=227,
              stratum = tibble::tibble(stratum = "All", p = 1),
              enroll_rate=tibble::tibble(duration=c(1,1,9),
                                         rate= c(5,5,47)),
              fail_rate=tibble::tibble(stratum="All",
                                       duration=c(100),
                                       fail_rate=log(2)/7,
                                       hr=0.65,
                                       dropout_rate=-log(1-0.05)/12),
              totalDuration=18,
              block=rep(c("experimental","control"),2),
              timing_type=1:5,
              rho_gamma=tibble::tibble(rho=0,gamma=0)
)

testthat::test_that("test for sim_fixed_n power comparing to gsDesign results with fixed duration in timing_type=1",{
  tt1test<-subset(test2,test2$cut=='Planned duration',select=c(Events, lnhr, z,Duration,Sim))
  expect_equal(object=sum(as.integer(tt1test$z<(-1.96)))/100, expected=0.93, tolerance=0.02)
})

testthat::test_that("test for sim_fixed_n power comparing to gsDesign results with target events in timing_type=2",{
  tt2test<-subset(test2,test2$cut=='Targeted events',select=c(Events, lnhr, z,Duration,Sim))
  expect_equal(object=sum(as.integer(tt2test$z<(-1.96)))/100, expected=0.90, tolerance=0.02)
})

testthat::test_that("test for events in the correct directions in timing_type=3 comparing to timing_type=2",{
  tt2test<-subset(test2,test2$cut=='Targeted events',select=c(Events, lnhr, z,Duration,Sim))
  tt3test<-subset(test2,test2$cut=='Minimum follow-up',select=c(Events, lnhr, z,Duration,Sim))
  ttvalue<-0
  for (i in 1:nrow(tt3test)){
    if ((tt3test$Duration[i]>tt2test$Duration[i])&(tt3test$Events[i]>=tt2test$Events[i])) {ttvalue[i]=1 }
    else if ((tt3test$Duration[i]<=tt2test$Duration[i])&(tt3test$Events[i]<=tt2test$Events[i])) {ttvalue[i]=1 }
    else {ttvalue[i]=0}
  }
  expect_equal(object=unique(ttvalue),expected=1)}
)

testthat::test_that("test for timing_type=4 outputs using timing_type 1 and 2 output",{
  tt1test<-subset(test2,test2$cut=='Planned duration',select=c(Events, lnhr, z,Duration,Sim))
  tt2test<-subset(test2,test2$cut=='Targeted events',select=c(Events, lnhr, z,Duration,Sim))
  tt4test<-subset(test2,test2$cut=='Max(planned duration, event cut)',select=c(Events, lnhr, z,Duration,Sim))
  tt4event=0
  for (i in 1:nrow(tt4test)){
    if (tt1test$Duration[i]<tt2test$Duration[i]) {tt4event[i]=tt2test$Events[i]}
    else {tt4event[i]=tt1test$Events[i]}
  }
  expect_equal(object=tt4event,expected=tt4test$Events)
})


testthat::test_that("test for timing_type=5 outputs using timing_type 2 and 3 output",{
  tt2test<-subset(test2,test2$cut=='Targeted events',select=c(Events, lnhr, z,Duration,Sim))
  tt3test<-subset(test2,test2$cut=='Minimum follow-up',select=c(Events, lnhr, z,Duration,Sim))
  tt5test<-subset(test2,test2$cut=='Max(min follow-up, event cut)',select=c(Events, lnhr, z,Duration,Sim))
  tt5event=0
  for (i in 1:nrow(tt5test)){
    if (tt2test$Duration[i]<tt3test$Duration[i]) {tt5event[i]=tt3test$Events[i]}
    else {tt5event[i]=tt2test$Events[i]}
  }
  expect_equal(object=tt5event,expected=tt5test$Events)
})

