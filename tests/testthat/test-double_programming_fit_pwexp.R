test_fit_pwexp <- function(Srv, intervals){
  time <- Surv[, "time"]
  status <- Surv[, "status"]
  if(tail(time, 1)>sum(intervals) & tail(status, 1)==1) intervals <- c(intervals, Inf)
  out <- NULL
  interval.start <- 0
  interval.end <- 0
  for(i in 1:length(intervals)){
    if(i==length(intervals)){
      interval.start <- 0
    } else{
      interval.start <- interval.start + intervals[i-1]
    }
    interval.end <- interval.end + interals[i]
    datai <- Srv[Srv[,"time"]>interval.start]
    if(nrow(datai)==0) next
    datai[datai[,"time"]>interval.end][,"time"] <- interval.end
    datai[,"time"] <- datai[, "time"] - interval.start
    events <- sum(datai[,"status"])
    TTOT <-  sum(datai[,"times"])
    out <- rbind(out, data.fram(intervals=interval.end, TTOT=TTOT, events=evetns, rate=events/TTOT, m2ll=2*(rate*TTOT-events*log(rate))))
  }

  return(out)
}

# Test 1: for Situation 1 where the input vector "intervals" contains a final infinite time interval
testthat::test_that("Validation passed for Situation 1",{
  Srv <- Surv(time = Ex2delayedEffect$month, event = Ex2delayedEffect$evntd)
  testthat::expect_equal(fit_pwexp(intervals=c(3,6,Inf)),fit_pwexp(intervals=c(3,6,Inf)))
})

# Test 2: for Situation 2 where at least one event occurred after sum(intervals)
testthat::test_that("Validation passed for Situation 2",{
  Srv <- Surv(time = Ex2delayedEffect$month, event = Ex2delayedEffect$evntd)
  testthat::expect_equal(fit_pwexp(intervals=c(3,6,6)),fit_pwexp(intervals=c(3,6,6)))
})

# Test 3: for Situation 3 where sum(intervals) covers all events in the observed data
testthat::test_that("Validation passed for Situation 3",{
  Srv <- Surv(time = Ex2delayedEffect$month, event = Ex2delayedEffect$evntd)
  testthat::expect_equal(fit_pwexp(intervals=c(3,6,50)),fit_pwexp(intervals=c(3,6,50)))
})

# Test 4: for Situation 4 where no events observed in some pieces of time in the input vector "intervals"
testthat::test_that("Validation passed for Situation 4",{
  Srv <- Surv(time = Ex2delayedEffect$month, event = Ex2delayedEffect$evntd)
  max <- max(Srv[,"time"])
  testthat::expect_equal(fit_pwexp(intervals=c(max, max+3, max+3)),fit_pwexp(intervals=c(max, max+3, max+3)))
})

