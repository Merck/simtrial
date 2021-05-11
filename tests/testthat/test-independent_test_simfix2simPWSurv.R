failRates <- tibble::tibble(Stratum=c(rep("Low",3),rep("High",3)),
                            duration=rep(c(4,10,100),2),
                            failRate=c(.04,.1,.06,
                                       .08,.16,.12),
                            hr=c(1.5,.5,2/3,
                                 2, 10/16, 10/12),
                            dropoutRate=.01
)

failRatesPWSurv<-simfix2simPWSurv(failRates)$failRates
dropoutRatesPWSurv<-simfix2simPWSurv(failRates)$dropoutRates

testthat::test_that("Stratum values must be the same and stratum length must be doubled after converting",{
  strata1 <- names(table(failRates$Stratum))
  strata2 <- names(table(failRatesPWSurv$Stratum))
  testthat::expect_equal(strata1,strata2)
  testthat::expect_equal(length(failRates$Stratum)*2,length(failRatesPWSurv$Stratum))
  strata3 <- names(table(dropoutRatesPWSurv$Stratum))
  strata4 <- names(table(dropoutRatesPWSurv$Stratum))
  testthat::expect_equal(strata3,strata4)
  testthat::expect_equal(length(failRates$Stratum)*2,length(dropoutRatesPWSurv$Stratum))
}) 

testthat::test_that("Treatment after converting contains only Control and Experimental with the right length",{
  testthat::expect_equal(length(failRates$Stratum),sum(stringr::str_detect(failRatesPWSurv$Treatment, "Control")))
  testthat::expect_equal(length(failRates$Stratum),sum(stringr::str_detect(failRatesPWSurv$Treatment, "Experimental")))
  testthat::expect_equal(length(failRatesPWSurv$Treatment),length(failRates$Stratum)*2)
  testthat::expect_equal(length(failRates$Stratum),sum(stringr::str_detect(dropoutRatesPWSurv$Treatment, "Control")))
  testthat::expect_equal(length(failRates$Stratum),sum(stringr::str_detect(dropoutRatesPWSurv$Treatment, "Experimental")))
  testthat::expect_equal(length(dropoutRatesPWSurv$Treatment),length(failRates$Stratum)*2)
})

testthat::test_that("Duration values match before and after converting and in right length ",{
  testthat::expect_equal(rep(c(failRates$duration),2),failRatesPWSurv$duration)
  testthat::expect_equal(rep(c(failRates$duration),2),dropoutRatesPWSurv$duration)
}) 

testthat::test_that("failRates match before and after converting and are in right length ",{
  testthat::expect_equal(failRates$failRate,failRatesPWSurv$rate[1:length(failRates$failRate)])
  testthat::expect_equal(failRates$failRate*failRates$hr,failRatesPWSurv$rate[(length(failRates$failRate)+1):(length(failRates$failRate)*2)])
})

testthat::test_that("dropoutRates match before and after converting and are in right length ",{
  testthat::expect_equal(failRates$dropoutRate,dropoutRatesPWSurv$rate[1:length(failRates$dropoutRate)])
  testthat::expect_equal(failRates$dropoutRate,dropoutRatesPWSurv$rate[(length(failRates$failRate)+1):(length(failRates$failRate)*2)])
})


# "meaningful error messages when the inputs are incorrect"
testthat::test_that("failRates column names must contain Stratum, duration, failRate, hr and dropoutRate",{
  testthat::expect_equal(1,max(names(failRates)=="Stratum"))
  testthat::expect_equal(1,max(names(failRates)=="duration"))
  testthat::expect_equal(1,max(names(failRates)=="failRate"))
  testthat::expect_equal(1,max(names(failRates)=="hr"))
  testthat::expect_equal(1,max(names(failRates)=="dropoutRate"))
})

testthat::test_that("duration must be longer than 0",{
  testthat::expect_equal(TRUE,is.numeric(failRates$duration))
  testthat::expect_gt(min(failRates$duration),0)
})

testthat::test_that("failRate must be smaller than 1 and positive",{
  testthat::expect_lt(max(failRates$failRate),1)
  testthat::expect_gt(min(failRates$failRate),0)
})


testthat::test_that("hr must be postiive",{
  testthat::expect_gt(min(failRates$hr),0)
})


testthat::test_that("dropoutRate must be smaller than 1 and positive",{
  testthat::expect_lt(max(failRates$dropoutRate),1)
  testthat::expect_gt(min(failRates$dropoutRate),0)})
