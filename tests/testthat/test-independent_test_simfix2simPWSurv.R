fail_rate <- tibble::tibble(Stratum=c(rep("Low",3),rep("High",3)),
                            duration=rep(c(4,10,100),2),
                            fail_rate=c(.04,.1,.06,
                                       .08,.16,.12),
                            hr=c(1.5,.5,2/3,
                                 2, 10/16, 10/12),
                            dropout_rate=.01
)

failRatesPWSurv<-simfix2simPWSurv(fail_rate)$fail_rate
dropoutRatesPWSurv<-simfix2simPWSurv(fail_rate)$dropout_rate

testthat::test_that("Stratum values must be the same and stratum length must be doubled after converting",{
  strata1 <- names(table(fail_rate$Stratum))
  strata2 <- names(table(failRatesPWSurv$Stratum))
  testthat::expect_equal(strata1,strata2)
  testthat::expect_equal(length(fail_rate$Stratum)*2,length(failRatesPWSurv$Stratum))
  strata3 <- names(table(dropoutRatesPWSurv$Stratum))
  strata4 <- names(table(dropoutRatesPWSurv$Stratum))
  testthat::expect_equal(strata3,strata4)
  testthat::expect_equal(length(fail_rate$Stratum)*2,length(dropoutRatesPWSurv$Stratum))
})

testthat::test_that("treatment after converting contains only control and experimental with the right length",{
  testthat::expect_equal(length(fail_rate$Stratum),sum(stringr::str_detect(failRatesPWSurv$treatment, "control")))
  testthat::expect_equal(length(fail_rate$Stratum),sum(stringr::str_detect(failRatesPWSurv$treatment, "experimental")))
  testthat::expect_equal(length(failRatesPWSurv$treatment),length(fail_rate$Stratum)*2)
  testthat::expect_equal(length(fail_rate$Stratum),sum(stringr::str_detect(dropoutRatesPWSurv$treatment, "control")))
  testthat::expect_equal(length(fail_rate$Stratum),sum(stringr::str_detect(dropoutRatesPWSurv$treatment, "experimental")))
  testthat::expect_equal(length(dropoutRatesPWSurv$treatment),length(fail_rate$Stratum)*2)
})

testthat::test_that("Duration values match before and after converting and in right length ",{
  testthat::expect_equal(rep(c(fail_rate$duration),2),failRatesPWSurv$duration)
  testthat::expect_equal(rep(c(fail_rate$duration),2),dropoutRatesPWSurv$duration)
})

testthat::test_that("fail_rate match before and after converting and are in right length ",{
  testthat::expect_equal(fail_rate$fail_rate,failRatesPWSurv$rate[1:length(fail_rate$fail_rate)])
  testthat::expect_equal(fail_rate$fail_rate*fail_rate$hr,failRatesPWSurv$rate[(length(fail_rate$fail_rate)+1):(length(fail_rate$fail_rate)*2)])
})

testthat::test_that("dropout_rate match before and after converting and are in right length ",{
  testthat::expect_equal(fail_rate$dropout_rate,dropoutRatesPWSurv$rate[1:length(fail_rate$dropout_rate)])
  testthat::expect_equal(fail_rate$dropout_rate,dropoutRatesPWSurv$rate[(length(fail_rate$fail_rate)+1):(length(fail_rate$fail_rate)*2)])
})


# "meaningful error messages when the inputs are incorrect"
testthat::test_that("fail_rate column names must contain Stratum, duration, fail_rate, hr and dropout_rate",{
  testthat::expect_equal(1,max(names(fail_rate)=="Stratum"))
  testthat::expect_equal(1,max(names(fail_rate)=="duration"))
  testthat::expect_equal(1,max(names(fail_rate)=="fail_rate"))
  testthat::expect_equal(1,max(names(fail_rate)=="hr"))
  testthat::expect_equal(1,max(names(fail_rate)=="dropout_rate"))
})

testthat::test_that("duration must be longer than 0",{
  testthat::expect_equal(TRUE,is.numeric(fail_rate$duration))
  testthat::expect_gt(min(fail_rate$duration),0)
})

testthat::test_that("fail_rate must be smaller than 1 and positive",{
  testthat::expect_lt(max(fail_rate$fail_rate),1)
  testthat::expect_gt(min(fail_rate$fail_rate),0)
})


testthat::test_that("hr must be postiive",{
  testthat::expect_gt(min(fail_rate$hr),0)
})


testthat::test_that("dropout_rate must be smaller than 1 and positive",{
  testthat::expect_lt(max(fail_rate$dropout_rate),1)
  testthat::expect_gt(min(fail_rate$dropout_rate),0)})
