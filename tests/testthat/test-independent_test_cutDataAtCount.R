TTEdata <- simPWSurv(n=200)
test_that("the input is a time-to-event data set", {

            testthat::expect_equal(1,max(names(TTEdata)=="Stratum"))
            testthat::expect_equal(1,max(names(TTEdata)=="enroll_time"))
            testthat::expect_equal(1,max(names(TTEdata)=="Treatment"))
            testthat::expect_equal(1,max(names(TTEdata)=="failTime"))
            testthat::expect_equal(1,max(names(TTEdata)=="dropoutTime"))
            testthat::expect_equal(1,max(names(TTEdata)=="fail"))
            testthat::expect_equal(1,max(names(TTEdata)=="cte"))

})

