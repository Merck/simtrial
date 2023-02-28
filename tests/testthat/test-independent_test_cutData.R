x <- simPWSurv(n=200)

test_that("x is a time-to-event data set", {
  testthat::expect_equal(1,max(names(x)=="Stratum"))
  testthat::expect_equal(1,max(names(x)=="enrollTime"))
  testthat::expect_equal(1,max(names(x)=="Treatment"))
  testthat::expect_equal(1,max(names(x)=="failTime"))
  testthat::expect_equal(1,max(names(x)=="dropoutTime"))
  testthat::expect_equal(1,max(names(x)=="fail"))
  testthat::expect_equal(1,max(names(x)=="cte"))
})

cutDate=10
xcut<-cut_data_by_date(x,cutDate)
test_that("only paitients recorded by cut_data_by_date are included", {
  Npts=dim(filter(x,enrollTime <= cutDate))[1]
  Nptscut=length(xcut$tte)
  testthat::expect_equal(Npts,Nptscut)
})

test_that("Time to event (tte) is cut off at the cutDate", {
  testthat::expect_lte(max(xcut$tte),cutDate)
})

test_that("the event variable is calculated correctly", {
  Nevent=sum(x$fail*(x$cte<=cutDate))
  Neventcut=dim(filter(xcut,event == 1))[1]
  testthat::expect_equal(Nevent,Neventcut)
})

