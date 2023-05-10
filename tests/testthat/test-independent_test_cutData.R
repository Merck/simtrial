x <- sim_pw_surv(n=200)

test_that("x is a time-to-event data set", {
  testthat::expect_equal(1,max(names(x)=="stratum"))
  testthat::expect_equal(1,max(names(x)=="enroll_time"))
  testthat::expect_equal(1,max(names(x)=="treatment"))
  testthat::expect_equal(1,max(names(x)=="fail_time"))
  testthat::expect_equal(1,max(names(x)=="dropout_time"))
  testthat::expect_equal(1,max(names(x)=="fail"))
  testthat::expect_equal(1,max(names(x)=="cte"))
})


cut_date=10
xcut<-cut_data_by_date(x,cut_date)
test_that("only paitients recorded by cut_data_by_date are included", {
  Npts=dim(filter(x,enroll_time <= cut_date))[1]
  Nptscut=length(xcut$tte)
  testthat::expect_equal(Npts,Nptscut)
})

test_that("Time to event (tte) is cut off at the cut_date", {
  testthat::expect_lte(max(xcut$tte),cut_date)
})

test_that("the event variable is calculated correctly", {
  Nevent=sum(x$fail*(x$cte<=cut_date))
  Neventcut=dim(filter(xcut,event == 1))[1]
  testthat::expect_equal(Nevent,Neventcut)
})

