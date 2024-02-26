test_that("x is a time-to-event data set", {
  x <- sim_pw_surv(n = 200)

  expect_equal(1, max(names(x) == "stratum"))
  expect_equal(1, max(names(x) == "enroll_time"))
  expect_equal(1, max(names(x) == "treatment"))
  expect_equal(1, max(names(x) == "fail_time"))
  expect_equal(1, max(names(x) == "dropout_time"))
  expect_equal(1, max(names(x) == "fail"))
  expect_equal(1, max(names(x) == "cte"))
})

test_that("only patients recorded by cut_data_by_date are included", {
  x <- sim_pw_surv(n = 200)
  cut_date <- 10
  xcut <- cut_data_by_date(x, cut_date)

  Npts <- dim(dplyr::filter(x, enroll_time <= cut_date))[1]
  Nptscut <- length(xcut$tte)

  expect_equal(Npts, Nptscut)
})

test_that("Time-to-event (TTE) is cut off at the cut_date", {
  x <- sim_pw_surv(n = 200)
  cut_date <- 10
  xcut <- cut_data_by_date(x, cut_date)

  expect_lte(max(xcut$tte), cut_date)
})

test_that("the event variable is calculated correctly", {
  x <- sim_pw_surv(n = 200)
  cut_date <- 10
  xcut <- cut_data_by_date(x, cut_date)

  Nevent <- sum(x$fail * (x$cte <= cut_date))
  Neventcut <- dim(dplyr::filter(xcut, event == 1))[1]

  expect_equal(Nevent, Neventcut)
})
