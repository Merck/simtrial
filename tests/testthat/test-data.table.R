test_that("functions that use data.table still return a data frame", {
  class_expected <- "data.frame"

  # counting_process()
  x <- data.frame(
    stratum = c(rep(1, 10), rep(2, 6)),
    treatment = rep(c(1, 1, 0, 0), 4),
    tte = 1:16,
    event = rep(c(0, 1), 8)
  )
  expect_identical(class(counting_process(x, arm = 1)), class_expected)

  # cut_data_by_date()
  x <- sim_pw_surv(n = 20)
  expect_identical(class(cut_data_by_date(x, 5)), class_expected)

  # fh_weight()
  expect_identical(class(fh_weight()), class_expected)

  # mb_weight()
  x <- sim_pw_surv()
  x <- cut_data_by_event(x, 125)
  x <- counting_process(x, arm = "experimental")
  expect_identical(class(mb_weight(x)), class_expected)

  # sim_fixed_n()
  expect_identical(class(sim_fixed_n(n_sim = 1)), class_expected)

  # sim_pw_surv()
  expect_identical(class(sim_pw_surv(n = 1)), class_expected)

  # simfix2simpwsurv()
  output <- simfix2simpwsurv()
  expect_identical(class(output$fail_rate), class_expected)
  expect_identical(class(output$dropout_rate), class_expected)
})

# simtrial functions accept any object that inherits "data.frame", eg tibble
# and data.table. These tests ensure that data.table-enabled functions make a
# copy instead of modifying the input object by reference
test_that("functions that use data.table do not modify input data table", {
  # confirm that tests can detect a data table that is modified by reference
  modify_by_reference <- function(x) {
    data.table::setDT(x)
    data.table::set(x, i = 1L, j = "second", value = 2)
    return(x[])
  }
  x <- data.table::data.table(first = 1)
  x_original <- data.table::copy(x)
  y <- modify_by_reference(x)
  expect_identical(x, y)
  expect_false(identical(x, x_original))

  # counting_process()
  x <- data.table::data.table(
    stratum = c(rep(1, 10), rep(2, 6)),
    treatment = rep(c(1, 1, 0, 0), 4),
    tte = 1:16,
    event = rep(c(0, 1), 8)
  )
  x_original <- data.table::copy(x)
  counting_process(x, arm = 1)
  expect_identical(x, x_original)

  # cut_data_by_date()
  x <- sim_pw_surv(n = 20)
  data.table::setDT(x)
  x_original <- data.table::copy(x)
  cut_data_by_date(x, 5)
  expect_identical(x, x_original)

  # fh_weight()
  x <- sim_pw_surv()
  x <- cut_data_by_event(x, 125)
  x <- counting_process(x, arm = "experimental")
  data.table::setDT(x)
  x_original <- data.table::copy(x)
  fh_weight(x)
  expect_identical(x, x_original)

  # get_analysis_date()
  x <- sim_pw_surv(n = 5)
  data.table::setDT(x)
  x_original <- data.table::copy(x)
  get_analysis_date(x, planned_calendar_time = 1)
  expect_identical(x, x_original)

  # get_cut_date_by_event()
  x <- sim_pw_surv(n = 5)
  data.table::setDT(x)
  x_original <- data.table::copy(x)
  get_cut_date_by_event(x, event = 1)
  expect_identical(x, x_original)

  # mb_weight()
  x <- sim_pw_surv()
  x <- cut_data_by_event(x, 125)
  x <- counting_process(x, arm = "experimental")
  data.table::setDT(x)
  x_original <- data.table::copy(x)
  mb_weight(x)
  expect_identical(x, x_original)

  # rpw_enroll()
  enroll_rate = data.table::data.table(
   rate = c(5, 15, 30),
   duration = c(100, 200, 100)
  )
  data.table::setDT(enroll_rate)
  enroll_rate_original <- data.table::copy(enroll_rate)
  rpw_enroll(n = 10, enroll_rate = enroll_rate)
  expect_identical(enroll_rate, enroll_rate_original)

  # sim_fixed_n()
  stratum <- data.table::data.table(stratum = "All", p = 1)
  enroll_rate <- data.table::data.table(duration = c(2, 2, 10), rate = c(3, 6, 9))
  fail_rate <- data.table::data.table(
    stratum = "All",
    duration = c(3, 100),
    fail_rate = log(2) / c(9, 18),
    hr = c(0.9, 0.6),
    dropout_rate = rep(0.001, 2)
  )
  data.table::setDT(stratum)
  stratum_original <- data.table::copy(stratum)
  data.table::setDT(enroll_rate)
  enroll_rate_original <- data.table::copy(enroll_rate)
  data.table::setDT(fail_rate)
  fail_rate_original <- data.table::copy(fail_rate)
  sim_fixed_n(
    n_sim = 1,
    stratum = stratum,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate
  )
  expect_identical(stratum, stratum_original)
  expect_identical(enroll_rate, enroll_rate_original)
  expect_identical(fail_rate, fail_rate_original)

  # sim_pw_surv()
  stratum <- data.table::data.table(stratum = "All", p = 1)
  enroll_rate <- data.table::data.table(rate = 9, duration = 1)
  fail_rate <- data.table::data.table(
    stratum = rep("All", 4),
    period = rep(1:2, 2),
    treatment = c(rep("control", 2), rep("experimental", 2)),
    duration = rep(c(3, 1), 2),
    rate = log(2) / c(9, 9, 9, 18)
  )
  dropout_rate <- data.table::data.table(
    stratum = rep("All", 2),
    period = rep(1, 2),
    treatment = c("control", "experimental"),
    duration = rep(100, 2),
    rate = rep(0.001, 2)
  )
  data.table::setDT(stratum)
  stratum_original <- data.table::copy(stratum)
  data.table::setDT(enroll_rate)
  enroll_rate_original <- data.table::copy(enroll_rate)
  data.table::setDT(fail_rate)
  fail_rate_original <- data.table::copy(fail_rate)
  data.table::setDT(dropout_rate)
  dropout_rate_original <- data.table::copy(dropout_rate)
  sim_pw_surv(
    n = 1,
    stratum = stratum,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    dropout_rate = dropout_rate
  )
  expect_identical(stratum, stratum_original)
  expect_identical(enroll_rate, enroll_rate_original)
  expect_identical(fail_rate, fail_rate_original)
  expect_identical(dropout_rate, dropout_rate_original)

  # simfix2simpwsurv()
  fail_rate <- data.table::data.table(
    stratum = "All",
    duration = c(3, 100),
    fail_rate = log(2) / c(9, 18),
    hr = c(0.9, 0.6),
    dropout_rate = rep(0.001, 2)
  )
  data.table::setDT(fail_rate)
  fail_rate_original <- data.table::copy(fail_rate)
  simfix2simpwsurv(fail_rate = fail_rate)
  expect_identical(fail_rate, fail_rate_original)
})
