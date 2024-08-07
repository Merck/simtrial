test_that("functions that use data.table still return a data frame", {
  class_expected <- "data.frame"

  # counting_process()
  x <- data.frame(
    stratum = c(rep(1, 10), rep(2, 6)),
    treatment = rep(c(1, 1, 0, 0), 4),
    tte = 1:16,
    event = rep(c(0, 1), 8)
  )
  expect_identical(class(counting_process(x, arm = 1)), c("counting_process", class_expected))

  # cut_data_by_date()
  x <- sim_pw_surv(n = 20)
  expect_identical(class(cut_data_by_date(x, 5)), class_expected)

  # early_zero_weight()
  x <- sim_pw_surv(n = 200)
  x <- cut_data_by_event(x, 125)
  x <- counting_process(x, arm = "experimental")
  expect_identical(class(early_zero_weight(x, early_period = 2)), class_expected)

  # fh_weight()
  expect_identical(class(fh_weight()), c("counting_process", class_expected))

  # mb_weight()
  x <- sim_pw_surv()
  x <- cut_data_by_event(x, 125)
  x <- counting_process(x, arm = "experimental")
  expect_identical(class(mb_weight(x)), class_expected)

  # sim_fixed_n()
  expect_identical(class(sim_fixed_n(n_sim = 1)), class_expected)

  # sim_pw_surv()
  expect_identical(class(sim_pw_surv(n = 1)), class_expected)

  # to_sim_pw_surv()
  output <- to_sim_pw_surv()
  expect_identical(class(output$fail_rate), class_expected)
  expect_identical(class(output$dropout_rate), class_expected)
})

# simtrial functions accept any object that inherits "data.frame", eg tibble
# and data.table. These tests ensure that data.table-enabled functions make a
# copy instead of modifying the input object by reference
test_that("functions that use data.table do not modify input data table", {
  skip_if_not_installed("gsDesign2")

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

  # early_zero_weight()
  x <- sim_pw_surv(n = 200)
  x <- cut_data_by_event(x, 125)
  x <- counting_process(x, arm = "experimental")
  data.table::setDT(x)
  x_original <- data.table::copy(x)
  early_zero_weight(x, early_period = 2)
  expect_identical(x, x_original)
  # stratified
  # Example 2: Stratified
  n <- 500
  # Two strata
  stratum <- c("Biomarker-positive", "Biomarker-negative")
  prevalence_ratio <- c(0.6, 0.4)
  # Enrollment rate
  enroll_rate <- gsDesign2::define_enroll_rate(
    stratum = rep(stratum, each = 2),
    duration = c(2, 10, 2, 10),
    rate = c(c(1, 4) * prevalence_ratio[1], c(1, 4) * prevalence_ratio[2])
  )
  enroll_rate$rate <- enroll_rate$rate * n / sum(enroll_rate$duration * enroll_rate$rate)
  # Failure rate
  med_pos <- 10 # Median of the biomarker positive population
  med_neg <- 8 # Median of the biomarker negative population
  hr_pos <- c(1, 0.7) # Hazard ratio of the biomarker positive population
  hr_neg <- c(1, 0.8) # Hazard ratio of the biomarker negative population
  fail_rate <- gsDesign2::define_fail_rate(
    stratum = rep(stratum, each = 2),
    duration = c(3, 1000, 4, 1000),
    fail_rate = c(log(2) / c(med_pos, med_pos, med_neg, med_neg)),
    hr = c(hr_pos, hr_neg),
    dropout_rate = 0.01
  )
  # Simulate data
  temp <- to_sim_pw_surv(fail_rate) # Convert the failure rate
  set.seed(2023)
  x <- sim_pw_surv(
    n = n, # Sample size
    # Stratified design with prevalence ratio of 6:4
    stratum = data.frame(stratum = stratum, p = prevalence_ratio),
    # Randomization ratio
    block = c("control", "control", "experimental", "experimental"),
    enroll_rate = enroll_rate, # Enrollment rate
    fail_rate = temp$fail_rate, # Failure rate
    dropout_rate = temp$dropout_rate # Dropout rate
  )
  x <- cut_data_by_event(x, 125)
  x <- counting_process(x, arm = "experimental")
  data.table::setDT(x)
  x_original <- data.table::copy(x)
  data.table::setDT(fail_rate)
  fail_rate_original <- data.table::copy(fail_rate)
  early_zero_weight(x, early_period = 2, fail_rate = fail_rate)
  expect_identical(x, x_original)
  expect_identical(fail_rate, fail_rate_original)

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

  # rpwexp_enroll()
  enroll_rate <- data.table::data.table(
    rate = c(5, 15, 30),
    duration = c(100, 200, 100)
  )
  data.table::setDT(enroll_rate)
  enroll_rate_original <- data.table::copy(enroll_rate)
  rpwexp_enroll(n = 10, enroll_rate = enroll_rate)
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

  # to_sim_pw_surv()
  fail_rate <- data.table::data.table(
    stratum = "All",
    duration = c(3, 100),
    fail_rate = log(2) / c(9, 18),
    hr = c(0.9, 0.6),
    dropout_rate = rep(0.001, 2)
  )
  data.table::setDT(fail_rate)
  fail_rate_original <- data.table::copy(fail_rate)
  to_sim_pw_surv(fail_rate = fail_rate)
  expect_identical(fail_rate, fail_rate_original)
})
