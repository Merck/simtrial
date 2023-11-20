library(dplyr)

test_that("early_zero_weight() with unstratified data", {
  # Example 1: Unstratified
  set.seed(123)
  input <- sim_pw_surv(n = 200)
  input <- cut_data_by_event(input, 125)
  input <- counting_process(input, arm = "experimental")
  output <- early_zero_weight(input, early_period = 2)

  observed <- output$weight
  expected <- rep(c(0, 1), c(15L, 110L))
  expect_equal(observed, expected)
})

test_that("early_zero_weight() with stratified data", {
  # Example 2: Stratified
  n <- 500
  # Two strata
  stratum <- c("Biomarker-positive", "Biomarker-negative")
  prevelance_ratio <- c(0.6, 0.4)
  # Enrollment rate
  enroll_rate <- gsDesign2::define_enroll_rate(
    stratum = rep(stratum, each = 2),
    duration = c(2, 10, 2, 10),
    rate = c(c(1, 4) * prevelance_ratio[1], c(1, 4) * prevelance_ratio[2])
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
  temp <- simfix2simpwsurv(fail_rate) # Convert the failure rate
  set.seed(2023)
  input <- sim_pw_surv(
    n = n, # Sample size
    # Stratified design with prevalence ratio of 6:4
    stratum = data.frame(stratum = stratum, p = prevelance_ratio),
    # Randomization ratio
    block = c("control", "control", "experimental", "experimental"),
    enroll_rate = enroll_rate, # Enrollment rate
    fail_rate = temp$fail_rate, # Failure rate
    dropout_rate = temp$dropout_rate # Dropout rate
  )
  input <- cut_data_by_event(input, 125)
  input <- counting_process(input, arm = "experimental")
  output <- early_zero_weight(input, early_period = 2, fail_rate = fail_rate)

  observed <- output$weight
  expected <- rep(c(0, 0.8, 0, 0.7), c(43L, 20L, 29L, 33L))
  expect_equal(observed, expected)
})

