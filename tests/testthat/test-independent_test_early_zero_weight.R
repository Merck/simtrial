library(dplyr)
library(gsDesign2)
library(tibble)


test_that("early_zero_weight() with unstratified data", {
  set.seed(123456)
  early_period = 2

  output <- sim_pw_surv(n = 200) |>
    cut_data_by_event(125) |>
    counting_process(arm = "experimental") |>
    early_zero_weight(early_period = early_period)

  observed <- output$weight
  expected <- if_else(output$tte < early_period,0,1)
  expect_equal(observed, expected)
})


test_that("early_zero_weight() with stratified data when fail_rate is not provided", {
  set.seed(123456)
  early_period = 2
  n <- 500
  # Two strata
  stratum <- c("Biomarker-positive", "Biomarker-negative")
  prevalence_ratio <- c(0.6, 0.4)
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
  temp <- to_sim_pw_surv(fail_rate) # Convert the failure rate
  x <- sim_pw_surv(
    n = n, # Sample size
    # Stratified design with prevalence ratio of 6:4
    stratum = tibble(stratum = stratum, p = prevalence_ratio),
    # Randomization ratio
    block = c("control", "control", "experimental", "experimental"),
    enroll_rate = enroll_rate, # Enrollment rate
    fail_rate = temp$fail_rate, # Failure rate
    dropout_rate = temp$dropout_rate # Dropout rate
  ) |>
    cut_data_by_event(125) |>
    counting_process(arm = "experimental")

  expect_error(early_zero_weight(x, early_period = early_period),"For stratified design to use `early_zero_weight\\(\\)`, `fail_rate` can't be `NULL`.")
})


test_that("early_zero_weight() with stratified data when fail_rate is not correctly provided", {
  set.seed(123456)
  early_period = 2
  n <- 500
  # Two strata
  stratum <- c("Biomarker-positive", "Biomarker-negative")
  prevalence_ratio <- c(0.6, 0.4)
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
  fail_rate_wrong <- tibble(
    stratum = c(rep('Biomarker-positive',3), rep('Biomarker-negative',2)),
    duration = c(3, 4, 1000, 4, 1000),
    fil_rate = c(log(2) / c(med_pos, med_pos, med_pos, med_neg, med_neg)),
    dropout_rate = rep(0.01,5),
    hr = c(1, 1, 0.7, 1, 0.8),
  )
  temp <- to_sim_pw_surv(fail_rate) # Convert the failure rate
  x <- sim_pw_surv(
    n = n, # Sample size
    # Stratified design with prevalence ratio of 6:4
    stratum = tibble(stratum = stratum, p = prevalence_ratio),
    # Randomization ratio
    block = c("control", "control", "experimental", "experimental"),
    enroll_rate = enroll_rate, # Enrollment rate
    fail_rate = temp$fail_rate, # Failure rate
    dropout_rate = temp$dropout_rate # Dropout rate
  ) |>
    cut_data_by_event(125) |>
    counting_process(arm = "experimental")

  expect_error(early_zero_weight(x, early_period = early_period, fail_rate = fail_rate_wrong),"`early_zero_weight\\(\\)` only allows delayed treatment effect, that is, 2 piece failure rate with HR = 1 at the first period.")
})


test_that("early_zero_weight() with stratified data when fail_rate is correctly provided", {
  set.seed(123456)
  early_period = 2
  n <- 500
  # Two strata
  stratum <- c("Biomarker-positive", "Biomarker-negative")
  prevalence_ratio <- c(0.6, 0.4)
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
  fail_rate_wrong <- tibble(
    stratum = c(rep('Biomarker-positive',3), rep('Biomarker-negative',2)),
    duration = c(3, 4, 1000, 4, 1000),
    fil_rate = c(log(2) / c(med_pos, med_pos, med_pos, med_neg, med_neg)),
    dropout_rate = rep(0.01,5),
    hr = c(1, 1, 0.7, 1, 0.8),
  )
  temp <- to_sim_pw_surv(fail_rate) # Convert the failure rate
  output <- sim_pw_surv(
    n = n, # Sample size
    # Stratified design with prevalence ratio of 6:4
    stratum = tibble(stratum = stratum, p = prevalence_ratio),
    # Randomization ratio
    block = c("control", "control", "experimental", "experimental"),
    enroll_rate = enroll_rate, # Enrollment rate
    fail_rate = temp$fail_rate, # Failure rate
    dropout_rate = temp$dropout_rate # Dropout rate
  ) |>
    cut_data_by_event(125) |>
    counting_process(arm = "experimental") |>
    early_zero_weight(early_period = early_period, fail_rate = fail_rate)

  observed <- output$weight
  expected <- if_else(output$stratum=='Biomarker-negative',if_else(output$tte<4,0,log(0.8)),if_else(output$tte<3,0,log(0.7)))
  expect_equal(observed, expected)
})
