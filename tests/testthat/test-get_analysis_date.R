# Simulate trial data
alpha <- 0.025
ratio <- 3
n <- 500
info_frac <- c(0.7, 1)
prevalence_ratio <- c(0.4, 0.6)
study_duration <- 48
# Two strata
stratum <- c("Biomarker-positive", "Biomarker-negative")
prevalence_ratio <- c(0.6, 0.4)
# enrollment rate
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
  duration = 1000,
  fail_rate = c(log(2) / c(med_pos, med_pos, med_neg, med_neg)),
  hr = c(hr_pos, hr_neg),
  dropout_rate = 0.01
)
# Simulate data
temp <- simfix2simpwsurv(fail_rate) # Convert the failure rate
set.seed(2023)
simulated_data <- sim_pw_surv(
  n = n, # Sample size
  # Stratified design with prevalence ratio of 6:4
  stratum = data.frame(stratum = stratum, p = prevalence_ratio),
  # Randomization ratio
  block = c("control", "control", "experimental", "experimental"),
  enroll_rate = enroll_rate, # Enrollment rate
  fail_rate = temp$fail_rate, # Failure rate
  dropout_rate = temp$dropout_rate # Dropout rate
)

test_that("planned_calendar_time", {
  observed <- get_analysis_date(simulated_data, planned_calendar_time = 24)
  expect_equal(observed, 24)
})

test_that("target_event_overall", {
  observed <- get_analysis_date(simulated_data, target_event_overall = 300)
  expect_equal(observed, 25.6150577826303)
})

test_that("planned_calendar_time + target_event_overall", {
  observed <- get_analysis_date(
    simulated_data,
    planned_calendar_time = 24,
    target_event_overall = 300
  )
  expect_equal(observed, 25.6150577826303)
})

test_that("target_event_per_stratum", {
  observed <- get_analysis_date(
    simulated_data,
    target_event_per_stratum = c(100, 200)
  )
  expect_equal(observed, 30.7886529927995)
})

test_that("target_event_overall + target_event_per_stratum", {
  observed <- get_analysis_date(
    simulated_data,
    target_event_overall = 150,
    target_event_per_stratum = c(100, NA)
  )
  expect_equal(observed, 18.3027216632238)
})

test_that("target_event_per_stratum + max_extension_for_target_event", {
  observed <- get_analysis_date(
    simulated_data,
    target_event_per_stratum = c(100, 200),
    max_extension_for_target_event = 30
  )
  expect_equal(observed, 30)
})

test_that("min_n_overall + min_followup", {
  observed <- get_analysis_date(
    simulated_data,
    min_n_overall = n * 0.8,
    min_followup = 12
  )
  expect_equal(observed, 28.8252059780061)
})

test_that("min_n_per_stratum + min_followup", {
  observed <- get_analysis_date(
    simulated_data,
    min_n_per_stratum = c(200, 160),
    min_followup = 12
  )
  expect_equal(observed, 27.3372780033387)
})

test_that("min_n_per_stratum + min_followup (requirement for only one stratum)", {
  observed <- get_analysis_date(
    simulated_data,
    min_n_per_stratum = c(200, NA),
    min_followup = 12
  )
  expect_equal(observed, 27.3372780033387)
})

test_that("min_n_overall + min_n_per_stratum + min_followup", {
  observed <- get_analysis_date(
    simulated_data,
    min_n_overall = n * 0.8,
    min_n_per_stratum = c(200, NA),
    min_followup = 12
  )
  expect_equal(observed, 28.8252059780061)
})

test_that("get_analysis_date() fails early with bad input", {
  # require non-negative number
  expect_error(get_analysis_date(simulated_data, planned_calendar_time = -1))
  expect_error(get_analysis_date(simulated_data, planned_calendar_time = -1))
  expect_error(get_analysis_date(simulated_data, planned_calendar_time = -1))
  expect_error(get_analysis_date(simulated_data, planned_calendar_time = -1))
  expect_error(get_analysis_date(simulated_data, planned_calendar_time = -1))
  expect_error(get_analysis_date(simulated_data, planned_calendar_time = -1))
  # require positive number
  expect_error(get_analysis_date(simulated_data, target_event_per_stratum = 0))
  expect_error(get_analysis_date(simulated_data, min_n_per_stratum = 0))
  # `min_n_overall` and `min_n_per_stratum` require `min_followup`.
  expect_error(
    get_analysis_date(simulated_data, min_n_overall = 1),
    "`min_followup` must be provided."
  )
  expect_error(
    get_analysis_date(simulated_data, min_n_per_stratum = 1),
    "`min_followup` must be provided."
  )
  # `min_n_overall` <= n
  expect_error(
    get_analysis_date(
      simulated_data,
      min_n_overall = n + 1,
      min_followup = 12
    ),
    "`min_n_overall` must be a positive number less than or equal to the total sample size."
  )
  # sum(min_n_per_stratum) <= n
  expect_error(
    get_analysis_date(
      simulated_data,
      min_n_per_stratum = c(n, 1),
      min_followup = 12
    ),
    "`min_n_per_stratum` must be a sum of positive numbers less than or equal to the total sample size."
  )
  # The following arguments require whole numbers: target_event_overall,
  # min_n_overall, target_event_per_stratum, min_n_per_stratum
  expect_error(
    get_analysis_date(simulated_data, target_event_overall = 300.1),
    "target_event_overall must be a single non-negative whole number \\(or NA\\)"
  )
  expect_error(
    get_analysis_date(
      simulated_data,
      min_n_overall = n * 0.8 + 0.1,
      min_followup = 12
    ),
    "min_n_overall must be a single non-negative whole number \\(or NA\\)"
  )
  expect_error(
    get_analysis_date(
      simulated_data,
      target_event_per_stratum = c(100.1, 200)
    ),
    "target_event_per_stratum must be a vector with only positive whole numbers and missing values"
  )
  expect_error(
    get_analysis_date(
      simulated_data,
      min_n_per_stratum = c(200.1, 160),
      min_followup = 12
    ),
    "min_n_per_stratum must be a vector with only positive whole numbers and missing values"
  )
})
