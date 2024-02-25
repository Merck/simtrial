test_that("planned_calendar_time", {
  simulated_data <- test_get_analysis_date()$simulated_data
  observed <- get_analysis_date(simulated_data, planned_calendar_time = 24)
  expect_equal(observed, 24)
})

test_that("target_event_overall", {
  simulated_data <- test_get_analysis_date()$simulated_data
  observed <- get_analysis_date(simulated_data, target_event_overall = 300)
  expect_equal(observed, 25.6150577826303)
})

test_that("planned_calendar_time + target_event_overall", {
  simulated_data <- test_get_analysis_date()$simulated_data
  observed <- get_analysis_date(
    simulated_data,
    planned_calendar_time = 24,
    target_event_overall = 300
  )
  expect_equal(observed, 25.6150577826303)
})

test_that("target_event_per_stratum", {
  simulated_data <- test_get_analysis_date()$simulated_data
  observed <- get_analysis_date(
    simulated_data,
    target_event_per_stratum = c(100, 200)
  )
  expect_equal(observed, 30.7886529927995)
})

test_that("target_event_overall + target_event_per_stratum", {
  simulated_data <- test_get_analysis_date()$simulated_data
  observed <- get_analysis_date(
    simulated_data,
    target_event_overall = 150,
    target_event_per_stratum = c(100, NA)
  )
  expect_equal(observed, 18.3027216632238)
})

test_that("target_event_per_stratum + max_extension_for_target_event", {
  simulated_data <- test_get_analysis_date()$simulated_data
  observed <- get_analysis_date(
    simulated_data,
    target_event_per_stratum = c(100, 200),
    max_extension_for_target_event = 30
  )
  expect_equal(observed, 30)
})

test_that("min_n_overall + min_followup", {
  res <- test_get_analysis_date()
  simulated_data <- res$simulated_data
  n <- res$n
  observed <- get_analysis_date(
    simulated_data,
    min_n_overall = n * 0.8,
    min_followup = 12
  )
  expect_equal(observed, 28.8252059780061)
})

test_that("min_n_per_stratum + min_followup", {
  simulated_data <- test_get_analysis_date()$simulated_data
  observed <- get_analysis_date(
    simulated_data,
    min_n_per_stratum = c(200, 160),
    min_followup = 12
  )
  expect_equal(observed, 27.3372780033387)
})

test_that("min_n_per_stratum + min_followup (requirement for only one stratum)", {
  simulated_data <- test_get_analysis_date()$simulated_data
  observed <- get_analysis_date(
    simulated_data,
    min_n_per_stratum = c(200, NA),
    min_followup = 12
  )
  expect_equal(observed, 27.3372780033387)
})

test_that("min_n_overall + min_n_per_stratum + min_followup", {
  res <- test_get_analysis_date()
  simulated_data <- res$simulated_data
  n <- res$n
  observed <- get_analysis_date(
    simulated_data,
    min_n_overall = n * 0.8,
    min_n_per_stratum = c(200, NA),
    min_followup = 12
  )
  expect_equal(observed, 28.8252059780061)
})

test_that("get_analysis_date() fails early with bad input", {
  res <- test_get_analysis_date()
  simulated_data <- res$simulated_data
  n <- res$n

  # Require non-negative number
  expect_error(get_analysis_date(simulated_data, planned_calendar_time = -1))
  expect_error(get_analysis_date(simulated_data, planned_calendar_time = -1))
  expect_error(get_analysis_date(simulated_data, planned_calendar_time = -1))
  expect_error(get_analysis_date(simulated_data, planned_calendar_time = -1))
  expect_error(get_analysis_date(simulated_data, planned_calendar_time = -1))
  expect_error(get_analysis_date(simulated_data, planned_calendar_time = -1))

  # Require positive number
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

  # The following arguments require whole numbers:
  # - target_event_overall
  # - min_n_overall
  # - target_event_per_stratum
  # - min_n_per_stratum
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
