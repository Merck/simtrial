test_that("rpwenroll handles 0 enrollment rate properly for 1st enrollment period (with duration of 1 time unit)", {
  n <- 50
  enrollRates <- tibble(duration = c(1, 2), rate = c(0, 5))
  x <- rpwenroll(n=n, enrollRates=enrollRates)
  
  expect_gt(x[1],enrollRates$duration[1])
})

test_that("rpwenroll handles 0 enrollment rate properly for final enrollment period", {
  n <- 5
  enrollRates <- tibble(duration = c(1, 2), rate = c(10, 0))
  expect_error(rpwenroll(n=n, enrollRates=enrollRates), NA)
 
  n <- 50
  enrollRates <- tibble(duration = c(1, 2), rate = c(10, 0))
  expect_error(rpwenroll(n=n, enrollRates=enrollRates), regexp = "Please specify > 0 enrollment rate for the last period; otherwise enrollment cannot finish.")
  
  n <- 5
  enrollRates <- tibble(duration = 1, rate = 0)
  expect_error(rpwenroll(n=n, enrollRates=enrollRates), regexp = "Please specify > 0 enrollment rate, otherwise enrollment cannot finish.")
  
})


