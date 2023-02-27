test_that("rpw_enroll handles 0 enrollment rate properly for 1st enrollment period (with duration of 1 time unit)", {
  n <- 50
  enrollRates <- tibble(duration = c(1, 2), rate = c(0, 5))
  x <- rpw_enroll(n=n, enrollRates=enrollRates)

  expect_gt(x[1],enrollRates$duration[1])
})

test_that("rpw_enroll handles 0 enrollment rate properly for final enrollment period", {

  set.seed(123)
  n <- 50
  enrollRates <- tibble(duration = c(1, 2), rate = c(10, 0))
  expect_error(rpw_enroll(n=n, enrollRates=enrollRates))

  n <- 5
  enrollRates <- tibble(duration = 1, rate = 0)
  expect_error(rpw_enroll(n=n, enrollRates=enrollRates))

})


