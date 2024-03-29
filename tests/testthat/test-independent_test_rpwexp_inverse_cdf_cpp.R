test_that("rpwexp_inverse_cdf_cpp handles 0 fail_rate for final period", {
  # 0 failure rate for last period
  s <- rpwexp_inverse_cdf_cpp(
    n = 100,
    fail_rate = data.frame(duration = c(0, 2), rate = c(5, 0))
  )
  # Expect all failure times = Inf if length of first time period is 0
  expect_equal(min(s), Inf)
})

test_that("rpwexp_inverse_cdf_cpp handles 0 fail rate properly for one period", {
  # 0 failure rate
  s <- rpwexp_inverse_cdf_cpp(
    n = 100,
    fail_rate = data.frame(duration = c(1), rate = c(0))
  )
  # Expect Inf fail time
  expect_equal(mean(is.infinite(s)), 1)
})

test_that("rpwexp_inverse_cdf_cpp handles 0 fail rate properly for multiple periods", {
  # 0 failure rate for 1st period (with duration of 1 time unit)
  s <- rpwexp_inverse_cdf_cpp(
    n = 100,
    fail_rate = data.frame(duration = c(1, 2), rate = c(0, 5))
  )
  # Expect all failure times > = 1 (no events in 0 failure rate period)
  expect_equal(mean(s >= 1), 1)
})
