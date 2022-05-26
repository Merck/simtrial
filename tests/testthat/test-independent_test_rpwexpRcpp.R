test_that("rpwexpRcpp handles 0 failrate for final period", {
  # 0 failure rate for last period
  s <- simtrial::rpwexp(n=100,
                        failRates=tibble(duration=c(0, 2),rate=c(5, 0)))
  sRcpp <- simtrial::rpwexpRcpp(n=100,
                                failRates=tibble(duration=c(0, 2),rate=c(5, 0)))
  # expect all failure times = Inf if length of first time period is 0
  testthat::expect_equal(sRcpp,s)
})

testthat::test_that("rpwexpRcpp handles 0 fail rate properly for one period",{
  # 0 failure rate
  s <- simtrial::rpwexp(n=100,
              failRates=tibble(duration=c(1),rate=c(0)))
  sRcpp <- simtrial::rpwexpRcpp(n=100,
                                failRates=tibble(duration=c(1),rate=c(0)))
  # expect Inf fail time
  testthat::expect_equal(sRcpp,s)
})


testthat::test_that("rpwexpRcpp handles 0 fail rate properly for multiple periods",{
  # 0 failure rate for 1st period (with duration of 1 time unit)
  set.seed(2022)
  s <- simtrial::rpwexp(n=100,
              failRates=tibble(duration=c(1, 2),rate=c(0, 5)))
  sRcpp <- simtrial::rpwexpRcpp(n=100,
                                failRates=tibble(duration=c(1, 2),rate=c(0, 5)))
  # expect all failure times > = 1 (no events in 0 failure rate period)
  testthat::expect_equal(sRcpp,s)
})

