# library(PWEALL)

# # test when set seed
# test_that("rpwexpinvRcpp handles 0 failrate for final period", {
#
#   n_rates <- 10
#   duration <- floor(runif(n_rates, 1, 5))
#   rate <- runif(n_rates, 0.0005, 0.095)
#   rate[n_rates] <- 0.0 # 0 failure rate for last period
#   tchange0 <- c(0, cumsum(duration))[1:(length(duration-1))]
#
#   n = 1000
#   set.seed(2022)
#   ref <- PWEALL::rpwe(n = n, rate = rate, tchange = tchange0)
#   s <- simtrial::rpwexpinvRcpp(n = n, rate = rate, duration = duration)
#
#   testthat::expect_equal(ref$r, s)
#
#   # expect Inf fail time exists
#   testthat::expect_equal(is.infinite(max(s)), TRUE)
# })
#
# test_that("rpwexpinvRcpp handles 0 failrate for final period", {
#   # 0 failure rate for last period
#   s <- simtrial::rpwexpinvRcpp(n = n, rate = c(5, 0), duration = c(0, 2))
#   # expect all failure times = Inf if length of first time period is 0
#   testthat::expect_equal(min(s),Inf)
# })
#
# testthat::test_that("rpwexpinvRcpp handles 0 fail rate properly for one period",{
#   # 0 failure rate
#   s <- simtrial::rpwexpinvRcpp(n = n, rate = c(0), duration = c(1))
#
#   # expect Inf fail time
#   testthat::expect_equal(mean(is.infinite(s)),1)
# })
#
# testthat::test_that("rpwexpinvRcpp handles 0 fail rate properly for multiple periods",{
#   # 0 failure rate for 1st period (with duration of 1 time unit)
#   s <- simtrial::rpwexpinvRcpp(n = n, rate = c(0, 5), duration = c(1, 2))
#
#   # expect all failure times > = 1 (no events in 0 failure rate period)
#   testthat::expect_equal(mean(s>=1),1)
# })
#
# test_that("rpwexpRcpp handles 0 failrate for final period", {
#   # 0 failure rate for last period
#   s <- simtrial::rpwexp(n=100,
#                         failRates=tibble(duration=c(0, 2),rate=c(5, 0)))
#   sRcpp <- simtrial::rpwexpRcpp(n=100,
#                                 failRates=tibble(duration=c(0, 2),rate=c(5, 0)))
#   # expect all failure times = Inf if length of first time period is 0
#   testthat::expect_equal(sRcpp,s)
# })
#
# testthat::test_that("rpwexpRcpp handles 0 fail rate properly for one period",{
#   # 0 failure rate
#   s <- simtrial::rpwexp(n=100,
#               failRates=tibble(duration=c(1),rate=c(0)))
#   sRcpp <- simtrial::rpwexpRcpp(n=100,
#                                 failRates=tibble(duration=c(1),rate=c(0)))
#   # expect Inf fail time
#   testthat::expect_equal(sRcpp,s)
# })
#
#
# testthat::test_that("rpwexpRcpp handles 0 fail rate properly for multiple periods",{
#   # 0 failure rate for 1st period (with duration of 1 time unit)
#   set.seed(2022)
#   s <- simtrial::rpwexp(n=100,
#               failRates=tibble(duration=c(1, 2),rate=c(0, 5)))
#   sRcpp <- simtrial::rpwexpRcpp(n=100,
#                                 failRates=tibble(duration=c(1, 2),rate=c(0, 5)))
#   # expect all failure times > = 1 (no events in 0 failure rate period)
#   testthat::expect_equal(sRcpp,s)
# })

