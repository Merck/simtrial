fail_rate <- tibble::tibble(
  stratum = c(rep("Low", 3), rep("High", 3)),
  duration = rep(c(4, 10, 100), 2),
  fail_rate = c(
    .04, .1, .06,
    .08, .16, .12
  ),
  hr = c(
    1.5, .5, 2 / 3,
    2, 10 / 16, 10 / 12
  ),
  dropout_rate = .01
)

failRatesPWSurv <- simfix2simpwsurv(fail_rate)$fail_rate
dropoutRatesPWSurv <- simfix2simpwsurv(fail_rate)$dropout_rate

testthat::test_that("stratgum values must be the same and stratum length must be doubled after converting", {
  stratum1 <- names(table(fail_rate$stratum))
  stratum2 <- names(table(failRatesPWSurv$stratum))
  testthat::expect_equal(stratum1, stratum2)
  testthat::expect_equal(length(fail_rate$stratum) * 2, length(failRatesPWSurv$stratum))
  stratum3 <- names(table(dropoutRatesPWSurv$stratum))
  stratum4 <- names(table(dropoutRatesPWSurv$stratum))
  testthat::expect_equal(stratum3, stratum4)
  testthat::expect_equal(length(fail_rate$stratum) * 2, length(dropoutRatesPWSurv$stratum))
})

testthat::test_that("treatment after converting contains only control and experimental with the right length", {
  testthat::expect_equal(length(fail_rate$stratum), sum(stringr::str_detect(failRatesPWSurv$treatment, "control")))
  testthat::expect_equal(length(fail_rate$stratum), sum(stringr::str_detect(failRatesPWSurv$treatment, "experimental")))
  testthat::expect_equal(length(failRatesPWSurv$treatment), length(fail_rate$stratum) * 2)
  testthat::expect_equal(length(fail_rate$stratum), sum(stringr::str_detect(dropoutRatesPWSurv$treatment, "control")))
  testthat::expect_equal(length(fail_rate$stratum), sum(stringr::str_detect(dropoutRatesPWSurv$treatment, "experimental")))
  testthat::expect_equal(length(dropoutRatesPWSurv$treatment), length(fail_rate$stratum) * 2)
})

testthat::test_that("Duration values match before and after converting and in right length ", {
  testthat::expect_equal(rep(c(fail_rate$duration), 2), failRatesPWSurv$duration)
  testthat::expect_equal(rep(c(fail_rate$duration), 2), dropoutRatesPWSurv$duration)
})

testthat::test_that("fail_rate match before and after converting and are in right length ", {
  testthat::expect_equal(fail_rate$fail_rate, failRatesPWSurv$rate[1:length(fail_rate$fail_rate)])
  testthat::expect_equal(fail_rate$fail_rate * fail_rate$hr, failRatesPWSurv$rate[(length(fail_rate$fail_rate) + 1):(length(fail_rate$fail_rate) * 2)])
})

testthat::test_that("dropout_rate match before and after converting and are in right length ", {
  testthat::expect_equal(fail_rate$dropout_rate, dropoutRatesPWSurv$rate[1:length(fail_rate$dropout_rate)])
  testthat::expect_equal(fail_rate$dropout_rate, dropoutRatesPWSurv$rate[(length(fail_rate$fail_rate) + 1):(length(fail_rate$fail_rate) * 2)])
})


# "meaningful error messages when the inputs are incorrect"
testthat::test_that("fail_rate column names must contain stratum, duration, fail_rate, hr and dropout_rate", {
  testthat::expect_equal(1, max(names(fail_rate) == "stratum"))
  testthat::expect_equal(1, max(names(fail_rate) == "duration"))
  testthat::expect_equal(1, max(names(fail_rate) == "fail_rate"))
  testthat::expect_equal(1, max(names(fail_rate) == "hr"))
  testthat::expect_equal(1, max(names(fail_rate) == "dropout_rate"))
})

testthat::test_that("duration must be longer than 0", {
  testthat::expect_equal(TRUE, is.numeric(fail_rate$duration))
  testthat::expect_gt(min(fail_rate$duration), 0)
})

testthat::test_that("fail_rate must be smaller than 1 and positive", {
  testthat::expect_lt(max(fail_rate$fail_rate), 1)
  testthat::expect_gt(min(fail_rate$fail_rate), 0)
})


testthat::test_that("hr must be postiive", {
  testthat::expect_gt(min(fail_rate$hr), 0)
})


testthat::test_that("dropout_rate must be smaller than 1 and positive", {
  testthat::expect_lt(max(fail_rate$dropout_rate), 1)
  testthat::expect_gt(min(fail_rate$dropout_rate), 0)
})
