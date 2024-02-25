test_that("stratum values must be the same and stratum length must be doubled after converting", {
  res <- test_simfix2simpwsurv()
  fail_rate <- res$fail_rate
  failRatesPWSurv <- res$failRatesPWSurv
  dropoutRatesPWSurv <- res$dropoutRatesPWSurv

  stratum1 <- names(table(fail_rate$stratum))
  stratum2 <- names(table(failRatesPWSurv$stratum))

  expect_equal(stratum1, stratum2)
  expect_equal(length(fail_rate$stratum) * 2, length(failRatesPWSurv$stratum))

  stratum3 <- names(table(dropoutRatesPWSurv$stratum))
  stratum4 <- names(table(dropoutRatesPWSurv$stratum))

  expect_equal(stratum3, stratum4)
  expect_equal(length(fail_rate$stratum) * 2, length(dropoutRatesPWSurv$stratum))
})

test_that("treatment after converting contains only control and experimental with the right length", {
  res <- test_simfix2simpwsurv()
  fail_rate <- res$fail_rate
  failRatesPWSurv <- res$failRatesPWSurv
  dropoutRatesPWSurv <- res$dropoutRatesPWSurv

  expect_equal(length(fail_rate$stratum), sum(str_detect(failRatesPWSurv$treatment, "control")))
  expect_equal(length(fail_rate$stratum), sum(str_detect(failRatesPWSurv$treatment, "experimental")))
  expect_equal(length(failRatesPWSurv$treatment), length(fail_rate$stratum) * 2)
  expect_equal(length(fail_rate$stratum), sum(str_detect(dropoutRatesPWSurv$treatment, "control")))
  expect_equal(length(fail_rate$stratum), sum(str_detect(dropoutRatesPWSurv$treatment, "experimental")))
  expect_equal(length(dropoutRatesPWSurv$treatment), length(fail_rate$stratum) * 2)
})

test_that("Duration values match before and after converting and in right length ", {
  res <- test_simfix2simpwsurv()
  fail_rate <- res$fail_rate
  failRatesPWSurv <- res$failRatesPWSurv
  dropoutRatesPWSurv <- res$dropoutRatesPWSurv

  expect_equal(rep(c(fail_rate$duration), 2), failRatesPWSurv$duration)
  expect_equal(rep(c(fail_rate$duration), 2), dropoutRatesPWSurv$duration)
})

test_that("fail_rate match before and after converting and are in right length ", {
  res <- test_simfix2simpwsurv()
  fail_rate <- res$fail_rate
  failRatesPWSurv <- res$failRatesPWSurv

  expect_equal(fail_rate$fail_rate, failRatesPWSurv$rate[seq_along(fail_rate$fail_rate)])
  expect_equal(fail_rate$fail_rate * fail_rate$hr, failRatesPWSurv$rate[(length(fail_rate$fail_rate) + 1):(length(fail_rate$fail_rate) * 2)])
})

test_that("dropout_rate match before and after converting and are in right length ", {
  res <- test_simfix2simpwsurv()
  fail_rate <- res$fail_rate
  dropoutRatesPWSurv <- res$dropoutRatesPWSurv

  expect_equal(fail_rate$dropout_rate, dropoutRatesPWSurv$rate[seq_along(fail_rate$dropout_rate)])
  expect_equal(fail_rate$dropout_rate, dropoutRatesPWSurv$rate[(length(fail_rate$fail_rate) + 1):(length(fail_rate$fail_rate) * 2)])
})

# "meaningful error messages when the inputs are incorrect"
test_that("fail_rate column names must contain stratum, duration, fail_rate, hr and dropout_rate", {
  res <- test_simfix2simpwsurv()
  fail_rate <- res$fail_rate

  expect_equal(1, max(names(fail_rate) == "stratum"))
  expect_equal(1, max(names(fail_rate) == "duration"))
  expect_equal(1, max(names(fail_rate) == "fail_rate"))
  expect_equal(1, max(names(fail_rate) == "hr"))
  expect_equal(1, max(names(fail_rate) == "dropout_rate"))
})

test_that("duration must be longer than 0", {
  res <- test_simfix2simpwsurv()
  fail_rate <- res$fail_rate

  expect_equal(TRUE, is.numeric(fail_rate$duration))
  expect_gt(min(fail_rate$duration), 0)
})

test_that("fail_rate must be smaller than 1 and positive", {
  res <- test_simfix2simpwsurv()
  fail_rate <- res$fail_rate

  expect_lt(max(fail_rate$fail_rate), 1)
  expect_gt(min(fail_rate$fail_rate), 0)
})

test_that("hr must be postiive", {
  res <- test_simfix2simpwsurv()
  fail_rate <- res$fail_rate

  expect_gt(min(fail_rate$hr), 0)
})

test_that("dropout_rate must be smaller than 1 and positive", {
  res <- test_simfix2simpwsurv()
  fail_rate <- res$fail_rate

  expect_lt(max(fail_rate$dropout_rate), 1)
  expect_gt(min(fail_rate$dropout_rate), 0)
})
