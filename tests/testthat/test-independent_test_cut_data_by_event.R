test_that("the input is a time-to-event data set", {
  TTEdata <- sim_pw_surv(n = 200)

  expect_equal(1, max(names(TTEdata) == "stratum"))
  expect_equal(1, max(names(TTEdata) == "enroll_time"))
  expect_equal(1, max(names(TTEdata) == "treatment"))
  expect_equal(1, max(names(TTEdata) == "fail_time"))
  expect_equal(1, max(names(TTEdata) == "dropout_time"))
  expect_equal(1, max(names(TTEdata) == "fail"))
  expect_equal(1, max(names(TTEdata) == "cte"))
})
