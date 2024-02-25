test_that("stratum percentage calculated from simulated dataset must be within
          the tolerance=0.002 of stratum in setup (0.4,0.6)", {
  res <- test_sim_pw_surv()
  expect_equal(
    object = c(
      sum(str_count(res$x$stratum, "Low")) / 400000,
      sum(str_count(res$x$stratum, "High")) / 400000
    ),
    expected = c(0.4, 0.6), tolerance = 0.002
  )
})

test_that("block calculated from simulated dataset equals size of 4 with 1:1
          randomization, which is 2 for each arm", {
  res <- test_sim_pw_surv()
  expect_equal(object = res$bktest1, expected = rep(2, length(res$bktest1)))
  expect_equal(object = res$bktest2, expected = rep(2, length(res$bktest2)))
})

test_that("fail_rate calculated from simulated dataset must be within the
          tolerance=0.1 of fail_rate in setting", {
  res <- test_sim_pw_surv()
  expect_equal(object = res$ratetest, expected = res$fail_rate$rate, tolerance = 0.1)
})

test_that("dropout_rate calculated from simulated dataset must be within
          the tolerance=0.0005 of dropout_rate=0.001 in setup", {
  res <- test_sim_pw_surv()
  duration <- 300
  drtest <- 0
  for (i in 1:duration) {
    drtest[i] <- sum(res$x$dropout_time <= i & res$x$dropout_time > (i - 1)) / 400000
  }
  expect_equal(object = drtest, expected = rep(0.001, 300), tolerance = 0.001)
})

test_that("enroll_rate calculated from simulated dataset must be within
          the relative tolerance=0.05 of enroll_rate in setup", {
  res <- test_sim_pw_surv()
  duration <- 300
  entest <- 0
  for (i in 1:duration) {
    entest[i] <- sum(res$x$enroll_time <= i & res$x$enroll_time > (i - 1))
  }
  entest1 <- entest[entest != 0]
  entestexp <- c(rep(100, 5), rep(3000, length(entest1) - 5))
  entest2 <- (entest1 - entestexp) / entestexp
  expect_equal(object = entest2, expected = rep(0, length(entest1)), tolerance = 0.05)
})

test_that("The actual number of events changes by changing total sample size", {
  res1 <- test_sim_pw_surv()
  res2 <- test_sim_pw_surv_2()
  expect_false(unique(res1$xevent$event == res2$zevent$event))
})

test_that("sim_pw_surv() fails early with mismatched treatment names", {
  block <- c(rep("x", 2), rep("y", 2))
  fail_rate <- data.frame(
    stratum = rep("All", 4),
    period = rep(1:2, 2),
    treatment = c(rep("x", 2), rep("y", 2)),
    duration = rep(c(3, 1), 2),
    rate = log(2) / c(9, 9, 9, 18)
  )
  dropout_rate <- data.frame(
    stratum = rep("All", 2),
    period = rep(1, 2),
    treatment = c("x", "y"),
    duration = rep(100, 2),
    rate = rep(0.001, 2)
  )

  expect_error(sim_pw_surv(block = block))
  expect_error(sim_pw_surv(fail_rate = fail_rate))
  expect_error(sim_pw_surv(dropout_rate = dropout_rate))
  # Works as long as treatment names are consistent
  expect_silent(
    xy <- sim_pw_surv(block = block, fail_rate = fail_rate, dropout_rate = dropout_rate)
  )
  expect_identical(sort(unique(xy$treatment)), c("x", "y"))
})
