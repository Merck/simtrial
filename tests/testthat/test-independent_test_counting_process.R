test_that("Counting Process Format without ties", {
  x <- data.frame(
    stratum = c(rep(1, 10), rep(2, 6)),
    treatment = rep(c(1, 1, 0, 0), 4),
    tte = 1:16,
    event = rep(c(0, 1), 8)
  )

  arm <- 1
  res_counting_process <- counting_process(x, arm)
  res_test <- surv_to_count(time = x$tte, status = x$event, trt = x$treatment, strats = x$stratum)

  res_test <- data.frame(subset(res_test, trt == 1)) |>
    subset(n.event > 0 & n.risk - tn.risk > 0 & tn.risk > 0)

  expect_equal(res_counting_process$o_minus_e, res_test$OminusE)
  expect_equal(res_counting_process$var_o_minus_e, res_test$Var)
})

test_that("Counting Process Format with ties", {
  x <- data.frame(
    stratum = c(rep(1, 10), rep(2, 6)),
    treatment = rep(c(1, 1, 0, 0), 4),
    tte = c(rep(1:4, each = 4)),
    event = rep(c(0, 1), 8)
  )
  arm <- 1
  res_counting_process <- counting_process(x, arm)
  res_test <- surv_to_count(time = x$tte, status = x$event, trt = x$treatment, strats = x$stratum)

  res_test <- data.frame(subset(res_test, trt == 1)) |>
    subset(n.event > 0 & n.risk - tn.risk > 0 & tn.risk > 0)

  expect_equal(res_counting_process$o_minus_e, res_test$OminusE)
  expect_equal(res_counting_process$var_o_minus_e, res_test$Var)
})
