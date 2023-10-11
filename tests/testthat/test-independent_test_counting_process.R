surv_to_count <- function(time, status, trt, strats) {
  db <- data.frame(time, status, trt, strats)

  # KM estimator by stratum
  tidy_survfit <- function(...) {
    .survfit <- summary(survival::survfit(...))
    n <- length(.survfit$time)
    data.frame(
      time = .survfit$time,
      n.risk = .survfit$n.risk,
      n.event = .survfit$n.event,
      surv = c(1, .survfit$surv[-n])
    ) # ensure left continuous
  }
  km <- db %>%
    group_by(strats) %>%
    dplyr::do(tidy_survfit(Surv(time, status) ~ 1, data = .))

  # KM estimator by stratum and treatment Group Predicted at Specified Time
  pred_survfit <- function(pred_time, ...) {
    .survfit <- survival::survfit(...)

    # At risk subjects at pred_time
    n.risk <- stepfun(.survfit$time, c(.survfit$n.risk, 0), right = TRUE)(pred_time)
    .x1 <- data.frame(time = pred_time, n.risk)

    # Number of Event
    .x2 <- data.frame(time = .survfit$time, n.event = .survfit$n.event) %>% subset(n.event > 0)

    merge(.x1, .x2, all = TRUE) %>% mutate(n.event = dplyr::if_else(is.na(n.event), 0, n.event))
  }

  km_by_trt <- db %>%
    group_by(strats, trt) %>%
    dplyr::do(pred_time = pred_survfit(km[km$strats == .$strats[1], ]$time,
      Surv(time, status) ~ 1,
      data = .
    )) %>%
    tidyr::unnest(cols = pred_time) %>%
    dplyr::rename(tn.risk = n.risk, tn.event = n.event)


  # Log Rank Expectation Difference and Variance
  res <- merge(km, km_by_trt, all = TRUE) %>%
    dplyr::arrange(trt, strats, time) %>%
    mutate(
      OminusE = tn.event - tn.risk / n.risk * n.event,
      Var = (n.risk - tn.risk) * tn.risk * n.event * (n.risk - n.event) / n.risk^2 / (n.risk - 1)
    )
}

testthat::test_that("Counting Process Format without ties", {
  x <- tibble(
    stratum = c(rep(1, 10), rep(2, 6)),
    treatment = rep(c(1, 1, 0, 0), 4),
    tte = 1:16,
    event = rep(c(0, 1), 8)
  )

  arm <- 1
  res_counting_process <- simtrial::counting_process(x, arm)
  res_test <- surv_to_count(time = x$tte, status = x$event, trt = x$treatment, strats = x$stratum)

  res_test <- tibble::as_tibble(subset(res_test, trt == 1)) %>%
    subset(n.event > 0 & n.risk - tn.risk > 0 & tn.risk > 0)

  testthat::expect_equal(res_counting_process$o_minus_e, res_test$OminusE)
  testthat::expect_equal(res_counting_process$var_o_minus_e, res_test$Var)
})


testthat::test_that("Counting Process Format with ties", {
  x <- tibble(
    stratum = c(rep(1, 10), rep(2, 6)),
    treatment = rep(c(1, 1, 0, 0), 4),
    tte = c(rep(1:4, each = 4)),
    event = rep(c(0, 1), 8)
  )
  arm <- 1
  res_counting_process <- counting_process(x, arm)
  res_test <- surv_to_count(time = x$tte, status = x$event, trt = x$treatment, strats = x$stratum)

  res_test <- tibble::as_tibble(subset(res_test, trt == 1)) %>%
    subset(n.event > 0 & n.risk - tn.risk > 0 & tn.risk > 0)

  testthat::expect_equal(res_counting_process$o_minus_e, res_test$OminusE)
  testthat::expect_equal(res_counting_process$var_o_minus_e, res_test$Var)
})
