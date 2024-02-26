# Helper functions used by test-independent_test_simfix2simpwsurv.R

test_simfix2simpwsurv <- function() {
  fail_rate <- data.frame(
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

  failRatesPWSurv <- to_sim_pw_surv(fail_rate)$fail_rate
  dropoutRatesPWSurv <- to_sim_pw_surv(fail_rate)$dropout_rate

  list(
    "fail_rate" = fail_rate,
    "failRatesPWSurv" = failRatesPWSurv,
    "dropoutRatesPWSurv" = dropoutRatesPWSurv
  )
}
