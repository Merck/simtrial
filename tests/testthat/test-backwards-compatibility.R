skip_if_offline()
skip_on_cran()

# see generate-backwards-compatible-test-data.R for how test data is generated
# with a previous version of simtrial
reference <- Sys.getenv("SIMTRIAL_TEST_BACKWARDS_COMPATIBILITY_REF")
if (nchar(reference) > 0) {
  source("generate-backwards-compatible-test-data.R")
  generate_test_data(reference = reference, outdir = "fixtures/backwards-compatibility")
  library("simtrial")
} else {
  skip(message = "Not testing backwards compatibility")
}

test_that("cut_data_by_date()", {
  set.seed(12345)
  observed <- cut_data_by_date(
    x = sim_pw_surv(n = 20),
    cut_date = 5
  )
  expected <- readRDS("fixtures/backwards-compatibility/cut_data_by_date_ex1.rds")
  expect_equivalent(as.data.frame(observed), as.data.frame(expected))
})

test_that("get_cut_date_by_event()", {
  set.seed(12345)
  x <- sim_pw_surv(
    n = 200,
    stratum = data.frame(
      stratum = c("Positive", "Negative"),
      p = c(.5, .5)
    ),
    fail_rate = data.frame(
      stratum = rep(c("Positive", "Negative"), 2),
      period = rep(1, 4),
      treatment = c(rep("control", 2), rep("experimental", 2)),
      duration = rep(1, 4),
      rate = log(2) / c(6, 9, 9, 12)
    ),
    dropout_rate = data.frame(
      stratum = rep(c("Positive", "Negative"), 2),
      period = rep(1, 4),
      treatment = c(rep("control", 2), rep("experimental", 2)),
      duration = rep(1, 4),
      rate = rep(.001, 4)
    )
  )
  observed <- get_cut_date_by_event(subset(x, stratum == "Positive"), event = 50)
  expected <- readRDS("fixtures/backwards-compatibility/get_cut_date_by_event_ex1.rds")
  expect_equivalent(as.data.frame(observed), as.data.frame(expected))
})

test_that("counting_process()", {
  # Example 1
  x <- data.frame(
    stratum = c(rep(1, 10), rep(2, 6)),
    treatment = rep(c(1, 1, 0, 0), 4),
    tte = 1:16,
    event = rep(c(0, 1), 8)
  )
  observed <- counting_process(x, arm = 1)
  expected <- readRDS("fixtures/backwards-compatibility/counting_process_ex1.rds")
  expect_equivalent(as.data.frame(observed), as.data.frame(expected))

  # Example 2
  set.seed(12345)
  x <- sim_pw_surv(n = 400)
  y <- cut_data_by_event(x, 150)
  observed <- counting_process(y, arm = "experimental")
  expected <- readRDS("fixtures/backwards-compatibility/counting_process_ex2.rds")
  expect_equivalent(as.data.frame(observed), as.data.frame(expected))

  # Example 3
  # Counting Process Format with ties
  x <- data.frame(
    stratum = c(rep(1, 10), rep(2, 6)),
    treatment = rep(c(1, 1, 0, 0), 4),
    tte = c(rep(1:4, each = 4)),
    event = rep(c(0, 1), 8)
  )
  arm <- 1
  observed <- counting_process(x, arm)
  expected <- readRDS("fixtures/backwards-compatibility/counting_process_ex3.rds")
  expect_equivalent(as.data.frame(observed), as.data.frame(expected))
})

test_that("fh_weight()", {
  # Example 1
  # Use default enrollment and event rates at cut at 100 events
  set.seed(12345)
  x <- sim_pw_surv(n = 200)
  x <- cut_data_by_event(x, 100)
  x <- counting_process(x, arm = "experimental")

  # Compute the corvariance between FH(0, 0), FH(0, 1) and FH(1, 0)
  observed <- fh_weight(x, rho_gamma = data.frame(rho = c(0, 0, 1), gamma = c(0, 1, 0)))
  expected <- readRDS("fixtures/backwards-compatibility/wlr_ex1.rds")
  expect_equivalent(as.data.frame(observed), as.data.frame(expected))
  observed <- fh_weight(x, rho_gamma = data.frame(rho = c(0, 0, 1), gamma = c(0, 1, 0)), return_variance = TRUE)
  expected <- readRDS("fixtures/backwards-compatibility/wlr_ex1_var.rds")
  expect_equivalent(as.data.frame(observed), as.data.frame(expected))
  observed <- fh_weight(x, rho_gamma = data.frame(rho = c(0, 0, 1), gamma = c(0, 1, 0)), return_corr = TRUE)
  expected <- readRDS("fixtures/backwards-compatibility/wlr_ex1_cor.rds")
  expect_equivalent(as.data.frame(observed), as.data.frame(expected))

  # Example 2
  # Use default enrollment and event rates at cut of 100 events
  set.seed(12345)
  x <- sim_pw_surv(n = 200)
  x <- cut_data_by_event(x, 100)
  x <- counting_process(x, arm = "experimental")
  observed <- fh_weight(x, rho_gamma = data.frame(rho = c(0, 0), gamma = c(0, 1)), return_corr = TRUE)
  expected <- readRDS("fixtures/backwards-compatibility/wlr_ex2.rds")
  expect_equivalent(as.data.frame(observed), as.data.frame(expected))
})

test_that("rpw_enroll()", {
  set.seed(12345)
  observed <- rpw_enroll(
    n = 1e5,
    enroll_rate = data.frame(
      rate = c(5, 15, 30),
      duration = c(100, 200, 100)
    )
  )
  expected <- readRDS("fixtures/backwards-compatibility/rpw_enroll_ex1.rds")
  expect_equal(observed, expected)

  # Example 2
  # Exponential enrollment
  set.seed(12345)
  observed <- rpw_enroll(
    n = 1e5,
    enroll_rate = data.frame(rate = .03, duration = 1)
  )
  expected <- readRDS("fixtures/backwards-compatibility/rpw_enroll_ex2.rds")
  expect_equal(observed, expected)
})

test_that("simfix2simpwsurv()", {
  # Example 1
  # Convert standard input
  observed <- simfix2simpwsurv()
  expected <- readRDS("fixtures/backwards-compatibility/simfix2simpwsurv_ex1.rds")
  expect_equivalent(
    as.data.frame(observed$fail_rate),
    as.data.frame(expected$fail_rate)
  )
  expect_equivalent(
    as.data.frame(observed$dropout_rate),
    as.data.frame(expected$dropout_rate)
  )

  # Example 2
  # Stratified example
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
  observed <- simfix2simpwsurv(fail_rate)
  expected <- readRDS("fixtures/backwards-compatibility/simfix2simpwsurv_ex2.rds")
  expect_equivalent(
    as.data.frame(observed$fail_rate),
    as.data.frame(expected$fail_rate)
  )
  expect_equivalent(
    as.data.frame(observed$dropout_rate),
    as.data.frame(expected$dropout_rate)
  )
})

test_that("sim_pw_surv()", {
  # Example 1
  set.seed(12345)
  observed <- sim_pw_surv(n = 20)
  expected <- readRDS("fixtures/backwards-compatibility/sim_pw_surv_ex1.rds")
  expect_equivalent(as.data.frame(observed), as.data.frame(expected))

  # Example 2
  # 3:1 randomization
  set.seed(12345)
  observed <- sim_pw_surv(
    n = 20,
    block = c(rep("experimental", 3), "control")
  )
  expected <- readRDS("fixtures/backwards-compatibility/sim_pw_surv_ex2.rds")
  expect_equivalent(as.data.frame(observed), as.data.frame(expected))

  # Example 3
  # Simulate 2 stratum; will use defaults for blocking and enrollRates
  set.seed(12345)
  observed <- sim_pw_surv(
    n = 20,
    # 2 stratum,30% and 70% prevalence
    stratum = data.frame(stratum = c("Low", "High"), p = c(.3, .7)),
    fail_rate = data.frame(
      stratum = c(rep("Low", 4), rep("High", 4)),
      period = rep(1:2, 4),
      treatment = rep(c(
        rep("control", 2),
        rep("experimental", 2)
      ), 2),
      duration = rep(c(3, 1), 4),
      rate = c(.03, .05, .03, .03, .05, .08, .07, .04)
    ),
    dropout_rate = data.frame(
      stratum = c(rep("Low", 2), rep("High", 2)),
      period = rep(1, 4),
      treatment = rep(c("control", "experimental"), 2),
      duration = rep(1, 4),
      rate = rep(.001, 4)
    )
  )
  expected <- readRDS("fixtures/backwards-compatibility/sim_pw_surv_ex3.rds")
  expect_equivalent(as.data.frame(observed), as.data.frame(expected))

  # Example 4
  # If you want a more rectangular entry for a data frame
  fail_rate <- list(
    data.frame(stratum = "Low", period = 1, treatment = "control", duration = 3, rate = .03),
    data.frame(stratum = "Low", period = 1, treatment = "experimental", duration = 3, rate = .03),
    data.frame(stratum = "Low", period = 2, treatment = "experimental", duration = 3, rate = .02),
    data.frame(stratum = "High", period = 1, treatment = "control", duration = 3, rate = .05),
    data.frame(stratum = "High", period = 1, treatment = "experimental", duration = 3, rate = .06),
    data.frame(stratum = "High", period = 2, treatment = "experimental", duration = 3, rate = .03)
  )
  fail_rate <- do.call(rbind, fail_rate)
  dropout_rate <- list(
    data.frame(stratum = "Low", period = 1, treatment = "control", duration = 3, rate = .001),
    data.frame(stratum = "Low", period = 1, treatment = "experimental", duration = 3, rate = .001),
    data.frame(stratum = "High", period = 1, treatment = "control", duration = 3, rate = .001),
    data.frame(stratum = "High", period = 1, treatment = "experimental", duration = 3, rate = .001)
  )
  dropout_rate <- do.call(rbind, dropout_rate)
  set.seed(12345)
  observed <- sim_pw_surv(
    n = 12,
    stratum = data.frame(stratum = c("Low", "High"), p = c(.3, .7)),
    fail_rate = fail_rate,
    dropout_rate = dropout_rate
  )
  expected <- readRDS("fixtures/backwards-compatibility/sim_pw_surv_ex4.rds")
  expect_equivalent(as.data.frame(observed), as.data.frame(expected))
})

test_that("sim_fixed_n()", {
  # Example 1
  # Show output structure
  set.seed(12345)
  observed <- sim_fixed_n(n = 2)
  expected <- readRDS("fixtures/backwards-compatibility/sim_fixed_n_ex1.rds")
  expect_equal(observed, expected)

  # Example 2
  # Example with 2 tests: logrank and FH(0,1)
  set.seed(12345)
  observed <- sim_fixed_n(n = 2, rho_gamma = data.frame(rho = 0, gamma = c(0, 1)))
  expected <- readRDS("fixtures/backwards-compatibility/sim_fixed_n_ex2.rds")
  expect_equal(observed, expected)

  # Example 3
  # Power by test
  # Only use cuts for events, events + min follow-up
  set.seed(12345)
  observed <- sim_fixed_n(
    n_sim = 2,
    timing_type = c(2, 5),
    rho_gamma = data.frame(rho = 0, gamma = c(0, 1))
  )
  expected <- readRDS("fixtures/backwards-compatibility/sim_fixed_n_ex3.rds")
  expect_equal(observed, expected)
})

test_that("mb_weight()", {
  # Use default enrollment and event rates at cut at 100 events
  # For transparency, may be good to set either `delay` or `w_max` to `Inf`
  set.seed(12345)
  x <- sim_pw_surv(n = 200)
  x <- cut_data_by_event(x, 125)
  x <- counting_process(x, arm = "experimental")

  # Example 1
  # Compute Magirr-Burman weights with `delay = 6`
  ZMB <- mb_weight(x, delay = 6, w_max = Inf)
  S <- with(ZMB, sum(o_minus_e * mb_weight))
  V <- with(ZMB, sum(var_o_minus_e * mb_weight^2))
  z <- S / sqrt(V)

  # Compute p-value of modestly weighted logrank of Magirr-Burman
  observed <- pnorm(z)
  expected <- readRDS("fixtures/backwards-compatibility/mb_weight_ex1.rds")
  expect_equal(observed, expected)

  # Example 2
  # Now compute with maximum weight of 2 as recommended in Magirr, 2021
  ZMB2 <- mb_weight(x, delay = Inf, w_max = 2)
  S <- with(ZMB2, sum(o_minus_e * mb_weight))
  V <- with(ZMB2, sum(var_o_minus_e * mb_weight^2))
  z <- S / sqrt(V)

  # Compute p-value of modestly weighted logrank of Magirr-Burman
  observed <- pnorm(z)
  expected <- readRDS("fixtures/backwards-compatibility/mb_weight_ex2.rds")
  expect_equal(observed, expected)
})

test_that("pvalue_maxcombo()", {
  # Example 1
  set.seed(12345)
  x <- sim_fixed_n(
    n_sim = 1,
    timing_type = 5,
    rho_gamma = data.frame(
      rho = c(0, 0, 1),
      gamma = c(0, 1, 1)
    )
  )
  observed <- pvalue_maxcombo(x)
  expected <- readRDS("fixtures/backwards-compatibility/pvalue_maxcombo_ex1.rds")
  expect_equal(observed, expected)

  # Example 2
  # Only use cuts for events, events + min follow-up
  set.seed(12345)
  xx <- sim_fixed_n(
    n_sim = 100,
    timing_type = 5,
    rho_gamma = data.frame(
      rho = c(0, 0, 1),
      gamma = c(0, 1, 1)
    )
  )

  # MaxCombo power estimate for cutoff at max of targeted events, minimum follow-up
  observed <- as.numeric(by(xx, xx$sim, pvalue_maxcombo))
  expected <- readRDS("fixtures/backwards-compatibility/pvalue_maxcombo_ex2.rds")
  expect_equal(observed, expected)
})

test_that("rpwexp()", {
  # Example 1
  # Exponential failure times
  observed <- rpwexp(
    n = 10000,
    fail_rate = data.frame(rate = 5, duration = 1)
  )
  expected <- readRDS("fixtures/backwards-compatibility/rpwexp_ex1.rds")
  expect_equal(observed, expected)

  # Example 2
  # Get 10k piecewise exponential failure times.
  # Failure rates are 1 for time 0 to 0.5, 3 for time 0.5 to 1, and 10 for > 1.
  # Intervals specifies duration of each failure rate interval
  # with the final interval running to infinity.
  observed <- rpwexp(
    n = 1e4,
    fail_rate = data.frame(rate = c(1, 3, 10), duration = c(.5, .5, 1))
  )
  expected <- readRDS("fixtures/backwards-compatibility/rpwexp_ex2.rds")
  expect_equal(observed, expected)
})
