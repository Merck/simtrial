# Helper functions used by test-double_programming_simfix.R

test_simfix <- function() {
  # Study design using gsDesign
  alpha <- 0.025
  gamma <- c(5, 5, 47)
  R <- c(1, 1, 9)
  median.c <- 7
  hr <- 0.65
  dropout <- 0.05
  # Set power = 0.9 with target events = 227
  PE <- gsDesign::nSurv(
    alpha = alpha, beta = c(1 - 0.9), sided = 1, lambdaC = log(2) / median.c, hr = hr,
    eta = -log(1 - dropout) / 12, gamma = gamma, R = R, T = 18
  )
  # Set power = 0.93 with duration = 18
  PE <- gsDesign::nSurv(
    alpha = alpha, beta = c(1 - 0.93), sided = 1, lambdaC = log(2) / median.c, hr = hr,
    eta = -log(1 - dropout) / 12, gamma = gamma, R = R, T = 18
  )

  # Test for power comparing sim_fixed_n results with simple study design
  set.seed(1234)
  test2 <- sim_fixed_n(
    n_sim = 100,
    sample_size = 434,
    target_event = 227,
    stratum = data.frame(stratum = "All", p = 1),
    enroll_rate = data.frame(
      duration = c(1, 1, 9),
      rate = c(5, 5, 47)
    ),
    fail_rate = data.frame(
      stratum = "All",
      duration = c(100),
      fail_rate = log(2) / 7,
      hr = 0.65,
      dropout_rate = -log(1 - 0.05) / 12
    ),
    total_duration = 18,
    block = rep(c("experimental", "control"), 2),
    timing_type = 1:5,
    rho_gamma = data.frame(rho = 0, gamma = 0)
  )

  list("test2" = test2)
}
