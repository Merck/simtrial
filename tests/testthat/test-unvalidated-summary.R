test_that("summary.simtrial_gs_wlr() returns consistent results for one-sided design", {
  # Test code adapted from example in ?summary.summary.simtrial_gs_wlr

  # Parameters for enrollment
  enroll_rampup_duration <- 4 # Duration for enrollment ramp up
  enroll_duration <- 16 # Total enrollment duration
  enroll_rate <- gsDesign2::define_enroll_rate(
    duration = c(
      enroll_rampup_duration, enroll_duration - enroll_rampup_duration),
    rate = c(10, 30))

  # Parameters for treatment effect
  delay_effect_duration <- 3 # Delay treatment effect in months
  median_ctrl <- 9 # Survival median of the control arm
  median_exp <- c(9, 14) # Survival median of the experimental arm
  dropout_rate <- 0.001
  fail_rate <- gsDesign2::define_fail_rate(
    duration = c(delay_effect_duration, 100),
    fail_rate = log(2) / median_ctrl,
    hr = median_ctrl / median_exp,
    dropout_rate = dropout_rate)

  # Other related parameters
  alpha <- 0.025 # Type I error
  beta <- 0.1 # Type II error
  ratio <- 1 # Randomization ratio (experimental:control)

  # Build a one-sided group sequential design
  design <- gsDesign2::gs_design_ahr(
    enroll_rate = enroll_rate, fail_rate = fail_rate,
    ratio = ratio, alpha = alpha, beta = beta,
    analysis_time = c(12, 24, 36),
    upper = gsDesign2::gs_spending_bound,
    upar = list(sf = gsDesign::sfLDOF, total_spend = alpha),
    lower = gsDesign2::gs_b,
    lpar = rep(-Inf, 3))

  # Define cuttings of 2 IAs and 1 FA
  ia1_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[1]))
  ia2_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[2]))
  fa_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[3]))

  # Run simulations
  set.seed(1)
  simulation <- sim_gs_n(
    n_sim = 3,
    sample_size = ceiling(design$analysis$n[3]),
    enroll_rate = design$enroll_rate,
    fail_rate = design$fail_rate,
    test = wlr,
    cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
    weight = fh(rho = 0, gamma = 0.5))

  # Summarize simulations
  observed <- simulation |>
    summary(bound = gsDesign::gsDesign(k = 3, test.type = 1, sfu = gsDesign::sfLDOF)$upper$bound)
  expected <- data.frame(
    analysis = c(1, 2, 3),
    sim_n = c(369.3333333333333, 505, 505),
    sim_event = c(97, 305, 405),
    sim_time = c(12.877359569828519, 24.990283397668506, 37.20491262038222),
    sim_upper_prob = c(NA, 1, NA)
  ) |>
    structure(
      class = c("simtrial_gs_wlr", "data.frame"),
      compare_with_design = "no",
      method = "FH(rho=0, gamma=0.5)"
    )
  expect_equal(observed, expected)

  # Summarize simulation and compare with the planned design
  observed <- simulation |> summary(design = design)
  expected <- data.frame(
    analysis = c(1, 2, 3),
    asy_upper_prob = c(0.00014865936645545522, 0.5723215057363614, 0.9000000002116888),
    sim_upper_prob = c(NA, 1, NA),
    sim_event = c(97, 305, 405),
    sim_n = c(369.3333333333333, 505, 505),
    sim_time = c(12.877359569828519, 24.990283397668506, 37.20491262038222),
    asy_time = c(12, 24, 36),
    asy_n = c(353.04671034431556, 504.3524433490222, 504.3524433490222),
    asy_event = c(96.77457617908364, 304.00996193840484, 404.14196474655887)
  ) |>
    structure(
      class = c("simtrial_gs_wlr", "data.frame"),
      compare_with_design = "yes",
      design_type = "one-sided",
      method = "FH(rho=0, gamma=0.5)"
    )
  expect_equal(observed, expected, tolerance = 1e-6)
})

test_that("summary.simtrial_gs_wlr() returns consistent results for two-sided design", {
  # Parameters for enrollment
  enroll_rampup_duration <- 4 # Duration for enrollment ramp up
  enroll_duration <- 16 # Total enrollment duration
  enroll_rate <- gsDesign2::define_enroll_rate(
    duration = c(
      enroll_rampup_duration, enroll_duration - enroll_rampup_duration),
    rate = c(10, 30))

  # Parameters for treatment effect
  delay_effect_duration <- 3 # Delay treatment effect in months
  median_ctrl <- 9 # Survival median of the control arm
  median_exp <- c(9, 14) # Survival median of the experimental arm
  dropout_rate <- 0.001
  fail_rate <- gsDesign2::define_fail_rate(
    duration = c(delay_effect_duration, 100),
    fail_rate = log(2) / median_ctrl,
    hr = median_ctrl / median_exp,
    dropout_rate = dropout_rate)

  # Other related parameters
  alpha <- 0.025 # Type I error
  beta <- 0.1 # Type II error
  ratio <- 1 # Randomization ratio (experimental:control)

  # Build a two-sided group sequential design
  design <- gsDesign2::gs_design_ahr(
    enroll_rate = enroll_rate, fail_rate = fail_rate,
    ratio = ratio, alpha = alpha, beta = beta,
    analysis_time = c(12, 24, 36),
    upper = gsDesign2::gs_spending_bound,
    upar = list(sf = gsDesign::sfLDOF, total_spend = alpha),
    lower = gsDesign2::gs_spending_bound,
    lpar = list(sf = gsDesign::sfLDOF, total_spend = beta))

  # Define cuttings of 2 IAs and 1 FA
  ia1_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[1]))
  ia2_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[2]))
  fa_cut <- create_cut(target_event_overall = ceiling(design$analysis$event[3]))

  # Run simulations
  set.seed(1)
  simulation <- sim_gs_n(
    n_sim = 3,
    sample_size = ceiling(design$analysis$n[3]),
    enroll_rate = design$enroll_rate,
    fail_rate = design$fail_rate,
    test = wlr,
    cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
    weight = fh(rho = 0, gamma = 0.5))

  # Summarize simulations
  observed <- simulation |>
    summary(bound = gsDesign::gsDesign(k = 3, test.type = 1, sfu = gsDesign::sfLDOF)$upper$bound)
  expected <- data.frame(
    analysis = c(1, 2, 3),
    sim_n = c(366.6666666666667, 535, 535),
    sim_event = c(103, 323, 429),
    sim_time = c(12.363838412468121, 24.374413483785986, 36.116791896100885),
    sim_upper_prob = c(NA, 0.6666666666666666, 1)
  ) |>
    structure(
      compare_with_design = "no",
      class = c("simtrial_gs_wlr", "data.frame"),
      method = "FH(rho=0, gamma=0.5)"
    )
  expect_equal(observed, expected)

  # Summarize simulation and compare with the planned design
  observed <- simulation |> summary(design = design)
  expected <- data.frame(
    analysis = c(1, 2, 3),
    asy_upper_prob = c(0.00016250401737420353, 0.6011019363189855, 0.9000000001924918),
    asy_lower_prob = c(0.0007883883873094952, 0.05707064419933058, 0.10004018006137042),
    sim_upper_prob = c(NA, 0.6666666666666666, 1),
    sim_lower_prob = rep(NA_real_, 3L),
    sim_event = c(103, 323, 429),
    sim_n = c(366.6666666666667, 535, 535),
    sim_time = c(12.363838412468121, 24.374413483785986, 36.116791896100885),
    asy_time = c(12, 24, 36),
    asy_n = c(374.08958620608826, 534.4136945801262, 534.4136945801262),
    asy_event = c(102.54269505243633, 322.13006815203613, 428.2303047466704)
  ) |>
    structure(
      compare_with_design = "yes",
      design_type = "two-sided",
      class = c("simtrial_gs_wlr", "data.frame"),
      method = "FH(rho=0, gamma=0.5)"
    )
  expect_equal(observed, expected, tolerance = 1e-6)
})
