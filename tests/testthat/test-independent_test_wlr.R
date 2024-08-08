# Unstratified, FH (Fleming-Harrington) ----
# Check value when Fleming-Harrington weight is used
test_that("wlr() with FH weight on unstratified data", {
  # Example 1: Unstratified
  set.seed(123456)

  base <- sim_pw_surv(n = 200) |> cut_data_by_event(125)
  basec <- base |> counting_process(arm = "experimental")

  rho <- c(0, 0, 1, 1)
  gamma <- c(0, 1, 0, 1)
  observed <- c()
  expected <- c()
  for (i in 1:length(rho)) {
    output <- base |>
      wlr(weight = fh(rho = rho[i], gamma = gamma[i]))
    observed[i] <- output$z

    basec <- basec |> dplyr::mutate(weight = s^(rho[i]) * (1 - s)^(gamma[i]))
    z <- -sum(basec$o_minus_e * basec$weight) / sqrt(sum(basec$weight^2 * basec$var_o_minus_e))
    expected[i] <- z
  }

  expect_equal(observed, expected)
})

# Stratified, FH (Fleming-Harrington) ----
# Check value when Fleming-Harrington weight is used
test_that("wlr() with FH weight on stratified data", {
  # Example 1: Stratified
  set.seed(123456)

  n <- 500
  # Two strata
  stratum <- c("Biomarker-positive", "Biomarker-negative")
  prevalence_ratio <- c(0.6, 0.4)
  enroll_rate <- gsDesign2::define_enroll_rate(
    stratum = rep(stratum, each = 2),
    duration = c(2, 10, 2, 10),
    rate = c(c(1, 4) * prevalence_ratio[1], c(1, 4) * prevalence_ratio[2])
  )
  enroll_rate$rate <- enroll_rate$rate * n / sum(enroll_rate$duration * enroll_rate$rate)
  # Failure rate
  med_pos <- 10 # Median of the biomarker positive population
  med_neg <- 8 # Median of the biomarker negative population
  hr_pos <- c(1, 0.7) # Hazard ratio of the biomarker positive population
  hr_neg <- c(1, 0.8) # Hazard ratio of the biomarker negative population
  fail_rate <- gsDesign2::define_fail_rate(
    stratum = rep(stratum, each = 2),
    duration = c(3, 1000, 4, 1000),
    fail_rate = c(log(2) / c(med_pos, med_pos, med_neg, med_neg)),
    hr = c(hr_pos, hr_neg),
    dropout_rate = 0.01
  )
  temp <- to_sim_pw_surv(fail_rate) # Convert the failure rate
  base <- sim_pw_surv(
    n = n, # Sample size
    # Stratified design with prevalence ratio of 6:4
    stratum = data.frame(stratum = stratum, p = prevalence_ratio),
    # Randomization ratio
    block = c("control", "control", "experimental", "experimental"),
    enroll_rate = enroll_rate, # Enrollment rate
    fail_rate = temp$fail_rate, # Failure rate
    dropout_rate = temp$dropout_rate # Dropout rate
  ) |> cut_data_by_event(125)
  basec <- base |> counting_process(arm = "experimental")

  rho <- c(0, 0, 1, 1)
  gamma <- c(0, 1, 0, 1)
  observed <- c()
  expected <- c()
  for (i in 1:length(rho)) {
    output <- base |> wlr(weight = fh(rho = rho[i], gamma = gamma[i]))
    observed[i] <- output$z

    basec <- basec |> dplyr::mutate(weight = s^(rho[i]) * (1 - s)^(gamma[i]))
    z <- -sum(basec$o_minus_e * basec$weight) / sqrt(sum(basec$weight^2 * basec$var_o_minus_e))
    expected[i] <- z
  }

  expect_equal(observed, expected)
})

# Unstratified, MB (Magirr and Burman) ----
# Check value when Magirr and Burman weight is used
test_that("wlr() with MB weight on unstratified data", {
  # Example 1: Unstratified
  set.seed(123456)

  base <- sim_pw_surv(n = 200) |> cut_data_by_event(125)
  basec <- base |> counting_process(arm = "experimental")

  delay <- c(4, 4, 7, 7)
  w_max <- c(2, 3, 2, 3)
  observed <- c()
  expected <- c()
  for (i in 1:length(delay)) {
    output <- base |> wlr(weight = mb(delay = delay[i], w_max = w_max[i]))
    observed[i] <- output$z

    wht <- basec |>
      dplyr::filter(tte <= delay[i]) |>
      dplyr::group_by(stratum) |>
      dplyr::summarise(mx = max(1 / s)) |>
      dplyr::mutate(mx = pmin(mx, w_max[i]))
    tmp <- basec |>
      dplyr::full_join(wht, by = c("stratum")) |>
      dplyr::mutate(weight = pmin(1 / s, mx))
    z <- -sum(tmp$o_minus_e * tmp$weight) / sqrt(sum(tmp$weight^2 * tmp$var_o_minus_e))
    expected[i] <- z
  }

  expect_equal(observed, expected)
})

# Stratified, MB (Magirr and Burman) ----
# Check value when Magirr and Burman weight is used
test_that("wlr() with MB weight on stratified data", {
  # Example 2: Stratified
  set.seed(123456)

  n <- 500
  # Two strata
  stratum <- c("Biomarker-positive", "Biomarker-negative")
  prevalence_ratio <- c(0.6, 0.4)
  enroll_rate <- gsDesign2::define_enroll_rate(
    stratum = rep(stratum, each = 2),
    duration = c(2, 10, 2, 10),
    rate = c(c(1, 4) * prevalence_ratio[1], c(1, 4) * prevalence_ratio[2])
  )
  enroll_rate$rate <- enroll_rate$rate * n / sum(enroll_rate$duration * enroll_rate$rate)
  # Failure rate
  med_pos <- 10 # Median of the biomarker positive population
  med_neg <- 8 # Median of the biomarker negative population
  hr_pos <- c(1, 0.7) # Hazard ratio of the biomarker positive population
  hr_neg <- c(1, 0.8) # Hazard ratio of the biomarker negative population
  fail_rate <- gsDesign2::define_fail_rate(
    stratum = rep(stratum, each = 2),
    duration = c(3, 1000, 4, 1000),
    fail_rate = c(log(2) / c(med_pos, med_pos, med_neg, med_neg)),
    hr = c(hr_pos, hr_neg),
    dropout_rate = 0.01
  )
  temp <- to_sim_pw_surv(fail_rate) # Convert the failure rate
  base <- sim_pw_surv(
    n = n, # Sample size
    # Stratified design with prevalence ratio of 6:4
    stratum = data.frame(stratum = stratum, p = prevalence_ratio),
    # Randomization ratio
    block = c("control", "control", "experimental", "experimental"),
    enroll_rate = enroll_rate, # Enrollment rate
    fail_rate = temp$fail_rate, # Failure rate
    dropout_rate = temp$dropout_rate # Dropout rate
  ) |> cut_data_by_event(125)
  basec <- base |> counting_process(arm = "experimental")

  delay <- c(4, 4, 7, 7)
  w_max <- c(2, 3, 2, 3)
  observed <- c()
  expected <- c()
  for (i in 1:length(delay)) {
    output <- base |> wlr(weight = mb(delay = delay[i], w_max = w_max[i]))
    observed[i] <- output$z

    wht <- basec |>
      dplyr::filter(tte <= delay[i]) |>
      dplyr::group_by(stratum) |>
      dplyr::summarise(mx = max(1 / s)) |>
      dplyr::mutate(mx = pmin(mx, w_max[i]))
    tmp <- basec |>
      dplyr::full_join(wht, by = c("stratum")) |>
      dplyr::mutate(weight = pmin(1 / s, mx))
    z <- -sum(tmp$o_minus_e * tmp$weight) / sqrt(sum(tmp$weight^2 * tmp$var_o_minus_e))
    expected[i] <- z
  }

  expect_equal(observed, expected)
})

# Unstratified, early_zero_weight ----
# Check value when early_zero_weight is used
test_that("wlr() with early_zero_weight on unstratified data", {
  # Example 1: Unstratified
  set.seed(123456)

  base <- sim_pw_surv(n = 200) |> cut_data_by_event(125)
  basec <- base |> counting_process(arm = "experimental")

  early_period <- c(2, 4, 6)
  observed <- c()
  expected <- c()
  for (i in 1:length(early_period)) {
    output <- base |> wlr(weight = early_zero(early_period = early_period[i]))
    observed[i] <- output$z

    # WLR using early_zero_weight yields the same results as directly removing the events happening earlier than `early_period`
    tmp <- basec |> dplyr::filter(tte >= early_period[i])
    # tmp <- basec |> mutate(weight=if_else(tte<early_period,0,1))
    z <- -sum(tmp$o_minus_e) / sqrt(sum(tmp$var_o_minus_e))
    expected <- c(expected, z)
  }
  expect_equal(observed, expected)
})

# Stratified, early_zero_weight ----
# Check value when early_zero_weight is used
test_that("wlr() with early_zero_weight on stratified data", {
  # Example 2: Stratified
  set.seed(123456)

  n <- 500
  # Two strata
  stratum <- c("Biomarker-positive", "Biomarker-negative")
  prevalence_ratio <- c(0.6, 0.4)
  enroll_rate <- gsDesign2::define_enroll_rate(
    stratum = rep(stratum, each = 2),
    duration = c(2, 10, 2, 10),
    rate = c(c(1, 4) * prevalence_ratio[1], c(1, 4) * prevalence_ratio[2])
  )
  enroll_rate$rate <- enroll_rate$rate * n / sum(enroll_rate$duration * enroll_rate$rate)
  # Failure rate
  med_pos <- 10 # Median of the biomarker positive population
  med_neg <- 8 # Median of the biomarker negative population
  hr_pos <- c(1, 0.7) # Hazard ratio of the biomarker positive population
  hr_neg <- c(1, 0.8) # Hazard ratio of the biomarker negative population
  fail_rate <- gsDesign2::define_fail_rate(
    stratum = rep(stratum, each = 2),
    duration = c(3, 1000, 4, 1000),
    fail_rate = c(log(2) / c(med_pos, med_pos, med_neg, med_neg)),
    hr = c(hr_pos, hr_neg),
    dropout_rate = 0.01
  )
  temp <- to_sim_pw_surv(fail_rate) # Convert the failure rate
  base <- sim_pw_surv(
    n = n, # Sample size
    # Stratified design with prevalence ratio of 6:4
    stratum = data.frame(stratum = stratum, p = prevalence_ratio),
    # Randomization ratio
    block = c("control", "control", "experimental", "experimental"),
    enroll_rate = enroll_rate, # Enrollment rate
    fail_rate = temp$fail_rate, # Failure rate
    dropout_rate = temp$dropout_rate # Dropout rate
  ) |> cut_data_by_event(125)
  basec <- base |> counting_process(arm = "experimental")

  early_period <- 2 # Except being the input, not actually used
  output <- base |> wlr(weight = early_zero(early_period = early_period, fail_rate = fail_rate))
  observed <- output$z

  tmp <- basec |> dplyr::mutate(
    weight = dplyr::if_else(
      stratum == "Biomarker-negative",
      dplyr::if_else(tte < 4, 0, log(0.8)),
      dplyr::if_else(tte < 3, 0, log(0.7))
    )
  )
  z <- -sum(tmp$o_minus_e * tmp$weight) / sqrt(sum(tmp$weight^2 * tmp$var_o_minus_e))
  expected <- z

  expect_equal(observed, expected)
})
