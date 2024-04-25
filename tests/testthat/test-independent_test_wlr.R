#### unstratified, FH (Fleming-Harrington) ----
# Check value when Fleming-Harrington weight is used
test_that("wlr() with FH weight on unstratified data", {
  # Example 1: Unstratified
  set.seed(123456)
  base <- sim_pw_surv(n = 200) |>
    cut_data_by_event(125) #|>
  output <- base |>
    wlr(weight = fh(rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1)))

  observed <- output$z
  base <- base |> counting_process(arm = "experimental")
  expected <- c()
  for (i in 1:length(observed)) {
    base <- base |> mutate(weight=s^(output$rho[i])*(1-s)^(output$gamma[i]))
    z <- sum(base$o_minus_e*base$weight)/sqrt(sum(base$weight^2*base$var_o_minus_e))
    expected <- c(expected,z)
  }
  expect_equal(observed, expected)
})


#### stratified, FH (Fleming-Harrington) ----
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
  enroll_rate$rate <- enroll_rate$rate * n / sum(enroll_rate$duration * enroll_rate$rate)  #??
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
  set.seed(123456)
  base <- sim_pw_surv(
    n = n, # Sample size
    # Stratified design with prevalence ratio of 6:4
    stratum = tibble(stratum = stratum, p = prevalence_ratio),
    # Randomization ratio
    block = c("control", "control", "experimental", "experimental"),
    enroll_rate = enroll_rate, # Enrollment rate
    fail_rate = temp$fail_rate, # Failure rate
    dropout_rate = temp$dropout_rate # Dropout rate
  ) |>
    cut_data_by_event(125)

  output <- base |>
    wlr(weight = fh(rho = c(0, 0, 1, 1), gamma = c(0, 1, 0, 1)))

  observed <- output$z
  base <- base |> counting_process(arm = "experimental")
  expected <- c()
  for (i in 1:length(observed)) {
    base <- base |> mutate(weight=s^(output$rho[i])*(1-s)^(output$gamma[i]))
    z <- sum(base$o_minus_e*base$weight)/sqrt(sum(base$weight^2*base$var_o_minus_e))
    expected <- c(expected,z)
  }
  expect_equal(observed, expected)
})


#### unstratified, MB (Magirr and Burman) ----
# Check value when Magirr and Burman weight is used
test_that("wlr() with MB weight on unstratified data", {
  # Example 1: Unstratified
  set.seed(123456)
  delay <- 4
  w_max <- 2
  base <- sim_pw_surv(n = 200) |>
    cut_data_by_event(125)
  output <- base |>
    wlr(weight = mb(delay = delay, w_max = w_max))

  observed <- output$z
  base <- base |> counting_process(arm = "experimental")
  base2 <- base |> filter(tte<=delay)
  expected <- c()
  for (i in 1:length(observed)) {
    wht <- base2 |> group_by(stratum) %>% summarise(mx = max(1/s)) |> mutate(mx = pmin(mx,w_max))
    base <- base |> full_join(wht, by=c('stratum')) |> mutate(weight=pmin(1/s,mx))
    z <- sum(base$o_minus_e*base$weight)/sqrt(sum(base$weight^2*base$var_o_minus_e))
    expected <- c(expected,z)
  }
  expect_equal(observed, expected)
  expect_equal(1+1, 2)
})


#### stratified, MB (Magirr and Burman) ----
# Check value when Magirr and Burman weight is used
test_that("wlr() with MB weight on stratified data", {
  # Example 2: Stratified
  set.seed(123456)
  delay <- 4
  w_max <- 2
  n <- 500
  # Two strata
  stratum <- c("Biomarker-positive", "Biomarker-negative")
  prevalence_ratio <- c(0.6, 0.4)
  enroll_rate <- gsDesign2::define_enroll_rate(
    stratum = rep(stratum, each = 2),
    duration = c(2, 10, 2, 10),
    rate = c(c(1, 4) * prevalence_ratio[1], c(1, 4) * prevalence_ratio[2])
  )
  enroll_rate$rate <- enroll_rate$rate * n / sum(enroll_rate$duration * enroll_rate$rate)  #??
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
  set.seed(123456)
  base <- sim_pw_surv(
    n = n, # Sample size
    # Stratified design with prevalence ratio of 6:4
    stratum = tibble(stratum = stratum, p = prevalence_ratio),
    # Randomization ratio
    block = c("control", "control", "experimental", "experimental"),
    enroll_rate = enroll_rate, # Enrollment rate
    fail_rate = temp$fail_rate, # Failure rate
    dropout_rate = temp$dropout_rate # Dropout rate
  ) |>
    cut_data_by_event(125)

  output <- base |>
    wlr(weight = mb(delay = delay, w_max = w_max))

  observed <- output$z
  base <- base |> counting_process(arm = "experimental")
  base2 <- base |> filter(tte<=delay)
  expected <- c()
  for (i in 1:length(observed)) {
    wht <- base2 |> group_by(stratum) %>% summarise(mx = max(1/s)) |> mutate(mx = pmin(mx,w_max))
    base <- base |> full_join(wht, by=c('stratum')) |> mutate(weight=pmin(1/s,mx))
    z <- sum(base$o_minus_e*base$weight)/sqrt(sum(base$weight^2*base$var_o_minus_e))
    expected <- c(expected,z)
  }
  expect_equal(observed, expected)
})


#### unstratified, early_zero_weight ----
# Check value when early_zero_weight is used
test_that("wlr() with early_zero_weight on unstratified data", {
  # Example 1: Unstratified
  set.seed(123456)
  early_period = 4
  base <- sim_pw_surv(n = 200) |>
    cut_data_by_event(125)
  output <- base |>
    wlr(weight = early_zero(early_period = early_period))

  observed <- output$z
  # WLR using early_zero_weight yields the same results as directly removing the events happening earlier than `early_period`
  base <- base |> counting_process(arm = "experimental") %>% filter(tte>=early_period)
  expected <- c()
  for (i in 1:length(observed)) {
    # base <- base |> mutate(weight=if_else(tte<early_period,0,1))
    z <- sum(base$o_minus_e)/sqrt(sum(base$var_o_minus_e))
    expected <- c(expected,z)
  }
  expect_equal(observed, expected)
})


#### stratified, early_zero_weight ----
# Check value when early_zero_weight is used
test_that("wlr() with early_zero_weight on stratified data", {
  # Example 2: Stratified
  set.seed(123456)
  early_period = 4
  n <- 500
  # Two strata
  stratum <- c("Biomarker-positive", "Biomarker-negative")
  prevalence_ratio <- c(0.6, 0.4)
  enroll_rate <- gsDesign2::define_enroll_rate(
    stratum = rep(stratum, each = 2),
    duration = c(2, 10, 2, 10),
    rate = c(c(1, 4) * prevalence_ratio[1], c(1, 4) * prevalence_ratio[2])
  )
  enroll_rate$rate <- enroll_rate$rate * n / sum(enroll_rate$duration * enroll_rate$rate)  #??
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
  set.seed(123456)
  x <- sim_pw_surv(
    n = n, # Sample size
    # Stratified design with prevalence ratio of 6:4
    stratum = tibble(stratum = stratum, p = prevalence_ratio),
    # Randomization ratio
    block = c("control", "control", "experimental", "experimental"),
    enroll_rate = enroll_rate, # Enrollment rate
    fail_rate = temp$fail_rate, # Failure rate
    dropout_rate = temp$dropout_rate # Dropout rate
  ) |>
    cut_data_by_event(125)

  output <- base |>
    wlr(weight = early_zero(early_period = early_period))

  observed <- output$z
  # WLR using early_zero_weight yields the same results as directly removing the events happening earlier than `early_period`
  base <- base |> counting_process(arm = "experimental") %>% filter(tte>=early_period)
  expected <- c()
  for (i in 1:length(observed)) {
    # base <- base |> mutate(weight=if_else(tte<early_period,0,1))
    z <- sum(base$o_minus_e)/sqrt(sum(base$var_o_minus_e))
    expected <- c(expected,z)
  }
  expect_equal(observed, expected)
})

