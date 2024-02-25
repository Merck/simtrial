# Helper functions used by test-unvalidated-get_analysis_date.R

test_get_analysis_date <- function() {
  # Simulate trial data
  alpha <- 0.025
  ratio <- 3
  n <- 500
  info_frac <- c(0.7, 1)
  prevalence_ratio <- c(0.4, 0.6)
  study_duration <- 48
  # Two strata
  stratum <- c("Biomarker-positive", "Biomarker-negative")
  prevalence_ratio <- c(0.6, 0.4)
  # enrollment rate
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
    duration = 1000,
    fail_rate = c(log(2) / c(med_pos, med_pos, med_neg, med_neg)),
    hr = c(hr_pos, hr_neg),
    dropout_rate = 0.01
  )
  # Simulate data
  temp <- to_sim_pw_surv(fail_rate) # Convert the failure rate
  set.seed(2023)
  simulated_data <- sim_pw_surv(
    n = n, # Sample size
    # Stratified design with prevalence ratio of 6:4
    stratum = data.frame(stratum = stratum, p = prevalence_ratio),
    # Randomization ratio
    block = c("control", "control", "experimental", "experimental"),
    enroll_rate = enroll_rate, # Enrollment rate
    fail_rate = temp$fail_rate, # Failure rate
    dropout_rate = temp$dropout_rate # Dropout rate
  )

  list("n" = n, "simulated_data" = simulated_data)
}
