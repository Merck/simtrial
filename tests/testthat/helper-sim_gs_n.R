# Helper functions used by test-unvalidated-sim_gs_n.R

test_enroll_rate <- function() {
  # parameters for enrollment
  enroll_rampup_duration <- 4 # duration for enrollment ramp up
  enroll_duration <- 16       # total enrollment duration
  enroll_rate <- gsDesign2::define_enroll_rate(
    duration = c(enroll_rampup_duration,
                 enroll_duration - enroll_rampup_duration),
    rate = c(10, 30)
  )
  return(enroll_rate)
}

test_fail_rate <- function() {
  # parameters for treatment effect
  delay_effect_duration <- 3  # delay treatment effect in months
  median_col <- 9             # survival median of the control arm
  median_exp <- c(9, 14)      # survival median of the experimental arm
  dropout_rate <- 0.001
  fail_rate <- gsDesign2::define_fail_rate(
    duration = c(delay_effect_duration, 100),
    fail_rate = log(2) /  median_col,
    hr = median_col / median_exp,
    dropout_rate = dropout_rate
  )
  return(fail_rate)
}

test_cutting <- function() {
  # other related parameters
  alpha <- 0.025              # type I error
  beta <- 0.1                 # type II error
  ratio <- 1                  # randomization ratio (exp:col)
  # Define cuttings of 2 IAs and 1 FA
  # IA1
  # The 1st interim analysis will occur at the later of the following 3 conditions:
  # - At least 20 months have passed since the start of the study
  # - At least 100 events have occurred
  # - At least 20 months have elapsed after enrolling 200/400 subjects, with a
  #   minimum of 20 months follow-up
  # However, if events accumulation is slow, we will wait for a maximum of 24 months.
  ia1 <- create_cutting(
    planned_calendar_time = 20,
    target_event_overall = 100,
    max_extension_for_target_event = 24,
    min_n_overall = 200,
    min_followup = 20
  )
  # IA2
  # The 2nd interim analysis will occur at the later of the following 3 conditions:
  # - At least 32 months have passed since the start of the study
  # - At least 250 events have occurred
  # - At least 10 months after IA1
  # However, if events accumulation is slow, we will wait for a maximum of 34 months.
  ia2 <- create_cutting(
    planned_calendar_time = 32,
    target_event_overall = 200,
    max_extension_for_target_event = 34,
    min_time_after_previous_analysis = 10
  )
  # FA
  # The final analysis will occur at the later of the following 2 conditions:
  # - At least 45 months have passed since the start of the study
  # - At least 300 events have occurred
  fa <- create_cutting(
    planned_calendar_time = 45,
    target_event_overall = 350
  )

  return(list(ia1 = ia1, ia2 = ia2, fa = fa))
}
